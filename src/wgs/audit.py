"""Live re-verification of a release against the public sources it cites (ADR-018).

Every finding in a release was computed from a dated snapshot of public knowledge. `wgs audit` asks each
source's live API whether it still says what we say it says, and writes a dated report next to the releases:

* ClinVar (E-utilities esummary): classification class and review stars for every cited record.
* GWAS Catalog (REST v2): the exact study + SNP association: effect allele, effect size, p-value, PubMed ID.
* PGS Catalog (REST): score exists, reported trait, publication PubMed ID.
* PubMed (E-utilities esummary): every cited PubMed ID exists; stored titles match.
* Links: every distinct non-per-variant URL, plus a sample of per-record URLs, answers without an error.

A difference is not automatically a bug: sources are revised. Each one is labelled `changed upstream` when
the live record was re-evaluated after our snapshot, otherwise `mismatch` (something we must explain or fix).
"""

from __future__ import annotations

import contextlib
import html
import json
import math
import random
import re
import time
import urllib.error
import urllib.request
from collections.abc import Callable, Iterable
from datetime import date
from pathlib import Path

import duckdb

from .config import Config
from .knowledge.clinvar import STARS
from .knowledge.http import UA
from .pipeline import console, now, write_json

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
GWAS = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"
PGS = "https://www.pgscatalog.org/rest/score/"


def _get(url: str, tries: int = 3, timeout: int = 30) -> bytes:
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            # NCBI E-utilities answers bursts with a transient 400 as well as 429
            retry = e.code >= 500 or e.code == 429 or (e.code == 400 and "eutils.ncbi" in url)
            if not retry or i == tries - 1:
                raise
        except (urllib.error.URLError, TimeoutError):
            if i == tries - 1:
                raise
        time.sleep(1.5 * (i + 1))
    raise RuntimeError("unreachable")


def _json(url: str) -> dict:
    return json.loads(_get(url))


def sig_class(s: str | None) -> str:
    """Python twin of knowledge.clinvar.sig_class_sql, applied to the first term of a ClinVar classification."""
    if not s:
        return "other"
    s = re.split(r"[|;,]", s.strip().lower().replace(" ", "_"))[0]
    for prefix, cls in (("conflicting", "conflicting"), ("pathogenic", "pathogenic"),
                        ("likely_pathogenic", "likely_pathogenic"), ("uncertain", "uncertain"), ("vus", "uncertain"),
                        ("benign", "benign"), ("likely_benign", "likely_benign")):
        if s.startswith(prefix):
            return cls
    if s == "drug_response":
        return "drug_response"
    if s in ("risk_factor", "likely_risk_allele", "established_risk_allele"):
        return "risk_factor"
    if s == "protective":
        return "protective"
    if s in ("association", "affects", "confers_sensitivity"):
        return "association"
    return "other"


def stars(review_status: str | None) -> int:
    return STARS.get((review_status or "").strip().lower().replace(" ", "_"), 0)


def _batches(xs: list, n: int) -> Iterable[list]:
    for i in range(0, len(xs), n):
        yield xs[i:i + n]


def _release(cfg: Config, release: str | None) -> tuple[dict, dict, duckdb.DuckDBPyConnection]:
    rd = cfg.releases_dir
    index = json.loads((rd / "index.json").read_text())
    rel = next((r for r in index["releases"] if r["id"] == release), None) if release else index["releases"][0]
    if rel is None:
        raise ValueError(f"no release {release}")
    m = json.loads((rd / rel["manifest"]).read_text())
    con = duckdb.connect()
    for name, t in m["tables"].items():
        con.execute(f"CREATE VIEW {name} AS SELECT * FROM read_parquet('{(rd / t['path']).as_posix()}')")
    return rel, m, con


def _rows(con, sql: str) -> list[dict]:
    rel = con.sql(sql)
    return [dict(zip(rel.columns, r, strict=True)) for r in rel.fetchall()]


def _sources(con, m: dict) -> list[dict]:
    """Every (table, finding, source) triple that cites a source record."""
    out = []
    for table, key in (("claims", "claim_id"), ("brain_variants", "gene || ':' || pos")):
        if table not in m["tables"]:
            continue
        cols = {c[0] for c in con.execute(f"DESCRIBE {table}").fetchall()}
        if "sources" not in cols:
            continue
        for r in _rows(con, f"SELECT {key} AS id, sources FROM {table}"):
            out.extend({"table": table, "id": r["id"], **s} for s in json.loads(r["sources"] or "[]"))
    return out


def check_clinvar(con, m: dict, snapshot: str | None) -> list[dict]:
    stored: dict[str, dict] = {}
    for table in ("claims", "brain_variants"):
        if table not in m["tables"]:
            continue
        cols = {c[0] for c in con.execute(f"DESCRIBE {table}").fetchall()}
        if "clinvar_id" in cols:
            q = (f"SELECT CAST(clinvar_id AS VARCHAR) AS id, clinvar_significance AS sig, clinvar_stars AS stars "
                 f"FROM {table} WHERE clinvar_id IS NOT NULL")
        elif "sources" in cols:
            q = (f"SELECT json_extract_string(s, '$.record') AS id, clinvar_significance AS sig, "
                 f"clinvar_stars AS stars "
                 f"FROM (SELECT unnest(from_json(sources, '[\"json\"]')) AS s, * FROM {table}) "
                 f"WHERE json_extract_string(s, '$.source') = 'clinvar' AND clinvar_significance IS NOT NULL")
        else:
            continue
        for r in _rows(con, q):
            stored.setdefault(r["id"], r)
    out = []
    ids = sorted(stored)
    for batch in _batches(ids, 150):
        res = _json(f"{EUTILS}?db=clinvar&retmode=json&id={','.join(batch)}").get("result", {})
        for vid in batch:
            mine, live = stored[vid], res.get(vid)
            row = {"source": "clinvar", "record": vid, "url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{vid}/",
                   "stored": f"{sig_class(mine['sig'])} {mine['stars']}★"}
            if not live or "error" in live:
                out.append({**row, "status": "missing", "live": None, "detail": "record not returned by ClinVar"})
                continue
            g = live.get("germline_classification") or {}
            lc, ls = sig_class(g.get("description")), stars(g.get("review_status"))
            row["live"] = f"{lc} {ls}★ ({g.get('description')}; evaluated {(g.get('last_evaluated') or '')[:10]})"
            if lc == sig_class(mine["sig"]) and ls == int(mine["stars"] or 0):
                out.append({**row, "status": "ok"})
            else:
                evaluated = (g.get("last_evaluated") or "")[:10].replace("/", "-")
                newer = bool(snapshot and evaluated and evaluated > snapshot)
                out.append({**row, "status": "changed upstream" if newer else "mismatch",
                            "detail": "ClinVar re-evaluated after our snapshot; `wgs knowledge refresh` picks it up"
                            if newer else "differs from the live record"})
        time.sleep(0.4)
    return out


def _num(s) -> float | None:
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", str(s or ""))
    return float(m.group()) if m else None


def check_gwas(con, m: dict) -> list[dict]:
    if "trait_snps" not in m["tables"]:
        return []
    out = []
    for r in _rows(con, "SELECT rsid, study, effect_allele, or_beta, mlog10p, pmid, effect_kind FROM trait_snps "
                        "WHERE study IS NOT NULL"):
        url = f"https://www.ebi.ac.uk/gwas/studies/{r['study']}"
        row = {"source": "gwas_catalog", "record": f"{r['study']} {r['rsid']}", "url": url,
               "stored": f"{r['rsid']}-{r['effect_allele']} {r['effect_kind']} {r['or_beta']} p=1e-{r['mlog10p']:.1f} "
                         f"PMID {r['pmid']}"}
        try:
            assoc = _json(f"{GWAS}?rs_id={r['rsid']}&accession_id={r['study']}&size=50")
        except urllib.error.HTTPError as e:
            out.append({**row, "status": "error", "live": None, "detail": f"HTTP {e.code}"})
            continue
        found = (assoc.get("_embedded") or {}).get("associations") or []
        best, why = None, "no association for this SNP in this study"
        for a in found:
            alleles = [x.split("-", 1)[-1] for x in a.get("snp_effect_allele") or []]
            size = (a.get("or_value") if r["effect_kind"] == "or" else a.get("beta") if r["effect_kind"] == "beta"
                    else a.get("or_value") or a.get("beta"))
            if r["effect_allele"] not in alleles and not (r["effect_allele"] is None and alleles in (["?"], [])):
                why = f"effect allele now {'/'.join(alleles)}"
                continue
            v = _num(size)
            if v is None or r["or_beta"] is None or not math.isclose(v, r["or_beta"], rel_tol=0.01):
                why = f"effect size now {size}"
                continue
            p = a.get("p_value")
            if p and r["mlog10p"] and p > 0 and abs(-math.log10(p) - r["mlog10p"]) > 0.5:
                why = f"p-value now {p}"
                continue
            best = a
            break
        if best is None:
            out.append({**row, "status": "mismatch", "live": f"{len(found)} association(s)", "detail": why})
        else:
            pm = str(best.get("pubmed_id") or "")
            ok = pm == str(r["pmid"])
            out.append({**row, "status": "ok" if ok else "mismatch",
                        "live": f"{best.get('snp_effect_allele')} {best.get('or_value') or best.get('beta')} "
                                f"p={best.get('p_value')} PMID {pm}",
                        **({} if ok else {"detail": "PubMed ID differs"})})
        time.sleep(0.2)
    return out


def check_pgs(con, m: dict) -> list[dict]:
    if "pgs_scores" not in m["tables"]:
        return []
    out = []
    for r in _rows(con, "SELECT DISTINCT pgs_id, trait_reported, pmid FROM pgs_scores ORDER BY pgs_id"):
        url = f"https://www.pgscatalog.org/score/{r['pgs_id']}/"
        row = {"source": "pgs_catalog", "record": r["pgs_id"], "url": url,
               "stored": f"{r['trait_reported']} · PMID {r['pmid']}"}
        try:
            live = _json(PGS + r["pgs_id"])
        except urllib.error.HTTPError as e:
            out.append({**row, "status": "missing" if e.code == 404 else "error", "live": None,
                        "detail": f"HTTP {e.code}"})
            continue
        lp = str((live.get("publication") or {}).get("PMID") or "")
        row["live"] = f"{live.get('trait_reported')} · PMID {lp}"
        bad = []
        if (live.get("trait_reported") or "").strip().lower() != (r["trait_reported"] or "").strip().lower():
            bad.append("reported trait differs")
        if r["pmid"] and lp != str(r["pmid"]).split(".")[0]:
            bad.append("PubMed ID differs")
        out.append({**row, "status": "mismatch" if bad else "ok", **({"detail": "; ".join(bad)} if bad else {})})
        time.sleep(0.15)
    return out


def _norm_title(t: str | None) -> str:
    return re.sub(r"[^a-z0-9]", "", (t or "").lower())


def check_pubmed(con, m: dict) -> list[dict]:
    cited: dict[str, str | None] = {}
    if "pgs_scores" in m["tables"]:
        for r in _rows(con, "SELECT DISTINCT CAST(pmid AS VARCHAR) AS pmid, title FROM pgs_scores "
                            "WHERE pmid IS NOT NULL UNION SELECT DISTINCT CAST(eval_pmid AS VARCHAR), NULL "
                            "FROM pgs_scores WHERE eval_pmid IS NOT NULL"):
            cited.setdefault(str(r["pmid"]).split(".")[0], r["title"])
    if "trait_snps" in m["tables"]:
        for r in _rows(con, "SELECT DISTINCT CAST(pmid AS VARCHAR) AS pmid FROM trait_snps WHERE pmid IS NOT NULL"):
            cited.setdefault(r["pmid"], None)
    out = []
    ids = sorted(p for p in cited if p.isdigit())
    for batch in _batches(ids, 150):
        res = _json(f"{EUTILS}?db=pubmed&retmode=json&id={','.join(batch)}").get("result", {})
        for p in batch:
            live = res.get(p) or {}
            row = {"source": "pubmed", "record": p, "url": f"https://pubmed.ncbi.nlm.nih.gov/{p}/",
                   "stored": cited[p] or "(cited by ID)"}
            if not live or "error" in live:
                out.append({**row, "status": "missing", "live": None, "detail": "PubMed has no such record"})
                continue
            row["live"] = f"{live.get('title')} ({live.get('source')}, {live.get('pubdate')})"
            t = cited[p]
            same = not t or _norm_title(t)[:60] == _norm_title(live.get("title"))[:60]
            out.append({**row, "status": "ok" if same else "mismatch", **({} if same else {"detail": "title differs"})})
        time.sleep(0.4)
    return out


# Per-record pages are numerous and already verified through the APIs above; only a sample is link-checked.
PER_RECORD = re.compile(r"clinvar/variation/|gnomad\.broadinstitute\.org/variant|ensembl\.org/.*/Variation|"
                        r"pubmed\.ncbi|gene\.sfari\.org/database/human-gene/")


def check_links(con, m: dict, sample: int = 40, seed: int = 7) -> list[dict]:
    urls: set[str] = {s["url"] for s in _sources(con, m) if s.get("url")}
    for table, col in (("pgx_drugs", "url"), ("pgs_scores", "trait_url")):
        if table in m["tables"]:
            urls |= {r["u"] for r in _rows(con, f"SELECT DISTINCT {col} AS u FROM {table} WHERE {col} LIKE 'http%'")}
    fixed = sorted(u for u in urls if not PER_RECORD.search(u))
    per = sorted(u for u in urls if PER_RECORD.search(u))
    picked = fixed + random.Random(seed).sample(per, min(sample, len(per)))
    out = []
    for u in picked:
        row = _probe(u)
        if row["status"] in ("broken", "unavailable"):  # often transient: try once more after a pause
            time.sleep(5)
            row = _probe(u)
        out.append(row)
        time.sleep(0.1)
    return out


def _probe(u: str) -> dict:
    try:
        req = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0 " + UA, "Range": "bytes=0-2047"})
        with urllib.request.urlopen(req, timeout=30) as r:
            return {"source": "link", "record": u, "url": u, "status": "ok", "live": f"HTTP {r.status}"}
    except urllib.error.HTTPError as e:
        # Some sites reject scripted clients (403/429) while serving browsers; that is not a dead link.
        # A server error (5xx) means the site is down right now, not that the page is gone.
        status = "blocked" if e.code in (401, 403, 429) else "unavailable" if e.code >= 500 else "broken"
        return {"source": "link", "record": u, "url": u, "status": status, "live": f"HTTP {e.code}"}
    except TimeoutError as e:
        return {"source": "link", "record": u, "url": u, "status": "unavailable", "live": type(e).__name__}
    except Exception as e:  # noqa: BLE001 - network errors of any kind are reported, not raised
        if isinstance(getattr(e, "reason", None), TimeoutError):
            return {"source": "link", "record": u, "url": u, "status": "unavailable", "live": "TimeoutError"}
        return {"source": "link", "record": u, "url": u, "status": "broken", "live": type(e).__name__}


def _page_text(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 " + UA})
    with urllib.request.urlopen(req, timeout=45) as r:
        html = r.read().decode("utf-8", "replace")
    return re.sub(r"<[^>]+>", " ", html)


def quote_numbers(quote: str) -> list[str]:
    """The numbers a curated quote relies on, without thousands separators: '1 in 8' -> ['1', '8']."""
    return [n.replace(",", "") for n in re.findall(r"\d[\d,]*(?:\.\d+)?", quote or "")]


def numbers_present(quote: str, text: str) -> list[str]:
    """Numbers from the quote that do NOT appear in the source text (empty = verified)."""
    text = html.unescape(text).replace("\u00b7", ".").replace("\u2009", " ").replace("\u00a0", " ")
    have = set(quote_numbers(re.sub(r"(?<=\d) (?=\d{3}\b)", "", text)))
    return [n for n in quote_numbers(quote) if n not in have]


def check_curated(con, m: dict) -> list[dict]:
    """Hand-curated baseline risks and per-SD effects: the source is live and still states the quoted numbers."""
    from . import absrisk

    t = absrisk.table()
    out = []
    for kind in ("baseline", "per_sd"):
        for r in t.get(kind, []):
            rec = f"{kind}: {r.get('disease') or r.get('pgs_id')} ({r.get('sex', '')})".replace(" ()", "")
            row = {"source": "curated", "record": rec, "url": r.get("url"), "stored": r.get("quote")}
            text = ""
            try:
                if r.get("pmid"):
                    text += _get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id="
                                 f"{r['pmid']}&rettype=abstract&retmode=text").decode("utf-8", "replace")
                if r.get("url"):
                    try:
                        text += _page_text(r["url"])
                    except urllib.error.HTTPError as e:
                        if not text:
                            raise
                        row["detail"] = f"page HTTP {e.code}; checked the PubMed abstract"
            except urllib.error.HTTPError as e:
                out.append({**row, "status": "blocked" if e.code in (401, 403, 429) else "broken",
                            "live": f"HTTP {e.code}"})
                continue
            except Exception as e:  # noqa: BLE001
                out.append({**row, "status": "broken", "live": type(e).__name__})
                continue
            missing = numbers_present(r.get("quote") or "", text)
            for _ in range(3):  # some sites (PMC) intermittently serve a short stub page instead of the article
                if not missing or not r.get("url"):
                    break
                time.sleep(3)
                with contextlib.suppress(Exception):
                    missing = numbers_present(r.get("quote") or "", text + _page_text(r["url"]))
            out.append({**row, "status": "mismatch" if missing else "ok",
                        "live": f"numbers not found on the page: {', '.join(missing)}" if missing else "quote found"})
            time.sleep(0.34)
    return out


CHECKS: dict[str, Callable] = {"clinvar": check_clinvar, "gwas": check_gwas, "pgs": check_pgs,
                               "pubmed": check_pubmed, "links": check_links, "curated": check_curated}


def audits_dir(cfg: Config) -> Path:
    return cfg.workspace / "audits"


def latest(cfg: Config) -> dict | None:
    d = audits_dir(cfg)
    reports = sorted(d.glob("*/summary.json")) if d.exists() else []
    return json.loads(reports[-1].read_text()) if reports else None


def run(cfg: Config, release: str | None = None, only: list[str] | None = None, link_sample: int = 40) -> dict:
    rel, m, con = _release(cfg, release)
    snapshot = (m.get("knowledge") or {}).get("clinvar")
    if isinstance(snapshot, dict):
        snapshot = snapshot.get("version")
    results: dict[str, list[dict]] = {}
    for name, fn in CHECKS.items():
        if only and name not in only:
            continue
        console.print(f"[bold]audit[/] {name} …")
        try:
            results[name] = (fn(con, m, snapshot) if name == "clinvar" else
                             fn(con, m, link_sample) if name == "links" else fn(con, m))
        except Exception as e:  # noqa: BLE001 - one unreachable source must not hide the others
            results[name] = [{"source": name, "status": "error", "detail": f"{type(e).__name__}: {e}"}]
    out = audits_dir(cfg) / f"{date.today().isoformat()}_release-{rel['id']}"
    if only and (out / "report.json").exists():  # a partial re-run updates today's full audit
        results = {**json.loads((out / "report.json").read_text()), **results}
    counts = {k: {s: sum(r["status"] == s for r in v) for s in sorted({r["status"] for r in v})}
              for k, v in results.items()}
    summary = {"release": rel["id"], "audited_at": now(), "clinvar_snapshot": snapshot, "counts": counts,
               "problems": [r for v in results.values() for r in v
                            if r["status"] in ("mismatch", "missing", "broken", "error")]}
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / "report.json", results)
    write_json(out / "summary.json", summary)
    # Beside the immutable release objects, so the dashboard can show how this release fared against live sources.
    (cfg.releases_dir / "audits").mkdir(exist_ok=True)
    write_json(cfg.releases_dir / "audits" / f"{rel['id']}.json", summary)
    console.print(f"Audit written to {out}")
    for k, c in counts.items():
        console.print(f"  {k}: " + ", ".join(f"{n} {s}" for s, n in c.items()))
    return summary
