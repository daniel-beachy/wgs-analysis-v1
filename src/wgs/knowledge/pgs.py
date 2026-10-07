"""Polygenic scores: the PGS Catalog (scores + independent evaluations) and its 1000 Genomes reference panel.

Which score represents a trait is decided by a rule, not by hand (ADR-017):

1. only scores mapped to exactly one trait and without interaction terms;
2. rank by the number of *independent* publications that evaluated the score (not the developers' own),
   then by how many people those independent evaluations included, then by the size of all evaluations,
   then by recency.

``data/pgs_featured.csv`` lists the traits the dashboard features, which section and body part each belongs to,
and — rarely — a pinned score with the reason the rule's pick was overridden. Every other trait whose best score
has at least one independent evaluation is offered in the "all well-tested scores" list.
"""

from __future__ import annotations

import csv
import hashlib
import re
import shutil
import tarfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from . import Built, Source, Upstream, http

FTP = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs"
META = f"{FTP}/metadata/pgs_all_metadata.tar.gz"
FEATURED = Path(__file__).parent / "data" / "pgs_featured.csv"
REF_URL = f"{FTP}/resources/pgsc_1000G_v1.tar.zst"
ALTERNATES = 2  # next-ranked scores kept per featured trait, in case the pick can't be matched


def scoring_url(pgs_id: str, build: str = "GRCh37") -> str:
    return f"{FTP}/scores/{pgs_id}/ScoringFiles/Harmonized/{pgs_id}_hmPOS_{build}.txt.gz"


def _num(s: str | None) -> float | None:
    """'1.71 [1.68, 1.74]' -> 1.71 ; '' -> None."""
    m = re.match(r"\s*(-?\d+(?:\.\d+)?(?:[eE]-?\d+)?)", s or "")
    return float(m.group(1)) if m else None


def _int(s: str | None) -> int:
    try:
        return int(float(s or 0))
    except ValueError:
        return 0


def _csv(p: Path) -> list[dict]:
    with open(p, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def select(meta: Path, featured: list[dict]) -> list[dict]:
    """Apply the selection rule. Returns one row per trait — featured traits first (each followed by up to
    ALTERNATES next-ranked stand-ins), then the well-tested rest."""
    scores = {r["Polygenic Score (PGS) ID"]: r for r in _csv(meta / "pgs_all_metadata_scores.csv")}
    sets = {r["PGS Sample Set (PSS)"]: r for r in _csv(meta / "pgs_all_metadata_evaluation_sample_sets.csv")}
    ind_pubs: dict[str, set] = {}
    ind_n: dict[str, int] = {}
    all_n: dict[str, int] = {}
    for p in _csv(meta / "pgs_all_metadata_performance_metrics.csv"):
        k = p["Evaluated Score"]
        s = scores.get(k)
        if not s:
            continue
        n = _int(sets.get(p["PGS Sample Set (PSS)"], {}).get("Number of Individuals"))
        all_n[k] = all_n.get(k, 0) + n
        if p["PGS Publication (PGP) ID"] != s["PGS Publication (PGP) ID"]:
            ind_pubs.setdefault(k, set()).add(p["PGS Publication (PGP) ID"])
            ind_n[k] = ind_n.get(k, 0) + n

    def key(k: str) -> tuple:
        return (len(ind_pubs.get(k, ())), ind_n.get(k, 0), all_n.get(k, 0), scores[k]["Release Date"], k)

    ranked: dict[str, list[str]] = {}
    for k, s in scores.items():
        t = s["Mapped Trait(s) (EFO ID)"]
        if not t or "," in t or "|" in t or s["Number of Interaction Terms"] not in ("", "0"):
            continue
        ranked.setdefault(t, []).append(k)
    for ks in ranked.values():
        ks.sort(key=key, reverse=True)
    best = {t: ks[0] for t, ks in ranked.items()}

    def row(t: str, k: str, f: dict | None) -> dict:
        s = scores[k]
        return {"trait_id": t, "trait_label": s["Mapped Trait(s) (EFO label)"], "pgs_id": k,
                "label": (f or {}).get("label") or s["Mapped Trait(s) (EFO label)"],
                "section": (f or {}).get("section") or "", "body_system": (f or {}).get("body_system") or "",
                "sex": (f or {}).get("sex") or "", "featured": f is not None, "alternate_for": "",
                "rule_pick": best.get(t), "pinned": bool(f and f.get("pin")), "why": (f or {}).get("why") or "",
                "independent_pubs": len(ind_pubs.get(k, ())), "independent_n": ind_n.get(k, 0),
                "evaluated_n": all_n.get(k, 0), "n_variants": _int(s["Number of Variants"])}

    out, seen = [], set()
    for f in featured:
        t = f["trait_id"]
        k = f.get("pin") or best.get(t)
        if not k or k not in scores:
            raise ValueError(f"featured trait {t} has no usable score in the PGS Catalog")
        out.append(row(t, k, f))
        seen.add(t)
        # Stand-ins, scored alongside in case the pick can't be matched to your data (the pipeline promotes the
        # best-ranked one that scores; the rest are not shown).
        out.extend({**row(t, alt, f), "featured": False, "pinned": False, "alternate_for": t}
                   for alt in [a for a in ranked.get(t, []) if a != k][:ALTERNATES])
    for t, k in sorted(best.items(), key=lambda kv: scores[kv[1]]["Mapped Trait(s) (EFO label)"].lower()):
        if t not in seen and ind_pubs.get(k):
            out.append(row(t, k, None))
    return out


class PgsCatalog(Source):
    id = "pgs_catalog"
    title = "PGS Catalog"
    homepage = "https://www.pgscatalog.org"
    licence = ("Catalog metadata: EMBL-EBI terms of use (open, cite Lambert et al. Nat Genet 2021); each scoring "
               "file carries its own licence, recorded per score")
    cadence = "releases most weeks; new scores and new independent evaluations"
    description = ("The open catalogue of published polygenic scores — recipes that add up thousands to millions of "
                   "small genetic effects — together with every published test of how well each score works.")
    schema = 2  # 2: featured traits carry next-ranked alternates

    def probe(self, pin: str | None = None) -> list[Upstream]:
        data = FEATURED.read_bytes()
        return [http.head(META), Upstream(url=f"vendored:{FEATURED.name}",
                                          etag=hashlib.sha256(data).hexdigest()[:16], size=len(data))]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        tgz = http.fetch(META, scratch / "meta.tar.gz", upstream[0].size)
        meta = scratch / "meta"
        meta.mkdir()
        with tarfile.open(tgz) as tf:
            tf.extractall(meta, filter="data")
        hits = list(meta.rglob("pgs_all_metadata_scores.csv"))
        meta = hits[0].parent
        featured = _csv(FEATURED)
        sel = select(meta, featured)

        con = duckdb.connect()
        con.execute("SET preserve_insertion_order=true")

        def csvr(name: str) -> str:
            return f"read_csv('{meta / name}', header=true, all_varchar=true)"

        tables = {}
        con.execute(f"""COPY (SELECT "Polygenic Score (PGS) ID" AS pgs_id, "PGS Name" AS name,
            "Reported Trait" AS trait_reported,
            "Mapped Trait(s) (EFO label)" trait_label, "Mapped Trait(s) (EFO ID)" trait_id,
            "PGS Development Method" AS method, "PGS Development Details/Relevant Parameters" AS method_details,
            TRY_CAST("Number of Variants" AS BIGINT) AS n_variants, "PGS Publication (PGP) ID" AS pgp_id,
            "Publication (PMID)" pmid, "Publication (doi)" doi,
            "Ancestry Distribution (%) - Source of Variant Associations (GWAS)" ancestry_gwas,
            "Ancestry Distribution (%) - Score Development/Training" AS ancestry_training,
            "Ancestry Distribution (%) - PGS Evaluation" AS ancestry_evaluation,
            "Release Date" AS released, "License/Terms of Use" AS licence
            FROM {csvr('pgs_all_metadata_scores.csv')}) TO '{out / 'scores.parquet'}' (FORMAT parquet)""")
        tables["scores"] = out / "scores.parquet"
        con.execute(f"""COPY (SELECT "PGS Publication/Study (PGP) ID" AS pgp_id, "First Author" AS first_author,
            "Title" AS title, "Journal Name" AS journal, "Publication Date" AS published, "PubMed ID (PMID)" pmid,
            "digital object identifier (doi)" doi
            FROM {csvr('pgs_all_metadata_publications.csv')}) TO '{out / 'publications.parquet'}' (FORMAT parquet)""")
        tables["publications"] = out / "publications.parquet"
        con.execute(f"""COPY (SELECT "Ontology Trait ID" AS trait_id, "Ontology Trait Label" AS label,
            "Ontology Trait Description" AS description, "Ontology URL" AS url
            FROM {csvr('pgs_all_metadata_efo_traits.csv')}) TO '{out / 'traits.parquet'}' (FORMAT parquet)""")
        tables["traits"] = out / "traits.parquet"
        con.execute(f"""COPY (
            SELECT p."PGS Performance Metric (PPM) ID" AS ppm_id, p."Evaluated Score" AS pgs_id,
                   p."PGS Sample Set (PSS)" pss_id, p."PGS Publication (PGP) ID" AS pgp_id,
                   p."PGS Publication (PGP) ID" <> s."PGS Publication (PGP) ID" AS independent,
                   p."Reported Trait" AS trait_reported, p."Covariates Included in the Model" AS covariates,
                   p."PGS Performance: Other Relevant Information" AS info,
                   p."Hazard Ratio (HR)" hr, p."Odds Ratio (OR)" "or", p."Beta" AS beta,
                   p."Area Under the Receiver-Operating Characteristic Curve (AUROC)" auroc,
                   p."Concordance Statistic (C-index)" cindex, p."Other Metric(s)" other,
                   TRY_CAST(TRY_CAST(e."Number of Individuals" AS DOUBLE) AS BIGINT) AS n,
                   TRY_CAST(TRY_CAST(e."Number of Cases" AS DOUBLE) AS BIGINT) AS cases,
                   e."Broad Ancestry Category" AS ancestry, e."Cohort(s)" cohorts
            FROM {csvr('pgs_all_metadata_performance_metrics.csv')} p
            JOIN {csvr('pgs_all_metadata_scores.csv')} s ON s."Polygenic Score (PGS) ID" = p."Evaluated Score"
            LEFT JOIN {csvr('pgs_all_metadata_evaluation_sample_sets.csv')} e
                   ON e."PGS Sample Set (PSS)" = p."PGS Sample Set (PSS)"
            WHERE p."Evaluated Score" IN (SELECT unnest(?))) TO '{out / 'evaluations.parquet'}' (FORMAT parquet)""",
                    [[r["pgs_id"] for r in sel]])
        tables["evaluations"] = out / "evaluations.parquet"

        pq.write_table(pa.Table.from_pylist(sel), out / "selected.parquet")
        tables["selected"] = out / "selected.parquet"

        sdir = out / "scoring"
        sdir.mkdir()

        def get(pid: str) -> None:
            http.fetch(scoring_url(pid), scratch / "scoring" / f"{pid}.txt.gz")
            shutil.move(scratch / "scoring" / f"{pid}.txt.gz", sdir / f"{pid}_hmPOS_GRCh37.txt.gz")

        with ThreadPoolExecutor(6) as ex:
            list(ex.map(get, sorted({r["pgs_id"] for r in sel})))
        version = (upstream[0].last_modified or "")[:10] or "unknown"
        return Built(version=version, tables=tables, upstream=upstream,
                     notes={"citation": "Lambert SA, et al. The Polygenic Score Catalog as an open database for "
                                        "reproducibility and systematic evaluation. Nat Genet. 2021;53:420-425.",
                            "selected": len(sel), "featured": sum(r["featured"] for r in sel),
                            "scoring_dir": "scoring"})


class PgsReference(Source):
    id = "pgs_reference"
    title = "PGS Catalog 1000 Genomes reference panel"
    homepage = "https://pgsc-calc.readthedocs.io/en/latest/how-to/prepare.html"
    licence = "1000 Genomes Project data (open; Fort Lauderdale principles), packaged by the PGS Catalog"
    cadence = "rarely (panel format versions)"
    description = ("2,504 people from 26 populations, sequenced by the 1000 Genomes Project. Your scores are placed "
                   "among the people whose genetic ancestry is most similar to yours, which is what makes a "
                   "percentile meaningful.")

    def probe(self, pin: str | None = None) -> list[Upstream]:
        return [http.head(REF_URL)]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        cache = self.workspace / "cache" / "pgs" / "pgsc_1000G_v1.tar.zst"
        tar = cache if cache.exists() and cache.stat().st_size == upstream[0].size else \
            http.fetch(REF_URL, cache, upstream[0].size)
        _untar_zstd(tar, out)
        psam = next(out.glob("GRCh37_*.psam"))
        con = duckdb.connect()
        con.execute(f"COPY (SELECT * FROM read_csv('{psam}', delim='\\t', header=true, all_varchar=true)) "
                    f"TO '{out / 'samples.parquet'}' (FORMAT parquet)")
        meta = (out / "meta.txt").read_text().strip() if (out / "meta.txt").exists() else ""
        return Built(version="v1", tables={"samples": out / "samples.parquet"}, upstream=upstream,
                     notes={"citation": "1000 Genomes Project Consortium. A global reference for human genetic "
                                        "variation. Nature. 2015;526:68-74.", "meta": meta})


def _untar_zstd(tar: Path, out: Path) -> None:
    """Extract the GRCh37 files (macOS and Windows tar lack zstd, so do it in Python)."""
    import zstandard

    with open(tar, "rb") as fh, zstandard.ZstdDecompressor().stream_reader(fh) as zr, \
            tarfile.open(fileobj=zr, mode="r|") as tf:
        for m in tf:
            if m.name.startswith("GRCh37") or m.name == "meta.txt":
                tf.extract(m, out, filter="data")
