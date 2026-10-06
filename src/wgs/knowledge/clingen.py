"""ClinGen: expert curation of which genes really cause which diseases, how dosage-sensitive they are,
and how *actionable* a finding is (can anything be done about it?).

A pathogenic-looking variant only matters if its gene is genuinely linked to the disease. ClinGen's
gene–disease validity scale (Definitive → Strong → Moderate → Limited → Disputed/Refuted) is the
evidence the tool uses for that link. The exports are regenerated daily, so the version is the
export date but a rebuild only happens when the curated content actually changes.

Actionability (adult and paediatric contexts) scores gene–condition pairs on severity, likelihood,
effectiveness and nature of intervention (each 0–3, summed to an "overall" score up to 12).
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import urllib.request
from pathlib import Path

import duckdb

from . import Built, Source, Upstream
from .http import UA

VALIDITY = "https://search.clinicalgenome.org/kb/gene-validity/download"
DOSAGE = "https://search.clinicalgenome.org/kb/gene-dosage/download"
ACTIONABILITY = {ctx: f"https://actionability.clinicalgenome.org/ac/{ctx}/api/summ?flavor=flat"
                 for ctx in ("Adult", "Pediatric")}


def _get(url: str) -> str:
    with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": UA}), timeout=120) as r:
        return r.read().decode("utf-8-sig")


def _rows(text: str) -> tuple[str, list[str], list[list[str]]]:
    rows = list(csv.reader(io.StringIO(text)))
    created = next((r[0].split(":", 1)[1].strip() for r in rows if r and r[0].startswith("FILE CREATED")), "")
    data = [r for r in rows if r and not r[0].startswith(("+", "CLINGEN", "FILE CREATED", "WEBPAGE"))]
    return created, data[0], data[1:]


class ClinGen(Source):
    id = "clingen"
    title = "ClinGen gene curation"
    homepage = "https://search.clinicalgenome.org/"
    licence = "CC0 1.0"
    cadence = "continuous (daily export)"
    schema = 2  # 2: adds actionability
    description = ("Expert panels' verdicts on whether a gene truly causes a disease (Definitive … Refuted), "
                   "whether losing or gaining a copy is harmful, and how actionable a finding is.")

    _cache: dict[str, str] = {}

    def probe(self, pin: str | None = None) -> list[Upstream]:
        out = []
        for url in (VALIDITY, DOSAGE):
            text = self._cache[url] = _get(url)
            _, head, rows = _rows(text)
            digest = hashlib.sha256(repr(sorted(map(tuple, rows))).encode()).hexdigest()[:16]
            out.append(Upstream(url=url, etag=digest, size=len(rows)))
        for url in ACTIONABILITY.values():
            text = self._cache[url] = _get(url)
            rows = json.loads(text)["rows"]
            digest = hashlib.sha256(repr(sorted(map(tuple, rows))).encode()).hexdigest()[:16]
            out.append(Upstream(url=url, etag=digest, size=len(rows)))
        return out

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        con = duckdb.connect()
        tables, created = {}, ""
        for url, name in ((VALIDITY, "gene_validity"), (DOSAGE, "dosage")):
            created, head, rows = _rows(self._cache.get(url) or _get(url))
            cols = [h.lower().replace(" (hgnc)", "").replace(" (mondo)", "").replace(" ", "_") for h in head]
            p = scratch / f"{name}.csv"
            with open(p, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(cols)
                w.writerows(rows)
            dst = out / f"{name}.parquet"
            con.execute(f"COPY (SELECT * FROM read_csv('{p}', header=true, all_varchar=true)) "
                        f"TO '{dst}' (FORMAT parquet)")
            tables[name] = dst
        rows = []
        for ctx, url in ACTIONABILITY.items():
            d = json.loads(self._cache.get(url) or _get(url))
            cols = d["columns"]
            for r in d["rows"]:
                x = dict(zip(cols, r, strict=False))
                rows.append({"context": ctx, "doc": x["docId"], "gene": x["geneOrVariant"], "disease": x["disease"],
                             "omim": x["omim"], "status": x["status-overall"], "outcome": x["outcome"],
                             "intervention": x["intervention"], "severity": x["severity"],
                             "likelihood": x["likelihood"], "effectiveness": x["effectiveness"],
                             "nature_of_intervention": x["natureOfIntervention"], "overall": x["overall"],
                             "url": f"https://actionability.clinicalgenome.org/ac/{ctx}/ui/stg2SummaryRpt?doc={x['docId']}"})
        p = scratch / "actionability.json"
        p.write_text(json.dumps(rows))
        dst = out / "actionability.parquet"
        con.execute(f"COPY (SELECT * FROM read_json('{p}', format='array', columns={{"
                    + ", ".join(f"'{c}': 'VARCHAR'" for c in rows[0]) + f"}})) TO '{dst}' (FORMAT parquet)")
        tables["actionability"] = dst
        return Built(version=created or "unknown", tables=tables, upstream=upstream)
