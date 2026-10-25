"""Reference-population metadata for the Ancestry section (ADR-020).

The 1000 Genomes reference panel labels each person with a population code (GBR, FIN, …). The International
Genome Sample Resource (IGSR) publishes the official name, description, sampling location and sample count for
each code; this source stores them so the dashboard can say "British in England and Scotland" with a link to the
population's IGSR page instead of a hand-written label.
"""

from __future__ import annotations

import datetime
import hashlib
import json
import urllib.request
from pathlib import Path

import duckdb

from . import Built, Source, Upstream

API = "https://www.internationalgenome.org/api/beta/population/_search"
PAGE = "https://www.internationalgenome.org/data-portal/population/{code}"


def _get() -> list[dict]:
    body = json.dumps({"size": 500, "query": {"match_all": {}}}).encode()
    req = urllib.request.Request(API, data=body, headers={"Content-Type": "application/json",
                                                          "User-Agent": "wgs-analysis"})
    with urllib.request.urlopen(req, timeout=60) as r:
        hits = json.load(r)["hits"]["hits"]
    return [h["_source"] for h in hits]


def rows(records: list[dict]) -> list[dict]:
    """One row per population in the 1000 Genomes phase 3 panel (the 26 populations the reference uses) and in the
    Human Genome Diversity Project (54, used for finer ancestry). HGDP populations have no IGSR code; their IGSR
    page ID ("YorubaHGDP") is the population name in the HGDP sample file plus "HGDP", so it serves as the code."""
    out = []
    for p in records:
        colls = [c.get("title", "") for c in p.get("dataCollections") or []]
        hgdp = "Human Genome Diversity Project" in colls
        if "1000 Genomes phase 3 release" not in colls and not hgdp:
            continue
        code = p.get("code") or (p.get("elasticId") if hgdp else None)
        if not code:
            continue
        sp = p.get("superpopulation") or {}
        out.append({"code": code, "name": p.get("name"), "description": p.get("description"),
                    "superpop": sp.get("code"), "superpop_name": sp.get("name"),
                    "latitude": p.get("latitude"), "longitude": p.get("longitude"),
                    "samples": (p.get("samples") or {}).get("count"), "url": PAGE.format(code=code),
                    "panel": "hgdp" if hgdp and not p.get("code") else "1000genomes"})
    return sorted(out, key=lambda r: (r["panel"], r["superpop"] or "", r["code"]))


class IgsrPopulations(Source):
    id = "igsr_populations"
    title = "1000 Genomes populations (IGSR)"
    homepage = "https://www.internationalgenome.org/data-portal/population"
    licence = "IGSR data: free to use without restriction (Fort Lauderdale principles)"
    cadence = "rarely"
    schema = 2
    description = ("The official names and sampling locations of the 26 populations in the 1000 Genomes reference "
                   "panel and the 54 of the Human Genome Diversity Project, used to describe which reference group "
                   "your DNA is most similar to.")

    def __init__(self) -> None:
        self._records: list[dict] | None = None

    def probe(self, pin: str | None = None) -> list[Upstream]:
        self._records = _get()
        rs = rows(self._records)
        digest = hashlib.sha256(json.dumps(rs, sort_keys=True).encode()).hexdigest()[:16]
        return [Upstream(url=API, etag=digest, size=len(rs))]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        rs = rows(self._records if self._records is not None else _get())
        if sum(r["panel"] == "1000genomes" for r in rs) < 26:
            raise RuntimeError(f"IGSR returned only {len(rs)} phase 3 populations (expected 26)")
        src = scratch / "populations.json"
        src.write_text(json.dumps(rs))
        dst = out / "populations.parquet"
        duckdb.connect().execute(f"COPY (SELECT * FROM read_json('{src}')) TO '{dst}' (FORMAT parquet)")
        return Built(version=datetime.date.today().isoformat(), tables={"populations": dst}, upstream=upstream,
                     notes={"citation": "Fairley S, et al. The International Genome Sample Resource (IGSR) "
                                        "collection of open human genomic variation resources. Nucleic Acids Res. "
                                        "2020;48:D941-D947."})
