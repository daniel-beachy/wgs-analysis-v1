"""STRchive: repeat-expansion disease loci (ADR-021).

Short tandem repeats are stretches where a few letters repeat ("CAGCAGCAG…"). At about 80 known places, a repeat
that grows past a threshold causes disease (Huntington's, fragile X, myotonic dystrophy, …). STRchive (Dashnow lab)
curates each locus: normal, intermediate and disease-causing repeat counts, inheritance, age of onset, gene–disease
evidence and PubMed references, plus a ready-made ExpansionHunter catalogue for each genome build.

Hiatt L, et al. STRchive: a dynamic resource detailing population-level and locus-specific insights at tandem
repeat disease loci. Genome Med. 2025;17:29. PMID 40140942.
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import duckdb

from . import Built, Source, Upstream

REPO = "dashnowlab/STRchive"
LATEST = f"https://api.github.com/repos/{REPO}/releases/latest"
RAW = "https://raw.githubusercontent.com/" + REPO + "/{tag}/data/{path}"
LOCI = "STRchive-loci.json"
CATALOG = "catalogs/STRchive-disease-loci.hg19.expansionhunter.json"
PAGE = "https://strchive.org/loci/{id}/"

KEEP = ["id", "disease_id", "gene", "chrom", "start_hg19", "stop_hg19", "disease", "age_onset", "age_onset_min",
        "typ_age_onset_min", "detection", "location_in_gene", "benign_min", "benign_max", "intermediate_min",
        "intermediate_max", "pathogenic_min", "pathogenic_max", "motif_len", "ref_copies"]
LISTS = ["evidence", "inheritance", "association_type", "reference_motif_reference_orientation",
         "pathogenic_motif_reference_orientation", "references", "genereviews", "omim", "mondo"]


def _get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "wgs-analysis", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read()


def rows(loci: list[dict]) -> list[dict]:
    out = []
    for r in loci:
        row = {k: r.get(k) for k in KEEP}
        for k in LISTS:
            v = r.get(k) or []
            row[k] = [str(x) for x in (v if isinstance(v, list) else [v])]
        row["url"] = PAGE.format(id=r["id"])
        row["pmids"] = [x.split(":", 1)[1] for x in row["references"] if x.startswith("pmid:")]
        out.append(row)
    return out


class Strchive(Source):
    id = "strchive"
    title = "STRchive (repeat-expansion disease loci)"
    homepage = "https://strchive.org"
    licence = "MIT"
    cadence = "monthly"
    description = ("Every known place where a repeated stretch of DNA causes disease when it grows too long, with the "
                   "normal and disease-causing repeat counts and the studies behind them.")

    def probe(self, pin: str | None = None) -> list[Upstream]:
        tag = pin or json.loads(_get(LATEST))["tag_name"]
        return [Upstream(url=RAW.format(tag=tag, path=LOCI), etag=tag),
                Upstream(url=RAW.format(tag=tag, path=CATALOG), etag=tag)]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        loci = json.loads(_get(upstream[0].url))
        catalog = json.loads(_get(upstream[1].url))
        if len(loci) < 50 or len(catalog) < 50:
            raise RuntimeError(f"STRchive looks incomplete: {len(loci)} loci, {len(catalog)} catalogue entries")
        src = scratch / "loci.json"
        src.write_text(json.dumps(rows(loci)))
        dst = out / "loci.parquet"
        duckdb.connect().execute(f"COPY (SELECT * FROM read_json('{src}')) TO '{dst}' (FORMAT parquet)")
        (out / "catalog_hg19.json").write_text(json.dumps(catalog, indent=1))
        return Built(version=upstream[0].etag.lstrip("v"), tables={"loci": dst}, upstream=upstream,
                     notes={"citation": "Hiatt L, et al. STRchive: a dynamic resource detailing population-level "
                                        "and locus-specific insights at tandem repeat disease loci. Genome Med. 2025."})
