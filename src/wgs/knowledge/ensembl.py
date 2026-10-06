"""Ensembl gene models for GRCh37 (release 87, frozen in 2017): where genes and their exons are.

Used for two things: `bcftools csq` predicts what each variant does to a protein, and the
`genes` table maps genes to coordinates for look-ups and region summaries. GRCh37 annotation no
longer changes upstream, so this source is effectively static — still versioned like the others.
"""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path

import duckdb

from . import Built, Source, Upstream
from .http import fetch, head

URL = "https://ftp.ensembl.org/pub/grch37/current/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz"


class Ensembl(Source):
    id = "ensembl"
    title = "Ensembl gene models (GRCh37)"
    homepage = "https://grch37.ensembl.org/"
    licence = "No restrictions (Ensembl data use policy)"
    cadence = "frozen (GRCh37 release 87)"
    description = ("Where each gene, transcript and exon sits on the genome; lets the tool predict what a "
                   "variant does to a protein (missense, stop-gained, splice…).")

    def probe(self, pin: str | None = None) -> list[Upstream]:
        return [head(URL)]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        gff = fetch(upstream[0].url, scratch / "genes.gff3.gz", upstream[0].size)
        shutil.copyfile(gff, out / "annotation.gff3.gz")
        tsv = scratch / "genes.tsv"
        with gzip.open(gff, "rt") as fh, open(tsv, "w") as o:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or "gene_id=" not in f[8] or not f[8].startswith("ID=gene:"):
                    continue
                a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                o.write("\t".join([f[0], f[3], f[4], f[6], a.get("gene_id", ""), a.get("Name", ""),
                                   a.get("biotype", ""), a.get("description", "").replace("\t", " ")]) + "\n")
        dst = out / "genes.parquet"
        duckdb.connect().execute(f"""
          COPY (SELECT chrom, CAST(start AS INTEGER) AS start, CAST("end" AS INTEGER) AS "end", strand, gene_id,
                       NULLIF(symbol, '') AS symbol, biotype,
                       NULLIF(regexp_replace(url_decode(description), ' \\[Source:.*\\]$', ''), '') AS description
                FROM read_csv('{tsv}', delim='\t', header=false, quote='', all_varchar=true,
                  columns={{'chrom':'VARCHAR','start':'VARCHAR','end':'VARCHAR','strand':'VARCHAR',
                           'gene_id':'VARCHAR','symbol':'VARCHAR','biotype':'VARCHAR','description':'VARCHAR'}})
                ORDER BY chrom, start)
          TO '{dst}' (FORMAT parquet, COMPRESSION zstd)""")
        return Built(version="GRCh37.87", tables={"genes": dst}, upstream=upstream,
                     notes={"files": ["annotation.gff3.gz"]})
