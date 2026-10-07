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


MANE_DIR = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"


class Mane(Source):
    id = "mane"
    title = "MANE Select transcripts (NCBI & EMBL-EBI)"
    homepage = "https://www.ncbi.nlm.nih.gov/refseq/MANE/"
    licence = "Public domain (NCBI) / no restrictions (EMBL-EBI)"
    cadence = "a few releases a year"
    description = ("The one agreed 'standard' transcript per gene that clinical labs use to name protein changes "
                   "(e.g. GALT p.Asn314Asp), so names here match papers and lab reports.")

    def probe(self, pin: str | None = None) -> list[Upstream]:
        import re

        from .http import listing
        names = sorted(set(re.findall(r'href="(MANE\.GRCh38\.v[\d.]+\.summary\.txt\.gz)"', listing(MANE_DIR))))
        if not names:
            raise RuntimeError(f"no MANE summary file found at {MANE_DIR}")
        return [head(MANE_DIR + names[-1])]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        import re

        src = fetch(upstream[0].url, scratch / "mane.summary.txt.gz", upstream[0].size)
        dst = out / "transcripts.parquet"
        duckdb.connect().execute(f"""
          COPY (SELECT split_part("Ensembl_nuc", '.', 1) AS transcript, symbol AS gene, "RefSeq_nuc" AS refseq,
                       "Ensembl_prot" AS ensembl_prot, "RefSeq_prot" AS refseq_prot, "MANE_status" AS status
                FROM read_csv('{src}', delim='\t', header=true, all_varchar=true)
                ORDER BY transcript)
          TO '{dst}' (FORMAT parquet)""")
        version = re.search(r"v([\d.]+)\.summary", upstream[0].url).group(1)
        return Built(version=version, tables={"transcripts": dst}, upstream=upstream,
                     notes={"citation": "Morales J, et al. A joint NCBI and EMBL-EBI transcript set for clinical "
                                        "genomics and research. Nature. 2022;604:310-315."})
