"""Population allele frequencies: how common each variant is in other people.

* 1000 Genomes phase 3 (2,504 people, five continental groups): global and per-continent
  frequencies for every site, including ones where you match the reference. Powers the
  "how you compare" views and, later, ancestry.
* gnomAD v2.1.1 (≈141,000 people, exomes + genomes), as packaged by slivar's "gnotate" archive:
  the highest frequency in any population (popmax) and homozygote counts. This is what tells a
  truly rare variant apart from one that is merely absent from 2,504 people — a key piece of
  ACMG/AMP evidence (BA1/BS1/PM2).
"""

from __future__ import annotations

import subprocess
import zipfile
from pathlib import Path

import duckdb

from . import Built, Source, Upstream
from .http import fetch, head

KG = ("https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/"
      "ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz")
GNOTATE = "https://s3.amazonaws.com/slivar/gnomad.hg37.zip"


class ThousandGenomes(Source):
    id = "1000genomes"
    title = "1000 Genomes phase 3"
    homepage = "https://www.internationalgenome.org/"
    licence = "Fort Lauderdale / open (IGSR data reuse policy)"
    cadence = "frozen (2013 call set, v5b; NCBI mirror)"
    description = ("Allele frequencies from 2,504 people in 26 populations, grouped as African, American, "
                   "East Asian, European and South Asian.")

    def probe(self, pin: str | None = None) -> list[Upstream]:
        return [head(KG)]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        vcf = fetch(upstream[0].url, scratch / "kg.sites.vcf.gz", upstream[0].size)
        tsv = scratch / "kg.tsv"
        fmt = "%CHROM\t%POS\t%REF\t%ALT\t%AF\t%AFR_AF\t%AMR_AF\t%EAS_AF\t%EUR_AF\t%SAS_AF\n"
        with open(tsv, "w") as fh:
            norm = subprocess.Popen(["bcftools", "norm", "-m", "-any", "--no-version", "-Ou", str(vcf)],
                                    stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            subprocess.run(["bcftools", "query", "-e", 'ALT~"<"', "-f", fmt, "-"], stdin=norm.stdout, stdout=fh,
                           stderr=subprocess.DEVNULL, check=True)
            if norm.wait() != 0:
                raise RuntimeError("bcftools norm failed on 1000 Genomes sites")
        dst = out / "af.parquet"
        cols = ", ".join(f"TRY_CAST(NULLIF({c}, '.') AS FLOAT) AS {c}" for c in
                         ("af", "afr", "amr", "eas", "eur", "sas"))
        duckdb.connect().execute(f"""
          SET preserve_insertion_order=false;
          COPY (SELECT chrom, CAST(pos AS INTEGER) AS pos, ref, alt, {cols}
                FROM read_csv('{tsv}', delim='\t', header=false, quote='', all_varchar=true,
                  columns={{'chrom':'VARCHAR','pos':'VARCHAR','ref':'VARCHAR','alt':'VARCHAR','af':'VARCHAR',
                           'afr':'VARCHAR','amr':'VARCHAR','eas':'VARCHAR','eur':'VARCHAR','sas':'VARCHAR'}})
                ORDER BY chrom, pos)
          TO '{dst}' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 1000000)""")
        return Built(version="phase3-v5b", tables={"af": dst}, upstream=upstream,
                     notes={"people": 2504, "assembly": "GRCh37"})


class GnomAD(Source):
    id = "gnomad"
    title = "gnomAD v2.1.1 (via slivar)"
    homepage = "https://gnomad.broadinstitute.org/"
    licence = "CC0 (gnomAD); archive by slivar, MIT"
    cadence = "frozen (v2.1.1 is the last GRCh37-native release)"
    description = ("Frequencies from ≈141,000 people (125,748 exomes + 15,708 genomes). Used to judge whether a "
                   "variant is rare enough to plausibly cause a rare disease.")

    def probe(self, pin: str | None = None) -> list[Upstream]:
        return [head(GNOTATE)]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        dst = fetch(upstream[0].url, out / "gnomad.hg37.zip", upstream[0].size)
        with zipfile.ZipFile(dst) as z:
            fields = sorted({b.removeprefix("gnotate-").removesuffix(".bin")
                             for b in (n.split("/")[-1] for n in z.namelist())
                             if b.startswith("gnotate-") and b != "gnotate-variant.bin"})
            chroms = sorted({n.split("/")[1] for n in z.namelist() if n.count("/") >= 2})
        # slivar archives are not a table DuckDB can read; record a tiny parquet describing the fields instead.
        meta = out / "fields.parquet"
        con = duckdb.connect()
        con.execute("CREATE TABLE f(field VARCHAR)")
        con.executemany("INSERT INTO f VALUES (?)", [[f] for f in fields])
        con.execute(f"COPY f TO '{meta}' (FORMAT parquet)")
        return Built(version="v2.1.1", tables={"fields": meta}, upstream=upstream,
                     notes={"files": ["gnomad.hg37.zip"], "fields": fields, "chroms": chroms})
