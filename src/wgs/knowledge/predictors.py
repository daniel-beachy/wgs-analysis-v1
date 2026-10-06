"""Computational predictions: what a variant *might* do when no expert has classified it.

These are labelled "predicted" everywhere and can never raise a claim above **Limited** evidence
(ADR-014). Predictions are useful for exploration and for prioritising questions to ask a
professional, not for conclusions.

* **AlphaMissense** (Google DeepMind, 2023): a protein-structure-aware model scoring every possible
  single-letter protein change. The developers' classes: likely benign < 0.34 ≤ ambiguous ≤ 0.564 <
  likely pathogenic.
* **REVEL** (Ioannidis et al., 2016): an ensemble of 13 older predictors for rare missense variants.
  ClinGen calibrated it (Pejaver et al., 2022): ≥ 0.644 supporting, ≥ 0.773 moderate, ≥ 0.932 strong
  evidence of pathogenicity (ACMG PP3); ≤ 0.290 supports benign (BP4).
* **gnomAD gene constraint** (v4.1): how strongly natural selection removes broken copies of each gene
  from the healthy population. A low LOEUF (< 0.6 in v4) means losing one copy is rarely tolerated.
  Gene-level, so independent of genome build.
"""

from __future__ import annotations

import zipfile
from pathlib import Path

import duckdb

from . import Built, Source, Upstream
from .http import fetch, head

AM_URL = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg19.tsv.gz"
REVEL_URL = "https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip"
CONSTRAINT_URL = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/"
                  "gnomad.v4.1.constraint_metrics.tsv")


def _con(scratch: Path):
    con = duckdb.connect()
    con.execute("SET preserve_insertion_order=false")
    con.execute(f"SET temp_directory='{scratch / 'duckdb_tmp'}'")
    return con


class AlphaMissense(Source):
    id = "alphamissense"
    title = "AlphaMissense (GRCh37)"
    homepage = "https://github.com/google-deepmind/alphamissense"
    licence = "CC BY-NC-SA 4.0 (non-commercial)"
    cadence = "static (2023 release)"
    description = ("AI-predicted effect of every possible single-letter protein change (0 = likely harmless, "
                   "1 = likely damaging). A prediction, not a clinical classification.")

    def probe(self, pin: str | None = None) -> list[Upstream]:
        return [head(AM_URL)]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        src = fetch(AM_URL, scratch / "am.tsv.gz", upstream[0].size)
        dst = out / "scores.parquet"
        con = _con(scratch)
        con.execute(f"""
          COPY (
            SELECT replace(chrom, 'chr', '') AS chrom, pos, ref, alt, max(score)::FLOAT AS am_score,
                   arg_max(protein_variant, score) AS am_protein_variant
            FROM read_csv('{src}', delim='\t', header=false, skip=4, quote='',
                          columns={{'chrom':'VARCHAR','pos':'INTEGER','ref':'VARCHAR','alt':'VARCHAR',
                                   'genome':'VARCHAR','uniprot':'VARCHAR','transcript':'VARCHAR',
                                   'protein_variant':'VARCHAR','score':'DOUBLE','cls':'VARCHAR'}})
            GROUP BY ALL ORDER BY chrom, pos
          ) TO '{dst}' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 500000)""")
        n = con.execute(f"SELECT count(*) FROM '{dst}'").fetchone()[0]
        return Built(version="2023-hg19", tables={"scores": dst}, upstream=upstream,
                     notes={"variants": n, "citation": "Cheng J, et al. Science 2023;381:eadg7492"})


class Revel(Source):
    id = "revel"
    title = "REVEL v1.3"
    homepage = "https://sites.google.com/site/revelgenomics/"
    licence = "Free for non-commercial use (REVEL terms)"
    cadence = "static (v1.3, 2021)"
    description = ("Ensemble score for rare missense variants, calibrated by ClinGen into supporting / moderate "
                   "/ strong levels of computational evidence.")

    def probe(self, pin: str | None = None) -> list[Upstream]:
        return [head(REVEL_URL)]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        z = fetch(REVEL_URL, scratch / "revel.zip", upstream[0].size)
        with zipfile.ZipFile(z) as zf:
            name = zf.namelist()[0]
            zf.extract(name, scratch)
        z.unlink()
        csv = scratch / name
        dst = out / "scores.parquet"
        con = _con(scratch)
        con.execute(f"""
          COPY (
            SELECT chr AS chrom, TRY_CAST(hg19_pos AS INTEGER) AS pos, ref, alt, max(REVEL)::FLOAT AS revel
            FROM read_csv('{csv}', header=true, all_varchar=false,
                          types={{'chr':'VARCHAR','hg19_pos':'VARCHAR','REVEL':'DOUBLE'}})
            WHERE TRY_CAST(hg19_pos AS INTEGER) IS NOT NULL
            GROUP BY ALL ORDER BY chrom, pos
          ) TO '{dst}' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 500000)""")
        n = con.execute(f"SELECT count(*) FROM '{dst}'").fetchone()[0]
        return Built(version="1.3", tables={"scores": dst}, upstream=upstream,
                     notes={"variants": n, "citation": "Ioannidis NM, et al. Am J Hum Genet 2016;99:877-885"})


class GeneConstraint(Source):
    id = "constraint"
    title = "gnomAD gene constraint v4.1"
    homepage = "https://gnomad.broadinstitute.org/downloads#v4-constraint"
    licence = "CC0 1.0 (gnomAD)"
    cadence = "with gnomAD major releases"
    description = ("How intolerant each gene is to being broken (LOEUF, pLI) or changed (missense Z), measured "
                   "from which variants are missing in ~800,000 healthy people.")

    def probe(self, pin: str | None = None) -> list[Upstream]:
        return [head(CONSTRAINT_URL)]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        src = fetch(CONSTRAINT_URL, scratch / "constraint.tsv", upstream[0].size)
        dst = out / "genes.parquet"
        con = _con(scratch)
        con.execute(f"""
          COPY (
            SELECT gene, arg_max(transcript, (mane_select = 'true')::INT * 2 + (canonical = 'true')::INT) AS transcript,
                   arg_max(TRY_CAST("lof.oe_ci.upper" AS FLOAT),
                           (mane_select = 'true')::INT * 2 + (canonical = 'true')::INT) AS loeuf,
                   arg_max(TRY_CAST("lof.pLI" AS FLOAT),
                           (mane_select = 'true')::INT * 2 + (canonical = 'true')::INT) AS pli,
                   arg_max(TRY_CAST("mis.z_score" AS FLOAT),
                           (mane_select = 'true')::INT * 2 + (canonical = 'true')::INT) AS mis_z,
                   arg_max(TRY_CAST("lof.obs" AS INTEGER),
                           (mane_select = 'true')::INT * 2 + (canonical = 'true')::INT) AS lof_obs,
                   arg_max(TRY_CAST("lof.exp" AS FLOAT),
                           (mane_select = 'true')::INT * 2 + (canonical = 'true')::INT) AS lof_exp
            FROM read_csv('{src}', delim='\t', header=true, all_varchar=true, quote='')
            WHERE mane_select = 'true' OR canonical = 'true'
            GROUP BY gene ORDER BY gene
          ) TO '{dst}' (FORMAT parquet)""")
        n = con.execute(f"SELECT count(*) FROM '{dst}'").fetchone()[0]
        return Built(version="4.1", tables={"genes": dst}, upstream=upstream, notes={"genes": n})
