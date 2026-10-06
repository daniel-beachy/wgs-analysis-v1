"""ClinVar (NCBI): curated assertions about what individual variants do. Weekly GRCh37 VCF.

Review status → stars follows NCBI's own scale
(https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/). Stars measure how much review an
assertion received, not how severe the variant is.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import duckdb

from . import Built, Source, Upstream
from .http import fetch, head, listing

BASE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37"

STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
    "criteria_provided,_single_submitter": 1,
}

# ClinVar significance → a small fixed vocabulary the rest of the tool reasons about
SIG_CLASS_SQL = """
CASE
  WHEN s IS NULL THEN 'other'
  WHEN s LIKE 'conflicting%' THEN 'conflicting'
  WHEN s LIKE 'pathogenic%' THEN 'pathogenic'
  WHEN s LIKE 'likely_pathogenic%' THEN 'likely_pathogenic'
  WHEN s LIKE 'uncertain%' OR s LIKE 'vus%' THEN 'uncertain'
  WHEN s LIKE 'benign%' THEN 'benign'
  WHEN s LIKE 'likely_benign%' THEN 'likely_benign'
  WHEN s = 'drug_response' THEN 'drug_response'
  WHEN s IN ('risk_factor', 'likely_risk_allele', 'established_risk_allele') THEN 'risk_factor'
  WHEN s = 'protective' THEN 'protective'
  WHEN s IN ('association', 'affects', 'confers_sensitivity') THEN 'association'
  ELSE 'other'
END
"""


class ClinVar(Source):
    id = "clinvar"
    title = "ClinVar"
    homepage = "https://www.ncbi.nlm.nih.gov/clinvar/"
    licence = "Public domain (US Government work); see NCBI policies"
    cadence = "weekly"
    description = ("NCBI's archive of what labs and expert panels have concluded about specific variants "
                   "(pathogenic, benign, drug response…), with a review status shown as 0–4 stars.")

    def probe(self, pin: str | None = None) -> list[Upstream]:
        if pin:  # e.g. "20250106" → reproduce an old analysis exactly
            year = pin[:4]
            url = f"{BASE}/archive_2.0/{year}/clinvar_{pin}.vcf.gz"
            return [head(url), head(url + ".tbi")]
        return [head(f"{BASE}/clinvar.vcf.gz"), head(f"{BASE}/clinvar.vcf.gz.tbi")]

    @staticmethod
    def archive_dates() -> list[str]:
        out = []
        for y in ("2026", "2025"):
            out += re.findall(r"clinvar_(\d{8})\.vcf\.gz\"", listing(f"{BASE}/archive_2.0/{y}/"))
        return sorted(set(out), reverse=True)

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        vcf = fetch(upstream[0].url, scratch / "clinvar.vcf.gz", upstream[0].size)
        header = subprocess.run(["bcftools", "view", "-h", str(vcf)], capture_output=True, text=True,
                                check=True).stdout
        m = re.search(r"##fileDate=(\d{4})-?(\d{2})-?(\d{2})", header)
        version = f"{m.group(1)}-{m.group(2)}-{m.group(3)}" if m else upstream[0].last_modified[:10]
        tsv = scratch / "clinvar.tsv"
        fmt = ("%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/ALLELEID\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNDN\t"
               "%INFO/GENEINFO\t%INFO/MC\t%INFO/CLNVC\t%INFO/ORIGIN\t%INFO/CLNSIGCONF\t%INFO/RS\n")
        with open(tsv, "w") as fh:
            subprocess.run(["bcftools", "query", "-f", fmt, str(vcf)], stdout=fh, stderr=subprocess.DEVNULL, check=True)
        stars = " ".join(f"WHEN '{k}' THEN {v}" for k, v in STARS.items())
        dst = out / "variants.parquet"
        con = duckdb.connect()
        con.execute(f"""
        COPY (
          WITH raw AS (
            SELECT * FROM read_csv('{tsv}', delim='\t', header=false, quote='', escape='', all_varchar=true,
              columns={{'chrom':'VARCHAR','pos':'VARCHAR','id':'VARCHAR','ref':'VARCHAR','alt':'VARCHAR',
                       'allele_id':'VARCHAR','clnsig':'VARCHAR','revstat':'VARCHAR','disease':'VARCHAR',
                       'geneinfo':'VARCHAR','mc':'VARCHAR','vtype':'VARCHAR','origin':'VARCHAR',
                       'sigconf':'VARCHAR','rs':'VARCHAR'}})
          ), n AS (
            SELECT
              chrom, CAST(pos AS INTEGER) AS pos, ref, alt,
              CAST(id AS BIGINT) AS variation_id,
              CAST(NULLIF(allele_id, '.') AS BIGINT) AS allele_id,
              CASE WHEN rs = '.' THEN NULL ELSE 'rs' || rs END AS rsid,
              NULLIF(clnsig, '.') AS significance,
              lower(split_part(NULLIF(clnsig, '.'), '|', 1)) AS s,
              list_transform(string_split(lower(NULLIF(clnsig, '.')), '|')[2:], x -> x) AS sig_modifiers,
              NULLIF(revstat, '.') AS review_status,
              replace(NULLIF(disease, '.'), '_', ' ') AS conditions,
              list_transform(string_split(NULLIF(geneinfo, '.'), '|'), g -> split_part(g, ':', 1)) AS genes,
              list_distinct(list_transform(string_split(NULLIF(mc, '.'), ','), c -> split_part(c, '|', 2)))
                AS molecular_consequences,
              NULLIF(vtype, '.') AS variant_type,
              CAST(NULLIF(origin, '.') AS INTEGER) AS origin_bits,
              NULLIF(sigconf, '.') AS conflicting_detail
            FROM raw
            WHERE alt <> '.'
          )
          SELECT * EXCLUDE (s),
            {SIG_CLASS_SQL} AS sig_class,
            coalesce(lower(significance) LIKE '%low_penetrance%', false) AS low_penetrance,
            CASE review_status {stars} ELSE 0 END AS stars
          FROM n
          ORDER BY chrom, pos
        ) TO '{dst}' (FORMAT parquet, COMPRESSION zstd)
        """)
        n = con.execute(f"SELECT count(*) FROM '{dst}'").fetchone()[0]
        return Built(version=version, tables={"variants": dst}, upstream=upstream,
                     notes={"records": n, "assembly": "GRCh37"})
