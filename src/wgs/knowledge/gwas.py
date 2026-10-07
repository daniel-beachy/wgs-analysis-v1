"""GWAS Catalog (single-variant trait associations) and SFARI Gene (autism research gene list).

The GWAS Catalog records every published genome-wide association: a variant, a trait, the allele that goes
with the trait, the effect size and p-value, and the study. ``data/gwas_featured.csv`` names the few
well-known single-variant traits the dashboard features (variant + the catalog's mapped trait label); all
evidence — direction, effect, replication — comes from the catalog at build time.
"""

from __future__ import annotations

import csv
import datetime
import hashlib
import zipfile
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from . import Built, Source, Upstream
from .http import fetch, head
from .pharmcat import CHAIN, Chain

GWAS = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip"
GWAS_FEATURED = Path(__file__).parent / "data" / "gwas_featured.csv"
SFARI = "https://gene.sfari.org//wp-content/themes/sfari-gene/utilities/download-csv.php?api-endpoint=genes"
GWS = 7.30103  # -log10(5e-8), genome-wide significance


def _vendored(p: Path) -> Upstream:
    data = p.read_bytes()
    return Upstream(url=f"vendored:{p.name}", etag=hashlib.sha256(data).hexdigest()[:16], size=len(data))


class GwasCatalog(Source):
    id = "gwas_catalog"
    title = "GWAS Catalog"
    homepage = "https://www.ebi.ac.uk/gwas/"
    licence = "EMBL-EBI terms of use (open; cite Cerezo et al. Nucleic Acids Res 2025)"
    cadence = "weekly"
    description = ("Every published genome-wide association study, curated by EMBL-EBI and NHGRI: which variant "
                   "goes with which trait, in which direction, how strongly, and in how many people.")

    def probe(self, pin: str | None = None) -> list[Upstream]:
        return [head(GWAS), head(CHAIN), _vendored(GWAS_FEATURED)]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        z = fetch(GWAS, scratch / "gwas.zip", upstream[0].size)
        with zipfile.ZipFile(z) as zf:
            name = next(n for n in zf.namelist() if n.endswith(".tsv"))
            zf.extract(name, scratch)
        tsv = scratch / name
        with open(GWAS_FEATURED, newline="", encoding="utf-8") as fh:
            feat = list(csv.DictReader(fh))
        rsids = sorted({f["rsid"] for f in feat})
        con = duckdb.connect()
        con.execute(f"""CREATE TABLE a AS SELECT "STUDY ACCESSION" AS study, PUBMEDID AS pmid,
              "FIRST AUTHOR" AS first_author, "DATE" AS published, JOURNAL AS journal,
              "DISEASE/TRAIT" AS trait_reported, MAPPED_TRAIT AS mapped_trait, MAPPED_TRAIT_URI AS mapped_trait_uri,
              SNPS AS rsid, "STRONGEST SNP-RISK ALLELE" AS snp_allele,
              split_part("STRONGEST SNP-RISK ALLELE", '-', 2) AS effect_allele,
              TRY_CAST("RISK ALLELE FREQUENCY" AS DOUBLE) AS effect_allele_freq,
              TRY_CAST(PVALUE_MLOG AS DOUBLE) AS mlog10p, "P-VALUE (TEXT)" AS p_text,
              TRY_CAST("OR or BETA" AS DOUBLE) AS or_beta, "95% CI (TEXT)" AS ci_text,
              "INITIAL SAMPLE SIZE" AS sample, "REPLICATION SAMPLE SIZE" AS replication, MAPPED_GENE AS gene,
              CHR_ID AS chrom38, TRY_CAST(CHR_POS AS BIGINT) AS pos38
            FROM read_csv('{tsv}', delim='\t', header=true, all_varchar=true, quote='', ignore_errors=true)
            WHERE SNPS IN (SELECT unnest(?::VARCHAR[]))""", [rsids])
        dst = out / "associations.parquet"
        con.execute(f"COPY (SELECT * FROM a ORDER BY rsid, mlog10p DESC) TO '{dst}' (FORMAT parquet)")
        total = con.execute(f"SELECT count(*) FROM read_csv('{tsv}', delim='\t', header=true, all_varchar=true, "
                            "quote='', ignore_errors=true)").fetchone()[0]

        chain = Chain(fetch(CHAIN, scratch / "chain.gz", upstream[1].size))
        rows = []
        for f in feat:
            hits = con.execute("""SELECT * FROM a WHERE rsid = ? AND mapped_trait = ? AND mlog10p >= ?
                                  ORDER BY (effect_allele IN ('A','C','G','T')) DESC, mlog10p DESC""",
                               [f["rsid"], f["mapped_trait"], GWS]).fetchdf().to_dict("records")
            if not hits:
                raise ValueError(f"{f['rsid']} / {f['mapped_trait']}: no genome-wide significant association")
            b = hits[0]
            lifted = chain.lift(f"chr{b['chrom38']}", int(b["pos38"])) if b["pos38"] else None
            rows.append({**{k: f[k] for k in ("rsid", "mapped_trait", "label", "section", "body_system")},
                         "mapped_trait_uri": b["mapped_trait_uri"], "gene": b["gene"],
                         "chrom38": b["chrom38"], "pos38": int(b["pos38"]) if b["pos38"] else None,
                         "chrom37": lifted[0].removeprefix("chr") if lifted else None,
                         "pos37": lifted[1] if lifted else None,
                         "effect_allele": b["effect_allele"], "effect_allele_freq": b["effect_allele_freq"],
                         "or_beta": b["or_beta"], "ci_text": b["ci_text"], "mlog10p": b["mlog10p"],
                         "best_study": b["study"], "best_pmid": b["pmid"], "best_first_author": b["first_author"],
                         "best_published": b["published"], "best_trait_reported": b["trait_reported"],
                         "best_sample": b["sample"],
                         "publications": len({h["pmid"] for h in hits}), "studies": len({h["study"] for h in hits})})
        fdst = out / "featured.parquet"
        pq.write_table(pa.Table.from_pylist(rows), fdst)
        version = (upstream[0].last_modified or "")[:10] or "unknown"
        return Built(version=version, tables={"associations": dst, "featured": fdst}, upstream=upstream,
                     notes={"citation": "Cerezo M, et al. The NHGRI-EBI GWAS Catalog: standards for reusability, "
                                        "sustainability and diversity. Nucleic Acids Res. 2025;53:D998-D1005.",
                            "catalog_associations": total, "featured": len(rows)})


class Sfari(Source):
    id = "sfari"
    title = "SFARI Gene"
    homepage = "https://gene.sfari.org"
    licence = "Free reuse with citation (SFARI Gene terms; cite Abrahams et al. Mol Autism 2013)"
    cadence = "quarterly"
    description = ("The Simons Foundation's expert-scored list of genes implicated in autism research. Score 1 = "
                   "high confidence (repeatedly shown), 2 = strong candidate, 3 = suggestive evidence; "
                   "'syndromic' genes cause a broader condition that often includes autism.")

    def probe(self, pin: str | None = None) -> list[Upstream]:
        # The CSV endpoint is generated on request and has no Last-Modified, so fingerprint by content.
        p = Path(fetch(SFARI, Path(self.workspace) / "cache" / "sfari" / "probe.csv"))
        data = p.read_bytes()
        p.unlink()
        return [Upstream(url=SFARI, etag=hashlib.sha256(data).hexdigest()[:16], size=len(data))]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        src = fetch(SFARI, scratch / "sfari.csv")
        dst = out / "genes.parquet"
        con = duckdb.connect()
        con.execute(f"""COPY (SELECT "gene-symbol" AS gene, "gene-name" AS name, "ensembl-id" AS ensembl_id,
              chromosome, "genetic-category" AS category, TRY_CAST("gene-score" AS INTEGER) AS score,
              "syndromic" = '1' AS syndromic, TRY_CAST("number-of-reports" AS INTEGER) AS reports, eagle
            FROM read_csv('{src}', header=true, all_varchar=true)) TO '{dst}' (FORMAT parquet)""")
        n = con.execute(f"SELECT count(*) FROM '{dst}'").fetchone()[0]
        if n < 1000:
            raise RuntimeError(f"SFARI list looks truncated ({n} genes)")
        version = datetime.date.today().isoformat()  # upstream has no date; content changes are caught by etag
        return Built(version=version, tables={"genes": dst}, upstream=upstream,
                     notes={"citation": "Abrahams BS, et al. SFARI Gene 2.0: a community-driven knowledgebase for "
                                        "the autism spectrum disorders (ASDs). Mol Autism. 2013;4:36.",
                            "genes": n})
