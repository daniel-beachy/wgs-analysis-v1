"""What genes (and genetic conditions) are, in words: so the dashboard never needs hand-written gene blurbs.

* **MedlinePlus Genetics** (US National Library of Medicine): plain-language, expert-reviewed summaries
  of ~1,500 genes and ~1,300 genetic conditions, written for the public. Bulk XML, updated continuously.
* **NCBI Gene / RefSeq summaries**: a technical one-paragraph summary for almost every human gene,
  used when MedlinePlus has no page (e.g. CYP2D6, NAT2).

Every gene card in the dashboard shows the MedlinePlus text if available, otherwise the RefSeq summary,
with a link and the review date, so the wording is sourced, refreshable and never typed in by hand.
"""

from __future__ import annotations

import gzip
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import duckdb

from . import Built, Source, Upstream
from .curated import _write_json_rows
from .http import fetch, head

MEDLINEPLUS_XML = "https://medlineplus.gov/download/ghr-summaries.xml"
NCBI_SUMMARY = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_summary.gz"
NCBI_INFO = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"


def _text(el: ET.Element | None) -> str:
    """Flatten an XHTML fragment to plain text, one paragraph per line."""
    if el is None:
        return ""
    paras = []
    for p in el.iter():
        if p.tag.endswith("}p") or p.tag.endswith("}li"):
            t = re.sub(r"\s+", " ", "".join(p.itertext())).strip()
            if t:
                paras.append(t)
    return "\n".join(paras) or re.sub(r"\s+", " ", "".join(el.itertext())).strip()


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def parse_medlineplus(path: Path) -> tuple[list[dict], list[dict]]:
    genes, conditions = [], []
    for _, el in ET.iterparse(path, events=("end",)):
        kind = _local(el.tag)
        if kind not in ("gene-summary", "health-condition-summary"):
            continue
        f = {_local(c.tag): c for c in el}
        texts = {}
        for t in (f.get("text-list") if f.get("text-list") is not None else []):
            role = next((c.text for c in t if _local(c.tag) == "text-role"), None)
            html = next((c for c in t if _local(c.tag) == "html"), None)
            texts[role] = _text(html)
        name = (f["name"].text or "").strip() if "name" in f else ""
        url = (f["ghr-page"].text or "").strip() if "ghr-page" in f else ""
        reviewed = (f["reviewed"].text or "").strip() if "reviewed" in f else None
        if kind == "gene-summary":
            ncbi = None
            for k in el.iter():
                if _local(k.tag) == "db-key":
                    kv = {_local(c.tag): c.text for c in k}
                    if kv.get("db") == "NCBI Gene":
                        ncbi = kv.get("key")
            conds = [n.text for rc in el.iter() if _local(rc.tag) == "related-health-condition"
                     for n in rc if _local(n.tag) == "name" and n.text]
            genes.append({"gene": (f["gene-symbol"].text or "").strip(), "name": name, "url": url,
                          "function": texts.get("function", ""), "reviewed": reviewed, "ncbi_id": ncbi,
                          "conditions": conds})
        else:
            conditions.append({"name": name, "url": url, "description": texts.get("description", ""),
                               "reviewed": reviewed})
        el.clear()
    return genes, conditions


class MedlinePlus(Source):
    id = "medlineplus"
    title = "MedlinePlus Genetics"
    homepage = "https://medlineplus.gov/genetics/"
    licence = "Public domain (US government work); attribution requested"
    cadence = "continuous"
    description = ("Plain-language, expert-reviewed explanations of what each gene does and what each genetic "
                   "condition is, from the US National Library of Medicine. Used for every gene card.")

    def probe(self, pin: str | None = None) -> list[Upstream]:
        return [head(MEDLINEPLUS_XML)]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        xml = fetch(upstream[0].url, scratch / "ghr-summaries.xml")
        genes, conditions = parse_medlineplus(xml)
        if len(genes) < 1000:
            raise RuntimeError(f"MedlinePlus XML looks truncated ({len(genes)} genes)")
        con = duckdb.connect()
        dg, dc = out / "medlineplus_genes.parquet", out / "medlineplus_conditions.parquet"
        _write_json_rows(con, genes, scratch / "g.json", dg,
                         {"gene": "VARCHAR", "name": "VARCHAR", "url": "VARCHAR", "function": "VARCHAR",
                          "reviewed": "VARCHAR", "ncbi_id": "VARCHAR", "conditions": "VARCHAR[]"})
        _write_json_rows(con, conditions, scratch / "c.json", dc,
                         {"name": "VARCHAR", "url": "VARCHAR", "description": "VARCHAR", "reviewed": "VARCHAR"})
        latest = max((g["reviewed"] or "" for g in genes), default="")
        version = (upstream[0].last_modified or latest or "unknown")
        return Built(version=_short_date(version), tables={"medlineplus_genes": dg, "medlineplus_conditions": dc},
                     upstream=upstream, notes={"genes": len(genes), "conditions": len(conditions)})


class NcbiGene(Source):
    id = "ncbi_gene"
    title = "NCBI Gene (RefSeq summaries)"
    homepage = "https://www.ncbi.nlm.nih.gov/gene/"
    licence = "Public domain (NCBI)"
    schema = 2
    cadence = "daily"
    description = ("Official human gene names and the RefSeq curators' one-paragraph summary of each gene's "
                   "function. Used when MedlinePlus has no plain-language page.")

    def probe(self, pin: str | None = None) -> list[Upstream]:
        return [head(NCBI_INFO), head(NCBI_SUMMARY)]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        info = fetch(upstream[0].url, scratch / "Homo_sapiens.gene_info.gz")
        summ = fetch(upstream[1].url, scratch / "gene_summary.gz")
        human = scratch / "gene_summary_9606.tsv"
        with gzip.open(summ, "rt", encoding="utf-8", errors="replace") as src, open(human, "w") as dst:
            dst.write(src.readline())
            for line in src:
                if line.startswith("9606\t"):
                    dst.write(line)
        dst_p = out / "ncbi_genes.parquet"
        duckdb.connect().execute(f"""
          COPY (
            WITH i AS (SELECT * FROM read_csv('{info}', delim='\t', header=true, quote='', all_varchar=true)),
                 s0 AS (SELECT * FROM read_csv('{human}', delim='\t', header=true, quote='', all_varchar=true)),
                 -- one summary per gene: prefer RefSeq curators, then OMIM, then machine-written (Alliance)
                 s AS (SELECT GeneID, arg_min(Summary, CASE Source WHEN 'RefSeq' THEN 0 WHEN 'OMIM' THEN 1 ELSE 2 END)
                              AS Summary,
                              arg_min(Source, CASE Source WHEN 'RefSeq' THEN 0 WHEN 'OMIM' THEN 1 ELSE 2 END) AS Source
                       FROM s0 GROUP BY 1)
            SELECT CASE WHEN i.Symbol_from_nomenclature_authority NOT IN ('-', '')
                        THEN i.Symbol_from_nomenclature_authority ELSE i.Symbol END AS gene,
                   i.GeneID AS ncbi_id, i.description AS name, i.type_of_gene AS type,
                   i.chromosome AS chrom, s.Summary AS summary, s.Source AS summary_source
            FROM i LEFT JOIN s ON s.GeneID = i.GeneID
            WHERE i.Symbol IS NOT NULL ORDER BY 1
          ) TO '{dst_p}' (FORMAT parquet, COMPRESSION zstd)""")
        return Built(version=_short_date(upstream[0].last_modified or "unknown"), tables={"ncbi_genes": dst_p},
                     upstream=upstream)


def _short_date(http_date: str) -> str:
    """'Tue, 06 Oct 2026 05:15:58 GMT' -> '2026-10-06' (anything else passes through)."""
    from email.utils import parsedate_to_datetime

    if re.match(r"\d{4}-\d{2}-\d{2}", http_date or ""):
        return http_date[:10]
    try:
        return parsedate_to_datetime(http_date).date().isoformat()
    except (TypeError, ValueError):
        return re.sub(r"[^\w.-]", "-", http_date)
