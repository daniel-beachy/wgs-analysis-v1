"""Curated disease knowledge: which diseases a gene causes, how they are inherited, how rare they are,
what they are, and which genes clinicians consider medically actionable.

* **HPO annotations** (Monarch / JAX): gene → disease links (OMIM, Orphanet) and each disease's mode
  of inheritance. Released roughly monthly.
* **Mondo** (Monarch): a unified disease ontology that cross-references OMIM, Orphanet and others,
  with plain-language definitions. Used to name and explain conditions and to translate ClinVar's
  condition IDs into the IDs HPO and Orphanet use.
* **Orphanet / Orphadata**: rare-disease inheritance, typical age of onset and prevalence.
* **ACMG SF v3.3**: the 84 genes the American College of Medical Genetics recommends reporting as
  secondary findings because something can be done about them (surveillance, prevention, treatment).
  The list changes about once a year by publication, so it is vendored as a small CSV with citation.

None of these say anything about a specific variant; they give the context needed to decide what a
variant ClinVar (or a predictor) flags would mean for you.
"""

from __future__ import annotations

import hashlib
import json
import re
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import duckdb

from .. import inheritance as inh
from . import Built, Source, Upstream
from .http import UA, fetch, head

HPO_BASE = "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download"
MONDO_OBO = "https://github.com/monarch-initiative/mondo/releases/latest/download/mondo.obo"
ORPHA = "https://www.orphadata.com/data/xml"
ACMG_CSV = Path(__file__).parent / "data" / "acmg_sf_v3.3.csv"


def _github_latest(url: str) -> Upstream:
    """GitHub 'latest' assets redirect to a signed URL; the release asset's Last-Modified/ETag identify it."""
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=60) as r:
        size = r.headers.get("Content-Length")
        return Upstream(url=url, last_modified=r.headers.get("Last-Modified"),
                        etag=(r.headers.get("ETag") or "").strip('"') or None, size=int(size) if size else None)


def _write_json_rows(con, rows: list[dict], path: Path, dst: Path, types: dict[str, str]) -> None:
    path.write_text(json.dumps(rows))
    cols = ", ".join(f"'{k}': '{v}'" for k, v in types.items())
    con.execute(f"COPY (SELECT * FROM read_json('{path}', format='array', columns={{{cols}}})) "
                f"TO '{dst}' (FORMAT parquet, COMPRESSION zstd)")


class HPO(Source):
    id = "hpo"
    title = "HPO disease annotations"
    homepage = "https://hpo.jax.org/"
    licence = "HPO licence (free to use with attribution; content must not be altered)"
    cadence = "monthly"
    description = ("Which diseases each gene causes (from OMIM and Orphanet via the Human Phenotype Ontology) "
                   "and how each disease is inherited.")

    def probe(self, pin: str | None = None) -> list[Upstream]:
        return [_github_latest(f"{HPO_BASE}/phenotype.hpoa"), _github_latest(f"{HPO_BASE}/genes_to_disease.txt")]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        hpoa = fetch(upstream[0].url, scratch / "phenotype.hpoa")
        g2d = fetch(upstream[1].url, scratch / "genes_to_disease.txt")
        version, skip = "unknown", 0
        with open(hpoa) as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                skip += 1
                if line.startswith("#version:"):
                    version = line.split(":", 1)[1].strip()
        con = duckdb.connect()
        case = "CASE hpo_id " + " ".join(f"WHEN '{k}' THEN '{v}'" for k, v in inh.HPO.items()) + " END"
        dst_inh = out / "disease_inheritance.parquet"
        con.execute(f"""
          COPY (
            WITH a AS (SELECT * FROM read_csv('{hpoa}', delim='\t', header=true, skip={skip}, quote='',
                                               all_varchar=true))
            SELECT database_id AS disease_id, any_value(disease_name) AS disease_name,
                   list_sort(list_distinct(list({case}) FILTER (WHERE aspect = 'I' AND {case} IS NOT NULL)))
                     AS modes
            FROM a GROUP BY 1 ORDER BY 1
          ) TO '{dst_inh}' (FORMAT parquet, COMPRESSION zstd)""")
        dst_gd = out / "gene_disease.parquet"
        con.execute(f"""
          COPY (
            SELECT gene_symbol AS gene, disease_id, association_type
            FROM read_csv('{g2d}', delim='\t', header=true, quote='', all_varchar=true)
            ORDER BY 1, 2
          ) TO '{dst_gd}' (FORMAT parquet, COMPRESSION zstd)""")
        return Built(version=version, tables={"disease_inheritance": dst_inh, "gene_disease": dst_gd},
                     upstream=upstream)


class Mondo(Source):
    id = "mondo"
    title = "Mondo disease ontology"
    homepage = "https://mondo.monarchinitiative.org/"
    licence = "CC BY 4.0"
    cadence = "monthly"
    description = ("A single disease vocabulary that links OMIM, Orphanet and other catalogues, with plain-"
                   "language definitions used to explain conditions.")

    def probe(self, pin: str | None = None) -> list[Upstream]:
        return [_github_latest(MONDO_OBO)]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        obo = fetch(upstream[0].url, scratch / "mondo.obo")
        version, terms, cur = "unknown", [], None
        with open(obo, encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("data-version:"):
                    version = line.split("/")[-1].strip() if "/" in line else line.split(":", 1)[1].strip()
                elif line == "[Term]":
                    cur = {"mondo_id": None, "name": None, "definition": None, "xrefs": [], "obsolete": False}
                    terms.append(cur)
                elif line.startswith("[") or not line:
                    cur = None if line.startswith("[") else cur
                elif cur is not None:
                    tag, _, val = line.partition(": ")
                    if tag == "id":
                        cur["mondo_id"] = val
                    elif tag == "name":
                        cur["name"] = val
                    elif tag == "def":
                        m = re.match(r'"((?:[^"\\]|\\.)*)"', val)
                        cur["definition"] = m.group(1).replace('\\"', '"') if m else None
                    elif tag == "xref":
                        x = val.split(" ", 1)[0]
                        if x.startswith(("OMIM:", "Orphanet:", "OMIMPS:")):
                            cur["xrefs"].append(x.replace("Orphanet:", "ORPHA:"))
                    elif tag == "is_obsolete" and val == "true":
                        cur["obsolete"] = True
        rows = [{k: v for k, v in t.items() if k != "obsolete"} for t in terms
                if t["mondo_id"] and t["mondo_id"].startswith("MONDO:") and not t["obsolete"]]
        con = duckdb.connect()
        dst = out / "diseases.parquet"
        _write_json_rows(con, rows, scratch / "mondo.json", dst,
                         {"mondo_id": "VARCHAR", "name": "VARCHAR", "definition": "VARCHAR", "xrefs": "VARCHAR[]"})
        return Built(version=version, tables={"diseases": dst}, upstream=upstream, notes={"terms": len(rows)})


class Orphanet(Source):
    id = "orphanet"
    title = "Orphanet rare diseases"
    homepage = "https://www.orphadata.com/"
    licence = "CC BY 4.0"
    cadence = "twice yearly"
    description = "Rare-disease genes, inheritance, typical age of onset and how common each disease is."

    FILES = {"genes": "en_product6.xml", "ages": "en_product9_ages.xml", "prev": "en_product9_prev.xml"}

    def probe(self, pin: str | None = None) -> list[Upstream]:
        return [head(f"{ORPHA}/{f}") for f in self.FILES.values()]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        paths = {k: fetch(f"{ORPHA}/{f}", scratch / f) for k, f in self.FILES.items()}
        dis: dict[str, dict] = {}
        genes: list[dict] = []

        def disorders(p: Path):
            for _, el in ET.iterparse(p, events=("end",)):
                if el.tag == "Disorder" and el.find("OrphaCode") is not None:
                    yield el
                    el.clear()

        def d(el) -> dict:
            code = f"ORPHA:{el.findtext('OrphaCode')}"
            return dis.setdefault(code, {"disease_id": code, "name": el.findtext("Name"), "modes": [],
                                         "inheritance": [], "onset": [], "prevalence": None})

        for el in disorders(paths["ages"]):
            x = d(el)
            x["onset"] = [n.findtext("Name") for n in el.iter("AverageAgeOfOnset")]
            names = [n.findtext("Name") for n in el.iter("TypeOfInheritance")]
            x["inheritance"] = names
            x["modes"] = sorted({m for m in map(inh.from_text, names) if m})
        for el in disorders(paths["prev"]):
            x = d(el)
            best = None
            for pv in el.iter("Prevalence"):
                cls = pv.findtext("PrevalenceClass/Name")
                typ = pv.findtext("PrevalenceType/Name") or ""
                geo = pv.findtext("PrevalenceGeographic/Name") or ""
                if not cls:
                    continue
                score = (not cls.lower().startswith(("unknown", "not yet"))) * 4 + (typ == "Point prevalence") * 2 \
                    + (geo == "Worldwide") * 1
                if best is None or score > best[0]:
                    best = (score, f"{cls} ({typ.lower()}, {geo})" if typ else cls)
            x["prevalence"] = best[1] if best else None
        for el in disorders(paths["genes"]):
            x = d(el)
            genes.extend({"gene": a.findtext("Gene/Symbol"), "disease_id": x["disease_id"],
                          "association_type": a.findtext("DisorderGeneAssociationType/Name"),
                          "status": a.findtext("DisorderGeneAssociationStatus/Name")}
                         for a in el.iter("DisorderGeneAssociation"))
        con = duckdb.connect()
        dst_d, dst_g = out / "disorders.parquet", out / "gene_disorder.parquet"
        _write_json_rows(con, list(dis.values()), scratch / "dis.json", dst_d,
                         {"disease_id": "VARCHAR", "name": "VARCHAR", "modes": "VARCHAR[]",
                          "inheritance": "VARCHAR[]", "onset": "VARCHAR[]", "prevalence": "VARCHAR"})
        _write_json_rows(con, genes, scratch / "genes.json", dst_g,
                         {"gene": "VARCHAR", "disease_id": "VARCHAR", "association_type": "VARCHAR",
                          "status": "VARCHAR"})
        lm = max((u.last_modified or "") for u in upstream)
        return Built(version=lm[:10] or "unknown", tables={"disorders": dst_d, "gene_disorder": dst_g},
                     upstream=upstream, notes={"disorders": len(dis), "gene_links": len(genes)})


class AcmgSF(Source):
    id = "acmg_sf"
    title = "ACMG secondary findings list"
    homepage = "https://www.acmg.net/ACMG/Medical-Genetics-Practice-Resources/Practice-Guidelines.aspx"
    licence = "Gene list from a published guideline (Lee et al., Genet Med 2025; ACMG SF v3.3), cited"
    cadence = "about yearly (by publication; vendored with the tool)"
    description = ("The genes clinical labs report even when you were not tested for them, because something "
                   "can be done (screening, prevention or treatment) if you carry a harmful variant.")

    def probe(self, pin: str | None = None) -> list[Upstream]:
        data = ACMG_CSV.read_bytes()
        return [Upstream(url=f"vendored:{ACMG_CSV.name}", etag=hashlib.sha256(data).hexdigest()[:16],
                         size=len(data))]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        dst = out / "genes.parquet"
        duckdb.connect().execute(f"COPY (SELECT * FROM read_csv('{ACMG_CSV}', header=true, all_varchar=true)) "
                                 f"TO '{dst}' (FORMAT parquet)")
        return Built(version="v3.3", tables={"genes": dst}, upstream=upstream,
                     notes={"citation": "Lee K, et al. ACMG SF v3.3 list for reporting of secondary findings in "
                                        "clinical exome and genome sequencing. Genet Med. 2025."})
