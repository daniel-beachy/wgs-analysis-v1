"""Annotate: join your variants with local knowledge, then turn matches into graded claims (ADR-012/013).

1. ``bcftools csq`` (Ensembl gene models) predicts each variant's effect on genes/proteins.
2. ``slivar`` adds gnomAD v2.1.1 population-max frequencies (when that source is installed).
3. DuckDB joins ClinVar assertions and 1000 Genomes frequencies by chrom/pos/ref/alt.
4. Each ClinVar match that says something about health, drugs or traits becomes a *claim* with
   an evidence grade, a call-confidence grade and the reasons for both (see ``wgs.evidence``).

Inputs include the current version of every knowledge table, so ``wgs knowledge refresh``
followed by ``wgs run`` re-annotates automatically and the next release records what changed.
Every source is optional: with none installed the stage is skipped; with some, it does what it can.
"""

from __future__ import annotations

import gzip
import json
import re
import shutil
import subprocess
from pathlib import Path

import duckdb

from .. import evidence as ev
from ..knowledge import Store
from ..pipeline import CANONICAL, Context, Stage, run, write_json
from . import genotypes, reference, variants

CSQ_RANK = {
    # HIGH
    "frameshift": 9, "stop_gained": 9, "splice_acceptor": 9, "splice_donor": 9, "start_lost": 8, "stop_lost": 8,
    # MODERATE
    "missense": 6, "inframe_deletion": 6, "inframe_insertion": 6, "inframe_altering": 6,
    # LOW
    "splice_region": 4, "synonymous": 3, "stop_retained": 3, "start_retained": 3, "coding_sequence": 3,
    "feature_elongation": 3,
    # MODIFIER
    "5_prime_utr": 2, "3_prime_utr": 2, "non_coding": 1, "intron": 1,
}
IMPACT_SQL = ("CASE WHEN r >= 8 THEN 'HIGH' WHEN r >= 6 THEN 'MODERATE' WHEN r >= 3 THEN 'LOW' "
              "ELSE 'MODIFIER' END")
CLAIM_CLASSES = ("pathogenic", "likely_pathogenic", "drug_response", "risk_factor", "protective", "association",
                 "conflicting")
SIG_LABEL = {"pathogenic": "Pathogenic", "likely_pathogenic": "Likely pathogenic", "drug_response": "Drug response",
             "risk_factor": "Risk factor", "protective": "Protective", "association": "Association",
             "conflicting": "Conflicting"}
# Genes whose results people may not want to see unprompted (ACMG SF-style or high-impact risk genes).
SENSITIVE_GENES = {"APOE", "BRCA1", "BRCA2", "HTT", "PSEN1", "PSEN2", "APP", "MLH1", "MSH2", "MSH6", "PMS2",
                   "TP53", "PRNP", "C9orf72"}


def out_dir(ctx: Context) -> Path:
    p = ctx.work / "annotate"
    p.mkdir(parents=True, exist_ok=True)
    return p


def annotations_path(ctx: Context) -> Path:
    return out_dir(ctx) / "annotations.parquet"


def claims_path(ctx: Context) -> Path:
    return out_dir(ctx) / "claims.parquet"


def summary_path(ctx: Context) -> Path:
    return out_dir(ctx) / "annotate.json"


def _knowledge(ctx: Context) -> dict[str, Path | None]:
    s = Store(ctx.cfg)
    gn = s.table("gnomad", "fields")
    return {
        "genes": s.table("ensembl", "genes"),
        "gff": (s.table("ensembl", "genes").parent / "annotation.gff3.gz") if s.table("ensembl", "genes") else None,
        "clinvar": s.table("clinvar", "variants"),
        "validity": s.table("clingen", "gene_validity"),
        "dosage": s.table("clingen", "dosage"),
        "kg": s.table("1000genomes", "af"),
        "gnomad": (gn.parent / "gnomad.hg37.zip") if gn else None,
    }


def _inputs(ctx: Context) -> list[Path]:
    k = _knowledge(ctx)
    return [variants.variants_parquet(ctx)] + [p for p in k.values() if p] + \
        ([genotypes.out_path(ctx)] if genotypes.out_path(ctx).exists() else [])


def _available(ctx: Context) -> str | None:
    if not variants.vcf_file(ctx) or not variants.variants_parquet(ctx).exists():
        return "no VCF / variant table"
    k = _knowledge(ctx)
    if not k["clinvar"] and not k["gff"]:
        return "no knowledge installed — run `wgs knowledge refresh`"
    return None


def _vcf_naming(ctx: Context) -> dict[str, str]:
    """canonical name → name used in the VCF (e.g. '1' → 'chr1', 'MT' → 'chrM')."""
    out = subprocess.run(["bcftools", "view", "-h", variants.vcf_with_index(ctx)], capture_output=True, text=True,
                         check=True).stdout
    names = {ln.split("ID=", 1)[1].split(",", 1)[0].rstrip(">") for ln in out.splitlines()
             if ln.startswith("##contig=<ID=")}
    m = {}
    for c in CANONICAL:
        for cand in (["chrM", "chrMT", "MT", "M"] if c == "MT" else [f"chr{c}", c]):
            if cand in names:
                m[c] = cand
                break
    return m


def _gff_for(ctx: Context, gff: Path, naming: dict[str, str]) -> Path:
    """Ensembl uses '1'…'MT'; rewrite contig names to match the VCF/reference when they differ."""
    if all(k == v for k, v in naming.items()):
        return gff
    dst = ctx.scratch / f"genes.{gff.stat().st_size}.gff3.gz"
    if dst.exists():
        return dst
    tmp = dst.with_suffix(".tmp")
    with gzip.open(gff, "rt") as fh, gzip.open(tmp, "wt", compresslevel=1) as o:
        for line in fh:
            if line.startswith("#"):
                o.write(line)
                continue
            c, rest = line.split("\t", 1)
            if c in naming:
                o.write(f"{naming[c]}\t{rest}")
    tmp.replace(dst)
    return dst


def _gnomad_fields(zip_path: Path) -> list[str]:
    meta = json.loads((zip_path.parent / "source.json").read_text())
    return [f for f in meta.get("notes", {}).get("fields", []) if f.startswith("gnomad_")]


def _consequences(ctx: Context, k: dict, log: Path) -> tuple[Path | None, list[str]]:
    """Run csq (+ slivar) over alt-carrying records; return a TSV of per-allele annotations."""
    fasta = reference.reference_fasta(ctx)
    naming = _vcf_naming(ctx)
    tmp = ctx.scratch / "annotate"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    bcf = tmp / "alt.bcf"
    vcf = variants.vcf_with_index(ctx)
    cmd = f"bcftools view -i 'GT=\"alt\"' -Ou '{vcf}' | bcftools norm -m -any -Ou"
    have_csq = bool(k["gff"] and fasta)
    if have_csq:
        gff = _gff_for(ctx, k["gff"], naming)
        cmd += f" | bcftools csq -l -f '{fasta}' -g '{gff}' -Ou"
    run(f"{cmd} | bcftools view -Ob -o '{bcf}' --write-index", shell=True, log=log)
    gfields: list[str] = []
    if k["gnomad"]:
        gfields = _gnomad_fields(k["gnomad"])
        out = tmp / "alt.gnomad.bcf"
        run(["slivar", "expr", "--vcf", str(bcf), "--gnotate", str(k["gnomad"]), "-o", str(out)], log=log)
        bcf = out
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t" + ("%INFO/BCSQ" if have_csq else ".") + \
        "".join(f"\t%INFO/{f}" for f in gfields) + "\n"
    tsv = tmp / "annot.tsv"
    run(f"bcftools query -f '{fmt}' '{bcf}' > '{tsv}'", shell=True, log=log)
    return tsv, gfields


def build(ctx: Context) -> dict:
    log = ctx.logs / "annotate.log"
    k = _knowledge(ctx)
    tsv, gfields = _consequences(ctx, k, log)
    rev = {v: c for c, v in _vcf_naming(ctx).items()}
    con = duckdb.connect()
    con.execute(f"SET threads={ctx.threads}")
    con.execute("SET preserve_insertion_order=false")
    cols = {"chrom_raw": "VARCHAR", "pos": "INTEGER", "ref": "VARCHAR", "alt": "VARCHAR", "bcsq": "VARCHAR"}
    cols.update({f: "VARCHAR" for f in gfields})
    con.execute(f"CREATE TEMP TABLE raw AS SELECT * FROM read_csv('{tsv}', delim='\t', header=false, quote='', "
                f"columns={json.dumps(cols)})")
    con.execute("CREATE TEMP TABLE chrmap(raw VARCHAR, chrom VARCHAR)")
    con.executemany("INSERT INTO chrmap VALUES (?, ?)", [[a, b] for a, b in rev.items()])
    rank = "CASE t " + " ".join(f"WHEN '{t}' THEN {r}" for t, r in CSQ_RANK.items()) + " ELSE 0 END"
    # One row per allele with the most severe consequence over all transcripts (protein-coding preferred).
    con.execute(f"""
      CREATE TEMP TABLE csq AS
      WITH e AS (
        SELECT m.chrom, r.pos, r.ref, r.alt, unnest(string_split(r.bcsq, ',')) AS entry
        FROM raw r JOIN chrmap m ON m.raw = r.chrom_raw WHERE r.bcsq IS NOT NULL AND r.bcsq <> '.'
      ), p AS (
        SELECT *, string_split(entry, '|') AS f FROM e WHERE entry NOT LIKE '@%'
      ), s AS (
        SELECT chrom, pos, ref, alt, ltrim(f[1], '*') AS csq, f[2] AS gene, f[3] AS transcript, f[4] AS biotype,
               f[6] AS aa_change, f[7] AS dna_change,
               list_max(list_transform(string_split(ltrim(f[1], '*'), '&'), t -> {rank})) AS r
        FROM p
      )
      SELECT chrom, pos, ref, alt,
             arg_max(csq, r * 10 + (biotype = 'protein_coding')::INT) AS consequence,
             arg_max(gene, r * 10 + (biotype = 'protein_coding')::INT) AS gene,
             arg_max(transcript, r * 10 + (biotype = 'protein_coding')::INT) AS transcript,
             arg_max(nullif(aa_change, ''), r * 10 + (biotype = 'protein_coding')::INT) AS aa_change,
             arg_max(nullif(dna_change, ''), r * 10 + (biotype = 'protein_coding')::INT) AS dna_change,
             max(r) AS r,
             list_distinct(list(gene) FILTER (WHERE gene <> '')) AS genes
      FROM s GROUP BY ALL""")
    # slivar writes -1 when gnomAD has no value for an allele: that is "unknown", not "ultra-rare"
    gsel = ", ".join(f"CASE WHEN TRY_CAST(NULLIF(r.{f}, '.') AS DOUBLE) >= 0 "
                     f"THEN TRY_CAST(NULLIF(r.{f}, '.') AS DOUBLE) END AS {f}" for f in gfields)
    con.execute(f"""
      CREATE TEMP TABLE gn AS SELECT m.chrom, r.pos, r.ref, r.alt{', ' + gsel if gsel else ''}
      FROM raw r JOIN chrmap m ON m.raw = r.chrom_raw""")
    vp = variants.variants_parquet(ctx)
    gp = genotypes.out_path(ctx)
    rs_join = (f"LEFT JOIN (SELECT chrom, pos, min(rsid) rsid, min(genotype) provider_gt FROM '{gp}' "
               f"WHERE rsid LIKE 'rs%' GROUP BY ALL) g USING (chrom, pos)") if gp.exists() else ""
    rs_cols = "g.rsid AS provider_rsid, g.provider_gt" if gp.exists() else \
        "NULL::VARCHAR AS provider_rsid, NULL::VARCHAR AS provider_gt"
    cv = k["clinvar"]
    # Records with no germline classification (e.g. somatic/oncogenicity-only entries) say nothing about you
    cv_join = (f"LEFT JOIN (SELECT * FROM '{cv}' WHERE significance IS NOT NULL) cv USING (chrom, pos, ref, alt)"
               if cv else "")
    cv_cols = ("cv.variation_id AS clinvar_id, cv.significance AS clinvar_significance, cv.sig_class AS "
               "clinvar_class, cv.stars AS clinvar_stars, cv.review_status AS clinvar_review, cv.conditions AS "
               "clinvar_conditions, cv.low_penetrance AS clinvar_low_penetrance, cv.conflicting_detail AS "
               "clinvar_conflicts, cv.rsid AS clinvar_rsid") if cv else \
        ("NULL::BIGINT clinvar_id, NULL::VARCHAR clinvar_significance, NULL::VARCHAR clinvar_class, "
         "NULL::INT clinvar_stars, NULL::VARCHAR clinvar_review, NULL::VARCHAR clinvar_conditions, "
         "NULL::BOOLEAN clinvar_low_penetrance, NULL::VARCHAR clinvar_conflicts, NULL::VARCHAR clinvar_rsid")
    kg = k["kg"]
    kg_join = f"LEFT JOIN '{kg}' kg USING (chrom, pos, ref, alt)" if kg else ""
    kg_cols = ("kg.af AS af_1kg, kg.afr AS af_1kg_afr, kg.amr AS af_1kg_amr, kg.eas AS af_1kg_eas, "
               "kg.eur AS af_1kg_eur, kg.sas AS af_1kg_sas") if kg else \
        ", ".join(f"NULL::FLOAT {c}" for c in ("af_1kg", "af_1kg_afr", "af_1kg_amr", "af_1kg_eas", "af_1kg_eur",
                                               "af_1kg_sas"))
    gn_cols = ", ".join([f"gn.{f}" for f in gfields] +
                        ([] if "gnomad_popmax_af" in gfields else ["NULL::DOUBLE AS gnomad_popmax_af"]))
    # csq reports nothing for alleles outside every gene: call those intergenic (only if csq actually ran)
    csq_col = "coalesce(c.consequence, 'intergenic')" if k["gff"] and reference.reference_fasta(ctx) \
        else "c.consequence"
    out = annotations_path(ctx)
    tmp = out.with_suffix(".tmp")
    con.execute(f"""
      COPY (
        SELECT v.chrom, v.chrom_order, v.pos, coalesce({'cv.rsid, ' if cv else ''}{'g.rsid, ' if gp.exists() else ''}
               nullif(v.vcf_id, '.')) AS rsid,
               v.ref, v.alt, v.vtype, v.filter, v.gt, v.zygosity, v.gq, v.dp, v.vaf,
               c.gene, c.genes, {csq_col} AS consequence, {IMPACT_SQL.replace('r >=', 'c.r >=')} AS impact,
               c.transcript,
               c.aa_change, c.dna_change,
               {cv_cols}, {kg_cols}, {gn_cols}, {rs_cols}
        FROM '{vp}' v
        LEFT JOIN csq c USING (chrom, pos, ref, alt)
        LEFT JOIN gn USING (chrom, pos, ref, alt)
        {cv_join} {kg_join} {rs_join}
        WHERE v.zygosity IN ('het', 'hom', 'hemi')
        ORDER BY v.chrom_order, v.pos, v.alt
      ) TO '{tmp}' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 100000)""")
    tmp.replace(out)
    claims = _claims(ctx, con, k)
    stats = con.execute(f"""
      SELECT count(*), count(consequence), count(*) FILTER (WHERE impact = 'HIGH'),
             count(*) FILTER (WHERE impact = 'MODERATE'), count(clinvar_id), count(af_1kg),
             count(*) FILTER (WHERE af_1kg IS NULL AND filter = 'PASS')
      FROM '{out}'""").fetchone()
    by_class = dict(con.execute(f"SELECT clinvar_class, count(*) FROM '{out}' WHERE clinvar_class IS NOT NULL "
                                f"GROUP BY 1 ORDER BY 2 DESC").fetchall())
    summary = {
        "evidence_model": ev.MODEL_VERSION,
        "knowledge": Store(ctx.cfg).versions(),
        "alleles": stats[0], "with_consequence": stats[1], "high_impact": stats[2], "moderate_impact": stats[3],
        "in_clinvar": stats[4], "in_1000genomes": stats[5], "pass_not_in_1000genomes": stats[6],
        "clinvar_classes": by_class, "claims": claims,
    }
    write_json(summary_path(ctx), summary)
    shutil.rmtree(ctx.scratch / "annotate", ignore_errors=True)
    return {k2: v for k2, v in summary.items() if not isinstance(v, dict)}


def _provider_agrees(row: dict) -> str | None:
    pg = row.get("provider_gt")
    ref, alt = row["ref"], row["alt"]
    if not pg or len(ref) != 1 or len(alt) != 1 or len(pg) != 2 or not set(pg) <= set("ACGT"):
        return None
    expect = ref + alt if row["zygosity"] == "het" else alt * 2
    return "agree" if sorted(expect) == sorted(pg) else "disagree"


def _claims(ctx: Context, con, k: dict) -> int:
    """ClinVar-backed claims about your alleles. Later steps add curated PGx/trait/PRS claims alongside."""
    out = claims_path(ctx)
    if not k["clinvar"]:
        con.execute(f"COPY (SELECT NULL::VARCHAR claim_id WHERE false) TO '{out}' (FORMAT parquet)")
        return 0
    validity, moi = {}, {}
    if k["validity"]:
        for gene, cls, m in con.execute(f"SELECT gene_symbol, classification, moi FROM '{k['validity']}'").fetchall():
            r = ev.GENE_VALIDITY_RANK.get(cls, -1)
            if r > ev.GENE_VALIDITY_RANK.get(validity.get(gene), -1):
                validity[gene] = cls
            moi.setdefault(gene, set()).add(m)
    classes = ", ".join(f"'{c}'" for c in CLAIM_CLASSES)
    rows = con.execute(f"""
      SELECT * FROM '{annotations_path(ctx)}'
      WHERE clinvar_class IN ({classes})
        AND (clinvar_class <> 'conflicting' OR clinvar_conflicts ILIKE '%pathogenic%')
      ORDER BY chrom_order, pos""").fetchdf().to_dict("records")
    versions = Store(ctx.cfg).versions()
    claims = []
    for r in rows:
        r = {key: (None if _isnan(val) else val) for key, val in r.items()}
        gene = r["gene"] or (r["genes"][0] if r.get("genes") is not None and len(r["genes"]) else None)
        cls = r["clinvar_class"]
        popmax = r.get("gnomad_popmax_af")
        e_level, e_reasons = ev.clinvar_evidence(r["clinvar_stars"], cls, popmax_af=popmax, kg_af=r["af_1kg"],
                                                 gene_validity=validity.get(gene),
                                                 low_penetrance=bool(r["clinvar_low_penetrance"]))
        c_level, c_reasons = ev.call_confidence(filter=r["filter"], gq=r["gq"], dp=r["dp"], vaf=r["vaf"],
                                                zygosity=r["zygosity"], provider=_provider_agrees(r))
        copies = 2 if r["zygosity"] == "hom" else 1
        gene_moi = moi.get(gene, set())
        recessive = "AR" in gene_moi and not gene_moi & {"AD", "SD", "XL"}
        if cls == "drug_response":
            section = "pgx"
        elif cls == "association":
            section = "traits"
        elif cls in ("pathogenic", "likely_pathogenic", "conflicting") and copies == 1 and recessive:
            section = "carrier"
        else:
            section = "health"
        conditions = _conditions(r["clinvar_conditions"])
        statement = _statement(cls, conditions, copies, recessive, r)
        major = _alt_is_major(r)
        if major:
            statement += _major_note(major, copies)
        sources = [{"source": "clinvar", "version": versions.get("clinvar"), "record": str(r["clinvar_id"]),
                    "url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{r['clinvar_id']}/"}]
        if popmax is not None:
            sources.append({"source": "gnomad", "version": versions.get("gnomad"), "record": "popmax_af", "url": ""})
        if r["af_1kg"] is not None:
            sources.append({"source": "1000genomes", "version": versions.get("1000genomes"), "record": "AF",
                            "url": ""})
        if gene in validity:
            sources.append({"source": "clingen", "version": versions.get("clingen"), "record": gene,
                            "url": f"https://search.clinicalgenome.org/kb/genes?search={gene}"})
        claims.append({
            "claim_id": f"variant:{r['chrom']}-{r['pos']}-{r['ref']}-{r['alt']}:clinvar",
            "section": section, "kind": "clinvar_variant", "category": cls, "category_label": SIG_LABEL[cls],
            "subject": gene or r["rsid"] or f"{r['chrom']}:{r['pos']}", "gene": gene, "rsid": r["rsid"],
            "chrom": r["chrom"], "pos": int(r["pos"]), "ref": r["ref"], "alt": r["alt"],
            "genotype": _genotype(r), "copies": copies, "zygosity": r["zygosity"],
            "consequence": r["consequence"], "impact": r["impact"], "aa_change": r["aa_change"],
            "conditions": conditions, "statement": statement,
            "clinvar_significance": r["clinvar_significance"], "clinvar_stars": int(r["clinvar_stars"] or 0),
            "gene_validity": validity.get(gene), "inheritance": "/".join(sorted(m for m in gene_moi if m)) or None,
            "af_1kg": r["af_1kg"], "af_1kg_eur": r["af_1kg_eur"], "gnomad_popmax_af": popmax,
            "evidence_level": e_level, "evidence_reasons": e_reasons,
            "call_level": c_level, "call_reasons": c_reasons,
            "overall": ev.overall(e_level, c_level),
            "alt_is_major": major is not None,
            "sensitive": cls in ("pathogenic", "likely_pathogenic") or gene in SENSITIVE_GENES,
            "sources": json.dumps(sources),
        })
    import pandas as pd

    df = pd.DataFrame(claims)
    tmp = out.with_suffix(".tmp")
    if df.empty:
        con.execute(f"COPY (SELECT NULL::VARCHAR claim_id WHERE false) TO '{tmp}' (FORMAT parquet)")
    else:
        con.register("claims_df", df)
        con.execute(f"""COPY (SELECT * REPLACE (pos::INTEGER AS pos, copies::TINYINT AS copies,
                                 clinvar_stars::TINYINT AS clinvar_stars, af_1kg::DOUBLE AS af_1kg,
                                 af_1kg_eur::DOUBLE AS af_1kg_eur, gnomad_popmax_af::DOUBLE AS gnomad_popmax_af)
                        FROM claims_df ORDER BY section, claim_id) TO '{tmp}' (FORMAT parquet)""")
    tmp.replace(out)
    return len(claims)


FILLER = {"not provided", "not specified", "see cases", "other", "none provided"}


def _conflicts(raw: str | None) -> str:
    """'Pathogenic(2)|Benign_(4)' → 'Pathogenic 2, Benign 4' (stable across ClinVar's formatting changes)."""
    out = []
    for part in (raw or "").split("|"):
        m = re.match(r"\s*(.+?)\s*\((\d+)\)\s*$", part.replace("_", " "))
        if m:
            out.append(f"{m.group(1).capitalize()} {m.group(2)}")
        elif part.strip():
            out.append(part.replace("_", " ").strip())
    return ", ".join(out) or "details not given"


def _conditions(raw: str | None) -> str:
    out = []
    for c in (raw or "").split("|"):
        c = c.strip()
        if not c or c.lower() in FILLER:
            continue
        c = c.capitalize() if c.isupper() else c
        if c not in out:
            out.append(c)
    # ClinVar's condition order is arbitrary between releases; sort so reordering is not reported as a change.
    return "; ".join(sorted(out, key=str.casefold)) or "an unspecified condition"


def _isnan(v) -> bool:
    return isinstance(v, float) and v != v


def _genotype(r: dict) -> str:
    a, b = (r["ref"], r["alt"]) if r["zygosity"] == "het" else (r["alt"], r["alt"])
    if r["zygosity"] == "hemi":
        return r["alt"]
    return f"{_short(a)}/{_short(b)}"


def _short(s: str) -> str:
    return s if len(s) <= 6 else f"{s[:3]}…({len(s)})"


def _alt_is_major(r: dict) -> tuple[float, str] | None:
    """Population frequency of the ALT allele when it is the *majority* allele, else None.

    GRCh37 was assembled from a handful of people, so at some sites the reference carries the rare allele
    (Factor V Leiden, rs6025, is the classic example). There the "variant" ClinVar describes is the version
    most people have, and "you carry two copies" means "you are typical"."""
    for af, label in ((r.get("af_1kg"), "1000 Genomes, all populations"),
                      (r.get("gnomad_popmax_af"), "gnomAD, highest population")):
        if af is not None:
            return (float(af), label) if af >= 0.5 else None
    return None


def _major_note(major: tuple[float, str], copies: int) -> str:
    af, label = major
    s = (f" Note: this is the majority allele, on about {af:.0%} of chromosomes ({label}); the reference genome "
         "happens to carry the rarer one.")
    if copies == 2:
        s += (" Two copies means your genotype here is the typical one; the effect usually discussed for this "
              "site belongs to the rarer allele, which you do not carry.")
    else:
        s += " With one copy you also carry the rarer reference allele on your other chromosome."
    return s


def _statement(cls: str, conditions: str, copies: int, recessive: bool, r: dict) -> str:
    copy_txt = "two copies (both chromosomes)" if copies == 2 else "one copy"
    if cls in ("pathogenic", "likely_pathogenic"):
        label = "pathogenic" if cls == "pathogenic" else "likely pathogenic"
        s = f"ClinVar lists this variant as {label} for {conditions}. You carry {copy_txt}."
        if copies == 1 and recessive:
            s += " The condition is recessive, so one copy usually means you are a carrier, not affected."
        elif copies == 1:
            s += " Whether one copy matters depends on how the condition is inherited and its penetrance."
        return s
    if cls == "conflicting":
        return (f"Labs disagree about this variant ({_conflicts(r.get('clinvar_conflicts'))}) "
                f"for {conditions}. You carry {copy_txt}.")
    if cls == "drug_response":
        return f"ClinVar links this variant to drug response: {conditions}. You carry {copy_txt}."
    if cls == "risk_factor":
        return f"Reported as a risk factor for {conditions}. You carry {copy_txt}."
    if cls == "protective":
        return f"Reported as protective for {conditions}. You carry {copy_txt}."
    return f"Reported association with {conditions}. You carry {copy_txt}."


STAGE = Stage(
    name="annotate", version="1", title="Annotate with knowledge sources and grade claims", fn=build,
    inputs=_inputs,
    outputs=lambda ctx: [annotations_path(ctx), claims_path(ctx), summary_path(ctx)],
    available=_available,
    params=lambda ctx: {"evidence_model": ev.MODEL_VERSION},
)
