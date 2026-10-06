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
from .. import inheritance as inh
from ..interpret import Kb
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
SIG_LABEL = {"predicted": "Predicted damaging", "pathogenic": "Pathogenic",
             "likely_pathogenic": "Likely pathogenic", "drug_response": "Drug response",
             "risk_factor": "Risk factor", "protective": "Protective", "association": "Association",
             "conflicting": "Conflicting"}
BLOOD_GROUP = re.compile(r"blood group", re.I)
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
        "clinvar_subs": s.table("clinvar", "submissions"),
        "validity": s.table("clingen", "gene_validity"),
        "dosage": s.table("clingen", "dosage"),
        "kg": s.table("1000genomes", "af"),
        "gnomad": (gn.parent / "gnomad.hg37.zip") if gn else None,
        "actionability": s.table("clingen", "actionability"),
        "hpo_inh": s.table("hpo", "disease_inheritance"),
        "hpo_gd": s.table("hpo", "gene_disease"),
        "mondo": s.table("mondo", "diseases"),
        "orpha_dis": s.table("orphanet", "disorders"),
        "orpha_gd": s.table("orphanet", "gene_disorder"),
        "acmg": s.table("acmg_sf", "genes"),
        "constraint": s.table("constraint", "genes"),
        "am": s.table("alphamissense", "scores"),
        "revel": s.table("revel", "scores"),
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


def _columns(p: Path) -> set[str]:
    return {r[0] for r in duckdb.connect().execute(f"DESCRIBE SELECT * FROM '{p}'").fetchall()}


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
               "clinvar_conflicts, cv.rsid AS clinvar_rsid, " +
               ("cv.condition_ids AS clinvar_condition_ids" if "condition_ids" in _columns(cv)
                else "NULL::VARCHAR AS clinvar_condition_ids")) if cv else \
        ("NULL::BIGINT clinvar_id, NULL::VARCHAR clinvar_significance, NULL::VARCHAR clinvar_class, "
         "NULL::INT clinvar_stars, NULL::VARCHAR clinvar_review, NULL::VARCHAR clinvar_conditions, "
         "NULL::BOOLEAN clinvar_low_penetrance, NULL::VARCHAR clinvar_conflicts, NULL::VARCHAR clinvar_rsid, "
         "NULL::VARCHAR clinvar_condition_ids")
    kg = k["kg"]
    kg_join = f"LEFT JOIN '{kg}' kg USING (chrom, pos, ref, alt)" if kg else ""
    kg_cols = ("kg.af AS af_1kg, kg.afr AS af_1kg_afr, kg.amr AS af_1kg_amr, kg.eas AS af_1kg_eas, "
               "kg.eur AS af_1kg_eur, kg.sas AS af_1kg_sas") if kg else \
        ", ".join(f"NULL::FLOAT {c}" for c in ("af_1kg", "af_1kg_afr", "af_1kg_amr", "af_1kg_eas", "af_1kg_eur",
                                               "af_1kg_sas"))
    pred_join, pred_cols = "", []
    for key, col in (("am", "am_score"), ("revel", "revel")):
        if k[key]:
            pred_join += (f" LEFT JOIN (SELECT chrom, pos, ref, alt, {col} FROM '{k[key]}') {key}"
                          " USING (chrom, pos, ref, alt)")
            pred_cols.append(f"{key}.{col}")
        else:
            pred_cols.append(f"NULL::FLOAT AS {col}")
    if k["constraint"]:
        pred_join += f" LEFT JOIN '{k['constraint']}' cs ON cs.gene = c.gene"
        pred_cols += ["cs.loeuf", "cs.pli", "cs.mis_z"]
    else:
        pred_cols += ["NULL::FLOAT AS loeuf", "NULL::FLOAT AS pli", "NULL::FLOAT AS mis_z"]
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
               {cv_cols}, {kg_cols}, {gn_cols}, {rs_cols}, {", ".join(pred_cols)}
        FROM '{vp}' v
        LEFT JOIN csq c USING (chrom, pos, ref, alt)
        LEFT JOIN gn USING (chrom, pos, ref, alt)
        {cv_join} {kg_join} {rs_join} {pred_join}
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
    """Graded claims about your alleles: ClinVar assertions, plus computational predictions for rare variants
    nobody has classified. Later steps add curated PGx/trait/PRS claims alongside."""
    out = claims_path(ctx)
    if not k["clinvar"]:
        con.execute(f"COPY (SELECT NULL::VARCHAR claim_id WHERE false) TO '{out}' (FORMAT parquet)")
        return 0
    kb = Kb(k)
    classes = ", ".join(f"'{c}'" for c in CLAIM_CLASSES)
    ann = annotations_path(ctx)
    rows = con.execute(f"""
      SELECT * FROM '{ann}'
      WHERE clinvar_class IN ({classes})
        AND (clinvar_class <> 'conflicting' OR clinvar_conflicts ILIKE '%pathogenic%')
      ORDER BY chrom_order, pos""").fetchdf().to_dict("records")
    lof = " OR ".join(f"consequence LIKE '%{c}%'" for c in ev.LOF)
    predicted_rows = con.execute(f"""
      SELECT * FROM '{ann}'
      WHERE coalesce(clinvar_class, '') NOT IN ({classes}, 'benign', 'likely_benign')
        AND filter = 'PASS'
        AND (am_score > {ev.AM_LIKELY_PATHOGENIC} OR revel >= {ev.REVEL_SUPPORTING} OR {lof})
      ORDER BY chrom_order, pos""").fetchdf().to_dict("records")
    versions = Store(ctx.cfg).versions()
    subs = _submissions(con, k["clinvar_subs"], {r["clinvar_id"] for r in rows if r["clinvar_id"] is not None})
    claims = []
    for r, predicted in [(r, False) for r in rows] + [(r, True) for r in predicted_rows]:
        r = {key: (None if _isnan(val) else val) for key, val in r.items()}
        gene = r["gene"] or (r["genes"][0] if r.get("genes") is not None and len(r["genes"]) else None)
        copies = 2 if r["zygosity"] == "hom" else 1
        popmax = r.get("gnomad_popmax_af")
        if predicted:
            if not kb.disease_gene(gene):
                continue
            conds = kb.gene_conditions(gene)
            modes, basis = kb.modes_for(gene, [])
            e_level, e_reasons = ev.predicted_evidence(
                consequence=r["consequence"], am_score=r.get("am_score"), revel=r.get("revel"), loeuf=r.get("loeuf"),
                popmax_af=popmax, kg_af=r["af_1kg"], recessive_gene="AR" in modes)
            if e_level is None:
                continue
            cls = "predicted"
            lab_subs, also_listed = [], []
        else:
            cls = r["clinvar_class"]
            conds = kb.conditions(r["clinvar_conditions"], r.get("clinvar_condition_ids"))
            lab_subs = subs.get(int(r["clinvar_id"]), []) if r["clinvar_id"] is not None else []
            conds, also_listed = _asserted(conds, lab_subs, cls)
            modes, basis = kb.modes_for(gene, conds)
            e_level, e_reasons = ev.clinvar_evidence(r["clinvar_stars"], cls, popmax_af=popmax, kg_af=r["af_1kg"],
                                                     gene_validity=kb.validity.get(gene),
                                                     low_penetrance=bool(r["clinvar_low_penetrance"]))
        c_level, c_reasons = ev.call_confidence(filter=r["filter"], gq=r["gq"], dp=r["dp"], vaf=r["vaf"],
                                                zygosity=r["zygosity"], provider=_provider_agrees(r))
        role, role_reason = ev.role(modes, r["zygosity"], r["chrom"])
        condition_roles = [{"name": c.name, "role": ev.role(set(c.modes), r["zygosity"], r["chrom"])[0]}
                           for c in conds if c.modes]
        disease_like = cls in ("pathogenic", "likely_pathogenic", "conflicting", "predicted")
        if cls == "drug_response":
            section, group = "pgx", "drug_response"
        elif cls == "association" or (conds and all(BLOOD_GROUP.search(c.name) for c in conds)):
            # ClinVar files blood-group antigens as "pathogenic" phenotypes; they're traits, not disease.
            section, group = "traits", "association"
            disease_like = False
        elif disease_like and role == "carrier":
            section, group = "carrier", "carrier"
        elif cls in ("risk_factor", "protective") or r["clinvar_low_penetrance"]:
            section, group = "health", "risk"
        elif cls == "conflicting":
            section, group = "health", "uncertain"
        elif cls == "predicted":
            section, group = "health", "predicted"
        else:
            section, group = "health", "monogenic"
        if not disease_like:
            # Risk alleles, drug responses and associations aren't "inherited conditions": dosage says it all.
            role, role_reason, condition_roles = None, None, []
        acmg = kb.acmg.get(gene)
        reportable = bool(acmg) and cls in ("pathogenic", "likely_pathogenic") and ev.acmg_reportable(
            acmg["report"], copies=copies, zygosity=r["zygosity"], consequence=r["consequence"],
            aa_change=r["aa_change"])
        act = kb.actionability.get(gene)
        cond_text = "; ".join(c.name for c in conds) or "an unspecified condition"
        major = _alt_is_major(r)
        statement = _statement(cls, cond_text, copies, role, r, gene)
        by_role = {x: [c["name"] for c in condition_roles if c["role"] == x] for x in ("carrier", "possible")}
        if role == "possible" and by_role["carrier"] and by_role["possible"]:
            statement = statement.split(" One copy can matter")[0] + (
                f" Different conditions linked to this variant are inherited differently: for "
                f"{_names(by_role['carrier'])} (recessive) one copy makes you a carrier; for "
                f"{_names(by_role['possible'])} (dominant) one copy can matter, though many people never develop "
                "it.")
        if major:
            statement += _major_note(major, copies)
        if disease_like and basis != "unknown":
            e_reasons = e_reasons + [f"Inheritance: {inh.describe(modes)} (from the "
                                     f"{'named condition' if basis == 'condition' else 'gene'}). {role_reason}"]
        if also_listed:
            e_reasons = e_reasons + [f"Also named on the ClinVar record, but by fewer labs or only in gene-wide "
                                     f"submissions: {_names([c.name for c in also_listed], 4)}"]
        sources = []
        if not predicted:
            sources.append({"source": "clinvar", "version": versions.get("clinvar"), "record": str(r["clinvar_id"]),
                            "url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{r['clinvar_id']}/"})
        if r.get("am_score") is not None:
            sources.append({"source": "alphamissense", "version": versions.get("alphamissense"),
                            "record": f"{r['am_score']:.3f}", "url": ""})
        if r.get("revel") is not None:
            sources.append({"source": "revel", "version": versions.get("revel"), "record": f"{r['revel']:.3f}",
                            "url": ""})
        if popmax is not None:
            sources.append({"source": "gnomad", "version": versions.get("gnomad"), "record": "popmax_af", "url": ""})
        if r["af_1kg"] is not None:
            sources.append({"source": "1000genomes", "version": versions.get("1000genomes"), "record": "AF",
                            "url": ""})
        if gene in kb.validity:
            sources.append({"source": "clingen", "version": versions.get("clingen"), "record": gene,
                            "url": f"https://search.clinicalgenome.org/kb/genes?search={gene}"})
        if any(c.modes_from in ("HPO", "Orphanet") for c in conds) or basis == "gene":
            sources.extend({"source": sid, "version": versions[sid], "record": "inheritance", "url": ""}
                           for sid in ("hpo", "orphanet", "mondo") if versions.get(sid))
        if acmg:
            sources.append({"source": "acmg_sf", "version": versions.get("acmg_sf"), "record": gene, "url": ""})
        kind = "predicted_variant" if predicted else "clinvar_variant"
        tag = "predicted" if predicted else "clinvar"
        claims.append({
            "claim_id": f"variant:{r['chrom']}-{r['pos']}-{r['ref']}-{r['alt']}:{tag}",
            "section": section, "group": group, "kind": kind, "category": cls,
            "category_label": SIG_LABEL[cls],
            "subject": gene or r["rsid"] or f"{r['chrom']}:{r['pos']}", "gene": gene, "rsid": r["rsid"],
            "chrom": r["chrom"], "pos": int(r["pos"]), "ref": r["ref"], "alt": r["alt"],
            "genotype": _genotype(r), "copies": copies, "zygosity": r["zygosity"],
            "consequence": r["consequence"], "impact": r["impact"], "aa_change": r["aa_change"],
            "conditions": cond_text, "conditions_detail": json.dumps([c.as_dict() for c in conds]),
            "conditions_also_listed": "; ".join(c.name for c in also_listed) or None,
            "statement": statement,
            "clinvar_significance": r["clinvar_significance"], "clinvar_stars": int(r["clinvar_stars"] or 0),
            "gene_validity": kb.validity.get(gene),
            "inheritance": "/".join(sorted(modes)) or None, "inheritance_basis": basis, "role": role,
            "role_reason": role_reason, "condition_roles": json.dumps(condition_roles),
            "submitters": json.dumps([] if predicted else lab_subs[:25]),
            "acmg_sf": acmg["category"] if acmg else None, "acmg_rule": acmg["report"] if acmg else None,
            "acmg_reportable": reportable,
            "actionability": json.dumps(act) if act else None,
            "am_score": r.get("am_score"), "revel": r.get("revel"), "loeuf": r.get("loeuf"),
            "af_1kg": r["af_1kg"], "af_1kg_eur": r["af_1kg_eur"], "gnomad_popmax_af": popmax,
            "evidence_level": e_level, "evidence_reasons": e_reasons,
            "call_level": c_level, "call_reasons": c_reasons,
            "overall": ev.overall(e_level, c_level),
            "alt_is_major": major is not None,
            "sensitive": disease_like or gene in SENSITIVE_GENES or bool(acmg),
            "sources": json.dumps(sources),
        })
    _compound_hets(claims, kb)
    import pandas as pd

    df = pd.DataFrame(claims)
    tmp = out.with_suffix(".tmp")
    if df.empty:
        con.execute(f"COPY (SELECT NULL::VARCHAR claim_id WHERE false) TO '{tmp}' (FORMAT parquet)")
    else:
        con.register("claims_df", df)
        con.execute(f"""COPY (SELECT * REPLACE (pos::INTEGER AS pos, copies::TINYINT AS copies,
                                 clinvar_stars::TINYINT AS clinvar_stars, af_1kg::DOUBLE AS af_1kg,
                                 af_1kg_eur::DOUBLE AS af_1kg_eur, gnomad_popmax_af::DOUBLE AS gnomad_popmax_af,
                                 am_score::FLOAT AS am_score, revel::FLOAT AS revel, loeuf::FLOAT AS loeuf)
                        FROM claims_df ORDER BY section, claim_id) TO '{tmp}' (FORMAT parquet)""")
    tmp.replace(out)
    return len(claims)


def _names(xs: list[str], n: int = 3) -> str:
    return ", ".join(xs[:n]) + (f" and {len(xs) - n} more" if len(xs) > n else "")


def _submissions(con, path: Path | None, ids: set[int]) -> dict[int, list[dict]]:
    """Per-lab ClinVar submissions for the variants we're about to claim about (newest first)."""
    if not path or not ids:
        return {}
    con.register("want_ids", __import__("pandas").DataFrame({"variation_id": sorted(ids)}))
    out: dict[int, list[dict]] = {}
    for vid, sig, cls, conds, who, when, contributes in con.execute(f"""
        SELECT variation_id, significance, sig_class, conditions, submitter, evaluated, contributes
        FROM '{path}' JOIN want_ids USING (variation_id)
        ORDER BY variation_id, evaluated DESC NULLS LAST""").fetchall():
        out.setdefault(vid, []).append({"submitter": who, "significance": sig, "class": cls,
                                        "conditions": list(conds or []),
                                        "evaluated": when.isoformat() if when else None,
                                        "counted": bool(contributes)})
    return out


def _asserted(conds: list, lab_subs: list[dict], cls: str) -> tuple[list, list]:
    """Split conditions into (supported, also_listed) by how many labs called the variant (likely) pathogenic for them.

    ClinVar's merged record lists every condition any submitter mentioned. Some labs submit one record naming every
    disease of the gene (a GJB2 hearing-loss variant then appears linked to dominant skin syndromes), so submissions
    naming more than two conditions don't count as specific support. A condition is kept when at least 20% as many
    labs name it specifically as name the best-supported one."""
    if cls not in ("pathogenic", "likely_pathogenic", "conflicting") or not lab_subs:
        return conds, []
    support: dict[str, int] = {}
    for sub in lab_subs:
        names = [n.casefold() for n in sub["conditions"] if n.casefold() not in FILLER]
        if sub["class"] in ("pathogenic", "likely_pathogenic") and 0 < len(names) <= 2:
            for n in names:
                support[n] = support.get(n, 0) + 1
    if not support:
        return conds, []
    floor = max(support.values()) * 0.2
    keep = [c for c in conds if support.get(c.name.casefold(), 0) >= floor]
    return (keep, [c for c in conds if c not in keep]) if keep else (conds, [])


def _compound_hets(claims: list[dict], kb: Kb) -> None:
    """Two different disease variants in the same recessive gene may hit both copies (if on different
    chromosomes). Short reads can't tell which chromosome each is on, so flag them rather than call 'carrier'."""
    by_gene: dict[str, list[dict]] = {}
    for c in claims:
        if c["section"] == "carrier" and c["gene"]:
            by_gene.setdefault(c["gene"], []).append(c)
    for gene, cs in by_gene.items():
        if len(cs) < 2:
            continue
        if not any(c["category"] in ("pathogenic", "likely_pathogenic") for c in cs):
            # e.g. GALT's Duarte allele: two 'conflicting' variants that normally travel together on one chromosome
            for c in cs:
                others = ", ".join(o["aa_change"] or f"{o['chrom']}:{o['pos']}" for o in cs if o is not c)
                c["statement"] += (f" You also carry another listed variant in {gene} ({others}). Labs disagree "
                                   "about both, and nearby variants like these are often inherited together on the "
                                   "same chromosome, so this is most likely still carrier status.")
            continue
        for c in cs:
            c["section"], c["group"], c["role"] = "health", "compound", "possible"
            c["statement"] += (f" You also carry another flagged variant in {gene}. If the two are on different copies "
                               "of the chromosome (\"compound heterozygous\"), both copies of the gene are affected; "
                               "only family testing or long-read sequencing can tell.")
            acmg = kb.acmg.get(gene)
            if acmg and c["category"] in ("pathogenic", "likely_pathogenic"):
                c["acmg_reportable"] = ev.acmg_reportable(acmg["report"], copies=1, zygosity="het",
                                                          consequence=c["consequence"], aa_change=c["aa_change"],
                                                          compound=True)


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


def _statement(cls: str, conditions: str, copies: int, role: str | None, r: dict, gene: str | None) -> str:
    copy_txt = "two copies (both chromosomes)" if copies == 2 else "one copy"
    if r["zygosity"] == "hemi":
        copy_txt = "one copy (on your single X chromosome)"
    tail = {
        "carrier": " This condition is recessive, so one copy usually means you are a healthy carrier, not affected.",
        "affected": " With this inheritance pattern, your genotype is the one associated with the condition.",
        "possible": " One copy can matter for dominant conditions, though many people who carry such variants never "
                    "develop the condition (incomplete penetrance).",
        "unknown": " Whether one copy matters depends on how the condition is inherited and its penetrance.",
    }.get(role or "", "") if copies == 1 or role == "affected" else ""
    if cls in ("pathogenic", "likely_pathogenic"):
        label = "pathogenic" if cls == "pathogenic" else "likely pathogenic"
        return f"ClinVar lists this variant as {label} for {conditions}. You carry {copy_txt}.{tail}"
    if cls == "predicted":
        what = ("predicted to stop the gene working" if ev.is_lof(r["consequence"]) else
                "predicted to damage the protein")
        return (f"This rare variant in {gene} is {what} by computer models; no lab has classified it. {gene} is "
                f"linked to {conditions}. You carry {copy_txt}.{tail}")
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
    name="annotate", version="4", title="Annotate with knowledge sources and grade claims", fn=build,
    inputs=_inputs,
    outputs=lambda ctx: [annotations_path(ctx), claims_path(ctx), summary_path(ctx)],
    available=_available,
    params=lambda ctx: {"evidence_model": ev.MODEL_VERSION},
)
