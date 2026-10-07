"""Traits & Brain (ADR-017): single well-studied variants, and rare changes in brain-development genes.

* **Single-variant traits** — the GWAS Catalog associations listed in ``knowledge/data/gwas_featured.csv``
  (eye colour, lactase, ALDH2 …). Each is genotyped with the same rule as PharmCAT positions (callable map,
  GQ ≥ 20, never assume reference where coverage was poor), then the effect allele from the catalog is counted
  (on either strand) and set against how common your genotype is in the 1000 Genomes populations.
* **Brain-development genes** — your rare protein-changing variants in SFARI Gene autism-research genes,
  flagged only by the same calibrated predictor rule used everywhere else (``evidence.predicted_evidence``),
  so a long list of harmless rare variants does not masquerade as a finding.
"""

from __future__ import annotations

from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pysam

from .. import evidence as ev
from ..knowledge import Store
from ..pipeline import Context, Stage, write_json
from . import annotate, coverage, reference, variants
from .pgx import _genotype_all

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
POPS = ("afr", "amr", "eas", "eur", "sas")


def traits_dir(ctx: Context) -> Path:
    p = ctx.work / "traits"
    p.mkdir(parents=True, exist_ok=True)
    return p


def out_snps(ctx: Context) -> Path:
    return traits_dir(ctx) / "snps.parquet"


def out_brain(ctx: Context) -> Path:
    return traits_dir(ctx) / "brain_variants.parquet"


def out_summary(ctx: Context) -> Path:
    return traits_dir(ctx) / "traits.json"


def effect_copies(alleles: list[str] | None, effect: str | None, site: set[str]) -> tuple[int | None, str | None]:
    """How many copies of the catalog's effect allele you carry, and on which strand it was read.

    `site` = the alleles seen at this position (reference + 1000 Genomes alternates). The catalog reports
    alleles on the forward strand of GRCh38 but older studies sometimes used the other strand, so if the
    effect allele is not one of the site's alleles, its complement is tried. Palindromic sites (A/T, C/G)
    cannot be resolved that way and are read as forward."""
    if alleles is None or effect not in COMP:
        return None, None
    if effect in site:
        return sum(a == effect for a in alleles), "forward"
    if COMP[effect] in site:
        return sum(a == COMP[effect] for a in alleles), "reverse"
    return None, None


CALL_REASON = {
    "variant": "A passing variant call with good genotype quality at this exact position",
    "reference": "No variant called here and the position is well covered, so you match the reference",
}


def gwas_effect(or_beta: float | None, ci_text: str | None, effect: str | None, trait: str | None) -> dict:
    """Read the GWAS Catalog's effect convention: a 'unit increase/decrease' note marks a beta (per copy of the
    effect allele, sign in the note); otherwise the number is an odds ratio."""
    if or_beta is None or or_beta != or_beta or effect is None:
        return {"effect_kind": None, "effect_direction": None, "effect_text": None}
    note = (ci_text or "").lower()
    if "increase" in note or "decrease" in note:
        up = "increase" in note
        return {"effect_kind": "beta", "effect_direction": 1 if up else -1,
                "effect_text": f"Each copy of {effect} moves the study’s measure of “{trait}” "
                               f"{'up' if up else 'down'} (β = {or_beta:.3g} per copy)"}
    if or_beta == 1:
        return {"effect_kind": "or", "effect_direction": 0,
                "effect_text": f"Each copy of {effect} leaves the odds of “{trait}” unchanged (odds ratio 1)"}
    up = or_beta > 1
    return {"effect_kind": "or", "effect_direction": 1 if up else -1,
            "effect_text": f"Each copy of {effect} {'raises' if up else 'lowers'} the odds of “{trait}” "
                           f"(odds ratio {or_beta:.3g} per copy)"}


def genotype_share(copies: int | None, p: float | None) -> float | None:
    """Share of people with the same number of effect-allele copies (Hardy–Weinberg from allele frequency p)."""
    if copies is None or p is None:
        return None
    return {0: (1 - p) ** 2, 1: 2 * p * (1 - p), 2: p * p}[copies]


def _snps(ctx: Context) -> list[dict]:
    st = Store(ctx.cfg)
    feat = st.table("gwas_catalog", "featured")
    af = st.table("1000genomes", "af")
    con = duckdb.connect()
    rows = con.execute(f"SELECT * FROM '{feat}' WHERE chrom37 IS NOT NULL").fetchdf().to_dict("records")
    # Positions in the shape the PharmCAT genotyper expects (any base may be the alternate).
    tmp = traits_dir(ctx) / "_positions.parquet"
    uniq = {(r["rsid"], r["chrom37"], int(r["pos37"])) for r in rows}
    fa_ref = {}
    fa = pysam.FastaFile(str(reference.reference_fasta(ctx)))
    prefixed = "chr1" in fa.references
    for rsid, c, p in uniq:
        name = f"chr{c}" if prefixed else c
        fa_ref[rsid] = fa.fetch(name, p - 1, p).upper()
    pq.write_table(pa.Table.from_pylist([
        {"rsid": rsid, "gene": None, "chrom38": c, "pos38": p, "chrom37": f"chr{c}", "pos37": p, "ref": fa_ref[rsid],
         "alts": [b for b in "ACGT" if b != fa_ref[rsid]]} for rsid, c, p in sorted(uniq)]), tmp)
    calls = {c["rsid"]: c for c in _genotype_all(ctx, tmp)}
    tmp.unlink()
    out = []
    for r in rows:
        c = calls[r["rsid"]]
        ref = fa_ref[r["rsid"]]
        alleles = None
        if c["gt"] != "./.":
            idx = [int(i) for i in c["gt"].split("/")]
            bases = [ref] + c["alts"]
            alleles = [bases[i] for i in idx]
        kg = con.execute(f"SELECT alt, af, {','.join(POPS)} FROM '{af}' WHERE chrom = ? AND pos = ? AND ref = ? "
                         "AND length(alt) = 1", [r["chrom37"], int(r["pos37"]), ref]).fetchall()
        site = {ref} | {k[0] for k in kg}
        effect = r["effect_allele"] if r["effect_allele"] in COMP else None
        copies, strand = effect_copies(alleles, effect, site)
        e_fwd = effect if strand != "reverse" else (COMP[effect] if effect else None)
        freq = {}
        if e_fwd and kg:
            for k in kg:
                vals = dict(zip(("af",) + POPS, k[1:], strict=True))
                for pop, v in vals.items():
                    if v is None:
                        continue
                    if e_fwd == k[0]:
                        freq[pop] = float(v)
                    elif e_fwd == ref:
                        freq[pop] = freq.get(pop, 1.0) - float(v)
        out.append({
            "rsid": r["rsid"], "label": r["label"], "mapped_trait": r["mapped_trait"], "section": r["section"],
            "body_system": r["body_system"], "gene": r["gene"], "chrom": r["chrom37"], "pos": int(r["pos37"]),
            "ref": ref, "genotype": "/".join(alleles) if alleles else None, "status": c["status"],
            "note": str(c["note"]) if c.get("note") else None,
            "effect_allele": effect, "effect_allele_strand": strand, "effect_copies": copies,
            "effect_af": freq.get("af"), "effect_af_eur": freq.get("eur"),
            "effect_af_pops": {k: round(v, 4) for k, v in freq.items()} or None,
            "genotype_share": genotype_share(copies, freq.get("af")),
            "genotype_share_eur": genotype_share(copies, freq.get("eur")),
            "or_beta": r["or_beta"], "ci_text": r["ci_text"], "mlog10p": r["mlog10p"],
            "catalog_effect_freq": r["effect_allele_freq"], "study": r["best_study"], "pmid": r["best_pmid"],
            "first_author": r["best_first_author"], "published": r["best_published"],
            "trait_reported": r["best_trait_reported"], "sample": r["best_sample"],
            "publications": r["publications"], "studies": r["studies"],
            "direction_known": effect is not None,
            **gwas_effect(r["or_beta"], r["ci_text"], effect, r["best_trait_reported"]),
        })
        level, reasons = ev.gwas_evidence(publications=r["publications"], mlog10p=r["mlog10p"],
                                          direction_known=effect is not None)
        call = "High" if alleles else "Not callable"
        out[-1].update({"evidence_level": level, "evidence_reasons": reasons, "call_confidence": call,
                        "call_reasons": [CALL_REASON.get(c["status"], c.get("note") or c["status"])]
                        if alleles else [c.get("note") or c["status"]],
                        "overall": ev.overall(level, call)})
    pq.write_table(pa.Table.from_pylist(out), out_snps(ctx))
    return out


def _brain(ctx: Context) -> list[dict]:
    st = Store(ctx.cfg)
    genes = st.table("sfari", "genes")
    ann = annotate.annotations_path(ctx)
    con = duckdb.connect()
    rows = con.execute(f"""SELECT a.chrom, a.pos, a.rsid, a.ref, a.alt, a.gene, a.zygosity, a.gt, a.gq,
            a.consequence, a.impact, a.aa_change, a.am_score, a.revel, a.loeuf, a.gnomad_popmax_af, a.af_1kg,
            a.clinvar_significance, a.clinvar_stars, s.score AS sfari_score, s.syndromic, s.category, s.reports
        FROM '{ann}' a JOIN '{genes}' s ON s.gene = a.gene
        WHERE a.filter = 'PASS' AND coalesce(a.gq, 0) >= 20 AND a.impact IN ('HIGH', 'MODERATE')
        ORDER BY a.chrom_order, a.pos""").fetchdf().to_dict("records")
    out = []
    for r in rows:
        r = {k: (None if isinstance(v, float) and v != v else v) for k, v in r.items()}
        level, reasons = ev.predicted_evidence(
            consequence=r["consequence"], am_score=r["am_score"], revel=r["revel"], loeuf=r["loeuf"],
            popmax_af=r["gnomad_popmax_af"], kg_af=r["af_1kg"], recessive_gene=False)
        afs = [a for a in (r["gnomad_popmax_af"], r["af_1kg"]) if a is not None]
        rare = not afs or max(afs) < ev.PREDICTED_MAX_AF
        if not rare:
            continue
        out.append({**r, "flagged": level is not None, "evidence_level": level, "evidence_reasons": reasons})
    pq.write_table(pa.Table.from_pylist(out) if out else pa.table({"gene": pa.array([], pa.string())}),
                   out_brain(ctx))
    return out


def _run(ctx: Context) -> dict:
    snps = _snps(ctx)
    brain = _brain(ctx) if Store(ctx.cfg).current("sfari") else []
    st = Store(ctx.cfg)
    used = {sid: (st.current(sid) or {}).get("version") for sid in ("gwas_catalog", "1000genomes", "sfari")}
    res = {"knowledge": {k: v for k, v in used.items() if v}, "evidence_model": ev.MODEL_VERSION,
           "snps": len(snps), "snps_called": sum(s["genotype"] is not None for s in snps),
           "brain_rare": len(brain), "brain_flagged": sum(b["flagged"] for b in brain)}
    write_json(out_summary(ctx), res)
    return res


def _available(ctx: Context) -> str | None:
    st = Store(ctx.cfg)
    if not st.current("gwas_catalog"):
        return "needs the gwas_catalog knowledge source"
    if not st.current("1000genomes"):
        return "needs the 1000genomes knowledge source"
    if not reference.reference_fasta(ctx):
        return "needs a verified reference"
    if not variants.variants_parquet(ctx).exists() or not (coverage.out_dir(ctx) / "callable.parquet").exists():
        return "needs variants and the callable map"
    return None


def _inputs(ctx: Context) -> list[Path]:
    st = Store(ctx.cfg)
    return [variants.variants_parquet(ctx), coverage.out_dir(ctx) / "callable.parquet",
            annotate.annotations_path(ctx),
            st.table("gwas_catalog", "featured"), st.table("sfari", "genes")]


STAGE = Stage(
    name="traits", version="3", title="Traits (GWAS single variants) and brain-development genes", fn=_run,
    inputs=_inputs, params=lambda ctx: {"evidence_model": ev.MODEL_VERSION},
    outputs=lambda ctx: [out_snps(ctx), out_brain(ctx), out_summary(ctx)], available=_available,
)
