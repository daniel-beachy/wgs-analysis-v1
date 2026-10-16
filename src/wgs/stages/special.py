"""Specialised callers (ADR-021): what a standard VCF can't show, read again from the aligned reads.

Short-read variant callers handle most of the genome, but three kinds of place need their own tools:

1. **Genes with a near-identical twin.** SMN1 (spinal muscular atrophy) sits next to SMN2, GBA1 (Gaucher disease,
   Parkinson's risk) next to its pseudogene GBAP1. Reads can't be placed on one copy, so callers skip them. Illumina's
   SMNCopyNumberCaller and Gauchian count copies and look for the disease changes using the differences between the
   twins (Chen et al. 2020; Toffoli et al. 2022).
2. **Repeat expansions.** ExpansionHunter (Dolzhenko et al. 2019) sizes the repeat at every STRchive disease locus.
   Counts are compared with STRchive's normal / intermediate / disease ranges *relative to the reference genome*:
   the catalogue's regions sometimes include interruptions that clinical counts leave out (AR: 22 pure CAG plus an
   interruption and 5 more), so a count is converted to the clinical scale as the reference's longest pure run plus
   the difference from the reference. The call grade comes from ExpansionHunter's confidence interval, so a sparse
   GC-rich locus reads "most likely typical, low confidence" rather than a false all-clear or a false alarm.

Findings outside the typical range become ordinary graded claims (Health / Carrier); every checked gene and locus is
also listed, with its grade, so "nothing found" is visible and auditable.
"""

from __future__ import annotations

import json
import re
import shutil
import sys
import tomllib
from functools import cache
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from .. import evidence as ev
from ..knowledge import Store
from ..pipeline import Context, Stage, run, write_json
from . import pgx, qc, reference, variants

PMID = "https://pubmed.ncbi.nlm.nih.gov/{}/"
CITE = {
    "smn": ("SMNCopyNumberCaller", "Chen X, et al. Genet Med. 2020;22:945-953", "32066871"),
    "gba": ("Gauchian", "Toffoli M, et al. Commun Biol. 2022;5:670", "35794204"),
    "eh": ("ExpansionHunter", "Dolzhenko E, et al. Bioinformatics. 2019;35:4754-4756", "31134279"),
    "eh_clinical": ("Repeat expansions from genomes", "Ibañez K, et al. Lancet Neurol. 2022;21:234-245", "35182509"),
    "acmg_carrier": ("ACMG carrier screening", "Gregg AR, et al. Genet Med. 2021;23:1793-1806", "34285390"),
}
STRCHIVE_EVIDENCE = {"Definitive": "Strong", "Strong": "Strong", "Moderate": "Moderate", "Limited": "Limited",
                     "Provisional": "Limited"}
ADULT_ONSET = 18  # STRchive age_onset_min at or above this → adult-onset → click-to-reveal (AGENTS.md §5.9)
HLA_TABLE = Path(__file__).parents[1] / "knowledge" / "data" / "hla_associations.toml"
HLA_MIN_QUALITY = 20  # T1K allele quality; lower calls are graded Medium
IUPAC = {"N": "[ACGT]", "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]",
         "D": "[AGT]", "H": "[ACT]", "V": "[ACG]"}


def out_dir(ctx: Context) -> Path:
    p = ctx.work / "special"
    p.mkdir(parents=True, exist_ok=True)
    return p


def out_genes(ctx: Context) -> Path:
    return out_dir(ctx) / "genes.parquet"


def out_repeats(ctx: Context) -> Path:
    return out_dir(ctx) / "repeats.parquet"


def out_claims(ctx: Context) -> Path:
    return out_dir(ctx) / "claims_special.parquet"


def out_hla(ctx: Context) -> Path:
    return out_dir(ctx) / "hla_associations.parquet"


def out_summary(ctx: Context) -> Path:
    return out_dir(ctx) / "special.json"


def _env() -> Path:
    return Path(sys.prefix).parent / "special"


def _tool(name: str) -> str | None:
    p = _env() / "bin" / name
    return str(p) if p.exists() else None


def _source(key: str) -> dict:
    tool, cite, pmid = CITE[key]
    return {"source": tool, "version": None, "record": cite, "url": PMID.format(pmid)}


def _manifest(ctx: Context, d: Path) -> Path:
    aln = ctx.inv.first("cram") or ctx.inv.first("bam")
    m = d / "manifest.txt"
    m.write_text(str(aln.path.resolve()) + "\n")
    return m


def _build(ctx: Context) -> str:
    return "19" if "chr1" in reference.cram_contigs(ctx) else "37"


def _uneven(log: Path) -> bool:
    return log.exists() and "uneven coverage" in log.read_text().lower()


# ---------------------------------------------------------------------------------------------------------------
# SMN1 / SMN2
# ---------------------------------------------------------------------------------------------------------------


def smn_row(res: dict, uneven: bool) -> dict:
    """Interpret SMNCopyNumberCaller's JSON for one sample (rules in ADR-021)."""
    smn1, smn2, silent = res.get("SMN1"), res.get("SMN2"), res.get("g27134TG_CN") or 0
    raw = res.get("SMN1_CN_raw") or []
    mean_raw = sum(raw) / len(raw) if raw else None
    if smn1 is None:
        status, text = "not_determined", (
            "Your SMN1 copy number could not be determined: the reads sit between 1 and 2 copies "
            f"(average estimate {mean_raw:.1f})" if mean_raw is not None else "Your SMN1 copy number could not be "
            "determined" ) + ", so carrier status for spinal muscular atrophy is unknown from this data"
    elif smn1 == 0:
        status, text = "affected", "No working SMN1 copies were found, the genotype of spinal muscular atrophy"
    elif smn1 == 1:
        status, text = "carrier", ("One SMN1 copy: a carrier of spinal muscular atrophy (usually no effect on you; "
                                   "relevant for children if a partner is also a carrier)")
    elif smn1 == 2 and silent:
        status, text = "possible_silent_carrier", (
            "Two SMN1 copies, but you carry g.27134T>G, a variant that often marks two copies on the same chromosome "
            "(a 'silent' carrier)")
    else:
        status, text = "typical", f"{smn1} SMN1 copies: not a carrier of spinal muscular atrophy"
    call = "Not callable" if smn1 is None else ("Medium" if uneven else "High")
    reasons = ([f"SMN1 copy-number estimates across 16 sites: {', '.join(f'{x:.2f}' for x in raw)}"] if raw else [])
    reasons += [f"SMNCopyNumberCaller result: {res.get('Info') or 'called'}; median depth {res.get('Median_depth')}"]
    if uneven:
        reasons.append("The tool warns that coverage is uneven, which makes copy-number calls less certain")
    return {"gene": "SMN1", "condition": "Spinal muscular atrophy", "inheritance": "AR", "status": status,
            "statement": text, "copies": smn1, "detail": json.dumps({"SMN1": smn1, "SMN2": smn2,
                                                                     "SMN2_delta7_8": res.get("SMN2delta78"),
                                                                     "g27134TG": silent, "SMN1_raw": raw}),
            "evidence_level": "Strong",
            "evidence_reasons": ["SMN1 loss is the established cause of spinal muscular atrophy; ACMG recommends "
                                 "offering SMA carrier screening to everyone planning a pregnancy"],
            "call_level": call, "call_reasons": reasons,
            "next_step": ("A clinical SMA carrier test (SMN1 copy number by qPCR or MLPA) settles this"
                          if status == "not_determined" else None),
            "url": "https://www.ncbi.nlm.nih.gov/books/NBK1352/",
            "sources": json.dumps([_source("smn"), _source("acmg_carrier"),
                                   {"source": "GeneReviews", "version": None, "record": "Spinal Muscular Atrophy",
                                    "url": "https://www.ncbi.nlm.nih.gov/books/NBK1352/"}])}


def _smn(ctx: Context, d: Path) -> dict | None:
    tool = _tool("smn_caller.py")
    if not tool:
        return None
    work = d / "smn"
    shutil.rmtree(work, ignore_errors=True)
    log = ctx.logs / "special-smn.log"
    log.unlink(missing_ok=True)
    run([str(_env() / "bin" / "python"), tool, "--manifest", str(_manifest(ctx, d)), "--genome", _build(ctx),
         "--outDir", str(work), "--prefix", "smn", "--threads", str(min(ctx.threads, 6)),
         "--reference", str(reference.reference_fasta(ctx))], log=log)
    (_, res), = json.loads((work / "smn.json").read_text()).items()
    return smn_row(res, _uneven(log))


# ---------------------------------------------------------------------------------------------------------------
# GBA1
# ---------------------------------------------------------------------------------------------------------------


def gba_row(res: dict, uneven: bool) -> dict:
    counts = res.get("variant_raw_count") or {}
    found = [v for v in str(res.get("GBAP1_like_variant_exon9_11") or "").split(",") if v]
    other = [v for v in str(res.get("other_variants") or "").split(",") if v]
    depth = [sum(int(x) for x in str(c).split(",")) for c in counts.values()]
    deletion = res.get("deletion_CN")
    if deletion is not None and deletion < 2:
        status, text = "carrier", "One copy of GBA1 carries a deletion from recombination with its pseudogene"
    elif found or other:
        status = "carrier"
        text = f"GBA1 change{'s' if len(found + other) > 1 else ''} found: {', '.join(found + other)}"
    else:
        status, text = "typical", (f"None of the {len(counts)} disease-linked GBA1 changes Gauchian tests were found, "
                                   "and there is no deletion or pseudogene swap")
    call = "Not callable" if deletion is None else ("Medium" if uneven or (depth and min(depth) < 15) else "High")
    reasons = [f"Each tested position was read {min(depth)}–{max(depth)} times" if depth else "No read counts",
               f"GBA1 + GBAP1 copy number {res.get('deletion_CN_raw')} (2 expected per gene)"]
    if uneven:
        reasons.append("The tool warns that coverage is uneven, which makes copy-number calls less certain")
    return {"gene": "GBA1", "condition": "Gaucher disease; GBA1 variants also raise Parkinson's disease risk",
            "inheritance": "AR", "status": status, "statement": text, "copies": deletion,
            "detail": json.dumps({"tested": sorted(counts), "found": found + other, "deletion_CN": deletion}),
            "evidence_level": "Strong",
            "evidence_reasons": ["Biallelic GBA1 variants cause Gaucher disease; carrying one raises Parkinson's "
                                 "disease risk (GeneReviews)"],
            "call_level": call, "call_reasons": reasons, "next_step": None,
            "url": "https://www.ncbi.nlm.nih.gov/books/NBK1269/",
            "sources": json.dumps([_source("gba"),
                                   {"source": "GeneReviews", "version": None, "record": "Gaucher Disease",
                                    "url": "https://www.ncbi.nlm.nih.gov/books/NBK1269/"}])}


def _gba(ctx: Context, d: Path) -> dict | None:
    site = next((_env() / "lib").glob("python3*/site-packages/gauchian"), None)
    if not site:
        return None
    work = d / "gba"
    shutil.rmtree(work, ignore_errors=True)
    log = ctx.logs / "special-gba.log"
    log.unlink(missing_ok=True)
    # Gauchian ships its depth_calling helpers inside the package but imports them as a top-level module.
    code = ("import sys; sys.path.append(sys.argv[1]); sys.argv = ['gauchian'] + sys.argv[2:]; "
            "from gauchian.gauchian import run; run()")
    run([str(_env() / "bin" / "python"), "-c", code, str(site), "-m", str(_manifest(ctx, d)), "-g", _build(ctx),
         "-o", str(work), "-p", "gba", "-t", str(min(ctx.threads, 6)), "--reference",
         str(reference.reference_fasta(ctx))], log=log)
    (_, res), = json.loads((work / "gba.json").read_text()).items()
    return gba_row(res, _uneven(log))


# ---------------------------------------------------------------------------------------------------------------
# Repeat expansions
# ---------------------------------------------------------------------------------------------------------------


def pure_run(seq: str, motif: str) -> int:
    """Longest uninterrupted run of `motif` (IUPAC codes allowed) in `seq`, in any frame."""
    pat = re.compile("".join(IUPAC.get(c, c) for c in motif.upper()))
    n, best = len(motif), 0
    for f in range(n):
        run_ = 0
        for i in range(f, len(seq) - n + 1, n):
            run_ = run_ + 1 if pat.fullmatch(seq[i:i + n].upper()) else 0
            best = max(best, run_)
    return best


def to_clinical(count: int, ref_count: int, pure: int, benign_max: int | None = None) -> int:
    """A catalogue count on the clinical scale: the reference genome's clinical count plus the difference from it.

    The reference's clinical count is its longest pure run of the motif, capped at the top of the typical range: the
    reference genome is a typical allele at these loci, and a loose motif such as GCN can absorb a flanking codon.
    Insertion-type loci (the disease motif is absent from the reference, pure == 0) are counted directly."""
    if pure == 0:
        return count
    ref_clinical = min(pure, benign_max) if benign_max is not None else pure
    return ref_clinical + (count - ref_count)


def classify_repeat(alleles: list[int], upper: list[int], ref_count: int, pure: int, t: dict,
                    vcf_same: bool = False) -> tuple[str, str]:
    """(status, call level) for one locus. Counts move to the clinical scale relative to the reference genome.

    `vcf_same`: the VCF has no length-changing call inside the repeat (confident or unsure), i.e. a second,
    independent caller says your repeat is the reference length. With ExpansionHunter within one unit of the
    reference too, that is a typical allele with a High call (ExpansionHunter can differ by a unit at regions that
    are not a whole number of repeats)."""
    bmax, imin, pmin = t.get("benign_max"), t.get("intermediate_min"), t.get("pathogenic_min")
    if not alleles or all(a == 0 for a in alleles) and upper and max(upper) > 100:
        return "not_readable", "Not callable"
    if vcf_same and pure > 0 and all(abs(a - ref_count) <= 1 for a in alleles):
        return "typical", "High"
    if pmin is None:
        return "no_thresholds", "Low"
    if bmax is not None and pmin <= bmax:
        return "no_thresholds", "Low"  # inconsistent ranges in the source: do not classify
    conv = [to_clinical(a, ref_count, pure, bmax) for a in alleles]
    conv_hi = [to_clinical(h, ref_count, pure, bmax) for h in upper]
    reference_like = pure > 0 and all(abs(a - ref_count) <= 1 for a in alleles)
    if not reference_like and any(c >= pmin for c in conv):
        status = "disease_range"
    elif not reference_like and imin is not None and any(c >= imin for c in conv):
        status = "intermediate"
    else:
        status = "typical"
    if status != "typical":
        call = "Medium"  # sizing of long expansions from short reads is approximate (Ibañez 2022)
    elif any(h >= pmin for h in conv_hi):
        call = "Low"
    elif imin is not None and any(h >= imin for h in conv_hi):
        call = "Medium"
    else:
        call = "High"
    return status, call


def _region_seq(fasta: Path, region: str) -> str:
    import pysam

    chrom, rng = region.split(":")
    start, end = (int(x) for x in rng.split("-"))
    with pysam.FastaFile(str(fasta)) as fa:
        name = chrom if chrom in fa.references else chrom.removeprefix("chr")
        return fa.fetch(name, start, end).upper()


def repeat_rows(eh: dict, loci: dict[str, dict], fasta: Path, indels=lambda region: None) -> list[dict]:
    """`indels(region)` returns the VCF's length-changing PASS calls inside a region (list), or None if unknown."""
    rows = []
    for lid, res in eh["LocusResults"].items():
        t = loci.get(lid)
        if not t:
            continue
        motifs = list(t.get("pathogenic_motif_reference_orientation") or [])
        variants = list((res.get("Variants") or {}).values())
        # Compound loci list several sub-repeats; only the one with the disease motif is classified.
        v = next((x for x in variants if not motifs or x.get("RepeatUnit") in motifs), None)
        if v is None:
            continue
        gt, ci = v.get("Genotype") or "", v.get("GenotypeConfidenceInterval") or ""
        alleles = [int(x) for x in gt.split("/") if x.isdigit()]
        upper = [int(x.split("-")[1]) for x in ci.split("/") if "-" in x]
        seq = _region_seq(fasta, v["ReferenceRegion"])
        motif = v["RepeatUnit"]
        ref_count, pure = round(len(seq) / len(motif)), pure_run(seq, motif)
        found = indels(v["ReferenceRegion"])
        vcf_same = found is not None and not found
        status, call = classify_repeat(alleles, upper, ref_count, pure, t, vcf_same)
        ev_level = STRCHIVE_EVIDENCE.get((list(t.get("evidence") or []) or [""])[0], "Limited")
        disputed = any(e in ("Disputed", "Refuted") for e in (t.get("evidence") or []))
        rows.append({
            "locus_id": lid, "gene": t["gene"], "disease": t["disease"],
            "inheritance": ",".join(t.get("inheritance") or []), "motif": motif,
            "alleles": gt, "confidence_interval": ci,
            # Shown only where the catalogue region is a clean repeat (or an insertion-type locus): elsewhere the
            # conversion is good enough to classify against wide ranges but not to quote as a number.
            "clinical_estimate": None if status == "not_readable" or not (pure == 0 or abs(pure - ref_count) <= 1)
            else "/".join(str(to_clinical(ref_count if vcf_same and pure > 0 and abs(a - ref_count) <= 1 else a,
                                          ref_count, pure, t.get("benign_max"))) for a in alleles),
            # Only meaningful where the repeat was read: an unreadable locus has no VCF calls either.
            "vcf_agrees": vcf_same and status not in ("not_readable",), "vcf_indels": json.dumps(found or []),
            "reference_count": ref_count, "reference_pure_run": pure,
            "benign_max": t.get("benign_max"), "intermediate_min": t.get("intermediate_min"),
            "pathogenic_min": t.get("pathogenic_min"),
            "status": "disputed_locus" if disputed else status,
            "evidence_level": "Limited" if disputed else ev_level,
            "strchive_evidence": ",".join(t.get("evidence") or []),
            "call_level": call,
            "overall": ev.overall("Limited" if disputed else ev_level, call),
            "coverage": res.get("Coverage"), "sensitive": (t.get("age_onset_min") or 0) >= ADULT_ONSET,
            "age_onset": t.get("age_onset"), "url": t.get("url"), "pmids": list(t.get("pmids") or []),
        })
    return rows


def _repeats(ctx: Context, d: Path, st: Store) -> list[dict] | None:
    tool = _tool("ExpansionHunter")
    catalog, loci_t = st.file("strchive", "catalog_hg19.json"), st.table("strchive", "loci")
    if not tool or not catalog or not loci_t:
        return None
    sex = (_sex(ctx) or "").upper()
    prefix = d / "expansionhunter"
    run([tool, "--reads", str((ctx.inv.first("cram") or ctx.inv.first("bam")).path), "--reference",
         str(reference.reference_fasta(ctx)), "--variant-catalog", str(catalog), "--output-prefix", str(prefix),
         "--sex", "male" if sex.startswith("XY") else "female", "--threads", str(min(ctx.threads, 4))],
        log=ctx.logs / "special-repeats.log")
    (d / "expansionhunter_realigned.bam").unlink(missing_ok=True)
    eh = json.loads(Path(str(prefix) + ".json").read_text())
    loci = {r["id"]: r for r in duckdb.sql(f"SELECT * FROM '{loci_t}'").fetchdf().to_dict("records")}
    for r in loci.values():
        for k, v in list(r.items()):
            if isinstance(v, float) and v != v:  # pandas reads empty fields as NaN
                r[k] = None
        for k in ("evidence", "inheritance", "pathogenic_motif_reference_orientation", "pmids"):
            r[k] = list(r.get(k) if r.get(k) is not None else [])
        for k in ("benign_max", "intermediate_min", "pathogenic_min", "age_onset_min"):
            r[k] = None if r.get(k) is None or r[k] != r[k] else int(r[k])
    vp = variants.variants_parquet(ctx)
    con = duckdb.connect()

    def indels(region: str) -> list[str] | None:
        if not vp.exists():
            return None
        chrom, rng = region.split(":")
        start, end = (int(x) for x in rng.split("-"))
        return [f"{p}:{r}>{a} {g}" for p, r, a, g in con.execute(
            f"""SELECT pos, ref, alt, gt FROM '{vp}' WHERE chrom = ? AND pos BETWEEN ? AND ?
                AND filter IN ('PASS', 'NoCall') AND length(ref) <> length(alt)""",
            [chrom.removeprefix("chr"), start - 5, end + 5]).fetchall()]

    return repeat_rows(eh, loci, reference.reference_fasta(ctx), indels)


# ---------------------------------------------------------------------------------------------------------------
# HLA associations
# ---------------------------------------------------------------------------------------------------------------


@cache
def hla_table() -> list[dict]:
    return tomllib.loads(HLA_TABLE.read_text())["association"] if HLA_TABLE.exists() else []


def carries(typed: dict[str, list[str]], prefix: str) -> list[str]:
    """Your alleles matching an allele prefix at a field boundary: 'DQB1*02' matches DQB1*02:01, not DQB1*020."""
    gene, _, fields = prefix.partition("*")
    want = fields.split(":")
    return [f"{gene}*{a}" for a in typed.get(gene, []) if a.split(":")[:len(want)] == want]


def _quality(typed: dict[str, list[str]], quality: dict[str, list[int]], allele: str) -> int:
    gene, _, name = allele.partition("*")
    alleles, q = typed.get(gene, []), quality.get(gene, [])
    return q[alleles.index(name)] if name in alleles and alleles.index(name) < len(q) else 0


def hla_rows(typed: dict[str, list[str]], quality: dict[str, list[int]]) -> list[dict]:
    """Match your HLA types against the curated table (rules in ADR-021)."""
    out = []
    for a in hla_table():
        found = sorted({x for p in a["alleles"] for x in carries(typed, p)})
        genes = sorted({p.split("*")[0] for p in a["alleles"]})
        kind = a["kind"]
        if kind in ("rule_out", "haplotype") and not all(g in typed for g in genes):
            continue  # a gene the rule needs was not typed: say nothing rather than guess
        if kind == "haplotype":
            hit, direction = all(carries(typed, p) for p in a["alleles"]), a.get("direction", "risk")
        elif kind == "rule_out":
            hit, direction = not found, "protective"
        else:
            hit, direction = bool(found), kind
        cond = a["condition"] if a["condition"][:2].isupper() else a["condition"][0].lower() + a["condition"][1:]
        if kind == "rule_out":
            text = (f"You don't carry {a['markers']}. {a['condition']} is very unlikely without them" if hit
                    else f"You carry {', '.join(found)}, one of the types needed for {cond}; most people "
                         "who carry it never develop it")
        elif hit:
            what = a["markers"] if kind == "haplotype" else ", ".join(found)
            text = f"You carry {what}, linked to {'lower' if direction == 'protective' else 'higher'} risk of {cond}"
        else:
            text = f"You do not carry {a['markers']}"
        basis = found or [f"{g}*{x}" for g in genes for x in typed.get(g, [])]
        low = min([_quality(typed, quality, x) for x in basis] or [0]) < HLA_MIN_QUALITY
        caveats = [c for c in a.get("caveats", []) if all(carries(typed, p) for p in c.get("if_present", []))]
        out.append({
            "association_id": a["id"], "condition": a["condition"], "kind": kind, "direction": direction,
            "markers": a["markers"], "carried": found, "notable": hit, "statement": text,
            "quotes": json.dumps(a.get("quotes", []) + caveats),
            "evidence_level": "Strong", "evidence_reasons": ["Large, replicated association studies (quoted)"],
            "call_level": "Medium" if low else "High",
            "call_reasons": [f"HLA types read by T1K from your reads, allele quality "
                             f"{'below' if low else 'at least'} {HLA_MIN_QUALITY}"],
            "overall": ev.overall("Strong", "Medium" if low else "High"),
        })
    return out


def hla_claims(rows: list[dict]) -> list[dict]:
    out = []
    for r in rows:
        if not r["notable"]:
            continue
        protective = r["direction"] == "protective"
        out.append({
            "claim_id": f"special:hla:{r['association_id']}", "section": "health", "group": "hla", "kind": "hla",
            "category": "protective" if protective else "risk_factor",
            "category_label": "Protective association" if protective else "Risk association",
            "subject": ", ".join(r["carried"]) or r["markers"], "gene": "HLA", "conditions": r["condition"],
            "statement": r["statement"], "evidence_level": r["evidence_level"],
            "evidence_reasons": r["evidence_reasons"], "call_level": r["call_level"], "call_reasons": r["call_reasons"],
            "overall": ev.overall(r["evidence_level"], r["call_level"]), "sensitive": False,
            "sources": json.dumps([{"source": q["source"], "version": None, "record": f"PMID {q['pmid']}",
                                    "url": q["url"], "quote": q["quote"]} for q in json.loads(r["quotes"])]),
        })
    return out


def _hla(ctx: Context) -> list[dict]:
    p = pgx.hla_out(ctx)
    if not p.exists():
        return []
    genes = json.loads(p.read_text()).get("genes") or {}
    typed = {g.removeprefix("HLA-"): [a.lstrip("*") for a in v.get("alleles") or []] for g, v in genes.items()}
    quality = {g.removeprefix("HLA-"): list(v.get("quality") or []) for g, v in genes.items()}
    return hla_rows(typed, quality)


# ---------------------------------------------------------------------------------------------------------------
# Claims
# ---------------------------------------------------------------------------------------------------------------


def gene_claims(rows: list[dict]) -> list[dict]:
    out = []
    for g in rows:
        if g["status"] in ("typical", "not_determined"):
            continue
        affected = g["status"] == "affected"
        out.append({
            "claim_id": f"special:{g['gene']}:{g['status']}", "section": "health" if affected else "carrier",
            "group": "monogenic" if affected else "carrier", "kind": "copy_number", "category": "pathogenic",
            "category_label": "Established cause", "subject": g["gene"], "gene": g["gene"],
            "conditions": g["condition"], "statement": g["statement"], "inheritance": g["inheritance"],
            "inheritance_basis": "GeneReviews", "role": "affected" if affected else "carrier",
            "role_reason": g["statement"], "copies": g["copies"],
            "evidence_level": g["evidence_level"], "evidence_reasons": g["evidence_reasons"],
            "call_level": g["call_level"], "call_reasons": g["call_reasons"],
            "overall": ev.overall(g["evidence_level"], g["call_level"]), "sensitive": False, "sources": g["sources"],
        })
    return out


def repeat_claims(rows: list[dict]) -> list[dict]:
    out = []
    for r in rows:
        if r["status"] not in ("intermediate", "disease_range"):
            continue
        ar = r["inheritance"] == "AR"
        two = r["clinical_estimate"] and sum(int(x) >= (r["pathogenic_min"] or 10**9)
                                             for x in r["clinical_estimate"].split("/")) == 2
        role = "possible" if r["status"] == "intermediate" else ("carrier" if ar and not two else "affected")
        label = "intermediate range" if r["status"] == "intermediate" else "disease-causing range"
        out.append({
            "claim_id": f"special:repeat:{r['locus_id']}", "section": "carrier" if role == "carrier" else "health",
            "group": "carrier" if role == "carrier" else "repeat", "kind": "repeat", "category": "pathogenic",
            "category_label": "Repeat expansion", "subject": r["gene"], "gene": r["gene"], "conditions": r["disease"],
            "statement": (f"Your {r['gene']} repeat is about {r['clinical_estimate']} copies of {r['motif']}, in the "
                          f"{label} for {r['disease']} (STRchive: typical ≤ {r['benign_max']}, disease ≥ "
                          f"{r['pathogenic_min']}). Repeat sizes from short reads are estimates: confirm with a "
                          "clinical repeat test"),
            "inheritance": r["inheritance"], "inheritance_basis": "STRchive", "role": role,
            "role_reason": f"{label} ({r['inheritance']})", "genotype": r["alleles"],
            "evidence_level": r["evidence_level"],
            "evidence_reasons": [f"STRchive gene–disease evidence: {r['strchive_evidence']}"],
            "call_level": r["call_level"],
            "call_reasons": [f"ExpansionHunter estimate {r['alleles']} (95% interval {r['confidence_interval']})"],
            "overall": r["overall"], "sensitive": r["sensitive"],
            "sources": json.dumps([{"source": "STRchive", "version": None, "record": r["locus_id"], "url": r["url"]},
                                   _source("eh"), _source("eh_clinical")]
                                  + [{"source": "PubMed", "version": None, "record": p, "url": PMID.format(p)}
                                     for p in r["pmids"][:3]]),
        })
    return out


# ---------------------------------------------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------------------------------------------


def _sex(ctx: Context) -> str | None:
    p = qc.out_path(ctx)
    return (json.loads(p.read_text()).get("inferred_sex") or {}).get("call") if p.exists() else None


def _write(rows: list[dict], path: Path, key: str) -> None:
    pq.write_table(pa.Table.from_pylist(rows) if rows else pa.table({key: pa.array([], pa.string())}), path)


def _special(ctx: Context) -> dict:
    d = out_dir(ctx)
    st = Store(ctx.cfg)
    summary: dict = {"tools": {k: {"citation": v[1], "pmid": v[2]} for k, v in CITE.items()}, "errors": []}
    genes: list[dict] = []
    for name, fn in (("smn", _smn), ("gba", _gba)):
        try:
            row = fn(ctx, d)
        except Exception as e:  # an optional caller failing must not sink the stage
            summary["errors"].append(f"{name}: {e}")
            row = None
        if row:
            row.update(overall=ev.overall(row["evidence_level"], row["call_level"]), tool=CITE[name][0])
            genes.append(row)
    try:
        repeats = _repeats(ctx, d, st) or []
    except Exception as e:
        summary["errors"].append(f"repeats: {e}")
        repeats = []
    hla = _hla(ctx)
    _write(hla, out_hla(ctx), "association_id")
    summary["repeats"] = {s: sum(r["status"] == s for r in repeats) for s in sorted({r["status"] for r in repeats})}
    summary["knowledge"] = {"strchive": (st.current("strchive") or {}).get("version")}
    _write(genes, out_genes(ctx), "gene")
    _write(repeats, out_repeats(ctx), "locus_id")
    _write(gene_claims(genes) + repeat_claims(repeats) + hla_claims(hla), out_claims(ctx), "claim_id")
    write_json(out_summary(ctx), summary)
    return {"genes": {g["gene"]: g["status"] for g in genes}, "repeats": summary["repeats"],
            "errors": summary["errors"]}


def _available(ctx: Context) -> str | None:
    if not (ctx.inv.first("cram") or ctx.inv.first("bam")):
        return "needs the aligned reads (CRAM/BAM)"
    if not reference.reference_fasta(ctx):
        return "needs the verified reference FASTA"
    if not _tool("ExpansionHunter"):
        return "specialised callers are not installed (pixi install -e special)"
    return None


def _inputs(ctx: Context) -> list[Path]:
    aln = ctx.inv.first("cram") or ctx.inv.first("bam")
    return [p for p in (aln.path if aln else None, qc.out_path(ctx), pgx.hla_out(ctx), HLA_TABLE)
            if p and p.exists()]


def _params(ctx: Context) -> dict:
    st = Store(ctx.cfg)
    return {"strchive": (st.current("strchive") or {}).get("version"), "evidence_model": ev.MODEL_VERSION,
            "adult_onset": ADULT_ONSET}


STAGE = Stage(
    name="special", version="2", title="Specialised callers (SMN1, GBA1, repeat expansions, HLA)", fn=_special,
    inputs=_inputs,
    outputs=lambda ctx: [out_genes(ctx), out_repeats(ctx), out_hla(ctx), out_claims(ctx), out_summary(ctx)],
    available=_available, params=_params,
)
