"""Publish a dated, immutable *release* that the dashboard reads.

Layout (inside the workspace, e.g. ``Genomics/wgs-data/releases/``)::

    index.json                      list of releases, newest first; ``latest`` points at one
    <release-id>/manifest.json      what this release contains (tables, documents, provenance)
    objects/<name>-<key>.parquet    content-keyed data files, shared between releases
    objects/<name>-<key>.json

Objects are keyed by a hash of their *inputs* and the transform version, so an
unchanged table is never rewritten or duplicated: a new release that only
changes one table costs one table of disk. Each object's own SHA-256 is also
recorded so a copy of the drive can be verified. Releases are never edited
after they are written; ``wgs`` only appends. This is what later lets the
dashboard show "what changed" between two releases (ADR-011).
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
from collections.abc import Callable
from datetime import UTC, datetime
from pathlib import Path

import duckdb

from .. import __version__
from .. import evidence as ev
from ..diff import diff, summarise
from ..knowledge import Store
from ..modules import evaluate
from ..pipeline import CANONICAL, Context, Stage, write_json
from . import annotate, coverage, genotypes, qc, reads, reference, variants

SCHEMA = 1
BIN_BP = 100_000
CHROM_ORDER_SQL = "CASE chrom " + " ".join(f"WHEN '{c}' THEN {i}" for i, c in enumerate(CANONICAL)) + " ELSE 99 END"


def releases_dir(ctx: Context) -> Path:
    p = ctx.cfg.releases_dir
    (p / "objects").mkdir(parents=True, exist_ok=True)
    return p


def index_path(ctx: Context) -> Path:
    return ctx.cfg.releases_dir / "index.json"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _fp(paths: list[Path]) -> str:
    return "|".join(f"{p.name}:{p.stat().st_size}:{int(p.stat().st_mtime)}" for p in paths)


def _copy_sql(sql: str, dest: Path, row_group: int = 65_536) -> None:
    con = duckdb.connect()
    con.execute("SET preserve_insertion_order = true")
    con.execute(f"COPY ({sql}) TO '{dest}' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE {row_group})")
    con.close()


# --- table builders: each returns SQL for a COPY, reading workspace Parquet ------------------------------


def _variants_sql(ctx: Context) -> str:
    vp = variants.variants_parquet(ctx)
    gp = genotypes.out_path(ctx)
    if gp.exists():
        rs = (f"(SELECT chrom, pos, min(rsid) AS rsid FROM '{gp}' WHERE rsid LIKE 'rs%' GROUP BY ALL)")
        return (f"SELECT v.chrom, v.chrom_order, v.pos, coalesce(r.rsid, nullif(v.vcf_id, '.')) AS rsid, "
                f"v.\"ref\", v.alt, v.vtype, v.qual, v.filter, v.gt, v.zygosity, v.gq, v.dp, v.ad_ref, v.ad_alt, "
                f"v.vaf, v.multiallelic FROM '{vp}' v LEFT JOIN {rs} r USING (chrom, pos) "
                f"ORDER BY v.chrom_order, v.pos, v.alt")
    return (f"SELECT chrom, chrom_order, pos, nullif(vcf_id, '.') AS rsid, \"ref\", alt, vtype, qual, filter, gt, "
            f"zygosity, gq, dp, ad_ref, ad_alt, vaf, multiallelic FROM '{vp}' ORDER BY chrom_order, pos, alt")


def _genotypes_sql(ctx: Context) -> str:
    gp = genotypes.out_path(ctx)
    return f"SELECT *, {CHROM_ORDER_SQL} AS chrom_order FROM '{gp}' ORDER BY chrom_order, pos"


def _coverage_bins_sql(ctx: Context) -> str:
    cw = coverage.out_dir(ctx) / "coverage_windows.parquet"
    return (f"SELECT chrom, {CHROM_ORDER_SQL} AS chrom_order, (start // {BIN_BP}) * {BIN_BP} AS start, "
            f"sum(depth * (\"end\" - start)) / sum(\"end\" - start) AS depth "
            f"FROM '{cw}' GROUP BY ALL ORDER BY chrom_order, start")


def _callable_sql(ctx: Context) -> str:
    cp = coverage.out_dir(ctx) / "callable.parquet"
    return f"SELECT chrom, {CHROM_ORDER_SQL} AS chrom_order, start, \"end\", state FROM '{cp}' " \
           f"ORDER BY chrom_order, start"


def _annotations_sql(ctx: Context) -> str:
    return f"SELECT * FROM '{annotate.annotations_path(ctx)}' ORDER BY chrom_order, pos, alt"


def _claims_sql(ctx: Context) -> str:
    return f"SELECT * FROM '{annotate.claims_path(ctx)}' ORDER BY section, claim_id"


def _genes_in(ctx: Context) -> list[Path]:
    if not annotate.annotations_path(ctx).exists():
        return []
    p = Store(ctx.cfg).table("ensembl", "genes")
    return [p] if p else []


def _genes_sql(ctx: Context) -> str:
    return (f"SELECT symbol, gene_id, chrom, start, \"end\", strand, biotype, description "
            f"FROM '{_genes_in(ctx)[0]}' WHERE symbol IS NOT NULL ORDER BY symbol")


def _acmg_in(ctx: Context) -> list[Path]:
    st = Store(ctx.cfg)
    paths = [st.table("acmg_sf", "genes"), st.table("ensembl", "genes"), coverage.out_dir(ctx) / "callable.parquet"]
    return paths if annotate.annotations_path(ctx).exists() and all(paths) else []


def _acmg_sql(ctx: Context) -> str:
    acmg, genes, callable_ = _acmg_in(ctx)
    return f"""
    WITH g AS (SELECT symbol, chrom, start, "end" FROM '{genes}' WHERE symbol IN (SELECT gene FROM '{acmg}')
               QUALIFY row_number() OVER (PARTITION BY symbol ORDER BY "end" - start DESC) = 1),
         cov AS (SELECT g.symbol, sum(least(c."end", g."end") - greatest(c.start, g.start))
                        FILTER (WHERE c.state = 'CALLABLE') AS callable_bp, max(g."end" - g.start) AS span
                 FROM g JOIN '{callable_}' c ON c.chrom = g.chrom AND c.start < g."end" AND c."end" > g.start
                 GROUP BY g.symbol)
    SELECT a.gene, a.category, string_agg(a.condition, '; ') AS conditions, any_value(a.inheritance) AS inheritance,
           any_value(a.report) AS report, any_value(cov.span) AS gene_bp,
           round(coalesce(any_value(cov.callable_bp), 0) / nullif(any_value(cov.span), 0), 4) AS callable_fraction
    FROM '{acmg}' a LEFT JOIN cov ON cov.symbol = a.gene
    GROUP BY a.gene, a.category ORDER BY a.category, a.gene"""


TABLES: list[tuple[str, str, Callable[[Context], list[Path]], Callable[[Context], str], str]] = [
    # name, transform version, inputs, sql, description
    ("variants", "1", lambda c: [variants.variants_parquet(c)] + ([genotypes.out_path(c)]
                                                                    if genotypes.out_path(c).exists() else []),
     _variants_sql, "Every VCF record, one row per ALT allele, with rsIDs from the provider genotype file"),
    ("raw_genotypes", "1", lambda c: [genotypes.out_path(c)], _genotypes_sql,
     "Provider genotype file (23andMe-style), one row per position"),
    ("coverage_bins", "1", lambda c: [coverage.out_dir(c) / "coverage_windows.parquet"], _coverage_bins_sql,
     f"Mean read depth in {BIN_BP // 1000} kb bins along every chromosome"),
    ("callable", "1", lambda c: [coverage.out_dir(c) / "callable.parquet"], _callable_sql,
     "Genome partitioned into CALLABLE / LOW / NO_COVERAGE intervals"),
    ("annotations", "1", lambda c: [annotate.annotations_path(c)], _annotations_sql,
     "Your variant alleles with gene effect (Ensembl/bcftools csq), ClinVar, 1000 Genomes and gnomAD frequencies"),
    ("claims", "1", lambda c: [annotate.claims_path(c)], _claims_sql,
     "Graded statements about your genotypes: evidence strength, call confidence, reasons and sources"),
    ("acmg_genes", "1", _acmg_in, _acmg_sql,
     "ACMG secondary-findings genes checked, with how much of each gene region was callable in your data"),
    ("genes", "1", _genes_in, _genes_sql, "Ensembl gene coordinates (GRCh37) for searching by gene name"),
]


def _documents(ctx: Context) -> dict[str, Path]:
    docs = {
        "qc": qc.out_path(ctx),
        "coverage": coverage.out_dir(ctx) / "coverage.json",
        "alignment_stats": reads.aln_out(ctx),
        "fastq_stats": reads.fq_out(ctx),
        "reference_full": reference.full_decode_file(ctx),
        "annotate": annotate.summary_path(ctx),
    }
    return {k: v for k, v in docs.items() if v.exists()}


def _inventory(ctx: Context) -> dict:
    files = [{"kind": f.kind, "name": f.path.name, "size": f.size,
              "folder": f.path.parent.name} for f in ctx.inv.files if f.kind not in ("index", "bed")]
    return {"sample": ctx.sample, "files": files, "warnings": ctx.inv.warnings}


def _knowledge(ctx: Context, annotated: bool) -> dict:
    """Versions of every knowledge source used (only if the release actually contains annotations)."""
    if not annotated:
        return {"versions": {}, "evidence_model": None, "sources": {}}
    summary = json.loads(annotate.summary_path(ctx).read_text())
    used = summary.get("knowledge", {})
    idx = Store(ctx.cfg).index()["sources"]
    sources = {}
    for sid, ver in used.items():
        s = idx.get(sid, {})
        v = next((x for x in s.get("versions", []) if x["version"] == ver), {})
        sources[sid] = {"title": s.get("title"), "homepage": s.get("homepage"), "licence": s.get("licence"),
                        "cadence": s.get("cadence"), "description": s.get("description"), "version": ver,
                        "built": v.get("built"),
                        "upstream": [{"url": u["url"], "last_modified": u.get("last_modified")}
                                     for u in v.get("upstream", [])]}
    return {"versions": used, "evidence_model": summary.get("evidence_model", ev.MODEL_VERSION),
            "sources": sources}


def _put_object(rd: Path, name: str, key: str, ext: str, make: Callable[[Path], None]) -> dict:
    rel = f"objects/{name}-{key}.{ext}"
    dest = rd / rel
    if not dest.exists():
        tmp = dest.with_suffix(dest.suffix + ".tmp")
        make(tmp)
        os.replace(tmp, dest)
    return {"path": rel, "bytes": dest.stat().st_size, "sha256": _sha256(dest)}


def _key(*parts: str) -> str:
    return hashlib.sha256("\n".join(parts).encode()).hexdigest()[:12]


def _input_paths(ctx: Context) -> list[Path]:
    paths: list[Path] = []
    for _, _, inputs, _, _ in TABLES:
        try:
            ins = inputs(ctx)
        except Exception:
            continue
        paths += [p for p in ins if p.exists()]
    return sorted(set(paths) | set(_documents(ctx).values()))


def build(ctx: Context) -> dict:
    rd = releases_dir(ctx)
    tables = {}
    for name, ver, inputs, sql, desc in TABLES:
        ins = inputs(ctx)
        if not ins or not all(p.exists() for p in ins):
            continue
        key = _key(name, ver, _fp(ins))
        obj = _put_object(rd, name, key, "parquet", lambda dest, s=sql: _copy_sql(s(ctx), dest))
        con = duckdb.connect()
        obj["rows"] = con.execute(f"SELECT count(*) FROM '{rd / obj['path']}'").fetchone()[0]
        con.close()
        tables[name] = {**obj, "description": desc}

    documents = {}
    for name, src in _documents(ctx).items():
        key = _key(name, hashlib.sha256(src.read_bytes()).hexdigest())
        documents[name] = _put_object(rd, name, key, "json", lambda dest, s=src: shutil.copyfile(s, dest))
    for name, payload in (("inventory", _inventory(ctx)), ("sections", evaluate(ctx.inv, ctx.cfg))):
        text = json.dumps(payload, indent=1, default=str)
        documents[name] = _put_object(rd, name, _key(name, text), "json", lambda dest, t=text: dest.write_text(t))

    knowledge = _knowledge(ctx, "annotations" in tables)
    if knowledge["sources"]:
        text = json.dumps(knowledge, indent=1)
        documents["knowledge"] = _put_object(rd, "knowledge", _key("knowledge", text), "json",
                                             lambda dest, t=text: dest.write_text(t))
    content = {"schema": SCHEMA, "sample": ctx.sample, "pipeline_version": __version__,
               "tables": tables, "documents": documents,
               "knowledge": {"versions": knowledge["versions"], "evidence_model": knowledge["evidence_model"]}}
    digest = _key(json.dumps(content, sort_keys=True))

    idx_file = index_path(ctx)
    index = json.loads(idx_file.read_text()) if idx_file.exists() else {"schema": SCHEMA, "releases": []}
    latest = next((r for r in index["releases"] if r["sample"] == ctx.sample), None)
    if latest and latest.get("digest") == digest:
        return {"release": latest["id"], "new": False}

    # "What changed" is computed against the previous release, after the digest so it never causes a release.
    prev = json.loads((rd / latest["manifest"]).read_text()) if latest else None
    changes = diff(prev, rd, {k: rd / v["path"] for k, v in tables.items()}, content["knowledge"]["versions"],
                   content["knowledge"]["evidence_model"])
    changes["summary"] = summarise(changes)
    text = json.dumps(changes, indent=1, default=str)
    documents["changes"] = _put_object(rd, "changes", _key("changes", text), "json",
                                       lambda dest, t=text: dest.write_text(t))

    created = datetime.now(UTC)
    rid = created.strftime("%Y-%m-%d_%H%M%S")
    while (rd / rid).exists():  # two releases in the same second (e.g. two samples)
        rid += "x"
    manifest = {**content, "id": rid, "created": created.strftime("%Y-%m-%dT%H:%M:%SZ"), "digest": digest,
                "previous": latest["id"] if latest else None}
    write_json(rd / rid / "manifest.json", manifest)
    index["releases"].insert(0, {"id": rid, "sample": ctx.sample, "created": manifest["created"],
                                 "digest": digest, "manifest": f"{rid}/manifest.json",
                                 "previous": manifest["previous"], "knowledge": content["knowledge"]["versions"],
                                 "summary": changes["summary"]})
    index["latest"] = rid
    write_json(idx_file, index)
    return {"release": rid, "new": True, "tables": {k: v["rows"] for k, v in tables.items()}}


def _available(ctx: Context) -> str | None:
    return None if variants.variants_parquet(ctx).exists() or genotypes.out_path(ctx).exists() \
        else "no variant or genotype tables yet"


STAGE = Stage(
    name="publish",
    version="2",
    title="Publish a dashboard release",
    fn=build,
    inputs=_input_paths,
    outputs=lambda ctx: [index_path(ctx)],
    available=_available,
)
