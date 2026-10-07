"""Polygenic scores (ADR-017): add up many small genetic effects, then say where you sit among similar people.

A raw polygenic score is just a number; it only means something next to other people's scores — and only
people of similar genetic ancestry, because allele frequencies (and so score distributions) differ between
populations. So this stage follows the PGS Catalog's own calculator (pgsc_calc) step for step, using its
tools, on your genome plus the 2,504 people of the 1000 Genomes reference panel:

1. **Reference QC** (cached): unrelated people, common biallelic SNPs, LD-thinned → the variants used for
   ancestry PCA; plus reference allele frequencies.
2. **Genotype you** at every reference-panel variant a selected score uses (and the PCA variants) from your
   VCF *and* the callable map — a site with no variant record counts as reference only where coverage was
   good; anything uncertain is left missing, never guessed.
3. **Ancestry PCA**: FRAPOSA projects you onto the reference's principal components.
4. **Match** each score's variants to the panel (pgscatalog-match: strand flips handled, ambiguous A/T and
   G/C SNPs dropped, scores with < 75 % of variants matched are rejected), then **score** you and the panel
   with plink2 (missing genotypes filled with the reference allele frequency, as pgsc_calc does).
5. **Ancestry-adjust** (pgscatalog-ancestry-adjust): your percentile among the reference population most
   similar to you, plus PCA-normalised Z-scores.

Heavy tools live in the ``pgs`` and ``fraposa`` pixi environments; results are small tables under work/pgs.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import shutil
from collections.abc import Iterator
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from .. import evidence as ev
from ..knowledge import Store
from ..pipeline import Context, Stage, run, write_json
from . import coverage, qc, variants
from .pgx import MIN_GQ, _env_tool

HIGH_LD = Path(__file__).parents[1] / "knowledge" / "data" / "high-LD-regions-hg19-GRCh37.txt"
MIN_OVERLAP = 0.75
REF_SCORE_VERSION = "1"  # bump if the reference scoring command changes
NEAR_PAD = 25  # bp: how far apart two spellings of the same indel can sit in short repeats
MATCH_BATCH_BYTES = 400_000_000  # compressed scoring-file bytes per pgscatalog-match call (~8 GB of RAM)
AUTOSOMES = [str(i) for i in range(1, 23)]
REF_QC = ["--max-alleles", "2", "--snps-only", "just-acgt", "--geno", "0.1",
          "--mind", "0.1", "--maf", "0.05", "--hwe", "1e-4", "--autosome"]


def _files(d: Path, pattern: str) -> Iterator[Path]:
    """Glob that skips macOS AppleDouble (`._*`) sidecars written onto exFAT drives."""
    return iter(sorted(p for p in d.glob(pattern) if not p.name.startswith("._")))


def pgs_dir(ctx: Context) -> Path:
    p = ctx.work / "pgs"
    p.mkdir(parents=True, exist_ok=True)
    return p


def out_scores(ctx: Context) -> Path:
    return pgs_dir(ctx) / "scores.parquet"


def out_summary(ctx: Context) -> Path:
    return pgs_dir(ctx) / "pgs.json"


def out_ancestry(ctx: Context) -> Path:
    return pgs_dir(ctx) / "ancestry.parquet"


def out_genotype_qc(ctx: Context) -> Path:
    return pgs_dir(ctx) / "genotype_qc.parquet"


def _tool(name: str) -> str | None:
    return _env_tool("pgs", name)


def _plink_memory_mb() -> int:
    """plink2's workspace: the full panel (78 M variants) needs ~8 GB; half of RAM (min 6 GB) avoids swapping."""
    try:
        total = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (AttributeError, ValueError, OSError):
        total = 16 * 2**30
    return max(6000, int(total / 2**20 * 0.5))


def _plink(ctx: Context, *args: str, log: str) -> None:
    run([_tool("plink2"), "--threads", str(ctx.threads), "--memory", str(_plink_memory_mb()), *args],
        log=ctx.logs / f"pgs-{log}.log")


def _panel(ctx: Context) -> dict:
    st = Store(ctx.cfg)
    cur = st.current("pgs_reference")
    root = st.root / "pgs_reference" / cur["version"]
    pgen = next(_files(root, "GRCh37_*.pgen"))
    stem = pgen.with_suffix("")
    king = next(_files(root, "*king.cutoff.out.id"), None)
    return {"stem": stem, "psam": stem.with_suffix(".psam"), "king": king, "version": cur["version"],
            "pvar": next(_files(root, f"{stem.name}.pvar*"))}


def _catalog(ctx: Context) -> tuple[Path, list[dict], str]:
    st = Store(ctx.cfg)
    cur = st.current("pgs_catalog")
    root = st.root / "pgs_catalog" / cur["version"]
    sel = duckdb.sql(f"SELECT * FROM '{root / 'selected.parquet'}'").fetchdf().to_dict("records")
    return root / "scoring", sel, cur["version"]


# ---------------------------------------------------------------------------------------------------------------
# 1. Reference QC (cached by panel version)
# ---------------------------------------------------------------------------------------------------------------


def _reference_qc(ctx: Context, panel: dict) -> Path:
    d = pgs_dir(ctx) / "ref" / panel["version"]
    done = d / "done"
    if done.exists():
        return d
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    vzs = "vzs" if panel["pvar"].name.endswith(".zst") else ""
    pf = ["--pfile", str(panel["stem"])] + ([vzs] if vzs else [])
    rm = ["--remove", str(panel["king"])] if panel["king"] else []
    _plink(ctx, *pf, *rm, *REF_QC, "--make-pgen", "--out", str(d / "qc"), log="ref-qc")
    _plink(ctx, "--pfile", str(d / "qc"), "--indep-pairwise", "1000", "50", "0.05",
           "--exclude", "range", str(HIGH_LD), "--out", str(d / "thin"), log="ref-thin")
    # Frequencies over the whole panel (for filling a missing genotype with its expected dosage).
    _plink(ctx, *pf, "--freq", "--out", str(d / "all"), log="ref-freq")
    done.write_text("ok\n")
    return d


# ---------------------------------------------------------------------------------------------------------------
# 2. Genotype the sample at the panel's score + PCA sites
# ---------------------------------------------------------------------------------------------------------------


def _score_positions(con: duckdb.DuckDBPyConnection, scoring: Path, ids: list[str]) -> None:
    files = [str(scoring / f"{i}_hmPOS_GRCh37.txt.gz") for i in ids]
    con.execute(f"""CREATE TABLE score_pos AS SELECT DISTINCT CAST(hm_chr AS VARCHAR) AS chrom,
        CAST(hm_pos AS INTEGER) AS pos FROM read_csv({files!r}, comment='#', delim='\t', header=true,
        union_by_name=true, all_varchar=true) WHERE hm_pos IS NOT NULL AND hm_pos <> ''
        AND CAST(hm_chr AS VARCHAR) IN ({",".join(repr(c) for c in AUTOSOMES)})""")


def _pvar_header(path: Path) -> tuple[int, list[str]]:
    """(lines before the data, column names from the #CHROM line)."""
    import io

    import zstandard

    with open(path, "rb") as fh:
        stream = zstandard.ZstdDecompressor().stream_reader(fh) if path.name.endswith(".zst") else fh
        for n, line in enumerate(io.TextIOWrapper(stream, encoding="utf-8"), 1):
            if line.startswith("#CHROM"):
                return n, line[1:].rstrip("\n").split("\t")
    raise ValueError(f"{path}: no #CHROM header")


def _read_pvar(path: Path) -> str:
    comp = "zstd" if path.name.endswith(".zst") else "none"
    types = {"POS": "INTEGER"}
    skip, names = _pvar_header(path)
    cols = ", ".join(f"'{c}': '{types.get(c, 'VARCHAR')}'" for c in names)
    return (f"read_csv('{path}', delim='\t', skip={skip}, header=false, compression='{comp}', "
            f"columns={{{cols}}}, quote='', escape='', auto_detect=false)")


def genotype_sql(sites: str, vp: Path, cp: Path) -> str:
    """SQL giving (chrom, pos, id, ref, alt, gt) for every site in table `sites`.

    The same rule as the PharmCAT genotyper, vectorised: an exact PASS call with GQ ≥ MIN_GQ gives a dosage;
    a different allele at the position, a no-call, a filtered or low-quality call, a site inside one of your
    deletions, a different indel within NEAR_PAD bp of an indel site, or an overlapping MNP is missing;
    otherwise a confident RefCall or a callable position is homozygous reference."""
    return f"""
    WITH v AS (SELECT chrom, pos, ref, alt, filter, replace(gt, '|', '/') AS gt, gq FROM '{vp}'
               WHERE chrom IN ({",".join(repr(c) for c in AUTOSOMES)})),
    atpos AS (SELECT s.chrom, s.pos,
                count(*) FILTER (WHERE v.gt NOT IN ('0/0', './.') AND v.filter <> 'RefCall') AS n_alt,
                bool_or(v.gt = './.' OR v.filter = 'NoCall') AS nocall,
                bool_or(v.filter = 'RefCall' AND v.gq >= {MIN_GQ}) AS refcall
              FROM (SELECT DISTINCT chrom, pos FROM {sites}) s JOIN v USING (chrom, pos) GROUP BY ALL),
    exact AS (SELECT s.id, v.gt, v.filter, v.gq FROM {sites} s JOIN v
              ON v.chrom = s.chrom AND v.pos = s.pos AND v.ref = s.ref AND v.alt = s.alt
              WHERE v.filter <> 'RefCall'),
    dels AS (SELECT chrom, pos + 1 AS a, pos + length(ref) - 1 AS b FROM v
             WHERE length(ref) > length(alt) AND gt NOT IN ('0/0', './.') AND filter <> 'RefCall'),
    indel AS (SELECT DISTINCT s.id FROM {sites} s JOIN dels d ON d.chrom = s.chrom AND s.pos BETWEEN d.a AND d.b),
    -- The same indel can be written at different positions (left/right alignment in repeats), and a multi-base
    -- substitution can hide a SNP: a nearby non-matching indel or an overlapping MNP makes the site unknown.
    va AS (SELECT chrom, pos AS a, pos + length(ref) - 1 AS b, length(ref) <> length(alt) AS is_indel, ref, alt
           FROM v WHERE gt NOT IN ('0/0', './.') AND filter <> 'RefCall'
             AND (length(ref) <> length(alt) OR length(ref) > 1)),
    vb AS (SELECT *, unnest(range((a - {NEAR_PAD}) // 1000, (b + {NEAR_PAD}) // 1000 + 1)) AS bin FROM va),
    sb AS (SELECT id, chrom, pos, ref, alt, pos + length(ref) - 1 AS b, length(ref) <> length(alt) AS is_indel,
             unnest(range(pos // 1000, (pos + length(ref) - 1) // 1000 + 1)) AS bin FROM {sites}),
    near AS (SELECT DISTINCT sb.id FROM sb JOIN vb ON vb.chrom = sb.chrom AND vb.bin = sb.bin
             WHERE NOT (vb.a = sb.pos AND vb.ref = sb.ref AND vb.alt = sb.alt)
               AND CASE WHEN sb.is_indel AND vb.is_indel THEN vb.a <= sb.b + {NEAR_PAD} AND vb.b >= sb.pos - {NEAR_PAD}
                        WHEN NOT vb.is_indel THEN vb.a <= sb.b AND vb.b >= sb.pos ELSE false END),
    callable AS (SELECT DISTINCT s.id FROM {sites} s JOIN '{cp}' c ON c.chrom = s.chrom
                 AND s.pos > c.start AND s.pos <= c."end" WHERE c.state = 'CALLABLE')
    SELECT s.chrom, s.pos, s.id, s.ref, s.alt,
      CASE
        WHEN e.id IS NOT NULL THEN
          CASE WHEN coalesce(a.n_alt, 0) > 1 OR e.filter NOT IN ('PASS', '.') OR coalesce(e.gq, 0) < {MIN_GQ}
                 THEN './.'
               WHEN e.gt = '1/1' THEN '1/1'
               WHEN e.gt IN ('0/1', '1/0') THEN '0/1'
               ELSE './.' END
        WHEN coalesce(a.n_alt, 0) > 0 OR coalesce(a.nocall, false) OR i.id IS NOT NULL OR n.id IS NOT NULL
          THEN './.'
        WHEN coalesce(a.refcall, false) OR c.id IS NOT NULL THEN '0/0'
        ELSE './.'
      END AS gt
    FROM {sites} s
    LEFT JOIN exact e ON e.id = s.id
    LEFT JOIN atpos a ON a.chrom = s.chrom AND a.pos = s.pos
    LEFT JOIN indel i ON i.id = s.id
    LEFT JOIN near n ON n.id = s.id
    LEFT JOIN callable c ON c.id = s.id
    """


def _target(ctx: Context, panel: dict, ref: Path, scoring: Path, ids: list[str]) -> dict:
    d = pgs_dir(ctx) / "target"
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    tmp = ctx.scratch / "pgs"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    con = duckdb.connect(str(tmp / "pgs.duckdb"))
    con.execute(f"SET temp_directory='{tmp}'; SET memory_limit='6GB'; SET threads={ctx.threads}")
    _score_positions(con, scoring, ids)
    con.execute(f"CREATE TABLE pruned AS SELECT column0 AS id FROM read_csv('{ref / 'thin.prune.in'}', "
                "header=false, columns={'column0':'VARCHAR'})")
    con.execute(f"""CREATE TABLE sites AS SELECT CHROM AS chrom, POS AS pos, ID AS id, REF AS ref, ALT AS alt
        FROM {_read_pvar(panel['pvar'])} p
        WHERE CHROM IN ({",".join(repr(c) for c in AUTOSOMES)}) AND ALT NOT LIKE '%,%'
          AND ((CHROM, POS) IN (SELECT (chrom, pos) FROM score_pos) OR ID IN (SELECT id FROM pruned))""")
    vp, cp = variants.variants_parquet(ctx), coverage.out_dir(ctx) / "callable.parquet"
    con.execute(f"CREATE TABLE geno AS {genotype_sql('sites', vp, cp)}")
    stats = dict(zip(["sites", "called", "alt_alleles"], con.execute(
        "SELECT count(*), count(*) FILTER (WHERE gt <> './.'), "
        "sum(CASE gt WHEN '0/1' THEN 1 WHEN '1/1' THEN 2 ELSE 0 END) FROM geno").fetchone(), strict=True))
    # Calibration: your genotypes vs how common each allele is in the panel. Being called reference where ~everyone
    # carries the other allele means the calls are wrong (an indel-spelling mismatch once did exactly this).
    con.execute(f"""CREATE TABLE af AS SELECT ID AS id, TRY_CAST(ALT_FREQS AS DOUBLE) AS af
        FROM read_csv('{ref / 'all.afreq'}', delim='\t', header=true, all_varchar=true)""")
    con.execute(f"""COPY (WITH g AS (SELECT CASE WHEN length(ref) = 1 AND length(alt) = 1 THEN 'SNP' ELSE 'indel' END
            AS kind, af, gt FROM geno JOIN af USING (id) WHERE af IS NOT NULL),
          b AS (SELECT *, CASE WHEN af < 0.01 THEN 1 WHEN af < 0.1 THEN 2 WHEN af < 0.5 THEN 3 WHEN af < 0.9 THEN 4
                               WHEN af < 0.99 THEN 5 ELSE 6 END AS band FROM g)
        SELECT kind, band, ['<1%', '1–10%', '10–50%', '50–90%', '90–99%', '>99%'][band] AS panel_af_band,
          count(*) AS sites, avg(af) AS panel_af,
          100 * avg((gt = '0/0')::INT) AS hom_ref_pct, 100 * avg((gt = '0/1')::INT) AS het_pct,
          100 * avg((gt = '1/1')::INT) AS hom_alt_pct, 100 * avg((gt = './.')::INT) AS missing_pct,
          sum(CASE gt WHEN '0/1' THEN 1 WHEN '1/1' THEN 2 ELSE 0 END)
            / nullif(2 * count(*) FILTER (WHERE gt <> './.'), 0)
            AS your_af
        FROM b GROUP BY ALL ORDER BY kind, band) TO '{out_genotype_qc(ctx)}' (FORMAT parquet)""")
    vcf = tmp / "target.vcf"
    with open(vcf, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        for c in AUTOSOMES:
            fh.write(f"##contig=<ID={c}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tself\n")
    body = tmp / "body.tsv"
    con.execute(f"""COPY (SELECT chrom, pos, id, ref, alt, '.', '.', '.', 'GT', gt FROM geno
        ORDER BY CAST(chrom AS INTEGER), pos, id) TO '{body}' (HEADER false, DELIMITER '\t', QUOTE '')""")
    with open(vcf, "a") as out, open(body) as src:
        shutil.copyfileobj(src, out)
    con.close()
    _plink(ctx, "--vcf", str(vcf), "--make-pgen", "--out", str(d / "self"), log="target")
    shutil.rmtree(tmp, ignore_errors=True)
    return stats


# ---------------------------------------------------------------------------------------------------------------
# 3. Ancestry PCA with FRAPOSA
# ---------------------------------------------------------------------------------------------------------------


def _pca(ctx: Context, panel: dict, ref: Path) -> dict:
    d = pgs_dir(ctx) / "pca"
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    tgt = pgs_dir(ctx) / "target" / "self"
    # PCA variants you were confidently genotyped at.
    _plink(ctx, "--pfile", str(tgt), "--extract", str(ref / "thin.prune.in"), "--geno", "0",
           "--write-snplist", "--out", str(d / "shared"), log="pca-shared")
    _plink(ctx, "--pfile", str(ref / "qc"), "--extract", str(d / "shared.snplist"), "--make-bed",
           "--out", str(d / "reference"), log="pca-ref")
    _plink(ctx, "--pfile", str(tgt), "--extract", str(d / "shared.snplist"),
           "--alt1-allele", str(d / "reference.bim"), "5", "2", "--make-bed", "--out", str(d / "self"),
           log="pca-target")
    fr = _env_tool("fraposa", "fraposa")
    run([fr, str(d / "reference"), "--method", "oadp", "--dim_ref", "10"], log=ctx.logs / "pgs-fraposa-ref.log",
        cwd=d)
    run([fr, str(d / "reference"), "--method", "oadp", "--dim_ref", "10", "--stu_filepref", str(d / "self"),
         "--out", str(d / "self")], log=ctx.logs / "pgs-fraposa-self.log", cwd=d)
    n = len((d / "shared.snplist").read_text().splitlines())
    return {"pca_variants": n}


# ---------------------------------------------------------------------------------------------------------------
# 4. Match + score
# ---------------------------------------------------------------------------------------------------------------


def _normalise(ctx: Context, scoring: Path, ids: list[str], version: str) -> list[str]:
    """pgscatalog-format each scoring file once per catalog version (the 6.6 M-variant scores take ~10 min)."""
    norm = pgs_dir(ctx) / "normalised" / version
    want = {f"normalised_{i}_hmPOS_GRCh37.txt.gz" for i in ids}
    have = {p.name for p in _files(norm, "normalised_*")} if (norm / "done").exists() else set()
    if not want <= have:
        shutil.rmtree(norm, ignore_errors=True)
        norm.mkdir(parents=True)
        files = [str(scoring / f"{i}_hmPOS_GRCh37.txt.gz") for i in ids]
        run([_tool("pgscatalog-format"), "-s", *files, "-t", "GRCh37", "-o", str(norm), "-l", "format_log.json",
             "--threads", str(min(ctx.threads, 4))], log=ctx.logs / "pgs-format.log")
        (norm / "done").write_text("ok\n")
    return sorted(str(norm / n) for n in want)


def _match(ctx: Context, scoring: Path, ids: list[str], version: str) -> dict:
    d = pgs_dir(ctx) / "match"
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    normed = _normalise(ctx, scoring, ids, version)
    pvar = pgs_dir(ctx) / "target" / "self.pvar"
    # One target file holds every chromosome, so one pass replaces pgsc_calc's per-chromosome match + merge.
    # --keep_ambiguous: A/T and C/G SNPs are dropped by default because array data may be on either strand.
    # Here both sides are on the GRCh37 forward strand (harmonised hm_ files; calls against the reference),
    # so they can be matched safely, as pgsc_calc's documentation allows.
    # Matching writes ~10 GB of temporary IPC files: do it on the fast scratch disk, keep only the results.
    # pgscatalog-match holds every scoring file in memory at once; 220+ scores need more than 16 GB. Each score is
    # matched on its own, so matching in batches gives identical results: outputs are renamed self_b{n}_* and the
    # summaries are concatenated.
    tmp = ctx.scratch / "pgs-match"
    batches = match_batches(normed)
    summaries = []
    try:
        for n, batch in enumerate(batches):
            shutil.rmtree(tmp, ignore_errors=True)
            tmp.mkdir(parents=True)
            run([_tool("pgscatalog-match"), "-d", "self", "-s", *batch, "-t", str(pvar), "--min_overlap",
                 str(MIN_OVERLAP), "--outdir", str(tmp), "--combined", "--keep_ambiguous"],
                log=ctx.logs / "pgs-match.log")
            for f in tmp.iterdir():
                if f.is_file():
                    name = f"self_b{n}_{f.name.removeprefix('self_')}"
                    shutil.copyfile(f, d / name)
                    if name.endswith("_summary.csv"):
                        summaries.append(d / name)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    with open(d / "self_summary.csv", "w") as out:
        for i, f in enumerate(summaries):
            lines = f.read_text().splitlines(keepends=True)
            out.writelines(lines if i == 0 else lines[1:])
    return {"match_batches": len(batches)}


def match_batches(files: list[str], cap: int = MATCH_BATCH_BYTES) -> list[list[str]]:
    """Greedy-pack scoring files (largest first) into batches whose compressed size stays under `cap`."""
    batches: list[list[str]] = []
    sizes: list[int] = []
    for f in sorted(files, key=lambda f: -os.path.getsize(f)):
        sz = os.path.getsize(f)
        for i, total in enumerate(sizes):
            if total + sz <= cap:
                batches[i].append(f)
                sizes[i] += sz
                break
        else:
            batches.append([f])
            sizes.append(sz)
    return [sorted(b) for b in batches]


def _score(ctx: Context, ref: Path, panel: dict) -> list[Path]:
    d = pgs_dir(ctx) / "score"
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    tgt = pgs_dir(ctx) / "target" / "self"
    out = []
    vzs = ["vzs"] if panel["pvar"].name.endswith(".zst") else []
    refsub: Path | None = None
    for sf in sorted(_files(pgs_dir(ctx) / "match", "*.scorefile.gz")):
        stem = sf.name.removesuffix(".scorefile.gz").removeprefix("self_")
        # Restrict to the scorefile's own IDs: the 1000G panel has long indels with ID "." that would otherwise
        # trip plink's duplicate-ID check under --read-freq (none of them can be in a matched scorefile).
        ids = ctx.scratch / f"pgs-{stem}.ids"
        with gzip.open(sf, "rt") as fh, open(ids, "w") as idf:
            ncol = len(fh.readline().rstrip("\n").split("\t"))
            for line in fh:
                idf.write(line.split("\t", 1)[0] + "\n")
        cols = f"3-{ncol}" if ncol > 3 else "3"
        common = ["--extract", str(ids), "--score", str(sf), "header-read", "cols=+scoresums,+denom", "list-variants",
                  "--score-col-nums", cols, "--read-freq", str(ref / "all.afreq")]
        _plink(ctx, "--pfile", str(tgt), *common, "--out", str(d / f"self_{stem}"), log=f"score-self-{stem}")
        # The panel's scores depend only on the scorefile and the panel, never on your genotypes: cache them
        # (re-scoring needs one pass over the full panel to build the subset, ~30 min from an external drive).
        h = hashlib.sha256(sf.read_bytes() + f"|{panel['version']}|{cols}|{REF_SCORE_VERSION}".encode())
        cache = pgs_dir(ctx) / "score-ref-cache" / f"{stem}-{h.hexdigest()[:16]}.sscore"
        if not cache.exists():
            if refsub is None:
                refsub = _ref_subset(ctx, panel, vzs)
            _plink(ctx, "--pfile", str(refsub), *common, "--out", str(d / f"reference_{stem}"),
                   log=f"score-ref-{stem}")
            cache.parent.mkdir(parents=True, exist_ok=True)
            for old in _files(cache.parent, f"{stem}-*.sscore"):  # keep only the latest per chunk
                old.unlink()
            shutil.copyfile(d / f"reference_{stem}.sscore", cache)
        else:
            shutil.copyfile(cache, d / f"reference_{stem}.sscore")
        ids.unlink()
        out += [d / f"self_{stem}.sscore", d / f"reference_{stem}.sscore"]
    if refsub is not None:
        for f in refsub.parent.glob(refsub.name + ".*"):
            f.unlink()
    return out


def _ref_subset(ctx: Context, panel: dict, vzs: list[str]) -> Path:
    """Copy just the panel variants that any matched scorefile uses to the internal scratch disk, once per run.

    The full panel is 78 M variants (7 GB) on the external drive; reading it once per scorefile chunk took
    30–75 min each. Scoring a few-million-variant subset on the internal disk takes minutes and gives identical
    sums, because --score only ever reads the variants it is given.
    """
    union = ctx.scratch / "pgs-ref-union.ids"
    seen: set[str] = set()
    with open(union, "w") as out:
        for sf in sorted(_files(pgs_dir(ctx) / "match", "*.scorefile.gz")):
            with gzip.open(sf, "rt") as fh:
                next(fh)
                for line in fh:
                    vid = line.split("\t", 1)[0]
                    if vid not in seen:
                        seen.add(vid)
                        out.write(vid + "\n")
    sub = ctx.scratch / "pgs-refsub"
    _plink(ctx, "--pfile", str(panel["stem"]), *vzs, "--extract", str(union), "--make-pgen", "--out", str(sub),
           log="score-ref-subset")
    union.unlink()
    return sub


# ---------------------------------------------------------------------------------------------------------------
# 5. Aggregate + ancestry adjustment
# ---------------------------------------------------------------------------------------------------------------


def _adjust(ctx: Context, panel: dict, sscores: list[Path]) -> Path:
    d = pgs_dir(ctx) / "adjust"
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    run([_tool("pgscatalog-aggregate"), "-s", *map(str, sscores), "-o", str(d), "--no-split"],
        log=ctx.logs / "pgs-aggregate.log")
    agg = next(_files(d, "aggregated_scores.txt.gz"))
    pca = pgs_dir(ctx) / "pca"
    args = [_tool("pgscatalog-ancestry-adjust"), "-d", "self", "-r", "reference", "--psam", str(panel["psam"]),
            "--ref_pcs", str(pca / "reference.pcs"), "--target_pcs", str(pca / "self.pcs"), "-p", "SuperPop",
            "-s", str(agg), "-a", "RandomForest", "--n_popcomp", "5", "-n", "empirical", "mean", "mean+var",
            "--n_normalization", "4", "--outdir", str(d)]
    if panel["king"]:
        args[args.index("-p"):args.index("-p")] = ["-x", str(panel["king"])]
    run(args, log=ctx.logs / "pgs-adjust.log")
    return d


def _clean(rows: list[dict]) -> list[dict]:
    return [{k: (None if isinstance(v, float) and v != v else v) for k, v in r.items()} for r in rows]


def best_evaluation(evals: list[dict]) -> dict | None:
    """The evaluation to show: independent before the developers' own, then one reporting an effect size,
    then the largest sample."""
    def key(e: dict) -> tuple:
        has_effect = any(e.get(k) for k in ("or", "hr", "beta", "auroc", "cindex"))
        return (bool(e.get("independent")), has_effect, e.get("n") or 0)
    return max(evals, key=key) if evals else None


def effect_text(e: dict | None) -> str | None:
    if not e:
        return None
    parts = [f"{lab} {e[k]}" for k, lab in (("or", "OR"), ("hr", "HR"), ("beta", "β"), ("auroc", "AUROC"),
                                            ("cindex", "C-index")) if e.get(k)]
    return "; ".join(parts) or e.get("other")


def promote(sel: list[dict], scored: dict[str, dict]) -> list[dict]:
    """For each featured trait, show the best-ranked score that could actually be scored on your data.

    If the selection rule's pick failed (too few variants matched), the first scorable stand-in takes its place
    and the pick moves to the browse list with its reason; stand-ins that aren't needed are dropped."""
    def ok(r: dict) -> bool:
        you = scored.get(r["pgs_id"]) or {}
        return bool(you.get("passed")) and you.get("percentile") is not None

    alts: dict[str, list[dict]] = {}
    for r in sel:
        if r.get("alternate_for"):
            alts.setdefault(r["alternate_for"], []).append(r)
    out = []
    for r in sel:
        if r.get("alternate_for"):
            continue
        sub = next((a for a in alts.get(r["trait_id"], []) if ok(a)), None) if r["featured"] and not ok(r) else None
        if sub:
            out.append({**sub, "featured": True, "alternate_for": "",
                        "why": f"Stands in for {r['pgs_id']}, the selection rule's pick, which could not be scored "
                               "on your data"})
            out.append({**r, "featured": False, "section": "", "why": f"Replaced by {sub['pgs_id']}"})
        else:
            out.append(r)
    return out


def grade(sel: list[dict], scored: dict[str, dict], root: Path) -> list[dict]:
    """One row per selected score: your result, catalog metadata, best evaluation and the Model 3 grades."""
    con = duckdb.connect()
    meta = {r["pgs_id"]: r for r in _clean(con.execute(f"""SELECT s.pgs_id, s.name, s.trait_reported, s.method,
        s.ancestry_gwas, s.ancestry_training, s.ancestry_evaluation, s.released, s.licence, s.pmid, s.doi,
        p.first_author, p.title, p.journal, p.published FROM '{root / "scores.parquet"}' s
        LEFT JOIN '{root / "publications.parquet"}' p USING (pgp_id)""").fetchdf().to_dict("records"))}
    about = {r["trait_id"]: r for r in _clean(con.execute(f"""SELECT trait_id, description AS trait_description,
        url AS trait_url FROM '{root / "traits.parquet"}'""").fetchdf().to_dict("records"))}
    evals: dict[str, list[dict]] = {}
    for e in _clean(con.execute(f"""SELECT e.*, p.first_author AS eval_author, p.published AS eval_published,
            p.pmid AS eval_pmid FROM '{root / "evaluations.parquet"}' e
            LEFT JOIN '{root / "publications.parquet"}' p USING (pgp_id)""").fetchdf().to_dict("records")):
        evals.setdefault(e["pgs_id"], []).append(e)
    out = []
    for r in sel:
        pid = r["pgs_id"]
        you = scored.get(pid, {})
        m = meta.get(pid, {})
        b = best_evaluation(evals.get(pid, []))
        e_level, e_reasons = ev.pgs_evidence(independent_pubs=r["independent_pubs"],
                                             independent_n=r["independent_n"], evaluated_n=r["evaluated_n"])
        c_level, c_reasons = ev.pgs_call(match_rate=you.get("match_rate"), passed=you.get("passed"))
        if c_level != "Not callable" and you.get("percentile") is None:
            c_level, c_reasons = "Not callable", c_reasons + ["No ancestry-adjusted percentile was produced"]
        out.append({
            **{k: r[k] for k in ("pgs_id", "trait_id", "trait_label", "label", "section", "body_system", "sex",
                                 "featured", "pinned", "why", "independent_pubs", "independent_n", "evaluated_n",
                                 "n_variants")},
            **{k: you.get(k) for k in ("raw_score", "percentile", "z", "z_norm1", "z_norm2", "compared_with",
                                       "matched", "total", "match_rate", "passed")},
            **{k: m.get(k) for k in ("name", "trait_reported", "method", "ancestry_gwas", "ancestry_training",
                                     "ancestry_evaluation", "released", "licence", "pmid", "doi", "first_author",
                                     "title", "journal", "published")},
            **{k: about.get(r["trait_id"], {}).get(k) for k in ("trait_description", "trait_url")},
            "eval_effect": effect_text(b), "eval_n": b and b.get("n"), "eval_cases": b and b.get("cases"),
            "eval_ancestry": b and b.get("ancestry"), "eval_independent": b and b.get("independent"),
            "eval_trait": b and b.get("trait_reported"), "eval_covariates": b and b.get("covariates"),
            "eval_author": b and b.get("eval_author"), "eval_published": b and b.get("eval_published"),
            "eval_pmid": b and b.get("eval_pmid"),
            "evidence_level": e_level, "evidence_reasons": e_reasons, "call_confidence": c_level,
            "call_reasons": c_reasons, "overall": ev.overall(e_level, c_level),
        })
    return out


def _collect(ctx: Context, d: Path, sel: list[dict], stats: dict, version: str, catalog_root: Path) -> dict:
    con = duckdb.connect()
    pgs = next(_files(d, "self_pgs.txt.gz"))
    pop = next(_files(d, "self_popsimilarity.txt.gz"))
    summary = pgs_dir(ctx) / "match" / "self_summary.csv"
    con.execute(f"CREATE TABLE sel AS SELECT * FROM (VALUES {','.join('(?,?)' for _ in sel)}) t(pgs_id, trait_id)",
                [x for r in sel for x in (r["pgs_id"], r["trait_id"])])
    con.execute(f"CREATE TABLE raw AS SELECT * FROM read_csv('{pgs}', delim='\t', header=true)")
    cols = {c for (c,) in con.execute("SELECT column_name FROM (DESCRIBE raw)").fetchall()}
    want = [c for c in ("PGS", "SUM", "percentile_MostSimilarPop", "Z_MostSimilarPop", "Z_norm1", "Z_norm2")
            if c in cols]
    # The comparison population lives in the population-similarity file, not the score file.
    con.execute(f"""CREATE TABLE s AS SELECT {",".join(f'r."{c}"' for c in want)}, p.MostSimilarPop AS pop
        FROM raw r LEFT JOIN read_csv('{pop}', delim='\t', header=true, all_varchar=true) p
          ON p.sampleset = r.sampleset AND p.IID = r.IID
        WHERE r.sampleset = 'self' AND r.IID = 'self'""")
    con.execute(f"""CREATE TABLE m AS SELECT split_part(accession, '_', 1) AS pgs_id,
        sum("count") FILTER (WHERE match_status = 'matched') AS matched,
        sum("count") AS total, bool_or(score_pass) AS passed
        FROM read_csv('{summary}', header=true) GROUP BY 1""") if summary.exists() else \
        con.execute("CREATE TABLE m (pgs_id VARCHAR, matched BIGINT, total BIGINT, passed BOOLEAN)")
    rows = con.execute("""SELECT sel.pgs_id, s.SUM AS raw_score,
        s.percentile_MostSimilarPop AS percentile, s.Z_MostSimilarPop AS z, s.Z_norm1 AS z_norm1,
        s.Z_norm2 AS z_norm2, s.pop AS compared_with, m.matched, m.total,
        m.matched / nullif(m.total, 0) AS match_rate, coalesce(m.passed, false) AS passed
        FROM sel LEFT JOIN s ON split_part(s.PGS, '_', 1) = sel.pgs_id
        LEFT JOIN m ON m.pgs_id = sel.pgs_id""").fetchdf()
    scored = {r["pgs_id"]: r for r in _clean(rows.to_dict("records"))}
    pq.write_table(pa.Table.from_pylist(grade(promote(sel, scored), scored, catalog_root)), out_scores(ctx))
    con.execute(f"""COPY (SELECT * FROM read_csv('{pop}', delim='\t', header=true))
        TO '{out_ancestry(ctx)}' (FORMAT parquet)""")
    me = con.execute(f"""SELECT * FROM read_csv('{pop}', delim='\t', header=true)
        WHERE sampleset = 'self' AND IID = 'self'""").fetchdf().to_dict("records")
    n_ok = con.execute(f"SELECT count(*) FROM '{out_scores(ctx)}' WHERE percentile IS NOT NULL").fetchone()[0]
    st = Store(ctx.cfg)
    res = {"catalog": version, "scored": n_ok,
           "knowledge": {"pgs_catalog": version, "pgs_reference": (st.current("pgs_reference") or {}).get("version")},
           "evidence_model": ev.MODEL_VERSION, "selected": len(sel), "genotyping": stats,
           "most_similar_population": me[0].get("MostSimilarPop") if me else None,
           "rf_probabilities": {k.removeprefix("RF_P_"): round(v, 4) for k, v in (me[0] if me else {}).items()
                                if k.startswith("RF_P_")}}
    write_json(out_summary(ctx), res)
    return res


def _pgs(ctx: Context) -> dict:
    panel = _panel(ctx)
    scoring, sel, version = _catalog(ctx)
    ids = sorted({r["pgs_id"] for r in sel})
    ref = _reference_qc(ctx, panel)
    stats = _target(ctx, panel, ref, scoring, ids)
    stats.update(_pca(ctx, panel, ref))
    _match(ctx, scoring, ids, version)
    sscores = _score(ctx, ref, panel)
    d = _adjust(ctx, panel, sscores)
    res = _collect(ctx, d, sel, stats, version, scoring.parent)
    return {k: res[k] for k in ("scored", "selected", "most_similar_population")}


def _available(ctx: Context) -> str | None:
    st = Store(ctx.cfg)
    if not st.current("pgs_catalog"):
        return "needs the pgs_catalog knowledge source (wgs knowledge refresh --source pgs_catalog)"
    if not st.current("pgs_reference"):
        return "needs the pgs_reference knowledge source (wgs knowledge refresh --source pgs_reference)"
    if not _tool("plink2") or not _env_tool("fraposa", "fraposa"):
        return "PGS tools are not installed (pixi install -e pgs -e fraposa)"
    if not variants.variants_parquet(ctx).exists() or not (coverage.out_dir(ctx) / "callable.parquet").exists():
        return "needs variants and the callable map"
    vcf = ctx.inv.first("vcf")
    if vcf and (vcf.meta or {}).get("build") not in (None, "GRCh37"):
        return "the PGS stage currently scores GRCh37 genomes"
    return None


def _inputs(ctx: Context) -> list[Path]:
    return [variants.variants_parquet(ctx), coverage.out_dir(ctx) / "callable.parquet", qc.out_path(ctx)]


def _params(ctx: Context) -> dict:
    st = Store(ctx.cfg)
    return {"catalog": (st.current("pgs_catalog") or {}).get("version"),
            "reference": (st.current("pgs_reference") or {}).get("version"), "min_overlap": MIN_OVERLAP,
            "min_gq": MIN_GQ, "evidence_model": ev.MODEL_VERSION}


STAGE = Stage(
    name="pgs", version="2", title="Polygenic scores (PGS Catalog + 1000 Genomes)", fn=_pgs, inputs=_inputs,
    outputs=lambda ctx: [out_scores(ctx), out_summary(ctx), out_ancestry(ctx), out_genotype_qc(ctx)],
    available=_available,
    params=_params,
)

