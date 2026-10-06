"""Reference genome matching the CRAM: find it, or rebuild it from a public assembly and *prove* it matches.

tellmeGen's CRAM was aligned to `GRCh37.primary_assembly.par_y_masked.fa`, which is not distributed. That name and
the contig names (`chr1..22, chrX, chrY, chrM`, `GL000191.1`...) match Ensembl's GRCh37 primary assembly with `chr`
added and the chrY pseudo-autosomal regions replaced by N. We rebuild it, then verify it two independent ways:

1. CRAM 3.0 slices carry an MD5 of the reference span they were encoded against. htslib checks it while decoding,
   so decoding reads from every contig with our FASTA fails loudly if even one base differs.
2. Every REF allele in the VCF (7M+ records, including RefCall sites) must match our FASTA.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import random
import urllib.request
from pathlib import Path

from ..pipeline import Context, Stage, console, run, write_json

ENSEMBL_URL = ("https://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/"
               "Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz")
# GRCh37 chrY pseudo-autosomal regions (1-based, inclusive).
PAR_Y = [(10001, 2649520), (59034050, 59363566)]
OUT_NAME = "GRCh37.primary_assembly.par_y_masked.fa"


def ref_dir(ctx: Context) -> Path:
    p = ctx.cfg.cache_dir / "reference"
    p.mkdir(parents=True, exist_ok=True)
    return p


def ready_file(ctx: Context) -> Path:
    return ref_dir(ctx) / "reference.ready"


def built_fasta(ctx: Context) -> Path:
    return ref_dir(ctx) / OUT_NAME


def cram_contigs(ctx: Context) -> dict[str, int]:
    import pysam

    cram = ctx.inv.first("cram") or ctx.inv.first("bam")
    with pysam.AlignmentFile(str(cram.path), "rc" if cram.kind == "cram" else "rb", check_sq=False) as af:
        return dict(zip(af.references, af.lengths, strict=True))


def fai_contigs(fasta: Path) -> dict[str, int]:
    fai = Path(str(fasta) + ".fai")
    if not fai.exists():
        return {}
    return {ln.split("\t")[0]: int(ln.split("\t")[1]) for ln in fai.read_text().splitlines() if ln}


def reference_fasta(ctx: Context) -> Path | None:
    """The verified reference, if any."""
    rf = ready_file(ctx)
    if rf.exists():
        p = Path(json.loads(rf.read_text())["fasta"])
        if p.exists():
            return p
    return None


def _download(url: str, dest: Path) -> None:
    if dest.exists():
        return
    tmp = dest.with_suffix(dest.suffix + ".part")
    console.print(f"Downloading {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "wgs-analysis"})
    with urllib.request.urlopen(req) as r, open(tmp, "ab") as fh:
        while chunk := r.read(4 << 20):
            fh.write(chunk)
    tmp.replace(dest)


def _rename(name: str) -> str:
    if name == "MT":
        return "chrM"
    if name.isdigit() or name in ("X", "Y"):
        return "chr" + name
    return name


def _build_from_ensembl(ctx: Context, want: dict[str, int]) -> Path:
    src = ref_dir(ctx) / "src" / ENSEMBL_URL.rsplit("/", 1)[1]
    src.parent.mkdir(parents=True, exist_ok=True)
    _download(ENSEMBL_URL, src)
    seqs: dict[str, bytearray] = {}
    name, buf = None, bytearray()
    console.print("Reading Ensembl FASTA and renaming contigs…")
    with gzip.open(src, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name:
                    seqs[name] = buf
                name, buf = _rename(line[1:].split()[0].decode()), bytearray()
            else:
                buf.extend(line.rstrip())
        if name:
            seqs[name] = buf
    missing = [c for c in want if c not in seqs]
    wrong = [c for c in want if c in seqs and len(seqs[c]) != want[c]]
    if missing or wrong:
        raise RuntimeError(f"Ensembl GRCh37 does not match CRAM contigs: missing={missing[:5]} wrong_len={wrong[:5]}")
    y = seqs["chrY"]
    for start, end in PAR_Y:
        y[start - 1:end] = b"N" * (end - start + 1)
    out = built_fasta(ctx)
    tmp = out.with_suffix(".fa.tmp")
    with open(tmp, "wb") as fh:
        for c in want:  # same order as the CRAM header
            seq = seqs[c]
            fh.write(f">{c}\n".encode())
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + b"\n")
    tmp.replace(out)
    Path(str(out) + ".fai").unlink(missing_ok=True)
    run(["samtools", "faidx", str(out)], log=ctx.logs / "reference.log")
    return out


def _verify_cram(ctx: Context, fasta: Path, contigs: dict[str, int]) -> dict:
    """Decode reads from every contig; htslib aborts on any slice MD5 mismatch."""
    cram = ctx.inv.first("cram")
    rng = random.Random(42)
    regions = []
    for c, ln in contigs.items():
        if c.startswith("chr"):
            for _ in range(3):
                s = rng.randint(1, max(1, ln - 200_000))
                regions.append(f"{c}:{s}-{s + 200_000}")
        else:
            regions.append(c)
    log = ctx.logs / "reference.log"
    region_file = ctx.scratch / "verify_regions.txt"
    region_file.write_text("\n".join(regions) + "\n")
    # Full decoding (SAM text) forces sequence reconstruction, which is where htslib checks slice MD5s.
    res = run(["bash", "-o", "pipefail", "-c",
               f"samtools view -@ 4 -T '{fasta}' '{cram.path}' $(cat '{region_file}') | wc -l"],
              capture=True, check=False)
    err = res.stderr or ""
    with open(log, "a") as fh:
        fh.write(f"verify_cram rc={res.returncode}\n{err[-5000:]}\n")
    if res.returncode != 0 or "md5" in err.lower() or "mismatch" in err.lower():
        raise RuntimeError("CRAM decoding with the rebuilt reference failed (reference MD5 mismatch?):\n"
                           + err[-1500:])
    return {"regions_decoded": len(regions), "reads_decoded": int(res.stdout.strip().split()[-1] or 0)}


def _verify_vcf(ctx: Context, fasta: Path) -> dict:
    from .variants import vcf_file

    vcf = vcf_file(ctx)
    if not vcf:
        return {}
    res = run(f"bcftools norm --check-ref w -f '{fasta}' -Ou '{vcf}' 2>&1 >/dev/null | "
              f"grep -c '^REF_MISMATCH' || true", shell=True, capture=True, check=False)
    mismatches = int((res.stdout or "0").strip().splitlines()[-1] or 0)
    return {"vcf_ref_mismatches": mismatches}


def _md5s(fasta: Path) -> dict[str, str]:
    import pysam

    out = {}
    with pysam.FastaFile(str(fasta)) as fa:
        for c in fa.references:
            out[c] = hashlib.md5(fa.fetch(c).upper().encode()).hexdigest()
    return out


def build(ctx: Context) -> dict:
    want = cram_contigs(ctx)
    candidate = None
    for f in ctx.inv.of("reference"):
        if not Path(str(f.path) + ".fai").exists() and not str(f.path).endswith(".gz"):
            run(["samtools", "faidx", str(f.path)], log=ctx.logs / "reference.log", check=False)
        if fai_contigs(f.path) == want:
            candidate, origin = f.path, "provided"
            break
    if candidate is None:
        candidate, origin = _build_from_ensembl(ctx, want), "rebuilt from Ensembl GRCh37 release 75"
    cram_check = _verify_cram(ctx, candidate, want)
    vcf_check = _verify_vcf(ctx, candidate)
    if vcf_check.get("vcf_ref_mismatches"):
        raise RuntimeError(f"{vcf_check['vcf_ref_mismatches']} VCF REF alleles disagree with the reference")
    md5 = _md5s(candidate)
    info = {"fasta": str(candidate), "origin": origin, "contigs": len(want), **cram_check, **vcf_check,
            "md5": md5}
    write_json(ready_file(ctx), info)
    return {k: v for k, v in info.items() if k != "md5"}


def _available(ctx: Context) -> str | None:
    if not (ctx.inv.first("cram") or ctx.inv.first("bam")):
        return "no CRAM/BAM found (the reference is only needed to decode reads)"
    return None


STAGE = Stage(
    name="reference", version="1", title="Reference genome (find / rebuild + verify against CRAM)", fn=build,
    inputs=lambda ctx: [f.path for f in ctx.inv.of("cram")],
    outputs=lambda ctx: [ready_file(ctx)],
    available=_available,
)


def full_decode_file(ctx: Context) -> Path:
    return ref_dir(ctx) / "full_decode.json"


def _full_decode(ctx: Context) -> dict:
    """Decode every read in the CRAM once. htslib checks every slice's reference MD5 while doing so, so passing
    proves the reference is identical wherever any read lies, not just in sampled windows."""
    cram, fasta = ctx.inv.first("cram"), reference_fasta(ctx)
    log = ctx.logs / "reference_full.log"
    # Uncompressed BAM output forces full sequence reconstruction at minimal encoding cost.
    res = run(["bash", "-o", "pipefail", "-c",
               f"samtools view -@ {ctx.threads} -T '{fasta}' -u '{cram.path}' | samtools view -c -@ 2 -"],
              capture=True, check=False)
    err = res.stderr or ""
    with open(log, "a") as fh:
        fh.write(f"full_decode rc={res.returncode}\n{err[-5000:]}\n")
    if res.returncode != 0 or "md5" in err.lower() or "mismatch" in err.lower():
        raise RuntimeError("Full CRAM decode failed — reference differs somewhere:\n" + err[-1500:])
    info = {"reads_decoded": int(res.stdout.strip().split()[-1]), "md5_mismatches": 0}
    write_json(full_decode_file(ctx), info)
    return info


FULL = Stage(
    name="reference_full", version="1", title="Reference: decode every CRAM read (one-time proof)",
    fn=_full_decode, inputs=lambda ctx: [f.path for f in ctx.inv.of("cram")][:1] + [ready_file(ctx)],
    outputs=lambda ctx: [full_decode_file(ctx)],
    available=lambda ctx: None if ctx.inv.first("cram") and reference_fasta(ctx)
    else "needs a CRAM and a verified reference",
)
