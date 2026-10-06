"""Content-based input discovery.

Files are classified by what they *contain* (magic bytes, headers), not by
where they sit, so any folder layout works. Filenames are only used as hints
(e.g. pairing FASTQ R1/R2). Everything here is cheap: we read headers only.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import os
import re
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

from .config import REPO_ROOT, Config

# Contig lengths that identify a reference build.
BUILD_BY_CHR1 = {249250621: "GRCh37", 248956422: "GRCh38", 248387328: "T2T-CHM13"}
SKIP_DIRS = {".git", ".pixi", "node_modules", "__pycache__", ".venv", ".Trashes", ".Spotlight-V100",
             ".fseventsd", "$RECYCLE.BIN", "System Volume Information"}
MAX_DEPTH = 6

# Optional outputs of external tools that we ingest when present (never required).
EXTERNAL_PATTERNS: dict[str, list[str]] = {
    "pharmcat": [r".*\.report\.json$", r".*\.phenotype\.json$", r".*\.match\.json$"],
    "cyrius": [r".*cyrius.*\.tsv$", r".*cyp2d6.*\.tsv$"],
    "t1k": [r".*_genotype\.tsv$"],
    "mosdepth": [r".*\.mosdepth\.summary\.txt$", r".*\.mosdepth\.global\.dist\.txt$"],
    "happy": [r".*\.summary\.csv$", r".*\.extended\.csv$"],
    "vcfeval": [r".*vcfeval.*summary\.txt$", r".*/summary\.txt$"],
    "pgsc_calc": [r".*_pgs\.txt\.gz$", r".*aggregated_scores\.txt\.gz$"],
    "haplogrep": [r".*haplogrep.*\.(txt|tsv)$"],
    "yhaplo": [r".*haplogroups\..*\.txt$"],
    "expansionhunter": [r".*expansionhunter.*\.json$", r".*\.eh\.json$"],
    "verifybamid": [r".*\.selfSM$"],
}


@dataclass
class InputFile:
    kind: str
    path: Path
    size: int
    meta: dict = field(default_factory=dict)
    source: str = "discovered"  # or "config"

    @property
    def name(self) -> str:
        return self.path.name

    def fingerprint(self) -> str:
        """Cheap, stable identity: size + sha256 of first and last MiB."""
        h = hashlib.sha256(str(self.size).encode())
        with open(self.path, "rb") as fh:
            h.update(fh.read(1 << 20))
            if self.size > 2 << 20:
                fh.seek(-(1 << 20), os.SEEK_END)
                h.update(fh.read())
        return h.hexdigest()[:16]


def _read_head(path: Path, n: int = 4096) -> bytes:
    try:
        with open(path, "rb") as fh:
            return fh.read(n)
    except OSError:
        return b""


def _gz_lines(path: Path, limit_bytes: int = 4 << 20):
    """Yield decoded lines from a (b)gzip or plain text file, bounded."""
    raw = _read_head(path, 2)
    opener = gzip.open if raw == b"\x1f\x8b" else open
    read = 0
    with opener(path, "rt", errors="replace") as fh:
        for line in fh:
            read += len(line)
            yield line.rstrip("\n")
            if read > limit_bytes:
                return


def parse_vcf_header(path: Path) -> dict:
    meta: dict = {"contigs": {}, "samples": [], "gvcf": False, "caller": None, "build": None}
    for line in _gz_lines(path, limit_bytes=64 << 20):
        if line.startswith("##contig="):
            m_id = re.search(r"ID=([^,>]+)", line)
            m_len = re.search(r"length=(\d+)", line)
            if m_id and m_len:
                meta["contigs"][m_id.group(1)] = int(m_len.group(1))
        elif line.startswith("##GVCFBlock") or "<NON_REF>" in line or "<*>" in line and "ALT" in line:
            meta["gvcf"] = True
        elif line.startswith("##DeepVariant_version="):
            meta["caller"] = "DeepVariant " + line.split("=", 1)[1]
        elif line.startswith("##source=") and not meta["caller"]:
            meta["caller"] = line.split("=", 1)[1]
        elif line.startswith("##reference="):
            meta["reference"] = line.split("=", 1)[1]
        elif line.startswith("#CHROM"):
            meta["samples"] = line.split("\t")[9:]
            break
    contigs = meta["contigs"]
    chr1 = contigs.get("chr1") or contigs.get("1")
    meta["build"] = BUILD_BY_CHR1.get(chr1)
    meta["chr_prefix"] = "chr1" in contigs
    meta["n_contigs"] = len(contigs)
    meta["has_chrM"] = any(c in contigs for c in ("chrM", "MT", "chrMT"))
    meta["has_chrY"] = any(c in contigs for c in ("chrY", "Y"))
    del meta["contigs"]
    for ext in (".tbi", ".csi"):
        if Path(str(path) + ext).exists():
            meta["index"] = str(path) + ext
    return meta


def parse_alignment_header(path: Path) -> dict:
    meta: dict = {}
    try:
        import pysam

        mode = "rc" if path.suffix == ".cram" else "rb"
        with pysam.AlignmentFile(str(path), mode, check_sq=False) as af:
            h = af.header.to_dict()
        sq = h.get("SQ", [])
        lengths = {s["SN"]: s["LN"] for s in sq}
        chr1 = lengths.get("chr1") or lengths.get("1")
        meta["build"] = BUILD_BY_CHR1.get(chr1)
        meta["n_contigs"] = len(sq)
        meta["has_md5"] = all("M5" in s for s in sq) if sq else False
        rg = (h.get("RG") or [{}])[0]
        meta["samples"] = sorted({r.get("SM", "") for r in h.get("RG", [])})
        meta["platform"] = rg.get("PL") or None
        meta["reference_hint"] = rg.get("DS") or next((s.get("UR") for s in sq if s.get("UR")), None)
        meta["programs"] = sorted({p.get("PN") or p.get("ID", "").split(".")[0] for p in h.get("PG", [])})
    except Exception as exc:  # pragma: no cover - depends on file
        meta["error"] = str(exc)
    for ext in (".crai", ".bai", ".csi"):
        for cand in (Path(str(path) + ext), path.with_suffix(ext)):
            if cand.exists():
                meta["index"] = str(cand)
    return meta


def parse_fastq(path: Path) -> dict:
    meta: dict = {}
    m = re.search(r"(?:_R?([12]))(?:_\d+)?\.f(?:ast)?q(?:\.gz)?$", path.name)
    meta["mate"] = int(m.group(1)) if m else None
    try:
        lines = []
        for line in _gz_lines(path, limit_bytes=8192):
            lines.append(line)
            if len(lines) == 4:
                break
        meta["read_length"] = len(lines[1])
        meta["first_read_id"] = lines[0].split()[0][:60]
    except Exception:
        pass
    return meta


RAW_HEADER = re.compile(r"^(rs\d+|i\d+|[A-Za-z0-9_:.-]+)\t(\w+)\t(\d+)\t([ACGTDI0-]{1,2})$")


def _looks_like_raw_genotypes(lines: list[str]) -> bool:
    data = [ln for ln in lines if ln and not ln.startswith("#") and not ln.lower().startswith("rsid")]
    return len(data) >= 3 and all(RAW_HEADER.match(ln.strip()) for ln in data[:20])


def parse_report_bundle(path: Path) -> dict | None:
    """A zip containing only PDFs (e.g. tellmeGen per-condition report downloads)."""
    try:
        with zipfile.ZipFile(path) as zf:
            names = [i.filename for i in zf.infolist() if not i.is_dir()
                     and not i.filename.split("/")[-1].startswith("._") and not i.filename.startswith("__MACOSX")]
    except Exception:
        return None
    if names and all(n.lower().endswith(".pdf") for n in names):
        return {"zipped": True, "pdfs": len(names)}
    return None


def parse_raw_genotypes(path: Path) -> dict | None:
    """23andMe/tellmeGen-style `rsid chrom pos genotype` (plain or zipped)."""
    try:
        if zipfile.is_zipfile(path):
            with zipfile.ZipFile(path) as zf:
                members = [i for i in zf.infolist() if not i.filename.split("/")[-1].startswith("._")
                           and not i.is_dir() and not i.filename.lower().endswith(".pdf")]
                for info in members:
                    with zf.open(info) as fh:
                        head = io.TextIOWrapper(fh, errors="replace").read(16384).splitlines()[:-1]
                    if _looks_like_raw_genotypes(head):
                        return {"member": info.filename, "uncompressed_size": info.file_size, "zipped": True}
            return None
        if path.suffix.lower() in {".txt", ".tsv", ".csv"} or path.name.endswith(".txt.gz"):
            head = list(_gz_lines(path, 16384))[:-1]
            if _looks_like_raw_genotypes(head):
                return {"zipped": False}
    except Exception:
        return None
    return None


def classify(path: Path) -> InputFile | None:
    name = path.name
    if name.startswith(("._", ".")):
        return None
    try:
        size = path.stat().st_size
    except OSError:
        return None
    lower = name.lower()
    head = _read_head(path, 16)

    if lower.endswith((".tbi", ".csi", ".crai", ".bai", ".fai", ".gzi", ".dict")):
        return InputFile("index", path, size)
    if lower.endswith((".bed", ".bed.gz")):
        return InputFile("bed", path, size)
    if head.startswith(b"CRAM"):
        return InputFile("cram", path, size, parse_alignment_header(path))
    if lower.endswith(".bam"):
        return InputFile("bam", path, size, parse_alignment_header(path))
    if re.search(r"\.(vcf|bcf)(\.gz|\.bgz)?$", lower):
        meta = parse_vcf_header(path)
        kind = "gvcf" if meta.get("gvcf") or ".g.vcf" in lower or ".gvcf" in lower else "vcf"
        return InputFile(kind, path, size, meta)
    if re.search(r"\.f(ast)?q(\.gz)?$", lower):
        return InputFile("fastq", path, size, parse_fastq(path))
    if re.search(r"\.(fa|fasta|fna)(\.gz)?$", lower):
        return InputFile("reference", path, size, {"indexed": Path(str(path) + ".fai").exists()})
    if head.startswith(b"%PDF"):
        provider = "tellmeGen" if "tellmegen" in str(path).lower() else None
        return InputFile("report", path, size, {"provider": provider})
    if head.startswith(b"PK") and (bundle := parse_report_bundle(path)):
        provider = "tellmeGen" if "tellmegen" in str(path).lower() else None
        return InputFile("report", path, size, {"provider": provider, **bundle})
    if head.startswith(b"PK") or lower.endswith((".txt", ".tsv", ".txt.gz")):
        meta = parse_raw_genotypes(path)
        if meta is not None:
            return InputFile("raw_genotypes", path, size, meta)
    if lower == "manifest.sha256" or lower.endswith((".sha256", ".md5")):
        return InputFile("checksums", path, size)
    for tool, patterns in EXTERNAL_PATTERNS.items():
        if any(re.match(p, str(path).replace(os.sep, "/"), re.I) for p in patterns):
            return InputFile(f"external:{tool}", path, size)
    return None


def walk(root: Path, exclude: list[Path], max_depth: int = MAX_DEPTH):
    root = root.resolve()
    if not root.exists():
        return
    if root.is_file():
        yield root
        return
    excl = {e.resolve() for e in exclude}
    base_depth = len(root.parts)
    for dirpath, dirnames, filenames in os.walk(root):
        dp = Path(dirpath)
        if len(dp.parts) - base_depth >= max_depth:
            dirnames[:] = []
        dirnames[:] = sorted(d for d in dirnames if d not in SKIP_DIRS and not d.startswith("._")
                             and (dp / d).resolve() not in excl)
        for f in sorted(filenames):
            yield dp / f


@dataclass
class Inventory:
    files: list[InputFile]
    searched: list[Path]
    sample: str | None
    samples_seen: list[str]
    warnings: list[str] = field(default_factory=list)

    def of(self, kind: str) -> list[InputFile]:
        return [f for f in self.files if f.kind == kind]

    def first(self, kind: str) -> InputFile | None:
        items = self.of(kind)
        return items[0] if items else None

    def external(self) -> dict[str, list[InputFile]]:
        out: dict[str, list[InputFile]] = {}
        for f in self.files:
            if f.kind.startswith("external:"):
                out.setdefault(f.kind.split(":", 1)[1], []).append(f)
        return out

    @property
    def fastq_pair(self) -> tuple[InputFile, InputFile] | None:
        fq = self.of("fastq")
        r1 = [f for f in fq if f.meta.get("mate") == 1]
        r2 = [f for f in fq if f.meta.get("mate") == 2]
        return (r1[0], r2[0]) if r1 and r2 else None


def _belongs(f: InputFile, sample: str) -> bool:
    s = sample.lower()
    if s in f.path.name.lower():
        return True
    return any(s in (x or "").lower() for x in f.meta.get("samples", []))


def discover(cfg: Config) -> Inventory:
    exclude = [REPO_ROOT, cfg.workspace]
    found: dict[Path, InputFile] = {}
    for paths in cfg.explicit.values():
        for p in paths:
            f = classify(p) if p.exists() else None
            if f:
                f.source = "config"
                found[p.resolve()] = f
    for root in cfg.search:
        for p in walk(root, exclude):
            rp = p.resolve()
            if rp not in found:
                f = classify(p)
                if f:
                    found[rp] = f
    for root in cfg.external_dirs:
        for p in walk(root, []):
            rp = p.resolve()
            if rp not in found:
                f = classify(p)
                if f and f.kind.startswith("external:"):
                    found[rp] = f

    files = list(found.values())
    samples = sorted({s for f in files if f.kind in ("vcf", "gvcf") for s in f.meta.get("samples", [])})
    warnings: list[str] = []
    sample = cfg.sample
    if not sample and len(samples) == 1:
        sample = samples[0]
    elif not sample and len(samples) > 1:
        warnings.append(f"Several samples found ({', '.join(samples)}); pass --sample to choose. "
                        f"Using {samples[0]}.")
        sample = samples[0]
    if sample and len(samples) > 1:
        files = [f for f in files if f.kind in ("index", "bed", "checksums", "reference", "report")
                 or f.kind.startswith("external:") or _belongs(f, sample) or f.source == "config"]
    return Inventory(files=files, searched=cfg.search, sample=sample, samples_seen=samples, warnings=warnings)
