"""PharmCAT: turns genotypes at known pharmacogene positions into star-allele diplotypes, metaboliser
phenotypes and the matching CPIC / DPWG / FDA prescribing guidance.

PharmCAT is software *and* knowledge: each release bundles the allele definitions (PharmVar / CPIC),
phenotype mappings and guideline recommendations current at the time (curated by ClinPGx, the merged
PharmGKB + CPIC + PharmVar resource). So refreshing this source is how the Medicines section learns about
new guidelines, new star alleles and changed recommendations — and a new version shows up in
"What changed" like any other knowledge update.

The release ships positions on GRCh38. Your data may be on another build (tellmeGen uses GRCh37), so the
build also lifts every position to GRCh37 with UCSC's hg38ToHg19 chain. Only the *positions* are lifted;
the genotypes are then read directly from your VCF and checked against your reference (ADR-015).
"""

from __future__ import annotations

import gzip
import json
import urllib.request
from pathlib import Path

import duckdb
import numpy as np

from . import Built, Source, Upstream
from .http import UA, fetch, head

RELEASES = "https://api.github.com/repos/PharmGKB/PharmCAT/releases"
CHAIN = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
JAR = "pharmcat.jar"
POSITIONS_VCF = "pharmcat_positions.vcf"


def _release(pin: str | None) -> dict:
    url = f"{RELEASES}/tags/v{pin.lstrip('v')}" if pin else f"{RELEASES}/latest"
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/vnd.github+json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read())


def _assets(rel: dict) -> tuple[dict, dict]:
    jar = next(a for a in rel["assets"] if a["name"].endswith("-all.jar"))
    pos = next(a for a in rel["assets"] if a["name"].startswith("pharmcat_positions") and a["name"].endswith(".vcf"))
    return jar, pos


class Chain:
    """Minimal UCSC chain-file lifter for single positions (enough for ~1,200 pharmacogene sites)."""

    def __init__(self, path: Path, chroms: set[str] | None = None):
        blocks: dict[str, list[tuple[int, int, int, str, str, int, int]]] = {}
        rank = 0
        with gzip.open(path, "rt") as fh:
            keep = False
            for line in fh:
                f = line.split()
                if not f:
                    continue
                if f[0] == "chain":
                    rank += 1
                    t_name, t_start, q_name, q_size, q_strand, q_start = f[2], int(f[5]), f[7], int(f[8]), f[9], \
                        int(f[10])
                    keep = chroms is None or t_name in chroms
                    t, q = t_start, q_start
                    continue
                if not keep:
                    continue
                size = int(f[0])
                blocks.setdefault(t_name, []).append((t, t + size, q, q_name, q_strand, q_size, rank))
                if len(f) == 3:
                    t += size + int(f[1])
                    q += size + int(f[2])
        self._blocks = {}
        for c, bl in blocks.items():
            arr = np.array([(b[0], b[1], b[2], b[5], b[6]) for b in bl], dtype=np.int64)
            self._blocks[c] = (arr, [b[3] for b in bl], [b[4] for b in bl])

    def lift(self, chrom: str, pos: int, length: int = 1) -> tuple[str, int, str] | None:
        """1-based pos (and the span of `length` bases) → (chrom, 1-based pos, strand) or None if unmapped/split."""
        if chrom not in self._blocks:
            return None
        arr, qnames, strands = self._blocks[chrom]
        x0, x1 = pos - 1, pos - 1 + length
        hit = np.nonzero((arr[:, 0] <= x0) & (arr[:, 1] >= x1))[0]
        if not len(hit):
            return None
        i = hit[np.argmin(arr[hit, 4])]  # highest-scoring chain (file order)
        t0, _, q0, q_size, _ = arr[i]
        q = int(q0 + (x0 - t0))
        if strands[i] == "-":
            return qnames[i], int(q_size - q - length + 1), "-"
        return qnames[i], q + 1, "+"


def read_positions(vcf: Path) -> list[dict]:
    rows = []
    for line in vcf.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        rows.append({"chrom38": f[0], "pos38": int(f[1]), "rsid": f[2] if f[2] != "." else None, "ref": f[3],
                     "alts": f[4].split(",") if f[4] != "." else [], "gene": info.get("PX")})
    return rows


class PharmCAT(Source):
    id = "pharmcat"
    title = "PharmCAT + CPIC / DPWG / FDA guidance (ClinPGx)"
    homepage = "https://pharmcat.clinpgx.org/"
    licence = "MPL-2.0 (software); CC BY-SA 4.0 (ClinPGx data)"
    cadence = "a few releases a year"
    schema = 1
    description = ("Pharmacogenomics Clinical Annotation Tool: star-allele definitions, metaboliser phenotypes and "
                   "the prescribing guidelines of CPIC (US) and DPWG (Netherlands) plus FDA label annotations, "
                   "curated by ClinPGx (PharmGKB + CPIC + PharmVar).")

    _rel: dict | None = None

    def probe(self, pin: str | None = None) -> list[Upstream]:
        self._rel = _release(pin)
        jar, pos = _assets(self._rel)
        tag = self._rel["tag_name"]
        return [Upstream(url=jar["browser_download_url"], last_modified=jar["updated_at"], etag=tag, size=jar["size"]),
                Upstream(url=pos["browser_download_url"], last_modified=pos["updated_at"], etag=tag,
                         size=pos["size"]),
                head(CHAIN)]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        jar_up, pos_up, chain_up = upstream
        version = (jar_up.etag or "").lstrip("v")
        fetch(jar_up.url, out / JAR, jar_up.size)
        fetch(pos_up.url, out / POSITIONS_VCF, pos_up.size)
        chain_path = fetch(CHAIN, scratch / "hg38ToHg19.over.chain.gz", chain_up.size)
        rows = read_positions(out / POSITIONS_VCF)
        chain = Chain(chain_path, {r["chrom38"] for r in rows})
        for r in rows:
            hit = chain.lift(r["chrom38"], r["pos38"], len(r["ref"]))
            r["chrom37"], r["pos37"], r["strand37"] = hit if hit else (None, None, None)
        con = duckdb.connect()
        con.execute("CREATE TABLE p (chrom38 VARCHAR, pos38 INTEGER, rsid VARCHAR, ref VARCHAR, alts VARCHAR[], "
                    "gene VARCHAR, chrom37 VARCHAR, pos37 INTEGER, strand37 VARCHAR)")
        con.executemany("INSERT INTO p VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                        [[r[k] for k in ("chrom38", "pos38", "rsid", "ref", "alts", "gene", "chrom37", "pos37",
                                         "strand37")] for r in rows])
        dst = out / "positions.parquet"
        con.execute(f"COPY (SELECT * FROM p ORDER BY chrom38, pos38) TO '{dst}' (FORMAT parquet)")
        lifted = sum(1 for r in rows if r["chrom37"])
        return Built(version=version, tables={"positions": dst}, upstream=upstream,
                     notes={"positions": len(rows), "lifted_to_grch37": lifted,
                            "genes": sorted({r["gene"] for r in rows if r["gene"]}),
                            "release_notes": (self._rel or {}).get("html_url")})
