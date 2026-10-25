"""High-depth reference genomes for finer ancestry (ADR-020 addendum).

The continental shares use the 1000 Genomes phase 3 panel, whose 26 populations include only five European, five
East Asian, ... groups. The Human Genome Diversity Project (HGDP, 929 genomes from 54 populations; Bergström et al.
2020) adds groups such as Basque, Sardinian, Mbuti, Russian, Bedouin, Uygur and Papuan.

Mixing panels is only fair when they were read the same way. The phase 3 genotypes come from low-depth sequencing;
HGDP's from 30× genomes. In a test on European groups that difference alone became a PCA axis: HGDP's Tuscans sat
apart from 1000 Genomes' Tuscans, and a 30× genome (like yours) was pulled towards the HGDP groups. So this source
pairs HGDP with the 30× re-sequencing of the same 1000 Genomes people (Byrska-Bishop et al. 2022). With both at 30×,
the two Tuscan samples are no longer told apart.

Only the positions the ancestry PCA uses are kept (the 1000 Genomes panel thinned by the same rule as the polygenic
score stage, about 59,000), so the build downloads about 6 GB one chromosome at a time and stores about 40 MB.
Files are the plink 2 conversions published on the plink 2 resources page; both releases are on GRCh38, and
positions are lifted from GRCh37 with UCSC's hg19ToHg38 chain.
"""

from __future__ import annotations

import hashlib
import re
import shutil
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from ..pipeline import run
from . import Built, Source, Upstream, http
from .pharmcat import Chain

PAGE = "https://www.cog-genomics.org/plink/2.0/resources"
CHAIN = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
# File names on the resources page. The 30x 1000 Genomes files appear twice (2022 release first, then the 2021
# preliminary one), so the first link for each name is used.
PANELS = {
    "1kg30x": {"pgen": "chr{c}_hg38.pgen.zst", "pvar": "chr{c}_hg38_noannot.pvar.zst", "psam": "hg38_corrected.psam",
               "label": "1000 Genomes (30×, 2022)"},
    "hgdp": {"pgen": "hgdp_chr{c}.pgen.zst", "pvar": "hgdp_chr{c}_noannot.pvar.zst", "psam": "hgdp.psam",
             "label": "HGDP (30×, 2019)"},
}
# HGDP's regions mapped to the 1000 Genomes continental groups they are compared within. The Middle East and
# Oceania have no 1000 Genomes counterpart, so they are kept but never used as a comparison group.
HGDP_REGION = {"AFRICA": "AFR", "EUROPE": "EUR", "EAST_ASIA": "EAS", "CENTRAL_SOUTH_ASIA": "SAS", "AMERICA": "AMR"}
COMPLEMENT = str.maketrans("ACGT", "TGCA")
MISSING = 3


def _links() -> dict[str, str]:
    html = http.listing(PAGE)
    out: dict[str, str] = {}
    for m in re.finditer(r'href="(https://www\.dropbox\.com/[^"]+?/([\w.]+?\.(?:pgen\.zst|pvar\.zst|psam))\?[^"]*)"',
                         html):
        out.setdefault(m.group(2), m.group(1).replace("&amp;", "&"))
    return out


def _wanted() -> list[str]:
    names = []
    for p in PANELS.values():
        names.append(p["psam"])
        names += [p[k].format(c=c) for c in range(1, 23) for k in ("pgen", "pvar")]
    return names


def read_bed(stem: Path) -> tuple[list[list[str]], list[str], np.ndarray]:
    """PLINK 1 fileset → (bim rows, sample IDs, copies of each site's first allele; MISSING where not called)."""
    bim = [line.split() for line in Path(f"{stem}.bim").read_text().splitlines()]
    fam = [line.split()[1] for line in Path(f"{stem}.fam").read_text().splitlines()]
    n = len(fam)
    raw = np.fromfile(f"{stem}.bed", dtype=np.uint8)[3:].reshape(len(bim), (n + 3) // 4)
    codes = np.stack([(raw >> s) & 3 for s in (0, 2, 4, 6)], axis=2).reshape(len(bim), -1)[:, :n]
    # PLINK 1 codes: 00 = two copies of A1, 01 = missing, 10 = one copy, 11 = none
    return bim, fam, np.choose(codes, [2, MISSING, 1, 0]).astype(np.uint8)


def align(site: tuple[str, str, str], a1: str, a2: str) -> int:
    """How a panel's alleles (a1, a2 on GRCh38) relate to a GRCh37 site (ref, alt, strand of the lift):
    +1 = a1 is the alt allele, -1 = a1 is the ref allele, 0 = no match or strand-ambiguous (A/T, C/G)."""
    ref, alt, strand = site
    if {ref, alt} in ({"A", "T"}, {"C", "G"}):
        return 0
    if strand == "-":
        a1, a2 = a1.translate(COMPLEMENT), a2.translate(COMPLEMENT)
    return 1 if (a1, a2) == (alt, ref) else -1 if (a1, a2) == (ref, alt) else 0


def alt_copies(g: np.ndarray, sign: int) -> np.ndarray:
    """Copies of the GRCh37 alt allele, from copies of the panel's a1 (MISSING stays MISSING)."""
    return g if sign == 1 else np.where(g == MISSING, MISSING, 2 - g).astype(np.uint8)


class AncestryPanels(Source):
    id = "ancestry_panels"
    title = "High-depth reference genomes (HGDP + 1000 Genomes 30×)"
    homepage = PAGE
    licence = ("HGDP and 1000 Genomes data: open, free to use (Fort Lauderdale principles); plink 2 conversions by "
               "Christopher Chang")
    cadence = "never (2019 and 2022 releases)"
    description = ("929 genomes from 54 Human Genome Diversity Project populations plus the 1000 Genomes people, all "
                   "sequenced at 30×, kept at the ancestry positions only. Used to tell finer reference populations "
                   "apart inside your main continental group.")

    def __init__(self) -> None:
        self._links: dict[str, str] | None = None

    def probe(self, pin: str | None = None) -> list[Upstream]:
        self._links = _links()
        missing = [n for n in _wanted() if n not in self._links]
        if missing:
            raise RuntimeError(f"plink 2 resources page no longer lists {', '.join(missing[:3])}")
        digest = hashlib.sha256("\n".join(self._links[n] for n in _wanted()).encode()).hexdigest()[:16]
        ref = (self.store.current("pgs_reference") or {}).get("version") if self.store else None
        return [Upstream(url=PAGE, etag=f"{digest}|pgs_reference={ref}", size=len(_wanted()))]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        from ..stages import pgs  # the stage module owns the PCA position rule

        links = self._links or _links()
        if not self.store or not self.store.current("pgs_reference"):
            raise RuntimeError("needs the pgs_reference source first (wgs knowledge refresh --source pgs_reference)")
        plink = pgs._tool("plink2")
        if not plink:
            raise RuntimeError("plink2 is not installed (pixi install -e pgs)")
        panel = pgs.panel_files(self.store)

        def plink_run(*args: str, log: str) -> None:
            run([plink, "--memory", str(pgs._plink_memory_mb()), *args], log=scratch / f"{log}.log")

        pgs.reference_thin(panel, scratch, plink_run)
        ids = (scratch / "thin.prune.in").read_text().split()
        for f in scratch.glob("qc.*"):
            f.unlink()
        chain = Chain(http.fetch(CHAIN, scratch / "hg19ToHg38.over.chain.gz"), {f"chr{c}" for c in range(1, 23)})
        sites: dict[str, tuple[str, str, str]] = {}   # "chr:pos38" -> (GRCh37 id, ...)
        lifted: dict[str, list[tuple[int, str, str, str, str]]] = {}
        for sid in ids:
            c, pos, ref, alt = sid.split(":")
            hit = chain.lift(f"chr{c}", int(pos))
            if hit and hit[0] == f"chr{c}":
                lifted.setdefault(c, []).append((hit[1], sid, ref, alt, hit[2]))
        geno: dict[str, dict[str, np.ndarray]] = {k: {} for k in PANELS}
        samples: dict[str, list[str]] = {}
        for key, p in PANELS.items():
            psam = http.fetch(links[p["psam"]], scratch / f"{key}.psam")
            for c in map(str, range(1, 23)):
                bed1 = scratch / f"r{c}.bed1"
                bed1.write_text("".join(f"{c}\t{x[0]}\t{x[0]}\n" for x in sorted(lifted.get(c, []))))
                zst = http.fetch(links[p["pgen"].format(c=c)], scratch / f"{key}{c}.pgen.zst")
                pvar = http.fetch(links[p["pvar"].format(c=c)], scratch / f"{key}{c}.pvar.zst")
                run([plink, "--zst-decompress", str(zst), str(scratch / f"{key}{c}.pgen")])
                zst.unlink()
                shutil.copy(psam, scratch / f"{key}{c}.psam")
                plink_run("--pfile", str(scratch / f"{key}{c}"), "vzs", "--extract", "bed1", str(bed1),
                          "--snps-only", "just-acgt", "--max-alleles", "2", "--make-bed",
                          "--out", str(scratch / f"{key}-sub{c}"), log=f"{key}-extract")
                for f in (scratch / f"{key}{c}.pgen", pvar):
                    f.unlink()
                bim, fam, g = read_bed(scratch / f"{key}-sub{c}")
                samples[key] = fam
                by_pos = {x[0]: x for x in lifted.get(c, [])}
                seen: dict[int, int] = {}
                for row in bim:
                    seen[int(row[3])] = seen.get(int(row[3]), 0) + 1
                for i, row in enumerate(bim):
                    pos38 = int(row[3])
                    x = by_pos.get(pos38)
                    if not x or seen[pos38] > 1:
                        continue
                    sign = align((x[2], x[3], x[4]), row[4], row[5])
                    if sign:
                        geno[key][x[1]] = alt_copies(g[i], sign)
                        sites[x[1]] = (c, x[0], x[4])
        keep = [s for s in ids if all(s in geno[k] for k in PANELS)]
        if len(keep) < 20_000:
            raise RuntimeError(f"only {len(keep)} positions matched in both panels (expected about 45,000)")
        order = [(k, iid) for k in PANELS for iid in samples[k]]
        mat = np.vstack([np.concatenate([geno[k][s] for k in PANELS]) for s in keep])
        pq.write_table(pa.table({
            "id": keep, "chrom": [s.split(":")[0] for s in keep], "pos": [int(s.split(":")[1]) for s in keep],
            "pos38": [sites[s][1] for s in keep], "strand": [sites[s][2] for s in keep],
            "genotypes": pa.array([r.tobytes() for r in mat], pa.binary())}),
            out / "genotypes.parquet", compression="zstd")
        phase3 = {ln.split("\t")[0] for ln in panel["psam"].read_text().splitlines()[1:]}
        meta = {**_psam(scratch / "1kg30x.psam", "1kg30x", phase3), **_psam(scratch / "hgdp.psam", "hgdp", phase3)}
        pq.write_table(pa.Table.from_pylist([{"col": i, "panel": k, "panel_label": PANELS[k]["label"], "iid": iid,
                                              **meta[(k, iid)]} for i, (k, iid) in enumerate(order)]),
                       out / "samples.parquet")
        for f in scratch.glob("*-sub*"):
            f.unlink()
        return Built(version=f"2022-2019-{panel['version']}", upstream=upstream,
                     tables={"genotypes": out / "genotypes.parquet", "samples": out / "samples.parquet"},
                     notes={"positions": len(keep), "thinned_positions": len(ids),
                            "citation": "Bergström A, et al. Insights into human genetic variation and population "
                                        "history from 929 diverse genomes. Science. 2020;367:eaay5012. "
                                        "Byrska-Bishop M, et al. High-coverage whole-genome sequencing of the "
                                        "expanded 1000 Genomes Project cohort including 602 trios. Cell. "
                                        "2022;185:3426-3440."})


def _psam(path: Path, panel: str, phase3: set[str]) -> dict[tuple[str, str], dict]:
    lines = [ln.split("\t") for ln in path.read_text().splitlines()]
    head = [h.lstrip("#") for h in lines[0]]
    out = {}
    for f in lines[1:]:
        r = dict(zip(head, f, strict=False))
        if panel == "hgdp":
            out[(panel, r["IID"])] = {"population": r["population"], "superpop": HGDP_REGION.get(r["region"]),
                                      "region": r["region"], "related": False}
        else:
            # The 2022 release adds 698 relatives of the phase 3 people. Only the 2,504 phase 3 people are used
            # (relatives would find each other as nearest neighbours)
            out[(panel, r["IID"])] = {"population": r["Population"], "superpop": r["SuperPop"],
                                      "region": r["SuperPop"], "related": r["IID"] not in phase3}
    return out

