"""IPD-IMGT/HLA: the official catalogue of HLA alleles (the immune-system genes behind several serious drug
reactions, e.g. HLA-B*57:01 and abacavir). Released quarterly; each release adds newly discovered alleles.

Built into a T1K allele index (T1K is the HLA typer, ADR-015). T1K's ``t1k-build.pl`` must be on PATH, which
the pixi environment provides on Apple Silicon; elsewhere this source reports a clear failure and HLA typing
is simply skipped (or an HLA result from another machine is ingested from ``external/``).
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import duckdb

from . import Built, Source, Upstream
from .http import fetch, head, listing

BASE = "https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/"
DAT = BASE + "hla.dat.zip"
INDEX = "hla_dna_seq.fa"


def t1k_tool(name: str) -> str | None:
    """A T1K executable from PATH or the active environment's bin folder."""
    found = shutil.which(name) or shutil.which(name, path=str(Path(sys.prefix) / "bin"))
    return found


class ImgtHla(Source):
    id = "imgt_hla"
    title = "IPD-IMGT/HLA allele catalogue"
    homepage = "https://www.ebi.ac.uk/ipd/imgt/hla/"
    licence = "CC BY-ND 4.0"
    cadence = "quarterly"
    schema = 1
    description = ("Every known allele of the HLA immune genes. Used to type your HLA-A and HLA-B alleles, which "
                   "CPIC guidelines use to flag serious reactions to abacavir, carbamazepine, allopurinol and others.")

    def probe(self, pin: str | None = None) -> list[Upstream]:
        if pin:
            raise ValueError("imgt_hla cannot be pinned (EBI only serves the latest release)")
        version = re.search(r"version:\s*IPD-IMGT/HLA\s+([\d.]+)", listing(BASE + "release_version.txt"))
        up = head(DAT)
        up.etag = version.group(1) if version else up.etag
        return [up]

    def build(self, upstream: list[Upstream], scratch: Path, out: Path) -> Built:
        build_pl = t1k_tool("t1k-build.pl")
        if not build_pl:
            raise RuntimeError("T1K is not installed on this platform (pixi provides it on Apple Silicon only)")
        zp = fetch(DAT, scratch / "hla.dat.zip", upstream[0].size)
        with zipfile.ZipFile(zp) as z:
            z.extract("hla.dat", scratch)
        subprocess.run(["perl", build_pl, "-d", str(scratch / "hla.dat"), "--target", "HLA",
                        "-o", str(scratch / "idx")],
                       check=True, capture_output=True)
        shutil.copyfile(scratch / "idx" / INDEX, out / INDEX)
        names = [ln[1:].split()[0] for ln in (out / INDEX).read_text().splitlines() if ln.startswith(">")]
        con = duckdb.connect()
        con.execute("CREATE TABLE a (allele VARCHAR)")
        con.executemany("INSERT INTO a VALUES (?)", [[n] for n in names])
        dst = out / "alleles.parquet"
        con.execute(f"COPY (SELECT allele, split_part(allele, '*', 1) AS gene FROM a ORDER BY allele) "
                    f"TO '{dst}' (FORMAT parquet)")
        return Built(version=upstream[0].etag or "latest", tables={"alleles": dst}, upstream=upstream,
                     notes={"alleles": len(names)})
