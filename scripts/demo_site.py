"""Publish the HG002 demo dashboard on GitHub Pages (ADR-028).

  pixi run demo-site                export the demo release, upload it, start the Pages deploy
  pixi run demo-site --no-upload    export and pack only, to look at the files first

The repository holds code only, so the demo data travels as the asset `demo-site.tar.gz` of a GitHub release tagged
`demo-data-<release id>`. `.github/workflows/pages.yml` builds the viewer and adds the newest such asset. Run
`pixi run demo` first when the demo release itself needs rebuilding.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ASSET = "demo-site.tar.gz"


def gh(*args: str, check: bool = True) -> subprocess.CompletedProcess:
    print("  $ gh", " ".join(args), flush=True)
    return subprocess.run(["gh", *args], cwd=REPO, check=check)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-upload", action="store_true", help="export and pack only")
    a = ap.parse_args()

    out = Path(os.environ.get("TMPDIR", tempfile.gettempdir())) / "wgs-build" / "demo-site"
    subprocess.run([sys.executable, "-m", "wgs", "--config", str(REPO / "config" / "demo.toml"), "export-site",
                    "--out", str(out)], cwd=REPO, check=True)
    release = json.loads((out / "data" / "index.json").read_text())["latest"]
    tarball = out / ASSET
    with tarfile.open(tarball, "w:gz") as t:
        t.add(out / "data", arcname="data", filter=lambda i: None if Path(i.name).name.startswith("._") else i)
    print(f"Packed {tarball} ({tarball.stat().st_size / 1e6:.1f} MB)")
    if a.no_upload:
        return

    tag = f"demo-data-{release}"
    if subprocess.run(["gh", "release", "view", tag], cwd=REPO, capture_output=True).returncode == 0:
        gh("release", "upload", tag, str(tarball), "--clobber")
    else:
        gh("release", "create", tag, str(tarball), "--latest=false", "--title", f"Demo data {release}", "--notes",
           "Data for the public HG002 demo dashboard on GitHub Pages (Genome in a Bottle's public reference genome, "
           "not a personal genome). Not a software release; see ADR-028.")
    if gh("workflow", "run", "pages.yml", check=False).returncode:
        print("Couldn't start the deploy: push .github/workflows/pages.yml to main first (pushing it deploys).")


if __name__ == "__main__":
    main()
