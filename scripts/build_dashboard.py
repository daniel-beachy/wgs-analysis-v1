"""Build the Genome Dashboard: viewer (Svelte + DuckDB-WASM) embedded in a tiny Go launcher.

Outputs, placed next to the data so the drive is self-contained:
  <Genomics>/Genome Dashboard.app   macOS, universal (Apple Silicon + Intel)
  <Genomics>/Genome Dashboard.exe   Windows x64

Work happens in an internal scratch directory because exFAT drives cannot hold the symlinks
that npm needs. Usage:  pixi run build-dashboard [--out DIR] [--skip-npm]
"""

from __future__ import annotations

import argparse
import os
import plistlib
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
APP_NAME = "Genome Dashboard"
BUNDLE_ID = "io.github.daniel-beachy.genome-dashboard"


def run(cmd: list[str], cwd: Path, env: dict | None = None) -> None:
    print("  $", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True, env={**os.environ, **(env or {})})


def sync(src: Path, dst: Path, *exclude: str) -> None:
    args = ["rsync", "-a", "--delete", "--exclude", "._*", "--exclude", ".DS_Store"]
    for e in exclude:
        args += ["--exclude", e]
    run([*args, f"{src}/", f"{dst}/"], cwd=REPO)


def version() -> str:
    try:
        out = subprocess.run(["git", "describe", "--always", "--dirty"], cwd=REPO, capture_output=True, text=True)
        return out.stdout.strip() or "dev"
    except OSError:
        return "dev"


def icon_png(scratch: Path) -> Path:
    """Render the app icon (padded favicon) to a 1024 px PNG using macOS sips."""
    svg = (REPO / "viewer/public/favicon.svg").read_text().replace('viewBox="0 0 32 32"', 'viewBox="-3.5 -3.5 39 39"')
    src = scratch / "icon.svg"
    src.write_text(svg)
    png = scratch / "icon.png"
    subprocess.run(
        ["sips", "-s", "format", "png", "-Z", "1024", str(src), "--out", str(png)], check=True, capture_output=True
    )
    return png


def icns(png: Path, scratch: Path) -> Path:
    iconset = scratch / "AppIcon.iconset"
    shutil.rmtree(iconset, ignore_errors=True)
    iconset.mkdir()
    for size in (16, 32, 128, 256, 512):
        for scale in (1, 2):
            px = size * scale
            name = f"icon_{size}x{size}{'@2x' if scale == 2 else ''}.png"
            subprocess.run(
                ["sips", "-z", str(px), str(px), str(png), "--out", str(iconset / name)],
                check=True,
                capture_output=True,
            )
    out = scratch / "AppIcon.icns"
    run(["iconutil", "-c", "icns", str(iconset), "-o", str(out)], cwd=scratch)
    return out


def build(out_dir: Path, skip_npm: bool) -> None:
    scratch = Path(tempfile.gettempdir()) / "wgs-build"
    scratch.mkdir(exist_ok=True)
    viewer, launcher = scratch / "viewer", scratch / "launcher"
    print("1/4 Sync sources to", scratch)
    sync(REPO / "viewer", viewer, "node_modules", "dist")
    sync(REPO / "docs", scratch / "docs")
    sync(REPO / "launcher", launcher, "app", "bin", "*.syso")

    print("2/4 Build the viewer")
    if not skip_npm or not (viewer / "node_modules").exists():
        run(
            ["npm", "ci" if (viewer / "package-lock.json").exists() else "install", "--no-audit", "--no-fund"],
            cwd=viewer,
        )
        shutil.copy2(viewer / "package-lock.json", REPO / "viewer/package-lock.json")
    run(["npx", "svelte-check", "--threshold", "error"], cwd=viewer)
    app_dir = launcher / "app"
    shutil.rmtree(app_dir, ignore_errors=True)
    run(
        ["npx", "vite", "build", "--emptyOutDir", "--logLevel", "warn"],
        cwd=viewer,
        env={"WGS_VIEWER_OUT": str(app_dir)},
    )

    print("3/4 Compile the launcher")
    ver = version()
    ldflags = f"-s -w -X main.version={ver}"
    bins = launcher / "bin"
    shutil.rmtree(bins, ignore_errors=True)
    bins.mkdir()
    png = icon_png(scratch)
    # Windows icon + manifest via a .syso resource (best effort: needs the go-winres module).
    try:
        run(
            [
                "go",
                "run",
                "github.com/tc-hib/go-winres@v0.3.3",
                "simply",
                "--icon",
                str(png),
                "--manifest",
                "gui",
                "--product-name",
                APP_NAME,
                "--file-description",
                APP_NAME,
                "--arch",
                "amd64",
            ],
            cwd=launcher,
        )
    except subprocess.CalledProcessError:
        print("  (no Windows icon: go-winres unavailable)")
    gobuild = ["go", "build", "-trimpath", "-ldflags"]
    run(
        [*gobuild, ldflags + " -H windowsgui", "-o", str(bins / "win-amd64.exe")],
        cwd=launcher,
        env={"GOOS": "windows", "GOARCH": "amd64", "CGO_ENABLED": "0"},
    )
    for f in launcher.glob("*.syso"):
        f.unlink()
    for arch in ("arm64", "amd64"):
        run(
            [*gobuild, ldflags, "-o", str(bins / f"mac-{arch}")],
            cwd=launcher,
            env={"GOOS": "darwin", "GOARCH": arch, "CGO_ENABLED": "0"},
        )
    run(
        ["lipo", "-create", "-output", str(bins / "mac-universal"), str(bins / "mac-arm64"), str(bins / "mac-amd64")],
        cwd=bins,
    )

    print("4/4 Package into", out_dir)
    bundle = scratch / f"{APP_NAME}.app"
    shutil.rmtree(bundle, ignore_errors=True)
    (bundle / "Contents/MacOS").mkdir(parents=True)
    (bundle / "Contents/Resources").mkdir()
    shutil.copy2(bins / "mac-universal", bundle / "Contents/MacOS/genome-dashboard")
    shutil.copy2(icns(png, scratch), bundle / "Contents/Resources/AppIcon.icns")
    with open(bundle / "Contents/Info.plist", "wb") as fh:
        plistlib.dump(
            {
                "CFBundleName": APP_NAME,
                "CFBundleDisplayName": APP_NAME,
                "CFBundleIdentifier": BUNDLE_ID,
                "CFBundleExecutable": "genome-dashboard",
                "CFBundleIconFile": "AppIcon",
                "CFBundlePackageType": "APPL",
                "CFBundleShortVersionString": "0.1.0",
                "CFBundleVersion": ver,
                "LSMinimumSystemVersion": "11.0",
                "LSUIElement": True,  # no Dock icon: it is just a local web server for the browser tab
                "NSHighResolutionCapable": True,
            },
            fh,
        )
    run(
        ["codesign", "--force", "--deep", "--sign", "-", str(bundle)], cwd=scratch
    )  # ad-hoc signature (required on Apple Silicon)

    out_dir.mkdir(parents=True, exist_ok=True)
    dest_app = out_dir / f"{APP_NAME}.app"
    shutil.rmtree(dest_app, ignore_errors=True)
    # exFAT can't store macOS extended attributes; ditto writes them as ._ files, which is fine.
    run(["ditto", "--norsrc", str(bundle), str(dest_app)], cwd=scratch)
    shutil.copy2(bins / "win-amd64.exe", out_dir / f"{APP_NAME}.exe")
    for p in (dest_app, out_dir / f"{APP_NAME}.exe"):
        size = sum(f.stat().st_size for f in p.rglob("*") if f.is_file()) if p.is_dir() else p.stat().st_size
        print(f"  {p}  ({size / 1e6:.1f} MB)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--out",
        type=Path,
        default=REPO.parent,
        help="where to put the app + exe (default: the folder containing this repo)",
    )
    ap.add_argument("--skip-npm", action="store_true", help="reuse installed node_modules")
    a = ap.parse_args()
    if sys.platform != "darwin":
        sys.exit("The build script currently runs on macOS (it packages the .app with lipo/codesign/iconutil).")
    build(a.out.resolve(), a.skip_npm)


if __name__ == "__main__":
    main()
