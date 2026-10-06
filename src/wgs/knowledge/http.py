"""Downloads via curl (present on macOS, Linux and Windows 10+): resumable, follows redirects, shows progress."""

from __future__ import annotations

import email.utils
import subprocess
import sys
import urllib.request
from pathlib import Path

from . import Upstream

UA = "wgs-analysis (personal genome tool; https://github.com/daniel-beachy)"


def head(url: str) -> Upstream:
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=60) as r:
        lm = r.headers.get("Last-Modified")
        iso = email.utils.parsedate_to_datetime(lm).strftime("%Y-%m-%dT%H:%M:%SZ") if lm else None
        size = r.headers.get("Content-Length")
        return Upstream(url=url, last_modified=iso, etag=(r.headers.get("ETag") or "").strip('"') or None,
                        size=int(size) if size else None)


def listing(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode("utf-8", "replace")


def fetch(url: str, dest: Path, expected_size: int | None = None) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")
    subprocess.run(["curl", "-fL", "--retry", "5", "--retry-delay", "5", "-C", "-", "-A", UA,
                    "--progress-bar" if sys.stderr.isatty() else "--no-progress-meter",
                    "-o", str(part), url], check=True)
    if expected_size and part.stat().st_size != expected_size:
        raise RuntimeError(f"{url}: got {part.stat().st_size} bytes, expected {expected_size}")
    part.rename(dest)
    return dest
