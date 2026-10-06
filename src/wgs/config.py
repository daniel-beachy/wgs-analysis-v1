"""Configuration resolution.

Nothing about file locations is hard-coded. Resolution order (first wins):

1. CLI flags (``--data``, ``--workspace``, ``--config``)
2. Environment variables ``WGS_CONFIG``, ``WGS_DATA``, ``WGS_WORKSPACE``
3. A TOML config file (``WGS_CONFIG`` or ``<repo>/wgs.local.toml``), including
   ``[host."<hostname>"]`` overrides so one file can serve several machines
4. Defaults relative to the repository: search ``<repo>/..`` for inputs and use
   ``<repo>/../wgs-data`` as the workspace

Relative paths in a config file are resolved against the config file's own
directory, so the whole ``Genomics/`` folder can move between machines and
mount points without edits.
"""

from __future__ import annotations

import os
import socket
import tomllib
from dataclasses import dataclass, field
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_NAME = "wgs.local.toml"

INPUT_KEYS = ("vcf", "cram", "reference", "fastq_r1", "fastq_r2", "raw_genotypes", "reports", "truth_vcf",
              "truth_bed")


@dataclass
class Config:
    workspace: Path
    search: list[Path]
    sample: str | None = None
    explicit: dict[str, list[Path]] = field(default_factory=dict)
    external_dirs: list[Path] = field(default_factory=list)
    forbidden_strings: list[str] = field(default_factory=list)
    config_file: Path | None = None
    host: str = ""

    @property
    def cache_dir(self) -> Path:
        return self.workspace / "cache"

    @property
    def work_dir(self) -> Path:
        return self.workspace / "work"

    @property
    def releases_dir(self) -> Path:
        return self.workspace / "releases"

    @property
    def dashboard_dir(self) -> Path:
        return self.workspace / "dashboard"


def _resolve(p: str | os.PathLike, base: Path) -> Path:
    path = Path(os.path.expandvars(os.path.expanduser(str(p))))
    return path if path.is_absolute() else (base / path).resolve()


def _as_list(v) -> list:
    if v is None or v == "":
        return []
    return list(v) if isinstance(v, (list, tuple)) else [v]


def find_config_file(explicit: str | os.PathLike | None = None) -> Path | None:
    for candidate in (explicit, os.environ.get("WGS_CONFIG")):
        if candidate:
            p = Path(candidate).expanduser()
            if not p.is_file():
                raise FileNotFoundError(f"Config file not found: {p}")
            return p.resolve()
    default = REPO_ROOT / DEFAULT_CONFIG_NAME
    return default if default.is_file() else None


def load_config(
    config_path: str | os.PathLike | None = None,
    data: list[str] | None = None,
    workspace: str | None = None,
    sample: str | None = None,
) -> Config:
    host = socket.gethostname().split(".")[0]
    cfg_file = find_config_file(config_path)
    raw: dict = {}
    base = REPO_ROOT
    if cfg_file:
        raw = tomllib.loads(cfg_file.read_text())
        base = cfg_file.parent

    paths = dict(raw.get("paths", {}))
    for name, override in raw.get("host", {}).items():
        if name.lower() == host.lower():
            paths.update(override.get("paths", override))

    env_data = os.environ.get("WGS_DATA")
    search_raw = data or ([s for s in env_data.split(os.pathsep) if s] if env_data else None) or \
        _as_list(paths.get("search")) or [".."]
    search_base = Path.cwd() if (data or env_data) else base
    search = [_resolve(s, search_base) for s in search_raw]

    ws_raw = workspace or os.environ.get("WGS_WORKSPACE") or paths.get("workspace") or "../wgs-data"
    ws_base = Path.cwd() if (workspace or os.environ.get("WGS_WORKSPACE")) else base
    ws = _resolve(ws_raw, ws_base)

    explicit = {}
    for key in INPUT_KEYS:
        vals = _as_list(raw.get("inputs", {}).get(key))
        if vals:
            explicit[key] = [_resolve(v, base) for v in vals]

    external = [_resolve(d, base) for d in _as_list(raw.get("external", {}).get("dirs"))]
    external.append(ws / "external")

    return Config(
        workspace=ws,
        search=search,
        sample=sample or raw.get("sample", {}).get("id") or None,
        explicit=explicit,
        external_dirs=external,
        forbidden_strings=_as_list(raw.get("guard", {}).get("forbidden_strings")),
        config_file=cfg_file,
        host=host,
    )
