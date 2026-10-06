"""Small, dependency-free stage runner.

A stage is skipped when its *stamp* (stage version + input fingerprints +
parameters) matches the previous run and all outputs still exist. This makes
`wgs run` safe to repeat: only stages whose inputs or code changed recompute.
"""

from __future__ import annotations

import json
import os
import shlex
import subprocess
import time
from collections.abc import Callable
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from rich.console import Console

from .config import Config
from .discovery import Inventory

console = Console()

CANONICAL = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]


def canon_chrom(c: str) -> str:
    c = c[3:] if c.lower().startswith("chr") else c
    return "MT" if c in ("M", "MT") else c


def chrom_sort_key(c: str) -> tuple[int, str]:
    c = canon_chrom(c)
    return (CANONICAL.index(c), "") if c in CANONICAL else (99, c)


def now() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


def write_json(path: Path, obj) -> None:
    atomic_write_text(path, json.dumps(obj, indent=2, default=str))


def run(cmd: list[str] | str, *, log: Path | None = None, check: bool = True, capture: bool = False,
        shell: bool = False, env: dict | None = None) -> subprocess.CompletedProcess:
    pretty = cmd if isinstance(cmd, str) else shlex.join(str(c) for c in cmd)
    console.print(f"[dim]$ {pretty}[/]")
    if log:
        log.parent.mkdir(parents=True, exist_ok=True)
        with open(log, "a") as fh:
            fh.write(f"\n[{now()}] $ {pretty}\n")
    full_env = {**os.environ, **(env or {})}
    if capture:
        res = subprocess.run(cmd, shell=shell, capture_output=True, text=True, env=full_env)
    else:
        with open(log, "a") if log else open(os.devnull, "w") as fh:
            res = subprocess.run(cmd, shell=shell, stdout=subprocess.PIPE if capture else fh, stderr=fh, text=True,
                                 env=full_env)
    if check and res.returncode != 0:
        tail = ""
        if log and log.exists():
            tail = "\n".join(log.read_text().splitlines()[-15:])
        raise RuntimeError(f"Command failed ({res.returncode}): {pretty}\n{tail or (res.stderr or '')[-2000:]}")
    return res


@dataclass
class Context:
    cfg: Config
    inv: Inventory
    force: set[str] = field(default_factory=set)
    threads: int = max(2, (os.cpu_count() or 4) - 2)

    @property
    def sample(self) -> str:
        return self.inv.sample or "sample"

    @property
    def work(self) -> Path:
        p = self.cfg.work_dir / self.sample
        p.mkdir(parents=True, exist_ok=True)
        return p

    @property
    def logs(self) -> Path:
        p = self.work / "logs"
        p.mkdir(parents=True, exist_ok=True)
        return p

    @property
    def scratch(self) -> Path:
        """Fast local scratch for large temporary files (not the exFAT drive)."""
        base = Path(os.environ.get("WGS_SCRATCH") or os.environ.get("TMPDIR") or "/tmp")
        p = base / "wgs-scratch" / self.sample
        p.mkdir(parents=True, exist_ok=True)
        return p


@dataclass
class Stage:
    name: str
    version: str
    title: str
    fn: Callable[[Context], dict]
    inputs: Callable[[Context], list[Path]]
    outputs: Callable[[Context], list[Path]]
    available: Callable[[Context], str | None] = lambda ctx: None  # returns a reason when it cannot run
    params: Callable[[Context], dict] = lambda ctx: {}


def _fingerprint(p: Path) -> str:
    st = p.stat()
    return f"{st.st_size}:{int(st.st_mtime)}"


def stamp_path(ctx: Context, stage: Stage) -> Path:
    return ctx.work / "stamps" / f"{stage.name}.json"


def run_stage(ctx: Context, stage: Stage) -> dict:
    reason = stage.available(ctx)
    if reason:
        console.print(f"[yellow]↷ {stage.title}: skipped — {reason}[/]")
        return {"status": "skipped", "reason": reason}
    inputs = [p for p in stage.inputs(ctx) if p]
    stamp = {
        "stage": stage.name,
        "version": stage.version,
        "inputs": {str(p): _fingerprint(p) for p in inputs},
        "params": stage.params(ctx),
    }
    sp = stamp_path(ctx, stage)
    if sp.exists() and stage.name not in ctx.force and "all" not in ctx.force:
        old = json.loads(sp.read_text())
        same = {k: old.get(k) for k in stamp} == stamp
        if same and all(p.exists() for p in stage.outputs(ctx)):
            console.print(f"[green]✓ {stage.title}[/] [dim](up to date)[/]")
            return {"status": "cached", **old.get("result", {})}
    console.print(f"[bold cyan]▶ {stage.title}[/]")
    t0 = time.time()
    result = stage.fn(ctx) or {}
    stamp.update(result=result, finished=now(), seconds=round(time.time() - t0, 1))
    write_json(sp, stamp)
    console.print(f"[green]✓ {stage.title}[/] [dim]({stamp['seconds']} s)[/]")
    return {"status": "ran", **result}
