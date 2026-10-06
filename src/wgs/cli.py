"""Command-line interface: `pixi run wgs <command>`."""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from . import __version__
from .config import DEFAULT_CONFIG_NAME, REPO_ROOT, Config, load_config
from .discovery import Inventory, discover
from .modules import evaluate

app = typer.Typer(add_completion=False, no_args_is_help=True,
                  help="Personal whole-genome analysis. Run `wgs doctor` first.")
console = Console()
STATE: dict = {}


@app.callback()
def main(
    config: str = typer.Option(None, "--config", "-c", help="Path to a wgs TOML config file."),
    data: list[str] = typer.Option(None, "--data", "-d", help="Folder or file to search for inputs (repeatable)."),
    workspace: str = typer.Option(None, "--workspace", "-w", help="Where caches, results and the dashboard go."),
    sample: str = typer.Option(None, "--sample", "-s", help="Sample ID to analyse when several are found."),
):
    STATE["cfg"] = load_config(config, data, workspace, sample)


def _cfg() -> Config:
    return STATE["cfg"]


def human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return str(n)


def _tool_version(cmd: list[str]) -> str | None:
    exe = shutil.which(cmd[0])
    if not exe:
        return None
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
        text = (out.stdout or out.stderr).strip().splitlines()
        return text[0][:60] if text else "present"
    except Exception:
        return "present"


def _fs_type(path: Path) -> str:
    try:
        out = subprocess.run(["df", "-T" if os.uname().sysname == "Linux" else "-k", str(path)],
                             capture_output=True, text=True, timeout=10).stdout
        if os.uname().sysname == "Darwin":
            mount = out.strip().splitlines()[-1].split()[-1]
            info = subprocess.run(["diskutil", "info", mount], capture_output=True, text=True, timeout=10).stdout
            for line in info.splitlines():
                if "File System Personality" in line:
                    return line.split(":", 1)[1].strip()
        return out.strip().splitlines()[-1].split()[1]
    except Exception:
        return "unknown"


def describe(f) -> str:
    m = f.meta
    if f.kind in ("vcf", "gvcf"):
        bits = [m.get("build") or "build?", m.get("caller") or "caller?", f"samples={','.join(m['samples'])}"]
        bits.append("indexed" if m.get("index") else "no index")
        return ", ".join(bits)
    if f.kind in ("cram", "bam"):
        if m.get("error"):
            return "unreadable header: " + m["error"][:60]
        bits = [m.get("build") or "build?", f"SM={','.join(m.get('samples', []))}",
                "M5 tags" if m.get("has_md5") else "no M5 tags (exact reference FASTA required)",
                "indexed" if m.get("index") else "no index"]
        return ", ".join(bits)
    if f.kind == "fastq":
        return f"R{m.get('mate')}, {m.get('read_length')} bp reads"
    if f.kind == "raw_genotypes":
        return ("zip → " + m["member"]) if m.get("zipped") else "plain text"
    if f.kind == "report":
        return m.get("provider") or "PDF"
    if f.kind == "reference":
        return "indexed" if m.get("indexed") else "no .fai"
    return ""


def _short(path: Path, roots: list[Path]) -> str:
    for r in roots:
        try:
            return str(path.relative_to(r))
        except ValueError:
            continue
    return str(path)


@app.command()
def doctor():
    """Show what inputs were found, which sections can run, and environment health."""
    cfg = _cfg()
    console.rule(f"[bold]wgs-analysis {__version__}[/] — doctor")
    console.print(f"Host: [cyan]{cfg.host}[/]   Config: {cfg.config_file or '(none — using defaults)'}")
    console.print("Searching: " + ", ".join(str(p) for p in cfg.search))
    inv: Inventory = discover(cfg)

    t = Table(title=f"Inputs (sample: {inv.sample or '?'})", show_lines=False)
    for col in ("kind", "file", "size", "details"):
        t.add_column(col, overflow="fold")
    order = ["vcf", "gvcf", "cram", "bam", "fastq", "raw_genotypes", "reference", "report", "checksums"]
    for f in sorted(inv.files, key=lambda x: (order.index(x.kind) if x.kind in order else 99, str(x.path))):
        if f.kind in ("index", "bed"):
            continue
        t.add_row(f.kind, _short(f.path, cfg.search), human(f.size), describe(f))
    console.print(t)
    for w in inv.warnings:
        console.print(f"[yellow]⚠ {w}[/]")
    if not inv.files:
        console.print("[red]No inputs found.[/] Use --data PATH or set [paths].search in "
                      f"{DEFAULT_CONFIG_NAME} (run `wgs init`).")

    mt = Table(title="Dashboard sections")
    for col in ("section", "status", "details"):
        mt.add_column(col, overflow="fold")
    colors = {"ready": "green", "partial": "yellow", "unavailable": "red"}
    for row in evaluate(inv, cfg):
        detail = row["detail"] + (f"  [dim]({row['note']})[/]" if row["note"] else "")
        mt.add_row(row["title"], f"[{colors[row['status']]}]{row['status']}[/]", detail)
    console.print(mt)

    tt = Table(title="Tools")
    tt.add_column("tool")
    tt.add_column("version")
    for name, cmd in {"bcftools": ["bcftools", "--version"], "samtools": ["samtools", "--version"],
                      "tabix": ["tabix", "--version"], "java": ["java", "-version"],
                      "duckdb (python)": None}.items():
        if cmd is None:
            import duckdb
            v = duckdb.__version__
        else:
            v = _tool_version(cmd)
        tt.add_row(name, v or "[red]missing[/]")
    console.print(tt)

    ws = cfg.workspace
    ws.mkdir(parents=True, exist_ok=True)
    probe = ws / ".write-test"
    try:
        probe.write_text("ok")
        probe.unlink()
        writable = "[green]writable[/]"
    except OSError as exc:
        writable = f"[red]not writable: {exc}[/]"
    usage = shutil.disk_usage(ws)
    console.print(f"Workspace: {ws}  ({writable}, {human(usage.free)} free, filesystem {_fs_type(ws)})")


@app.command()
def init(force: bool = typer.Option(False, help="Overwrite an existing config.")):
    """Write wgs.local.toml with paths relative to the repo (portable across machines)."""
    target = REPO_ROOT / DEFAULT_CONFIG_NAME
    if target.exists() and not force:
        console.print(f"{target} already exists (use --force).")
        raise typer.Exit(1)
    cfg = _cfg()

    def rel(p: Path) -> str:
        try:
            return os.path.relpath(p, REPO_ROOT).replace(os.sep, "/")
        except ValueError:  # different drive on Windows
            return str(p)

    inv = discover(cfg)
    search = sorted({rel(s) for s in cfg.search})
    lines = [
        "# wgs-analysis local configuration (not committed).",
        "# Relative paths are resolved against this file's folder, so the drive can mount anywhere.",
        "[paths]",
        f"workspace = \"{rel(cfg.workspace)}\"",
        "search = [" + ", ".join(f'"{s}"' for s in search) + "]",
        "",
        "[sample]",
        f"id = \"{inv.sample or ''}\"",
        "",
        "# [inputs] explicit file overrides win over discovery, e.g.",
        "# vcf = \"../WGS/02_VCF/sample.vcf.gz\"",
        "",
        "# [host.\"other-laptop\"]   # per-machine overrides keyed by hostname",
        "# search = [\"D:/Genomics/WGS\"]",
        "",
        "[guard]",
        "# Strings that must never appear in committed files (checked by the pre-commit hook).",
        "forbidden_strings = [" + (f'"{inv.sample}"' if inv.sample else "") + "]",
        "",
    ]
    target.write_text("\n".join(lines))
    console.print(f"Wrote {target}")


@app.command("verify-inputs")
def verify_inputs():
    """Verify input files against MANIFEST.sha256 files found next to them (slow: reads everything)."""
    inv = discover(_cfg())
    manifests = inv.of("checksums")
    if not manifests:
        console.print("No checksum manifest found.")
        raise typer.Exit(1)
    bad = 0
    for man in manifests:
        for line in man.path.read_text().splitlines():
            if not line.strip():
                continue
            digest, rel = line.split(maxsplit=1)
            target = (man.path.parent / rel.strip().lstrip("*")).resolve()
            h = hashlib.sha256()
            with open(target, "rb") as fh:
                for chunk in iter(lambda: fh.read(8 << 20), b""):
                    h.update(chunk)
            ok = h.hexdigest() == digest
            bad += not ok
            console.print(f"{'[green]OK[/]' if ok else '[red]FAILED[/]'}  {target}")
    raise typer.Exit(1 if bad else 0)


@app.command()
def version():
    console.print(__version__)


@app.command("run")
def run_cmd(
    only: list[str] = typer.Option(None, "--only", help="Run only these stages (repeatable)."),
    force: list[str] = typer.Option(None, "--force", help="Recompute these stages ('all' for everything)."),
    threads: int = typer.Option(0, help="Worker threads (default: CPU count - 2)."),
):
    """Run the analysis pipeline. Stages whose inputs and code are unchanged are skipped."""
    from .pipeline import Context, run_stage
    from .stages.registry import STAGES

    cfg = _cfg()
    inv = discover(cfg)
    ctx = Context(cfg=cfg, inv=inv, force=set(force or []))
    if threads:
        ctx.threads = threads
    console.rule(f"wgs run — sample {ctx.sample}")
    names = {s.name for s in STAGES}
    for name in (only or []):
        if name not in names:
            console.print(f"[red]Unknown stage {name}. Known: {', '.join(sorted(names))}[/]")
            raise typer.Exit(2)
    failed = []
    for stage in STAGES:
        if only and stage.name not in only:
            continue
        try:
            run_stage(ctx, stage)
        except Exception as exc:
            import traceback

            failed.append(stage.name)
            with open(ctx.logs / f"{stage.name}.log", "a") as fh:
                fh.write(traceback.format_exc())
            console.print(f"[red]✗ {stage.title} failed:[/] {exc} [dim](traceback in logs/{stage.name}.log)[/]")
    if failed:
        console.print(f"[red]{len(failed)} stage(s) failed: {', '.join(failed)}. Logs: {ctx.logs}[/]")
        raise typer.Exit(1)


knowledge_app = typer.Typer(no_args_is_help=True, help="Public knowledge sources: status and refresh.")
app.add_typer(knowledge_app, name="knowledge")


def _knowledge_scratch() -> Path:
    import tempfile

    base = Path(os.environ.get("WGS_SCRATCH") or tempfile.gettempdir())
    p = base / "wgs-scratch"
    p.mkdir(parents=True, exist_ok=True)
    return p


@knowledge_app.command("status")
def knowledge_status():
    """Show which version of each knowledge source is in use and when it was last checked."""
    from .knowledge import Store
    from .knowledge.registry import SOURCES

    store = Store(_cfg())
    idx = store.index()["sources"]
    t = Table(title="Knowledge sources", show_lines=False)
    for c in ("source", "current version", "versions kept", "last checked", "tables", "cadence"):
        t.add_column(c)
    for src in SOURCES:
        e = idx.get(src.id, {})
        cur = store.current(src.id)
        tables = ", ".join(f"{k} ({v['rows']:,})" for k, v in (cur or {}).get("tables", {}).items())
        t.add_row(src.title, e.get("current", "[yellow]not downloaded[/]"), str(len(e.get("versions", []))),
                  (e.get("checked") or "")[:16], tables, src.cadence)
    console.print(t)


@knowledge_app.command("refresh")
def knowledge_refresh(
    source: list[str] = typer.Option(None, "--source", help="Only these sources (repeatable)."),
    pin: list[str] = typer.Option(None, "--pin", help="source=version, e.g. clinvar=20250106 (repeatable)."),
    force: bool = typer.Option(False, help="Rebuild even if upstream is unchanged."),
    keep: int = typer.Option(4, help="Versions to keep per source."),
):
    """Download new versions of public knowledge. Then `wgs run` re-annotates and records what changed."""
    from .knowledge import refresh

    pins = dict(p.split("=", 1) for p in (pin or []))
    res = refresh(_cfg(), _knowledge_scratch(), only=source or None, pin=pins, force=force, keep=keep)
    if any(v.startswith("failed") for v in res.values()):
        raise typer.Exit(1)


@app.command("query")
def query_cmd(
    sql: str = typer.Argument(None, help="SQL to run; '-' reads it from stdin. Omit to list the tables."),
    release: str = typer.Option(None, "--release", "-r", help="Release ID (default: newest)."),
    fmt: str = typer.Option("table", "--format", "-f", help="table | csv | json"),
    limit: int = typer.Option(200, help="Maximum rows to print."),
):
    """Run read-only SQL against a dashboard release (tables: claims, annotations, variants, ...)."""
    import json
    import sys

    import duckdb

    if sql == "-":
        sql = sys.stdin.read()
    rd = _cfg().releases_dir
    index = json.loads((rd / "index.json").read_text())
    rel = next((r for r in index["releases"] if r["id"] == release), None) if release else index["releases"][0]
    if rel is None:
        raise typer.BadParameter(f"no release {release}")
    m = json.loads((rd / rel["manifest"]).read_text())
    con = duckdb.connect()
    for name, t in m["tables"].items():
        con.execute(f"CREATE VIEW {name} AS SELECT * FROM read_parquet('{(rd / t['path']).as_posix()}')")
    for name, d in m.get("documents", {}).items():
        path = rd / (d["path"] if isinstance(d, dict) else d)
        con.execute(f"CREATE VIEW doc_{name} AS SELECT * FROM read_json_auto('{path.as_posix()}')")
    # Read-only sandbox: only this release folder is readable, nothing can be written.
    con.execute(f"SET allowed_directories = ['{rd.as_posix()}/']")
    con.execute("SET enable_external_access = false")
    con.execute("SET lock_configuration = true")
    if not sql:
        console.print(f"Release [bold]{rel['id']}[/] · sample {m['sample']}")
        for name, t in m["tables"].items():
            console.print(f"  [cyan]{name}[/] ({t['rows']:,} rows) — {t.get('description', '')}")
        console.print("  Documents (JSON) as views: " + ", ".join(f"doc_{d}" for d in m.get("documents", {})))
        return
    stmts = con.extract_statements(sql)
    if len(stmts) != 1 or stmts[0].type != duckdb.StatementType.SELECT:
        raise typer.BadParameter("only a single SELECT (or WITH … SELECT) query is allowed")
    rel_ = con.sql(sql)
    cols, rows = rel_.columns, rel_.fetchmany(limit)
    if fmt == "json":
        print(json.dumps([dict(zip(cols, r, strict=True)) for r in rows], default=str, indent=1))
    elif fmt == "csv":
        import csv

        w = csv.writer(sys.stdout)
        w.writerow(cols)
        w.writerows(rows)
    else:
        t = Table(show_lines=False)
        for c in cols:
            t.add_column(c, overflow="fold")
        for r in rows:
            t.add_row(*("" if v is None else str(v) for v in r))
        console.print(t)
