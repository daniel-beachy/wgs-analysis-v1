"""Demo input: Genome in a Bottle HG002's truth set in place of a personal genome (Step 9c, ADR-023 addendum).

The demo release runs the normal pipeline on public data, so the dashboard can be shown to anyone. The input is
GIAB's HG002 v4.2.1 truth VCF (GRCh37, chromosomes 1–22) from the `giab` knowledge source, reshaped so the stages
treat it like a sequencing VCF:

* only records inside GIAB's confident regions are kept, because outside them GIAB doesn't vouch for anything;
* contigs are renamed `1` → `chr1` to match the verified GRCh37 reference the stages already use;
* INFO and the FORMAT fields other than GT and GQ are dropped: DP and AD are summed over a dozen sequencing runs
  (hundreds of reads), which would read as one 30× genome's depth.

The confident regions become the callable map (coverage stage), and call confidence follows the "truth-set
genotype" rule in `evidence.truth_call`. Nothing here reads or writes personal data.
"""

from __future__ import annotations

from pathlib import Path

from .config import Config
from .knowledge import Store
from .pipeline import run

AUTOSOMES = [str(i) for i in range(1, 23)]


def truth_inputs(cfg: Config) -> dict[str, Path | None]:
    """The truth VCF and its confident regions: from [truth] in the config, else the giab knowledge source."""
    st = Store(cfg)
    return {"vcf": cfg.truth.get("vcf") or st.file("giab", "truth.vcf.gz"),
            "regions": cfg.truth.get("regions") or st.file("giab", "confident.bed")}


def input_dir(cfg: Config) -> Path:
    return cfg.workspace / "input"


def demo_vcf(cfg: Config) -> Path:
    return input_dir(cfg) / f"{cfg.sample or 'HG002'}.truth.vcf.gz"


def prepare(cfg: Config) -> Path:
    """Write the demo VCF (and its index) into the demo workspace; a no-op when it is newer than its sources."""
    src = truth_inputs(cfg)
    missing = [k for k, v in src.items() if not v or not Path(v).exists()]
    if missing:
        raise FileNotFoundError(f"Truth set {', '.join(missing)} not found: run `wgs knowledge refresh --source giab` "
                                "with your main config (the demo borrows its knowledge store).")
    out = demo_vcf(cfg)
    newest = max(Path(v).stat().st_mtime for v in src.values())
    if out.exists() and Path(str(out) + ".tbi").exists() and out.stat().st_mtime >= newest:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    names = out.parent / "chr-names.txt"
    names.write_text("".join(f"{c}\tchr{c}\n" for c in AUTOSOMES))
    label = cfg.truth.get("label") or "Genome in a Bottle HG002 truth set"
    header = out.parent / "demo-header.txt"
    header.write_text(f"##source={label} (demo input: confident regions only, contigs renamed to chr, "
                      "INFO and FORMAT other than GT/GQ removed)\n")
    tmp = out.with_suffix(".tmp.gz")
    log = cfg.workspace / "logs" / "demo-prepare.log"
    run(f"bcftools view -T '{src['regions']}' -Ou '{src['vcf']}' | "
        f"bcftools annotate -x INFO,^FORMAT/GT,FORMAT/GQ --rename-chrs '{names}' -h '{header}' -Ou | "
        f"bcftools view -Oz -o '{tmp}'", shell=True, log=log)
    tmp.replace(out)
    run(["bcftools", "index", "-f", "-t", str(out)], log=log)
    for f in (names, header):
        f.unlink()
    return out
