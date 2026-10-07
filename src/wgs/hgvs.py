"""Protein changes in HGVS notation (p.Asn314Asp), the form used by ClinVar, papers and lab reports.

`bcftools csq` writes amino-acid changes as ``314N>314D`` (and whole remaining protein sequences for
frameshifts). This converts them to the standard HGVS protein syntax so a name on the dashboard can be
searched for and matched against the literature. Positions are on whichever transcript was chosen
(MANE Select when available — see ``annotate``).
"""

from __future__ import annotations

import re

THREE = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "E": "Glu", "Q": "Gln", "G": "Gly",
         "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
         "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "Ter", "X": "Xaa", "U": "Sec"}
_RX = re.compile(r"^(\d+)([A-Z*]+)(?:>(\d+)([A-Z*]+))?$")


def _aa(s: str) -> str:
    return "".join(THREE.get(c, c) for c in s)


def protein(raw: str | None, consequence: str | None = None) -> str | None:
    """bcftools csq amino-acid field → HGVS p. notation; None if there is no protein change to name."""
    if not raw:
        return None
    m = _RX.match(raw)
    if not m:
        return None
    pos, ref, _, alt = int(m.group(1)), m.group(2), m.group(3), m.group(4)
    csq = consequence or ""
    if alt is None or ref == alt:
        return f"p.{_aa(ref[0])}{pos}="
    if "frameshift" in csq:
        i = next((k for k in range(min(len(ref), len(alt))) if ref[k] != alt[k]), min(len(ref), len(alt)))
        first = ref[i] if i < len(ref) else "*"
        new = alt[i] if i < len(alt) else ""
        if new == "*":
            return f"p.{_aa(first)}{pos + i}Ter"
        ter = alt.find("*", i)
        tail = f"Ter{ter - i + 1}" if ter >= 0 else "Ter?"
        return f"p.{_aa(first)}{pos + i}{_aa(new)}fs{tail}"
    if len(ref) == 1 and len(alt) == 1:
        if pos == 1 and ref == "M":
            return "p.Met1?"
        if ref == "*":
            return f"p.Ter{pos}{_aa(alt)}ext*?"
        return f"p.{_aa(ref)}{pos}{_aa(alt)}"
    # In-frame insertions/deletions/longer substitutions: trim the shared start and name what changed.
    k = 0
    while k < min(len(ref), len(alt)) and ref[k] == alt[k]:
        k += 1
    r, a, start = ref[k:], alt[k:], pos + k
    if not r and a:
        # The residue after the insertion isn't in the csq field, so write it as a valid delins of the last one.
        return f"p.{_aa(ref[-1])}{start - 1}delins{_aa(ref[-1] + a)}"
    span = f"{_aa(r[0])}{start}" + (f"_{_aa(r[-1])}{start + len(r) - 1}" if len(r) > 1 else "")
    if not a:
        return f"p.{span}del"
    if a.endswith("*") and len(a) == 1:
        return f"p.{_aa(r[0])}{start}Ter"
    return f"p.{span}delins{_aa(a)}"
