"""Modes of inheritance: one small vocabulary shared by every knowledge source and the evidence model.

How a condition is inherited decides what a variant means for *you*: one copy of a recessive
variant usually makes you a carrier (healthy, but able to pass it on), whereas one copy of a
dominant variant can be enough to cause the condition. HPO, Orphanet, ClinGen and ACMG each
spell these differently; they are normalised to the short codes below.
"""

from __future__ import annotations

import re

LABEL = {
    "AD": "Autosomal dominant",
    "AR": "Autosomal recessive",
    "SD": "Semidominant",
    "XLR": "X-linked recessive",
    "XLD": "X-linked dominant",
    "XL": "X-linked",
    "YL": "Y-linked",
    "MT": "Mitochondrial",
    "MF": "Multifactorial / polygenic",
    "DG": "Digenic / oligogenic",
    "SP": "Sporadic / somatic",
}

# HPO "Mode of inheritance" terms (HP:0000005 subtree) seen in phenotype.hpoa
HPO = {
    "HP:0000006": "AD", "HP:0000007": "AR", "HP:0032113": "SD", "HP:0001417": "XL", "HP:0001419": "XLR",
    "HP:0001423": "XLD", "HP:0001450": "YL", "HP:0001427": "MT", "HP:0001426": "MF", "HP:0010982": "MF",
    "HP:0010984": "DG", "HP:0010983": "DG", "HP:0003745": "SP", "HP:0001428": "SP", "HP:0001442": "SP",
    "HP:0012274": "AD", "HP:0012275": "AD", "HP:0025352": "AD",
}

# Orphanet TypeOfInheritance names (and ClinGen / ACMG spellings)
_TEXT = [
    (r"semi-?dominant", "SD"),
    (r"autosomal dominant|^ad$", "AD"),
    (r"autosomal recessive|^ar$", "AR"),
    (r"x-linked recessive|^xlr$", "XLR"),
    (r"x-linked dominant|^xld$", "XLD"),
    (r"x-linked|^xl$", "XL"),
    (r"y-linked", "YL"),
    (r"mitochondrial", "MT"),
    (r"multigenic|multifactorial|polygenic", "MF"),
    (r"digenic|oligogenic", "DG"),
    (r"sporadic|somatic", "SP"),
]


def from_text(s: str | None) -> str | None:
    """'Autosomal recessive' → 'AR' (None when not a recognised mode, e.g. 'Unknown')."""
    t = (s or "").strip().lower()
    for pat, code in _TEXT:
        if re.search(pat, t):
            return code
    return None


def cues(condition_names: str | None) -> set[str]:
    """Inheritance spelled out in condition names, e.g. 'Autosomal recessive nonsyndromic hearing loss 1A'."""
    t = (condition_names or "").lower()
    out = {code for phrase, code in (("autosomal dominant", "AD"), ("autosomal recessive", "AR"),
                                     ("x-linked recessive", "XLR"), ("x-linked dominant", "XLD"))
           if phrase in t}
    if "x-linked" in t and not out & X_LINKED:
        out.add("XL")
    return out


RECESSIVE = {"AR"}
DOMINANT = {"AD", "SD", "XLD"}
X_LINKED = {"XL", "XLR", "XLD"}


def describe(modes: set[str] | list[str] | None) -> str:
    return " / ".join(LABEL.get(m, m) for m in sorted(modes or [])) or "unknown"
