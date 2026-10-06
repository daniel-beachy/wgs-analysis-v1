"""Interpretation context: everything the curated sources say about genes and conditions, in one place.

Loaded once per annotate run from the current knowledge tables (each optional). Answers:

* Which conditions does this ClinVar record name, what are they, how are they inherited, how rare?
  (ClinVar condition IDs → Mondo → OMIM / Orphanet → HPO and Orphanet inheritance.)
* Failing that, how are this gene's diseases inherited? (ClinGen, HPO, Orphanet, ACMG.)
* Is the gene on the ACMG secondary-findings list, and what does ClinGen say about actionability?
* How tolerant is the gene of being broken? (gnomAD constraint.)
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import duckdb

from . import inheritance as inh

FILLER = {"not provided", "not specified", "see cases", "other", "none provided"}


@dataclass
class Condition:
    name: str
    ids: list[str] = field(default_factory=list)          # OMIM:… / ORPHA:… / MONDO:…
    modes: list[str] = field(default_factory=list)
    modes_from: str | None = None                          # 'HPO', 'Orphanet', 'name'
    definition: str | None = None
    onset: list[str] = field(default_factory=list)
    prevalence: str | None = None

    def as_dict(self) -> dict:
        return {"name": self.name, "ids": self.ids, "modes": self.modes, "modes_from": self.modes_from,
                "definition": self.definition, "onset": self.onset, "prevalence": self.prevalence}


class Kb:
    def __init__(self, paths: dict[str, Path | None]):
        con = duckdb.connect()

        def rows(key: str, sql: str) -> list[tuple]:
            p = paths.get(key)
            return con.execute(sql.format(p=p)).fetchall() if p else []

        self.hpo_modes = {d: (n, list(m or [])) for d, n, m in
                          rows("hpo_inh", "SELECT disease_id, disease_name, modes FROM '{p}'")}
        self.orpha = {d: (n, list(m or []), list(o or []), pv) for d, n, m, o, pv in
                      rows("orpha_dis", "SELECT disease_id, name, modes, onset, prevalence FROM '{p}'")}
        self.mondo: dict[str, tuple[str, str | None, list[str]]] = {}
        self.xref_to_mondo: dict[str, str] = {}
        for mid, name, definition, xrefs in rows("mondo", "SELECT mondo_id, name, definition, xrefs FROM '{p}'"):
            self.mondo[mid] = (name, definition, list(xrefs or []))
            for x in xrefs or []:
                self.xref_to_mondo.setdefault(x, mid)
        # gene → modes, from every source that says how the gene's diseases are inherited
        self.gene_modes: dict[str, set[str]] = {}
        self.gene_diseases: dict[str, list[str]] = {}
        for g, d, t in rows("hpo_gd", "SELECT gene, disease_id, association_type FROM '{p}'"):
            if t == "MENDELIAN":
                self.gene_diseases.setdefault(g, []).append(d)
                self.gene_modes.setdefault(g, set()).update(self.hpo_modes.get(d, ("", []))[1])
        for g, d, t in rows("orpha_gd", "SELECT gene, disease_id, association_type FROM '{p}'"):
            if g and t and t.startswith("Disease-causing"):
                self.gene_diseases.setdefault(g, []).append(d)
                self.gene_modes.setdefault(g, set()).update(self.orpha.get(d, ("", []))[1])
        self.validity: dict[str, str] = {}
        from .evidence import GENE_VALIDITY_RANK
        for g, cls, m in rows("validity", "SELECT gene_symbol, classification, moi FROM '{p}'"):
            if GENE_VALIDITY_RANK.get(cls, -1) > GENE_VALIDITY_RANK.get(self.validity.get(g), -1):
                self.validity[g] = cls
            if GENE_VALIDITY_RANK.get(cls, -1) >= 3:
                code = inh.from_text(m)
                if code:
                    self.gene_modes.setdefault(g, set()).add(code)
        self.acmg: dict[str, dict] = {}
        for g, cat, cond, inh_, rep in rows("acmg", "SELECT gene, category, condition, inheritance, report FROM '{p}'"):
            e = self.acmg.setdefault(g, {"category": cat, "conditions": [], "inheritance": inh_, "report": rep})
            e["conditions"].append(cond)
        self.actionability: dict[str, dict] = {}
        for ctx, gene, disease, outcome, intervention, overall, url in rows(
                "actionability", "SELECT context, gene, disease, outcome, intervention, overall, url FROM '{p}' "
                                 "WHERE status = 'Released'"):
            m = re.match(r"(\d+)", overall or "")
            if not m:
                continue
            score = int(m.group(1))
            for g in (gene or "").split(","):
                g = g.strip()
                best = self.actionability.get(g)
                if best is None or (ctx == "Adult", score) > (best["context"] == "Adult", best["score"]):
                    self.actionability[g] = {"context": ctx, "score": score, "code": overall, "disease": disease,
                                             "outcome": outcome, "intervention": intervention, "url": url}
        self.constraint = {g: {"loeuf": lo, "pli": pli, "mis_z": mz} for g, lo, pli, mz in
                           rows("constraint", "SELECT gene, loeuf, pli, mis_z FROM '{p}'")}

    # ---------------------------------------------------------------------------------------------------------

    def conditions(self, names: str | None, ids: str | None) -> list[Condition]:
        """ClinVar CLNDN / CLNDISDB (both '|'-separated per condition) → resolved conditions, filler removed."""
        name_list = (names or "").split("|")
        id_list = (ids or "").split("|") if ids else []
        out: dict[str, Condition] = {}
        for i, raw in enumerate(name_list):
            name = raw.strip()
            if not name or name.lower() in FILLER:
                continue
            name = name.capitalize() if name.isupper() else name
            c = Condition(name=name)
            group = id_list[i] if i < len(id_list) and len(id_list) == len(name_list) else ""
            mondo = None
            for x in group.split(","):
                x = x.strip()
                if x.startswith("MONDO:MONDO:"):
                    mondo = x[6:]
                elif x.startswith("OMIM:") and x[5:].isdigit():
                    c.ids.append(x)
                elif x.startswith("Orphanet:"):
                    c.ids.append("ORPHA:" + x.split(":", 1)[1].replace("ORPHA", ""))
            if mondo is None:
                mondo = next((self.xref_to_mondo[x] for x in c.ids if x in self.xref_to_mondo), None)
            if mondo:
                c.ids.insert(0, mondo)
                if mondo in self.mondo:
                    _, c.definition, xrefs = self.mondo[mondo]
                    c.ids += [x for x in xrefs if x.startswith(("OMIM:", "ORPHA:")) and x not in c.ids]
            for x in c.ids:
                if x in self.hpo_modes and self.hpo_modes[x][1] and not c.modes:
                    c.modes, c.modes_from = self.hpo_modes[x][1], "HPO"
                if x in self.orpha:
                    _, modes, onset, prev = self.orpha[x]
                    if modes and not c.modes:
                        c.modes, c.modes_from = modes, "Orphanet"
                    c.onset = c.onset or onset
                    c.prevalence = c.prevalence or (prev if prev and not prev.lower().startswith("unknown") else None)
            if not c.modes:
                cue = inh.cues(name)
                if cue:
                    c.modes, c.modes_from = sorted(cue), "name"
            key = name.casefold()
            if key not in out:
                out[key] = c
        return sorted(out.values(), key=lambda c: c.name.casefold())

    def disease_gene(self, gene: str | None) -> bool:
        return bool(gene) and (gene in self.gene_diseases or
                               self.validity.get(gene) in ("Definitive", "Strong", "Moderate"))

    def modes_for(self, gene: str | None, conds: list[Condition]) -> tuple[set[str], str]:
        """Inheritance for a claim: from its named conditions if any are known, else from the gene."""
        m = {x for c in conds for x in c.modes}
        if m:
            return m, "condition"
        g = self.gene_modes.get(gene or "", set())
        if gene in self.acmg:
            g = g | {self.acmg[gene]["inheritance"]}
        return g, ("gene" if g else "unknown")

    def gene_conditions(self, gene: str | None) -> list[Condition]:
        """Diseases linked to a gene (for predicted variants, which have no ClinVar condition)."""
        seen: dict[str, Condition] = {}
        for d in self.gene_diseases.get(gene or "", []):
            name = (self.hpo_modes.get(d) or self.orpha.get(d) or (None,))[0]
            mondo = self.xref_to_mondo.get(d)
            if mondo and mondo in self.mondo:
                name = self.mondo[mondo][0]
            if not name:
                continue
            name = name[0].upper() + name[1:]
            ids = ",".join(x.replace("ORPHA:", "Orphanet:") for x in [d])
            for c in self.conditions(name, ids):
                key = mondo or c.name.casefold()
                if key not in seen:
                    seen[key] = c
        return sorted(seen.values(), key=lambda c: c.name.casefold())
