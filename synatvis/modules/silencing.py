"""silencing — QUARANTINED, unvalidated heuristic (CLAUDE.md §2, §3.6, §9, §10).

Nuclear transgene silencing is the dominant Cr failure mode, but it is driven
mainly by chromatin/strain, NOT by CDS sequence. This module is therefore NOT a
sequence verdict and must never be presented as validated (``VALIDATED = False``;
the report marks its flags HEURISTIC). Its honest output is:

  * re-report the composition module's GC troughs as sequence-visible risk;
  * always emit one INFO flag carrying strain/promoter/intron guidance.
"""
from __future__ import annotations

from typing import Dict, List

from ..flags import Flag, Severity
from .composition import find_gc_troughs

NAME = "silencing"
VALIDATED = False   # QUARANTINED
DEFAULT_ON = True
SUMMARY = "HEURISTIC: sequence-visible silencing risk (GC troughs) + strain/intron guidance"

EVIDENCE = "Schroda 2019; Neupert/Bock UVM strains — silencing is chromatin/strain-driven"


def run(tx, profile: Dict) -> List[Flag]:
    cfg = profile["silencing"]
    flags: List[Flag] = []

    if cfg.get("reuse_composition_troughs", True):
        for start, end, min_gc in find_gc_troughs(tx.full, profile["composition"]):
            flags.append(Flag(
                module=NAME, severity=Severity.LOW,
                start=start, end=end, region=tx.region_at(start),
                message=(f"HEURISTIC: GC trough (min {min_gc:.0%}) is a "
                         f"sequence-visible silencing risk. NOT a verdict — "
                         f"silencing is mainly chromatin/strain-driven."),
                evidence=EVIDENCE,
                suggested_edit="raise local GC (see composition/codon)",
                detail={"min_gc": round(min_gc, 3), "heuristic": True},
            ))

    flags.append(Flag(
        module=NAME, severity=Severity.INFO,
        start=0, end=len(tx.full), region="transcript",
        message="HEURISTIC guidance: " + str(cfg.get("guidance", "")),
        evidence=EVIDENCE, suggested_edit=None,
        detail={"heuristic": True, "validated": False},
    ))
    return flags
