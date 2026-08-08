"""Module registry (CLAUDE.md §4, §8).

Registers every detection module with the flags registry in the CLAUDE.md §8
build order: low-noise validated modules first, the noisy pair (polya, splice)
next, then the quarantined silencing module last. A module never hard-codes a
host constant — it reads the active profile — so adding the chloroplast profile
requires zero edits here or in any ``modules/*.py``.
"""
from __future__ import annotations

from ..flags import ModuleSpec, register
from . import (
    codon,
    composition,
    cloning,
    structure,
    uorf,
    instability,
    polya,
    splice,
    silencing,
)

_MODULES = [
    codon, composition, cloning,      # low-noise, validated
    structure, uorf, instability,     # carried over, re-tuned
    polya, splice,                    # the noisy pair (Leg-2 gated)
    silencing,                        # quarantined, unvalidated
]

for _m in _MODULES:
    register(ModuleSpec(
        name=_m.NAME,
        run=_m.run,
        validated=_m.VALIDATED,
        default_on=_m.DEFAULT_ON,
        summary=_m.SUMMARY,
    ))

__all__ = ["codon", "composition", "cloning", "structure", "uorf",
           "instability", "polya", "splice", "silencing"]
