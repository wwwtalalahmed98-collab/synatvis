"""SynAT.Vis — Synthetic Algal Transcript Visualiser / validator.

A transcript-level red-flag scanner for recombinant gene cassettes destined for
nuclear expression in *Chlamydomonas reinhardtii*. It reports *why a designed
cassette might fail to produce an intact, well-translated transcript* in the
host, with a suggested synonymous fix per flag.

It is diagnostic, not predictive: it emits a ranked flag list, never a composite
"expression score" (see CLAUDE.md §1, §9). It is silent on proteolysis / protein
stability, and its nuclear report says nothing about a chloroplast construct.
"""
from __future__ import annotations

__version__ = "0.1.0"

from .seqio import Transcript, read_fasta, read_record  # noqa: E402
from .flags import Flag, Severity  # noqa: E402
from .scanner import scan  # noqa: E402
from .profiles import load_profile  # noqa: E402

__all__ = [
    "Transcript",
    "read_fasta",
    "read_record",
    "Flag",
    "Severity",
    "scan",
    "load_profile",
    "__version__",
]
