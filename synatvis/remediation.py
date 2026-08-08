"""Minimal synonymous-edit engine (CLAUDE.md §6).

For a CDS-internal, pattern-based flag (a TGTAA NUE, a Type IIS site, an ARE, an
accidental splice motif), compute the *minimal set of synonymous codon changes*
that removes the motif AND does not recreate it in the local window, preserving
the amino-acid sequence exactly. Candidate codons are ranked by the Cr profile's
adaptiveness so a fix never pushes GC/codon usage the wrong way. If no clean fix
exists, say so explicitly — never fabricate an edit.
"""
from __future__ import annotations

import itertools
import re
from dataclasses import dataclass
from typing import List, Optional, Sequence

from .codon_tables import CodonTable
from .util import iupac_to_regex


@dataclass
class Edit:
    """A concrete synonymous edit, or an explicit no-fix result."""

    ok: bool
    reason: str = ""
    cds_start: int = 0                 # CDS-relative start of the edited window
    changes: Optional[List[dict]] = None  # per-codon {index, aa, old, new}
    old_window: str = ""
    new_window: str = ""
    delta_adaptiveness: float = 0.0

    def describe(self) -> str:
        if not self.ok:
            return f"no clean synonymous fix — {self.reason}"
        parts = [
            f"codon {c['index']} {c['old']}->{c['new']} ({c['aa']})"
            for c in (self.changes or [])
        ]
        arrow = f"{self.old_window} -> {self.new_window}"
        return f"synonymous: {', '.join(parts)}  [{arrow}]"


def _has_any(seq: str, patterns: Sequence[re.Pattern]) -> bool:
    return any(p.search(seq) for p in patterns)


def _compile(avoid: Sequence[str]) -> List[re.Pattern]:
    return [re.compile(iupac_to_regex(a)) for a in avoid if a]


def synonymous_fix(
    cds: str,
    table: CodonTable,
    cds_start: int,
    cds_end: int,
    avoid: Sequence[str],
    flank: int = 6,
    max_codons_changed: int = 3,
) -> Edit:
    """Find the minimal synonymous edit removing every *avoid* pattern.

    Parameters
    ----------
    cds : the coding sequence (in-frame from index 0).
    cds_start, cds_end : half-open CDS-relative span of the motif to remove.
    avoid : motifs (IUPAC ok) that must NOT appear in the local window afterward,
            e.g. the NUE variants, or a Type IIS site plus its reverse complement.
    flank : nt of context on each side that is also kept motif-free.
    """
    n = len(cds)
    if cds_start < 0 or cds_end > n or cds_start >= cds_end:
        return Edit(False, "motif span is not inside the CDS")

    patterns = _compile(avoid)
    first_codon = cds_start // 3
    last_codon = (cds_end - 1) // 3
    codon_indices = list(range(first_codon, last_codon + 1))

    # local window (codon-aligned outer bounds + flank) used for the recreation test
    win_lo = max(0, first_codon * 3 - flank)
    win_hi = min(n, (last_codon + 1) * 3 + flank)

    def window_of(seq: str) -> str:
        return seq[win_lo:win_hi]

    if not _has_any(window_of(cds), patterns):
        return Edit(False, "motif not present in the local window (nothing to fix)")

    # per-codon synonymous options (best-adaptiveness first), excluding the original
    options = {}
    base_weight = 0.0
    for ci in codon_indices:
        codon = cds[ci * 3:ci * 3 + 3]
        base_weight += table.weight.get(codon, 0.0)
        syns = [c for c in table.synonyms(codon) if c != codon]
        options[ci] = syns

    # try changing k codons, smallest k first; within k, prefer highest total adaptiveness
    for k in range(1, min(max_codons_changed, len(codon_indices)) + 1):
        best: Optional[Edit] = None
        best_score = None
        for combo in itertools.combinations(codon_indices, k):
            choice_lists = [options[ci] for ci in combo]
            if any(not lst for lst in choice_lists):
                continue
            for picks in itertools.product(*choice_lists):
                mutated = list(cds)
                for ci, new in zip(combo, picks):
                    mutated[ci * 3:ci * 3 + 3] = list(new)
                mseq = "".join(mutated)
                if _has_any(window_of(mseq), patterns):
                    continue
                new_weight = sum(table.weight.get(new, 0.0) for new in picks)
                old_weight = sum(table.weight.get(cds[ci * 3:ci * 3 + 3], 0.0)
                                 for ci in combo)
                score = new_weight - old_weight
                if best_score is None or score > best_score:
                    changes = [
                        {"index": ci, "aa": table.aa_of(cds[ci * 3:ci * 3 + 3]),
                         "old": cds[ci * 3:ci * 3 + 3], "new": new}
                        for ci, new in zip(combo, picks)
                    ]
                    best = Edit(
                        ok=True,
                        cds_start=win_lo,
                        changes=changes,
                        old_window=window_of(cds),
                        new_window=window_of(mseq),
                        delta_adaptiveness=round(score, 4),
                    )
                    best_score = score
        if best is not None:
            return best

    return Edit(False, "constraint conflict (no synonymous change removes the motif)")
