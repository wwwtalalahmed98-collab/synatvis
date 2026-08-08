"""uorf — upstream AUG / uORF scanner in the 5'UTR (CLAUDE.md §3.5).

An AUG in the 5'UTR can initiate an upstream ORF and reduce initiation at the
main start. Flags every uAUG (when the profile enables it), raising severity when
the uAUG sits in a strong initiation context or opens a uORF that overlaps the
main start out of frame.
"""
from __future__ import annotations

from typing import Dict, List

from ..codon_tables import STOPS
from ..flags import Flag, Severity

NAME = "uorf"
VALIDATED = True
DEFAULT_ON = True
SUMMARY = "Upstream AUGs / uORFs in the 5'UTR that compete with the main start"

EVIDENCE = "ribosome scanning: 5'UTR AUGs divert initiation from the main ORF"


def _strong_context(utr5: str, pos: int) -> bool:
    """Kozak-like: a purine at -3 (or the CDS G at +4) marks a strong context."""
    minus3 = utr5[pos - 3] if pos - 3 >= 0 else ""
    return minus3 in ("A", "G")


def run(tx, profile: Dict) -> List[Flag]:
    cfg = profile["uorf"]
    if not cfg.get("flag_any_uaug", True):
        return []
    utr5 = tx.utr5.upper()
    flags: List[Flag] = []

    i = utr5.find("ATG")
    while i != -1:
        # find the first in-frame stop within the 5'UTR from this uAUG
        stop_at = None
        for k in range(i, len(utr5) - 2, 3):
            if utr5[k:k + 3] in STOPS:
                stop_at = k
                break
        strong = _strong_context(utr5, i)
        if stop_at is not None:
            kind = "uORF (opens and closes within the 5'UTR)"
            severity = Severity.MEDIUM if strong else Severity.LOW
            end = stop_at + 3
        else:
            # no stop before CDS: overlaps the main start; frame vs main AUG matters
            in_frame_with_main = ((len(utr5) - i) % 3 == 0)
            kind = ("uAUG in frame with main ATG (N-terminal extension)"
                    if in_frame_with_main else
                    "uAUG out of frame, overlapping the main start")
            severity = Severity.MEDIUM if (strong or not in_frame_with_main) else Severity.LOW
            end = len(utr5)
        flags.append(Flag(
            module=NAME, severity=severity,
            start=i, end=end, region="5utr",
            message=(f"Upstream AUG at 5'UTR position {i}: {kind}. "
                     f"{'Strong' if strong else 'Weak'} initiation context."),
            evidence=EVIDENCE,
            suggested_edit=("mutate the uAUG (e.g. ATG->ACG) in the 5'UTR "
                            "if it is not functional"),
            detail={"uaug_pos": i, "strong_context": strong},
        ))
        i = utr5.find("ATG", i + 1)
    return flags
