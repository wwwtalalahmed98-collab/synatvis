"""composition — GC-trough and homopolymer scanner (CLAUDE.md §3.1, §9).

INVERTED relative to the plant tool: in Cr nuclear, **LOW GC is the hazard**
(unfavourable/low GC triggers heterochromatinisation and epigenetic transgene
silencing; Barahimipour 2015, PMID 26402748). High GC is NORMAL and is never
flagged. This module flags local GC *troughs* and long homopolymers (a synthesis
and processivity risk). It hard-codes no constant — all thresholds come from the
active profile.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

from ..flags import Flag, Severity
from ..util import gc_fraction, windows, longest_homopolymer

NAME = "composition"
VALIDATED = True
DEFAULT_ON = True
SUMMARY = "Local GC troughs (silencing risk; LOW GC is the hazard) and homopolymer runs"

EVIDENCE = "Barahimipour 2015 Plant J (PMID 26402748): low/unfavourable GC drives Cr transgene silencing"


def find_gc_troughs(seq: str, cfg: Dict) -> List[Tuple[int, int, float]]:
    """Return merged low-GC troughs as ``(start, end, min_gc)``.

    A window counts as low if its GC is below ``gc_trough_warn``; contiguous low
    windows are merged; a trough is kept only if it spans >= ``min_trough_len``.
    """
    size = int(cfg["window"])
    step = int(cfg["step"])
    warn = float(cfg["gc_trough_warn"])
    min_len = int(cfg["min_trough_len"])

    low_spans: List[Tuple[int, int, float]] = []
    for start, sub in windows(seq, size, step):
        gc = gc_fraction(sub)
        if gc < warn:
            low_spans.append((start, start + size, gc))

    if not low_spans:
        return []
    merged: List[List[float]] = []
    for s, e, gc in low_spans:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
            merged[-1][2] = min(merged[-1][2], gc)
        else:
            merged.append([s, e, gc])
    return [(int(s), int(e), gc) for s, e, gc in merged if (e - s) >= min_len]


def run(tx, profile: Dict) -> List[Flag]:
    cfg = profile["composition"]
    full = tx.full
    flags: List[Flag] = []

    low = float(cfg["gc_trough_low"])
    for start, end, min_gc in find_gc_troughs(full, cfg):
        severity = Severity.MEDIUM if min_gc < low else Severity.LOW
        flags.append(Flag(
            module=NAME,
            severity=severity,
            start=start,
            end=end,
            region=tx.region_at(start),
            message=(f"GC trough: local GC drops to {min_gc:.0%} over "
                     f"{end - start} nt (target ~{float(cfg['target_gc']):.0%}). "
                     f"Low GC risks heterochromatinisation / transgene silencing."),
            evidence=EVIDENCE,
            suggested_edit=("raise local GC by choosing higher-GC synonymous codons "
                            "in this window (see codon module); if UTR, redesign"),
            detail={"min_gc": round(min_gc, 3), "target_gc": cfg["target_gc"],
                    "hazard": "low_gc"},
        ))

    hp_min = int(cfg["homopolymer_min"])
    # scan each region separately so the run does not span a UTR/CDS boundary silently
    for region, seq, offset in (("5utr", tx.utr5, 0),
                                ("cds", tx.cds, tx.cds_start),
                                ("3utr", tx.utr3, tx.cds_end)):
        i = 0
        s = seq.upper()
        while i < len(s):
            j = i
            while j < len(s) and s[j] == s[i]:
                j += 1
            run_len = j - i
            if run_len >= hp_min:
                flags.append(Flag(
                    module=NAME,
                    severity=Severity.LOW,
                    start=offset + i,
                    end=offset + j,
                    region=region,
                    message=(f"Homopolymer run: {run_len}x '{s[i]}'. Long runs cause "
                             f"synthesis errors and polymerase slippage."),
                    evidence="synthesis/processivity heuristic",
                    suggested_edit=("break the run with a synonymous substitution "
                                    "(CDS) or redesign (UTR)"),
                    detail={"base": s[i], "length": run_len},
                ))
            i = j
    return flags
