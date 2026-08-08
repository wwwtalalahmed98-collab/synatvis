"""polya — premature Cr poly(A) signal scanner (CLAUDE.md §2, §3.3, §7).

The Cr nuclear near-upstream element (NUE) is **UGUAA (TGTAA)** — ~52% of sites,
~10-30 nt upstream of the cleavage site — in a **G/C-rich** context (Shen 2008;
Zhao 2014). This is a FLIP from the plant AAUAAA / AU-rich logic and is NOT
lineage-gated. A premature NUE-in-context inside the 5'UTR or CDS risks a
truncated transcript.

This is the NOISIEST module. Its operating point (which context strength counts)
must be set from Leg-1 specificity on native Cr transcripts, not from assumption
(CLAUDE.md §7, §10). The thresholds here are the profile defaults, pending Leg-1.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

from ..flags import Flag, Severity
from ..remediation import synonymous_fix
from ..util import find_all, gc_fraction

NAME = "polya"
VALIDATED = True   # detection logic is validated by Leg-2; operating point by Leg-1
DEFAULT_ON = True
SUMMARY = "Premature Cr poly(A) signals (TGTAA in a G/C-rich context) that truncate the transcript"

EVIDENCE = "Shen 2008 (PMID 18493049); Zhao 2014 (PMC4025486): Cr NUE is UGUAA in G/C-rich context"


def _collect(seq: str, motifs: List[str]) -> List[Tuple[int, str]]:
    hits: List[Tuple[int, str]] = []
    for m in motifs:
        for pos in find_all(seq, m):
            hits.append((pos, m))
    hits.sort()
    # drop a hit whose span is fully contained in a previously kept hit
    kept: List[Tuple[int, str]] = []
    for pos, m in hits:
        end = pos + len(m)
        if kept and pos >= kept[-1][0] and end <= kept[-1][0] + len(kept[-1][1]):
            continue
        kept.append((pos, m))
    return kept


def run(tx, profile: Dict) -> List[Flag]:
    cfg = profile["polya"]
    motifs = [m.upper() for m in cfg.get("nue_motifs", [])]
    if not motifs:
        return []
    win = int(cfg["gc_context_window"])
    gcmin = float(cfg["gc_context_min"])
    table = profile.get("_codon_table")

    flags: List[Flag] = []
    for region, seq, offset in (("5utr", tx.utr5, 0),
                                ("cds", tx.cds, tx.cds_start),
                                ("3utr", tx.utr3, tx.cds_end)):
        for pos, motif in _collect(seq, motifs):
            end = pos + len(motif)
            ctx = seq[end:end + win]
            ctx_gc = gc_fraction(ctx)
            context_ok = len(ctx) >= win // 2 and ctx_gc >= gcmin

            if region == "3utr":
                if not context_ok:
                    continue
                severity = Severity.INFO
                msg = (f"Poly(A) signal {motif} with G/C-rich context in 3'UTR "
                       f"(context GC {ctx_gc:.0%}) — expected here, not a hazard.")
            else:
                if context_ok:
                    severity = Severity.HIGH if region == "cds" else Severity.MEDIUM
                    where = "inside the CDS" if region == "cds" else "in the 5'UTR"
                    msg = (f"Premature poly(A) signal {motif} {where} in a G/C-rich "
                           f"context (GC {ctx_gc:.0%}). Risks 3' cleavage and a "
                           f"truncated transcript.")
                else:
                    severity = Severity.LOW
                    msg = (f"{motif} in {region} but context is not G/C-rich "
                           f"(GC {ctx_gc:.0%} < {gcmin:.0%}); less likely to be a "
                           f"functional Cr NUE (Leg-1 operating point).")

            suggestion = None
            if region == "cds" and table is not None and context_ok:
                edit = synonymous_fix(tx.cds, table,
                                      cds_start=pos, cds_end=end, avoid=motifs)
                suggestion = edit.describe()
            flags.append(Flag(
                module=NAME, severity=severity,
                start=offset + pos, end=offset + end, region=region,
                message=msg, evidence=EVIDENCE,
                suggested_edit=suggestion,
                detail={"motif": motif, "context_gc": round(ctx_gc, 3),
                        "context_ok": context_ok},
            ))
    return flags
