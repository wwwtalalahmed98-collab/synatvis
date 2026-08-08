"""instability — AU-rich element (ARE) clustering (CLAUDE.md §3.5).

Clustered AU-rich elements (ATTTA core) are classic mRNA-destabilising signals,
most relevant in the 3'UTR but scanned everywhere. Carried over from the plant
tool and re-tuned to Cr. CDS-internal AREs get a minimal synonymous fix.
"""
from __future__ import annotations

from typing import Dict, List

from ..flags import Flag, Severity
from ..remediation import synonymous_fix
from ..util import find_all

NAME = "instability"
VALIDATED = True
DEFAULT_ON = True
SUMMARY = "Clustered AU-rich elements (AREs) that destabilise the transcript"

EVIDENCE = "ARE-mediated mRNA decay (AUUUA core); re-tuned to Cr"


def run(tx, profile: Dict) -> List[Flag]:
    cfg = profile["instability"]
    motifs = [m.upper() for m in cfg.get("are_motifs", [])]
    if not motifs:
        return []
    cluster_min = int(cfg["are_cluster_min"])
    window = int(cfg["are_window"])
    table = profile.get("_codon_table")

    flags: List[Flag] = []
    for region, seq, offset in (("5utr", tx.utr5, 0),
                                ("cds", tx.cds, tx.cds_start),
                                ("3utr", tx.utr3, tx.cds_end)):
        hits: List[int] = []
        for m in motifs:
            hits.extend(find_all(seq, m))
        hits.sort()
        if len(hits) < cluster_min:
            continue
        # find any window of `window` nt containing >= cluster_min hits
        reported_spans = []
        for a in range(len(hits)):
            in_win = [h for h in hits if hits[a] <= h < hits[a] + window]
            if len(in_win) >= cluster_min:
                span = (in_win[0], in_win[-1] + 5)
                if reported_spans and span[0] <= reported_spans[-1][1]:
                    continue
                reported_spans.append(span)
                suggestion = None
                if region == "cds" and table is not None:
                    edit = synonymous_fix(tx.cds, table,
                                          cds_start=span[0], cds_end=span[1],
                                          avoid=motifs)
                    suggestion = edit.describe()
                flags.append(Flag(
                    module=NAME, severity=Severity.MEDIUM,
                    start=offset + span[0], end=offset + span[1], region=region,
                    message=(f"ARE cluster: {len(in_win)} AU-rich elements within "
                             f"{window} nt in {region}. Destabilises the mRNA."),
                    evidence=EVIDENCE,
                    suggested_edit=suggestion,
                    detail={"n_are": len(in_win)},
                ))
    return flags
