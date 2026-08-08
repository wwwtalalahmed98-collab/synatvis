"""splice — re-scoped Cr splice-site check (CLAUDE.md §3.4, §9).

Introns are ASSETS in Cr, not hazards: short, GC-rich, and deliberately inserted
(e.g. rbcS2 intron 1) to make transgenes express (Baier 2018). The dicot AU-rich
cryptic-intron detector and dicot/monocot clade gating from the plant tool are
REMOVED. This module, at LOW priority, only:

  (a) recognises an intron-like feature (Cr-consensus donor GTRAG + acceptor YAG,
      short and GC-rich) and reports it as INFO — consistent with a deliberately
      inserted intron; verify its boundaries;
  (b) optionally flags an accidental strong Cr-type splice pair whose removal
      would be out of frame and thus fragment/frameshift the CDS.

Conservative by design (strict consensus) so it stays quiet on clean sequence;
gated by its Leg-2 injection test before it earns any higher severity.
"""
from __future__ import annotations

from typing import Dict, List

from ..flags import Flag, Severity
from ..util import find_all, gc_fraction

NAME = "splice"
VALIDATED = True
DEFAULT_ON = True
SUMMARY = "Cr-type splice sites: intended-intron sanity check and accidental splice pairs"

EVIDENCE = "Baier 2018 (NAR, PMC6061784); Schroda 2019 — Cr introns are short, GC-rich assets"


def run(tx, profile: Dict) -> List[Flag]:
    cfg = profile["splice"]
    donor = cfg.get("donor_consensus", "")
    acceptor = cfg.get("acceptor_consensus", "")
    if not donor or not acceptor:
        return []  # inert (e.g. chloroplast profile)

    min_i = int(cfg["min_intron"])
    max_i = int(cfg["max_intron"])
    cds = tx.cds.upper()
    offset = tx.cds_start

    donors = find_all(cds, donor, iupac=True)
    acceptors = set(find_all(cds, acceptor, iupac=True))

    flags: List[Flag] = []
    seen_spans = []
    for d in donors:
        # nearest valid acceptor forming an intron of allowed length
        for a in range(d + min_i, min(len(cds), d + max_i)):
            if a in acceptors:
                intron_start, intron_end = d, a + len(acceptor)
                intron = cds[intron_start:intron_end]
                if len(intron) < min_i:
                    continue
                if seen_spans and intron_start <= seen_spans[-1]:
                    break
                seen_spans.append(intron_end)
                gc = gc_fraction(intron)
                in_frame_removal = (len(intron) % 3 == 0)
                if gc >= 0.60:
                    flags.append(Flag(
                        module=NAME, severity=Severity.INFO,
                        start=offset + intron_start, end=offset + intron_end,
                        region="cds",
                        message=(f"Cr-type intron-like feature: donor {donor} / "
                                 f"acceptor {acceptor}, {len(intron)} nt, GC {gc:.0%}. "
                                 f"Consistent with a deliberately inserted intron "
                                 f"(e.g. {cfg.get('inserted_intron_name','')}); "
                                 f"verify boundaries."),
                        evidence=EVIDENCE, suggested_edit=None,
                        detail={"intron_len": len(intron), "intron_gc": round(gc, 3),
                                "in_frame_removal": in_frame_removal},
                    ))
                elif not in_frame_removal:
                    flags.append(Flag(
                        module=NAME, severity=Severity.LOW,
                        start=offset + intron_start, end=offset + intron_end,
                        region="cds",
                        message=(f"Possible accidental Cr-type splice pair "
                                 f"({len(intron)} nt, GC {gc:.0%}, out of frame). "
                                 f"If spliced, would frameshift/fragment the CDS."),
                        evidence=EVIDENCE, suggested_edit=None,
                        detail={"intron_len": len(intron), "intron_gc": round(gc, 3),
                                "in_frame_removal": in_frame_removal},
                    ))
                break  # one intron per donor
    return flags
