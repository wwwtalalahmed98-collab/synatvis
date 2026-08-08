"""structure — start codon, 5'UTR length, and start-region accessibility (§3.5).

Checks that the CDS begins with ATG, that the 5'UTR is not pathologically long,
and that the start codon is not buried in strong secondary structure. If
ViennaRNA is installed it folds a window spanning the start codon and flags a
low-accessibility start; otherwise it falls back to a GC-pairing heuristic and
labels the flag HEURISTIC. Thresholds come from the profile.
"""
from __future__ import annotations

from typing import Dict, List

from ..flags import Flag, Severity
from ..structure_energy import fold, has_vienna

NAME = "structure"
VALIDATED = True
DEFAULT_ON = True
SUMMARY = "Start codon presence, 5'UTR length, and mRNA folding dG (start-region + global)"


def run(tx, profile: Dict) -> List[Flag]:
    cfg = profile["structure"]
    flags: List[Flag] = []

    # 1. start codon present
    if not tx.cds.upper().startswith("ATG"):
        flags.append(Flag(
            module=NAME, severity=Severity.HIGH,
            start=tx.cds_start, end=tx.cds_start + 3, region="cds",
            message=f"CDS does not begin with ATG (starts with {tx.cds[:3]!r}).",
            evidence="translation initiation requires an AUG start",
            suggested_edit=None, detail={"start_codon": tx.cds[:3]},
        ))

    # 2. very long 5'UTR
    if len(tx.utr5) > int(cfg["utr5_max_len"]):
        flags.append(Flag(
            module=NAME, severity=Severity.INFO,
            start=0, end=len(tx.utr5), region="5utr",
            message=(f"5'UTR is long ({len(tx.utr5)} nt > {cfg['utr5_max_len']}). "
                     f"Long 5'UTRs add scanning burden and uORF opportunity."),
            evidence="ribosome scanning heuristic",
            suggested_edit=None, detail={"utr5_len": len(tx.utr5)},
        ))

    # 3. mRNA folding dG — start-region structure impedes initiation
    #    (Kudla 2009; LinearDesign / Zhang 2020). Metric = paired fraction of the
    #    folded initiation window; thresholds are per-backend, set from the 5,000
    #    native genes.
    before = int(cfg.get("struct_window_before", 15))
    after = int(cfg.get("struct_window_after", 45))
    start_abs = tx.cds_start
    window = tx.full[max(0, start_abs - before): start_abs + after]
    if len(window) >= 20:
        mfe, paired, backend = fold(window)
        if backend == "ViennaRNA":
            thr = float(cfg.get("struct_paired_max", 0.70))
            if paired > thr:
                flags.append(Flag(
                    module=NAME, severity=Severity.MEDIUM,
                    start=max(0, start_abs - before), end=start_abs + after,
                    region="5utr" if tx.utr5 else "cds",
                    message=(f"Structured start: {paired:.0%} of the initiation window "
                             f"is base-paired (dG {mfe:.1f} kcal/mol). Strong 5' "
                             f"structure impedes ribosome initiation."),
                    evidence="ViennaRNA fold of the initiation window (Kudla 2009; Zhang 2020)",
                    suggested_edit="reduce pairing near the AUG (synonymous changes / 5'UTR redesign)",
                    detail={"mfe": round(mfe, 2), "paired_fraction": round(paired, 3),
                            "backend": backend},
                ))
        elif backend == "heuristic":
            thr = float(cfg.get("struct_paired_max_heuristic", 0.80))
            if paired > thr:
                flags.append(Flag(
                    module=NAME, severity=Severity.LOW,
                    start=max(0, start_abs - before), end=start_abs + after,
                    region="5utr" if tx.utr5 else "cds",
                    message=(f"HEURISTIC (ViennaRNA absent): the initiation window is "
                             f"highly base-paired ({paired:.0%}, Nussinov) and may "
                             f"impede initiation. Install ViennaRNA for a real dG."),
                    evidence="Nussinov base-pairing fallback",
                    suggested_edit=None,
                    detail={"paired_fraction": round(paired, 3), "backend": backend},
                ))

    # 4. global mRNA structure — informational (ViennaRNA only; folding stabilises
    #    the mRNA, LinearDesign/Zhang 2020, so this is context, not a hazard).
    #    Gated on ViennaRNA: the Nussinov fallback is O(n^3) and must not run here.
    if has_vienna() and len(tx.cds) >= 60:
        seg = tx.cds[:300]                       # 5' region dG (cheap; most relevant)
        mfe_g, paired_g, backend_g = fold(seg)
        if backend_g == "ViennaRNA" and mfe_g is not None:
            per_nt = mfe_g / len(seg)
            flags.append(Flag(
                module=NAME, severity=Severity.INFO,
                start=tx.cds_start, end=tx.cds_start + len(seg), region="cds",
                message=(f"5' mRNA folding dG {per_nt:.2f} kcal/mol/nt "
                         f"({paired_g:.0%} paired). More structure tends to stabilise "
                         f"the mRNA (context, not a hazard)."),
                evidence="ViennaRNA global fold (Zhang 2020, LinearDesign)",
                suggested_edit=None,
                detail={"mfe_per_nt": round(per_nt, 3), "paired_fraction": round(paired_g, 3)},
            ))
    return flags
