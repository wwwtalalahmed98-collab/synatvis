"""codon — rare-codon-cluster scan on the Cr nuclear GC-biased table (§3.2).

Uses the profile's GC-biased Cr codon table. Codon usage in Cr drives BOTH
translational efficiency AND mRNA stability (Barahimipour 2015), so this is more
load-bearing here than CAI was in the plant tool — but it is still ONE input,
never a verdict. Flags clusters of rare (low-adaptiveness) codons and notes a
low whole-CDS Relative Codon Adaptation. No composite score is emitted.
"""
from __future__ import annotations

from typing import Dict, List

from ..codon_tables import attach_advanced, gene_tai
from ..flags import Flag, Severity

NAME = "codon"
VALIDATED = True
DEFAULT_ON = True
SUMMARY = "Rare-codon clusters, codon optimality (mRNA stability), codon-pair bias, and tRNA adaptation (tAI)"

EVIDENCE = "Barahimipour 2015 Plant J: Cr codon usage sets translational efficiency and mRNA stability"
EVIDENCE_OPT = "Presnyak 2015 Cell; Hanson 2017: codon optimality is a major determinant of mRNA stability"
EVIDENCE_CPB = "Coleman 2008: codon-pair bias affects translation; multi-criteria design (Demissie 2025)"
EVIDENCE_TAI = "dos Reis 2004 tAI; Cr tRNA gene counts from GtRNAdb (Crein5)"


def run(tx, profile: Dict) -> List[Flag]:
    table = profile["_codon_table"]
    cfg = profile["codon"]
    thresh = float(cfg["rare_weight_threshold"])
    win = int(cfg["cluster_window_codons"])
    min_rare = int(cfg["cluster_min_rare"])

    codons = [c for _, c in tx.codons()]
    weights = [table.weight.get(c, 0.0) for c in codons]
    rare_mask = [w < thresh for w in weights]

    flags: List[Flag] = []

    # sliding window over codons; merge overlapping hot windows into one region
    hot: List[int] = []  # list of codon-window start indices that trip the threshold
    for i in range(0, max(0, len(codons) - win + 1)):
        if sum(rare_mask[i:i + win]) >= min_rare:
            hot.append(i)

    # merge consecutive hot starts into spans
    spans = []
    for i in hot:
        if spans and i <= spans[-1][1]:
            spans[-1][1] = i + win
        else:
            spans.append([i, i + win])

    for c0, c1 in spans:
        c1 = min(c1, len(codons))
        rare_here = [(c0 + k, codons[c0 + k]) for k in range(c1 - c0)
                     if rare_mask[c0 + k]]
        abs_start = tx.abs_from_cds(c0 * 3)
        abs_end = tx.abs_from_cds(c1 * 3)
        # suggestion: best synonym for each rare codon
        subs = []
        for idx, cod in rare_here[:6]:
            best = next((s for s in table.synonyms(cod) if s != cod), None)
            if best and table.weight.get(best, 0) > table.weight.get(cod, 0):
                subs.append(f"codon {idx} {cod}->{best}")
        flags.append(Flag(
            module=NAME,
            severity=Severity.MEDIUM,
            start=abs_start,
            end=abs_end,
            region="cds",
            message=(f"Rare-codon cluster: {len(rare_here)} low-adaptiveness codons "
                     f"within {win} codons. In Cr this slows elongation and can "
                     f"destabilise the mRNA."),
            evidence=EVIDENCE,
            suggested_edit=("; ".join(subs) if subs else None),
            detail={"n_rare": len(rare_here),
                    "rare_codons": [c for _, c in rare_here],
                    "threshold": thresh},
        ))

    if codons:
        rca = table.rca(tx.cds)
        if rca < float(cfg["low_rca_warn"]):
            flags.append(Flag(
                module=NAME,
                severity=Severity.INFO,
                start=tx.cds_start,
                end=tx.cds_end,
                region="cds",
                message=(f"Whole-CDS Relative Codon Adaptation is low (RCA={rca:.2f} "
                         f"vs warn {float(cfg['low_rca_warn']):.2f}). One input, not a "
                         f"verdict — no expression score is implied."),
                evidence=EVIDENCE,
                suggested_edit=None,
                detail={"rca": round(rca, 3)},
            ))

    # --- Tier-A frontier metrics (multi-criteria panel) ---
    flags.extend(_advanced_checks(tx, profile, cfg, codons, table))
    return flags


def _advanced_checks(tx, profile, cfg, codons, table) -> List[Flag]:
    adv = attach_advanced(profile, profile.get("_base_dir", "."))
    opt, pairs, tai_w = adv["optimality"], adv["pairs"], adv["tai"]
    flags: List[Flag] = []
    if not codons:
        return flags

    # (1) codon-OPTIMALITY: windowed mean optimality below cutoff = de-optimized
    #     region (mRNA-stability risk), separate axis from adaptation.
    if opt:
        win = int(cfg.get("optimality_window", 30))
        cutoff = float(cfg.get("optimality_cutoff", -0.20))
        ovals = [opt.get(c, 0.0) for c in codons]
        spans = []
        for i in range(0, max(0, len(ovals) - win + 1)):
            if sum(ovals[i:i + win]) / win < cutoff:
                if spans and i <= spans[-1][1]:
                    spans[-1][1] = i + win
                else:
                    spans.append([i, i + win])
        for c0, c1 in spans:
            c1 = min(c1, len(codons))
            mean_o = sum(ovals[c0:c1]) / (c1 - c0)
            # suggest highest-optimality synonyms for the worst codons in the span
            worst = sorted(range(c0, c1), key=lambda k: ovals[k])[:6]
            subs = []
            for k in worst:
                cod = codons[k]
                alt = max((s for s in table.synonyms(cod)), key=lambda s: opt.get(s, -9),
                          default=None)
                if alt and opt.get(alt, -9) > opt.get(cod, -9) + 0.3:
                    subs.append(f"codon {k} {cod}->{alt}")
            flags.append(Flag(
                module=NAME, severity=Severity.MEDIUM,
                start=tx.abs_from_cds(c0 * 3), end=tx.abs_from_cds(c1 * 3), region="cds",
                message=(f"De-optimized region: mean codon optimality {mean_o:+.2f} over "
                         f"{c1 - c0} codons (non-optimal codons destabilise the mRNA)."),
                evidence=EVIDENCE_OPT,
                suggested_edit=("; ".join(subs) if subs else
                                "raise optimality with higher-optimality synonyms"),
                detail={"mean_optimality": round(mean_o, 3)},
            ))

    # (2) codon-PAIR bias: a run of consecutive strongly under-represented pairs.
    if pairs:
        cpb_cut = float(cfg.get("cpb_cutoff", -0.5))
        run_len = int(cfg.get("cpb_run_len", 5))
        cur, best = 0, None
        runs = []
        for i in range(len(codons) - 1):
            cps = pairs.get(codons[i] + codons[i + 1], 0.0)
            if cps < cpb_cut:
                cur += 1
                if cur >= run_len:
                    best = (i - cur + 1, i + 1)
            else:
                if best:
                    runs.append(best)
                cur, best = 0, None
        if best:
            runs.append(best)
        for a, b in runs:
            flags.append(Flag(
                module=NAME, severity=Severity.LOW,
                start=tx.abs_from_cds(a * 3), end=tx.abs_from_cds((b + 1) * 3), region="cds",
                message=(f"Codon-pair-deoptimized stretch: {b - a + 1} consecutive "
                         f"under-represented codon pairs (CPS < {cpb_cut}); a known "
                         f"translational drag."),
                evidence=EVIDENCE_CPB, suggested_edit=None,
                detail={"pair_run": b - a + 1},
            ))

    # (3) tAI: whole-CDS tRNA adaptation below the native operating point.
    if tai_w:
        t = gene_tai(tx.cds, tai_w)
        warn = float(cfg.get("low_tai_warn", 0.34))
        if 0 < t < warn:
            flags.append(Flag(
                module=NAME, severity=Severity.INFO,
                start=tx.cds_start, end=tx.cds_end, region="cds",
                message=(f"Low tRNA adaptation (tAI={t:.2f} vs native operating point "
                         f"{warn:.2f}); codons are poorly matched to the Cr tRNA pool."),
                evidence=EVIDENCE_TAI, suggested_edit=None,
                detail={"tai": round(t, 3)},
            ))
    return flags
