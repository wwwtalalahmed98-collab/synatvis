"""Leg 2 (sweep) — sensitivity curves for the tunable modules (CLAUDE.md §7).

The plain Leg-2 test checks that a single injected motif is detected. This
scaffold turns that into a *detection curve*: it injects a signal at increasing
strength and reports the fraction detected at each level. Read together with the
Leg-1 false-positive rate, each curve gives the module a specificity/sensitivity
tradeoff — so a threshold can be chosen as an operating point on a real curve,
not from a false-positive number alone.

Signals swept:
  * polya       — TGTAA with downstream G/C context of increasing strength
                  (operating point: gc_context_min).
  * codon       — an 8-codon window with an increasing number of rare codons
                  (operating point: cluster_min_rare).
  * composition — a fixed-length stretch of decreasing local GC
                  (operating point: gc_trough_low / min_trough_len).

It tests DETECTION geometry, not biological consequence — label it as such. Each
level uses several replicate insertions (varied position and randomised fill) so
the reported number is a rate, not a single yes/no.
"""
from __future__ import annotations

import random
from typing import Dict, List, Tuple

from ..codon_tables import load_for_profile
from ..flags import Severity
from ..profiles import load_profile, PACKAGE_DIR
from ..scanner import scan
from ..seqio import Transcript

# a long, clean, GC-rich preferred-codon parent (many insertion positions, no flags)
CLEAN_PARENT = "ATG" + ("GCCGGCCTGGAGAAGGCC" * 12) + "TAA"


def _rand_seq(n: int, gc: float, rng: random.Random) -> str:
    return "".join(rng.choice("GC") if rng.random() < gc else rng.choice("AT")
                   for _ in range(n))


def _detect(cds: str, module: str, profile: Dict, min_sev=Severity.MEDIUM) -> bool:
    res = scan(Transcript(cds=cds), profile=profile, only=[module])
    return any(f.severity >= min_sev for f in res.flags)


def _positions(r: int) -> int:
    return 3 * (3 + (r % 20))  # codon-boundary insertion points across the parent


def sweep_polya(profile: Dict, R: int = 20) -> Tuple[float, List]:
    win = int(profile["polya"]["gc_context_window"])
    thr = float(profile["polya"]["gc_context_min"])
    rng = random.Random(1)
    rows = []
    for g in (0.40, 0.55, 0.70, 0.80, 0.85, 0.90, 1.00):
        hits = 0
        for r in range(R):
            ctx = _rand_seq(win, g, rng)
            pos = _positions(r)
            cds = CLEAN_PARENT[:pos] + "TGTAA" + ctx + CLEAN_PARENT[pos:]
            hits += _detect(cds, "polya", profile)
        rows.append((g, hits / R))
    return thr, rows


def sweep_codon(profile: Dict, R: int = 10) -> Tuple[int, List]:
    table = load_for_profile(profile, PACKAGE_DIR)
    thr_w = float(profile["codon"]["rare_weight_threshold"])
    win = int(profile["codon"]["cluster_window_codons"])
    cmin = int(profile["codon"]["cluster_min_rare"])
    rare = [c for c in table.weight if 0 < table.weight[c] < thr_w]
    pref = ["GCC", "CTG", "GGC", "GAG", "AAG", "GAC"]
    rng = random.Random(2)
    rows = []
    for k in range(0, win + 1):
        hits = 0
        for r in range(R):
            block = [rng.choice(rare) for _ in range(k)] + \
                    [rng.choice(pref) for _ in range(win - k)]
            rng.shuffle(block)
            pos = _positions(r)
            cds = CLEAN_PARENT[:pos] + "".join(block) + CLEAN_PARENT[pos:]
            hits += _detect(cds, "codon", profile)
        rows.append((k, hits / R))
    return cmin, rows


def sweep_composition(profile: Dict, R: int = 20) -> Tuple[float, List]:
    low = float(profile["composition"]["gc_trough_low"])
    mlen = int(profile["composition"]["min_trough_len"])
    L = mlen + 15
    rng = random.Random(3)
    rows = []
    for g in (0.20, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55):
        hits = 0
        for r in range(R):
            block = _rand_seq(L, g, rng)
            pos = _positions(r)
            cds = CLEAN_PARENT[:pos] + block + CLEAN_PARENT[pos:]
            hits += _detect(cds, "composition", profile)
        rows.append((g, hits / R))
    return low, rows


def _bar(frac: float) -> str:
    return "#" * int(round(frac * 20))


def main(profile_name: str = "cr_nuclear") -> int:
    profile = load_profile(profile_name)
    print(f"Leg 2 (sweep) — sensitivity/detection curves (profile: {profile_name})")
    print("  [tests DETECTION geometry, not biological consequence]")

    thr, rows = sweep_polya(profile)
    print(f"\n  polya — detection vs downstream context GC "
          f"(operating point gc_context_min = {thr:.2f}):")
    for g, f in rows:
        mark = "  <- threshold" if abs(g - thr) < 1e-6 else ""
        print(f"    context GC {g:>4.0%}   detected {f:>5.0%}  {_bar(f)}{mark}")

    cmin, rows = sweep_codon(profile)
    print(f"\n  codon — detection vs rare codons in an 8-codon window "
          f"(operating point cluster_min_rare = {cmin}):")
    for k, f in rows:
        mark = "  <- threshold" if k == cmin else ""
        print(f"    {k} rare / 8      detected {f:>5.0%}  {_bar(f)}{mark}")

    low, rows = sweep_composition(profile)
    print(f"\n  composition — detection vs local GC of a {int(profile['composition']['min_trough_len'])+15}-nt "
          f"stretch (operating point gc_trough_low = {low:.2f}):")
    for g, f in rows:
        mark = "  <- threshold" if abs(g - low) < 1e-6 else ""
        print(f"    local GC {g:>4.0%}     detected {f:>5.0%}  {_bar(f)}{mark}")

    print("\n  Read with Leg-1 FP: each curve is the sensitivity side of that module's")
    print("  operating point. Extend R / the level grid, or add position/length axes,")
    print("  to refine any threshold on its full specificity-sensitivity tradeoff.")
    return 0
