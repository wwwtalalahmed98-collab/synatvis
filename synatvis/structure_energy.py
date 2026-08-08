"""RNA folding for the mRNA-ΔG axis (CLAUDE.md §3.5 extension).

Structure over the start codon impedes ribosome initiation (Kudla 2009; the
LinearDesign line of work, Zhang 2020). This module quantifies local structure
two ways:

* if **ViennaRNA** is installed, the real minimum free energy (MFE, kcal/mol) and
  the paired fraction of the MFE structure;
* otherwise a dependency-free **Nussinov** base-pairing maximiser gives the paired
  fraction (labelled a heuristic — no kcal/mol).

The flagging metric is the *paired fraction* (0..1), which is comparable between
backends; ΔG is reported as extra detail when ViennaRNA is present.
"""
from __future__ import annotations

from typing import Optional, Tuple

_PAIRS = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G"), ("G", "T"), ("T", "G")}

_VIENNA = None
# Nussinov is O(n^3) in pure Python; never run it on long sequences.
_NUSSINOV_MAX = 160


def has_vienna() -> bool:
    """True if ViennaRNA is importable (cached). Gate global folds on this."""
    global _VIENNA
    if _VIENNA is None:
        try:
            import RNA  # type: ignore  # noqa: F401
            _VIENNA = True
        except Exception:
            _VIENNA = False
    return _VIENNA


def _nussinov_pairs(seq: str, min_loop: int = 3) -> int:
    """Maximum number of base pairs (Nussinov DP). O(n^3); use on short windows."""
    n = len(seq)
    if n < min_loop + 2:
        return 0
    dp = [[0] * n for _ in range(n)]
    for span in range(min_loop + 1, n):
        for i in range(0, n - span):
            j = i + span
            best = dp[i + 1][j]                      # i unpaired
            if dp[i][j - 1] > best:                  # j unpaired
                best = dp[i][j - 1]
            if (seq[i], seq[j]) in _PAIRS:           # i,j pair
                inner = dp[i + 1][j - 1] if j - 1 >= i + 1 else 0
                if 1 + inner > best:
                    best = 1 + inner
            for k in range(i + 1, j):                # bifurcation
                v = dp[i][k] + dp[k + 1][j]
                if v > best:
                    best = v
            dp[i][j] = best
    return dp[0][n - 1]


def fold(seq: str) -> Tuple[Optional[float], float, str]:
    """Return ``(mfe_or_None, paired_fraction, backend)`` for *seq* (DNA letters)."""
    seq = seq.upper().replace("U", "T")
    if len(seq) < 6:
        return None, 0.0, "none"
    if has_vienna():
        import RNA  # type: ignore

        struct, mfe = RNA.fold(seq.replace("T", "U"))
        paired = sum(1 for ch in struct if ch in "()") / len(struct)
        return float(mfe), paired, "ViennaRNA"
    if len(seq) > _NUSSINOV_MAX:
        # too long for the pure-Python fallback; caller should gate global folds
        return None, 0.0, "heuristic-skipped"
    pairs = _nussinov_pairs(seq)
    return None, (2 * pairs) / len(seq), "heuristic"
