"""Codon-table loading, the standard genetic code, and adaptiveness (CLAUDE.md §4).

The active profile names a TSV of Cr nuclear codon usage. :class:`CodonTable`
exposes:

* ``translate`` / ``aa_of`` — the standard genetic code (host-independent);
* ``fraction`` — per-codon usage fraction within its amino-acid family;
* ``weight`` — relative adaptiveness ``w_i = fraction_i / max(fraction in aa)``,
  the building block for a Relative-Codon-Adaptation metric that here reflects
  BOTH translational efficiency and mRNA stability (Barahimipour 2015);
* ``synonyms`` — synonymous codons ranked by adaptiveness, used by remediation.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Tuple

# Standard genetic code (DNA alphabet). Host-independent; only *usage* is host-specific.
_BASES = "TCAG"
_AA = (
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
)
STANDARD_CODE: Dict[str, str] = {}
_i = 0
for _a in _BASES:
    for _b in _BASES:
        for _c in _BASES:
            STANDARD_CODE[_a + _b + _c] = _AA[_i]
            _i += 1

STOPS = frozenset(c for c, aa in STANDARD_CODE.items() if aa == "*")


class CodonTable:
    """Codon usage for one host profile plus the standard genetic code."""

    def __init__(self, fractions: Dict[str, float], source: str = "") -> None:
        self.source = source
        self.fraction: Dict[str, str] = {}
        # normalise: ensure every sense codon has a fraction (default 0)
        self.fraction = {c: float(fractions.get(c, 0.0)) for c in STANDARD_CODE}
        # group by amino acid
        self._by_aa: Dict[str, List[str]] = {}
        for codon, aa in STANDARD_CODE.items():
            self._by_aa.setdefault(aa, []).append(codon)
        # relative adaptiveness weights
        self.weight: Dict[str, float] = {}
        for aa, codons in self._by_aa.items():
            mx = max((self.fraction[c] for c in codons), default=0.0)
            for c in codons:
                self.weight[c] = (self.fraction[c] / mx) if mx > 0 else 0.0

    # -- genetic code ---------------------------------------------------
    def aa_of(self, codon: str) -> str:
        return STANDARD_CODE[codon.upper()]

    def translate(self, cds: str) -> str:
        cds = cds.upper()
        return "".join(
            STANDARD_CODE.get(cds[k:k + 3], "X")
            for k in range(0, len(cds) - len(cds) % 3, 3)
        )

    # -- usage ----------------------------------------------------------
    def synonyms(self, codon: str, exclude_stop: bool = True) -> List[str]:
        """Synonymous codons for *codon*, ranked by adaptiveness (best first)."""
        aa = self.aa_of(codon)
        codons = list(self._by_aa[aa])
        if exclude_stop and aa == "*":
            pass  # stops are synonymous only to stops; caller decides
        return sorted(codons, key=lambda c: (-self.weight[c], c))

    def is_rare(self, codon: str, threshold: float) -> bool:
        """True if the codon's adaptiveness weight is below *threshold*."""
        return self.weight.get(codon.upper(), 0.0) < threshold

    def rca(self, cds: str) -> float:
        """Geometric-mean Relative Codon Adaptation over sense codons (0..1)."""
        import math

        logs, n = 0.0, 0
        for k in range(0, len(cds) - len(cds) % 3, 3):
            codon = cds[k:k + 3].upper()
            if codon in STOPS or codon not in self.weight:
                continue
            w = self.weight[codon]
            if w <= 0:
                w = 1e-3  # floor so a single very-rare codon does not zero the score
            logs += math.log(w)
            n += 1
        return math.exp(logs / n) if n else 0.0


def load_tsv(path: str) -> CodonTable:
    """Load a codon-usage TSV (``codon<TAB>aa<TAB>fraction``; ``#`` comments)."""
    fractions: Dict[str, float] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) < 3 or parts[0].lower() == "codon":
                continue
            codon = parts[0].upper().replace("U", "T")
            try:
                fractions[codon] = float(parts[2])
            except ValueError:
                continue
    return CodonTable(fractions, source=os.path.basename(path))


def load_for_profile(profile: Dict, base_dir: str) -> CodonTable:
    """Load the codon table referenced by ``profile['codon']['table']``."""
    rel = profile.get("codon", {}).get("table")
    if not rel:
        raise ValueError("profile has no codon.table entry")
    path = rel if os.path.isabs(rel) else os.path.join(base_dir, rel)
    return load_tsv(path)


# ---------------------------------------------------------------------------
# Tier-A frontier tables (optimality, codon-pair bias, tAI) — loaded lazily
# ---------------------------------------------------------------------------
def _load_kv(path: str, key_col: int, val_col: int) -> Dict[str, float]:
    out: Dict[str, float] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) <= max(key_col, val_col) or parts[0].lower() in ("codon", "pair", "anticodon"):
                continue
            try:
                out[parts[key_col].upper()] = float(parts[val_col])
            except ValueError:
                continue
    return out


def attach_advanced(profile: Dict, base_dir: str) -> Dict:
    """Load and cache the optimality / codon-pair / tAI tables on the profile.

    Returns a dict with keys ``optimality`` {codon->float}, ``pairs`` {AB->cps},
    ``tai`` {codon->w}. Missing tables load as empty (the checks then no-op).
    """
    if "_advanced" in profile:
        return profile["_advanced"]
    cfg = profile.get("codon", {})

    def _p(key):
        rel = cfg.get(key)
        if not rel:
            return None
        return rel if os.path.isabs(rel) else os.path.join(base_dir, rel)

    adv = {"optimality": {}, "pairs": {}, "tai": {}}
    if _p("optimality_table") and os.path.isfile(_p("optimality_table")):
        adv["optimality"] = _load_kv(_p("optimality_table"), 0, 2)   # codon, aa, optimality
    if _p("codon_pair_table") and os.path.isfile(_p("codon_pair_table")):
        adv["pairs"] = _load_kv(_p("codon_pair_table"), 0, 1)        # pair, cps
    if _p("tai_table") and os.path.isfile(_p("tai_table")):
        adv["tai"] = _load_kv(_p("tai_table"), 0, 2)                 # codon, aa, w
    profile["_advanced"] = adv
    return adv


def gene_tai(cds: str, tai_w: Dict[str, float]) -> float:
    """Geometric-mean tAI of a CDS over sense codons with a cognate-tRNA weight."""
    import math
    logs, n = 0.0, 0
    for k in range(0, len(cds) - len(cds) % 3, 3):
        w = tai_w.get(cds[k:k + 3].upper())
        if w and w > 0:
            logs += math.log(w)
            n += 1
    return math.exp(logs / n) if n else 0.0
