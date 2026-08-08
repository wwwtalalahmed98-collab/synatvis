"""Protein-fate and PTM prediction from the coding sequence (opt-in, protein-level).

The scanner is transcript-level; this module translates the CDS and reads
sequence-encoded *protein* features that determine where the protein goes and how
it is post-translationally modified. These govern protein fate more than
expression level, so they are reported as a complementary layer and drive the cell
visualisation's routing; at most one gentle, clearly-labelled modifier touches the
expression index.

Detected (all heuristic, sequence-based, and labelled as such):
  * signal peptide  — N-terminal hydrophobic h-region (von Heijne 1985; SignalP concept)
  * transmembrane   — Kyte-Doolittle hydropathy windows (Kyte & Doolittle 1982)
  * N-glycosylation — sequon N-X-S/T, X != P (Gavel & von Heijne 1990)
  * O-glyc / disulfide — Ser/Thr density and cysteine count (indicative only)
  * predicted localisation — secreted / membrane / cytosolic
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List

from .codon_tables import STANDARD_CODE

# Kyte-Doolittle hydropathy
_KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
       "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
       "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2, "X": 0.0}
_SEQUON = re.compile(r"N[^P][ST]")

# average residue masses (monoisotopic-free, Da) for MW; ExPASy ProtParam values
_AA_MW = {"A": 71.08, "R": 156.19, "N": 114.10, "D": 115.09, "C": 103.14,
          "E": 129.12, "Q": 128.13, "G": 57.05, "H": 137.14, "I": 113.16,
          "L": 113.16, "K": 128.17, "M": 131.19, "F": 147.18, "P": 97.12,
          "S": 87.08, "T": 101.10, "W": 186.21, "Y": 163.18, "V": 99.13, "X": 110.0}
# pKa set for theoretical pI (EMBOSS/ProtParam-style)
_PKA = {"Nterm": 9.69, "Cterm": 2.34, "C": 8.33, "D": 3.65, "E": 4.25,
        "H": 6.0, "K": 10.53, "R": 12.48, "Y": 10.07}


def gravy(prot: str) -> float:
    """Grand average of hydropathy (Kyte & Doolittle 1982)."""
    if not prot:
        return 0.0
    return sum(_KD.get(c, 0.0) for c in prot) / len(prot)


def molecular_weight(prot: str) -> float:
    """Approximate molecular weight in Da (residue masses + one water)."""
    return round(sum(_AA_MW.get(c, 110.0) for c in prot) + 18.02, 1) if prot else 0.0


def _net_charge(prot: str, pH: float) -> float:
    pos = 1.0 / (1.0 + 10 ** (pH - _PKA["Nterm"]))
    for aa, pk in (("K", _PKA["K"]), ("R", _PKA["R"]), ("H", _PKA["H"])):
        pos += prot.count(aa) / (1.0 + 10 ** (pH - pk))
    neg = 1.0 / (1.0 + 10 ** (_PKA["Cterm"] - pH))
    for aa, pk in (("D", _PKA["D"]), ("E", _PKA["E"]), ("C", _PKA["C"]), ("Y", _PKA["Y"])):
        neg += prot.count(aa) / (1.0 + 10 ** (pk - pH))
    return pos - neg


def theoretical_pi(prot: str) -> float:
    """Isoelectric point by bisection on net charge (ProtParam method)."""
    if not prot:
        return 7.0
    lo, hi = 0.0, 14.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if _net_charge(prot, mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 2)


def aromaticity(prot: str) -> float:
    """Fraction of aromatic residues (F, W, Y)."""
    if not prot:
        return 0.0
    return round(sum(prot.count(a) for a in "FWY") / len(prot), 3)


def biophysics(prot: str) -> Dict:
    """Sequence-derived protein biophysical parameters (all standard formulas)."""
    return {
        "length": len(prot),
        "mw_kda": round(molecular_weight(prot) / 1000.0, 1),
        "pI": theoretical_pi(prot),
        "gravy": round(gravy(prot), 3),
        "aromaticity": aromaticity(prot),
    }


@dataclass
class ProteinFate:
    protein_len: int
    localization: str
    signal_peptide: bool
    signal_region: List[int] = field(default_factory=list)  # [start, end] in aa
    tm_count: int = 0
    tm_segments: List[List[int]] = field(default_factory=list)
    n_glyc_sites: List[int] = field(default_factory=list)     # aa positions of the N
    st_fraction: float = 0.0
    cys_count: int = 0
    disulfide_potential: bool = False
    notes: List[str] = field(default_factory=list)
    modifier: float = 1.0            # gentle expression modifier (secretory burden)
    modifier_reason: str = ""

    def to_dict(self) -> Dict:
        return {"protein_len": self.protein_len, "localization": self.localization,
                "signal_peptide": self.signal_peptide, "tm_count": self.tm_count,
                "n_glyc_sites": len(self.n_glyc_sites), "cys_count": self.cys_count,
                "disulfide_potential": self.disulfide_potential,
                "modifier": self.modifier, "modifier_reason": self.modifier_reason,
                "notes": self.notes}


def translate(cds: str) -> str:
    cds = cds.upper()
    aa = []
    for k in range(0, len(cds) - len(cds) % 3, 3):
        a = STANDARD_CODE.get(cds[k:k + 3], "X")
        if a == "*":
            break
        aa.append(a)
    return "".join(aa)


def _win_mean(seq: str, i: int, w: int) -> float:
    return sum(_KD.get(c, 0.0) for c in seq[i:i + w]) / w


def _tm_segments(prot: str, w: int = 19, thr: float = 2.0, skip: int = 0) -> List[List[int]]:
    segs, i, n = [], skip, len(prot)
    while i <= n - w:
        if _win_mean(prot, i, w) >= thr:
            j = i
            while j <= n - w and _win_mean(prot, j, w) >= thr:
                j += 1
            segs.append([i, j + w - 1])
            i = j + w
        else:
            i += 1
    return segs


def predict_protein_fate(transcript) -> ProteinFate:
    prot = translate(transcript.cds)
    n = len(prot)
    if n == 0:
        return ProteinFate(0, "cytosolic", False, notes=["empty translation"])

    # signal peptide: a strong, N-terminal hydrophobic h-region followed by a
    # hydrophilic c-region. Conservative bar (specific, not sensitive) — a coarse
    # heuristic; use SignalP for real calls.
    signal, sig_region = False, []
    head = prot[:35]
    best, best_i = 0.0, 0
    for i in range(1, min(len(head) - 11, 14)):
        m = _win_mean(head, i, 11)
        if m > best:
            best, best_i = m, i
    if best >= 2.3 and best_i <= 13:
        c_region = _win_mean(prot, best_i + 11, 6) if best_i + 17 <= n else 0.0
        if c_region < 1.0:                 # hydrophilic c-region / cleavage boundary
            signal, sig_region = True, [1, best_i + 13]

    # transmembrane segments (skip the signal-peptide region so it isn't double-counted)
    skip = sig_region[1] if signal else 0
    tm = _tm_segments(prot, skip=skip)

    # N-glycosylation sequons (exclude the last two residues)
    glyc = [m.start() + 1 for m in _SEQUON.finditer(prot[:-2])] if n > 3 else []

    st_fraction = (prot.count("S") + prot.count("T")) / n
    cys = prot.count("C")

    if signal and not tm:
        loc = "secreted"
    elif tm:
        loc = f"membrane ({len(tm)} TM)"
    else:
        loc = "cytosolic"

    fate = ProteinFate(
        protein_len=n, localization=loc, signal_peptide=signal, signal_region=sig_region,
        tm_count=len(tm), tm_segments=tm, n_glyc_sites=glyc, st_fraction=round(st_fraction, 3),
        cys_count=cys, disulfide_potential=(cys >= 2 and cys % 2 == 0),
    )
    if signal:
        fate.notes.append("N-terminal signal peptide → secretory pathway (ER → Golgi).")
    if glyc:
        fate.notes.append(f"{len(glyc)} N-glycosylation sequon(s) (N-X-S/T); Cr N-glycosylation "
                          f"differs from mammalian.")
    if fate.disulfide_potential:
        fate.notes.append(f"{cys} cysteines → possible disulfide bonds (oxidising ER lumen).")

    # gentle, honest secretory-burden modifier: a secreted, heavily-glycosylated
    # protein carries extra folding/QC load that can lower net yield.
    if loc == "secreted" and len(glyc) >= 4:
        fate.modifier = 0.90
        fate.modifier_reason = ("secretory + heavy N-glycosylation load "
                                "(folding/QC burden; low-confidence protein-fate factor)")
    return fate
