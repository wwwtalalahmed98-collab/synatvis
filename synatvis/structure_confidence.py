"""Self-explanatory structure-confidence layer.

DeepMind-style structure predictors (AlphaFold2/3) output per-residue confidence
(pLDDT) and a pairwise error map (PAE) — powerful, but you need an expert to read
the plots. This module does the interpretation SynAT.Vis is meant to add: it turns
confidence into plain language ("this part folds into a solid, well-defined shape";
"this stretch is likely floppy/disordered").

Two modes, both honest:
  * if real pLDDT is supplied (from an AlphaFold/Boltz plugin), it interprets that.
  * otherwise it computes a clearly-labelled SEQUENCE PROXY of order/disorder using
    the TOP-IDP scale (Campen et al. 2008) — an intrinsic-disorder predictor, NOT
    AlphaFold — so nothing is fabricated as a structure prediction.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .ptm import translate

# TOP-IDP disorder propensity (Campen et al. 2008); higher = more disorder-promoting
_TOPIDP = {
    "W": -0.884, "F": -0.697, "Y": -0.510, "I": -0.486, "M": -0.397, "L": -0.326,
    "V": -0.121, "N": 0.007, "C": 0.023, "T": 0.059, "A": 0.060, "G": 0.166,
    "R": 0.180, "D": 0.192, "H": 0.303, "Q": 0.318, "K": 0.586, "S": 0.341,
    "E": 0.736, "P": 0.987, "X": 0.0,
}


@dataclass
class StructureConfidence:
    source: str                     # "pLDDT (model)" or "sequence proxy (TOP-IDP)"
    mean_confidence: float          # 0-100
    per_residue: List[float] = field(default_factory=list)   # 0-100
    segments: List[Dict] = field(default_factory=list)       # ordered/disordered runs
    plain_summary: str = ""
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict:
        return {"source": self.source, "mean_confidence": self.mean_confidence,
                "n_residues": len(self.per_residue), "segments": self.segments,
                "plain_summary": self.plain_summary, "notes": self.notes}


def _smooth(vals: List[float], w: int = 11) -> List[float]:
    n = len(vals)
    out = []
    h = w // 2
    for i in range(n):
        a, b = max(0, i - h), min(n, i + h + 1)
        out.append(sum(vals[a:b]) / (b - a))
    return out


def _segments(conf: List[float], thr: float = 55.0, min_len: int = 6) -> List[Dict]:
    segs, i, n = [], 0, len(conf)
    while i < n:
        ordered = conf[i] >= thr
        j = i
        while j < n and (conf[j] >= thr) == ordered:
            j += 1
        if j - i >= min_len:
            segs.append({"start": i + 1, "end": j,
                         "kind": "folded domain" if ordered else "flexible / disordered",
                         "mean": round(sum(conf[i:j]) / (j - i), 0)})
        i = j
    # merge tiny gaps by relabelling handled implicitly; keep as-is for clarity
    return segs


def _plain(mean: float, segs: List[Dict], source: str) -> str:
    folded = [s for s in segs if s["kind"] == "folded domain"]
    flex = [s for s in segs if s["kind"].startswith("flexible")]
    if mean >= 70:
        head = "Overall this protein is predicted to fold into a solid, well-defined shape."
    elif mean >= 55:
        head = ("Overall this protein is mostly foldable, with some floppy regions that may "
                "wobble.")
    else:
        head = ("Overall this protein looks largely flexible/disordered — it may not settle "
                "into one fixed shape, which can affect stability, solubility and activity.")
    bits = [head]
    if folded:
        rr = ", ".join(f"{s['start']}–{s['end']}" for s in folded[:4])
        bits.append(f"Confidently-folded region(s): {rr} — treat these as the structured core.")
    if flex:
        rr = ", ".join(f"{s['start']}–{s['end']}" for s in flex[:4])
        bits.append(f"Likely-flexible region(s): {rr} — good places for linkers/tags, but "
                    f"watch for protease-sensitive or aggregation-prone stretches.")
    if "proxy" in source:
        bits.append("(This is a sequence-based order/disorder estimate, not an AlphaFold "
                     "structure — wire an AlphaFold3/Boltz plugin for a real fold.)")
    return "  ".join(bits)


def interpret_structure(transcript, plddt: Optional[List[float]] = None) -> StructureConfidence:
    """Plain-language structure confidence, from real pLDDT if given else a proxy."""
    prot = translate(transcript.cds)
    if not prot:
        return StructureConfidence(source="none", mean_confidence=0.0,
                                   plain_summary="No protein sequence to assess.")
    if plddt:
        conf = [max(0.0, min(100.0, float(v))) for v in plddt][:len(prot)]
        source = "pLDDT (model)"
    else:
        # TOP-IDP -> order confidence: map disorder propensity to a 0-100 order score
        raw = [_TOPIDP.get(a, 0.0) for a in prot]
        conf = [max(0.0, min(100.0, 60.0 - 45.0 * d)) for d in _smooth(raw)]
        source = "sequence proxy (TOP-IDP)"
    mean = round(sum(conf) / len(conf), 0)
    segs = _segments(conf)
    return StructureConfidence(
        source=source, mean_confidence=mean, per_residue=[round(c, 0) for c in conf],
        segments=segs, plain_summary=_plain(mean, segs, source),
        notes=["pLDDT>70 ≈ confident, <50 ≈ likely disordered (AlphaFold convention)."]
        if plddt else ["Order/disorder proxy; higher = more likely to fold. TOP-IDP scale."])
