"""Expression-propensity prediction — a transparent ENSEMBLE of models (opt-in).

The validated diagnostic core emits no expression score. This module is the
separate, opt-in predictor the design contract reserves (CLAUDE.md §1). Rather than
trust one metric, it combines several published, mechanistically-distinct models,
each anchored to the distribution of 5,000 well-expressed native Cr genes and shown
separately so the estimate is never a black box:

  * codon adaptation (RCA)      — translation efficiency & mRNA stability in Cr
                                  (Sharp & Li 1987; Barahimipour 2015; Presnyak 2015)
  * tRNA-pool adaptation (tAI)  — elongation efficiency (dos Reis 2004)
  * silencing resistance (GC)   — low GC drives silencing in Cr (Barahimipour 2015)
  * structural stability (ΔG)   — folded mRNA is longer-lived; the LinearDesign
                                  objective (Zhang 2020) — active when ViennaRNA is present
  * ML readouts (opt-in)        — CodonBERT / Saluki / UTR-LM, shown when installed

The ensemble index is the weighted mean of the AVAILABLE anchored models, then
multiplied by hard hazard gates (premature poly(A), silencing GC trough, uORF,
missing start). It is a MODEL, not a measured yield — uncalibrated pending
experimental data, but directionally validated (ranks Cr-optimised rescue
partners above native; native above foreign).
"""
from __future__ import annotations

import bisect
import json
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .codon_tables import STANDARD_CODE, attach_advanced, gene_tai, load_for_profile
from .flags import Severity
from .structure_energy import fold, has_vienna
from .util import gc_fraction

# anchored models: key -> (label, default weight, reference)
_MODELS = {
    "rca": ("codon adaptation -> translation & mRNA stability", 0.30,
            "Sharp & Li 1987; Barahimipour 2015; Presnyak 2015"),
    "tai": ("tRNA-pool adaptation -> elongation", 0.30, "dos Reis 2004"),
    "gc":  ("GC content -> silencing resistance", 0.20, "Barahimipour 2015"),
    "dg":  ("structural stability (mRNA folding dG)", 0.20, "Zhang 2020 (LinearDesign)"),
}


@dataclass
class ExpressionResult:
    epi: float
    band: str
    models: List[Dict] = field(default_factory=list)     # anchored ensemble members
    ml_readouts: List[Dict] = field(default_factory=list)  # opt-in ML predictions (shown, not blended)
    hazards: List[Dict] = field(default_factory=list)
    landscape: List[float] = field(default_factory=list)   # per-codon adaptiveness
    confidence: str = ""
    fate: Optional[Dict] = None                            # protein fate / PTM layer

    def to_dict(self) -> Dict:
        return {"expression_propensity_index": round(self.epi, 1), "band": self.band,
                "models": self.models, "ml_readouts": self.ml_readouts,
                "hazards": self.hazards, "n_landscape": len(self.landscape),
                "protein_fate": self.fate, "confidence": self.confidence,
                "is_model_not_measured": True}


def load_reference(base_dir: str) -> Optional[Dict]:
    path = os.path.join(base_dir, "data", "cr_expression_reference.json")
    if not os.path.isfile(path):
        return None
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def _percentile(value: float, arr: List[float]) -> float:
    if not arr:
        return 50.0
    return float(min(100, max(0, bisect.bisect_right(arr, value))))


def predict_expression(transcript, profile: Dict, scan_result=None,
                       run_ml: bool = True) -> ExpressionResult:
    base = profile.get("_base_dir", ".")
    ref = load_reference(base)
    table = load_for_profile(profile, base)
    tai = attach_advanced(profile, base)["tai"]
    pct = ref["percentiles"] if ref else {}
    weights = dict((k, v[1]) for k, v in _MODELS.items())
    weights.update(profile.get("expression", {}).get("weights", {}) or {})

    cds = transcript.cds
    # raw per-model scores (0-100), only for available models
    scores: Dict[str, float] = {}
    scores["rca"] = _percentile(table.rca(cds), pct.get("rca", [])) if pct else 50.0
    scores["tai"] = _percentile(gene_tai(cds, tai), pct.get("tai", [])) if pct else 50.0
    scores["gc"] = _percentile(gc_fraction(cds), pct.get("gc", [])) if pct else 50.0
    dg_available = has_vienna() and pct.get("dg") and len(cds) >= 60
    if dg_available:
        mfe, _pf, _be = fold(cds[:300])
        # more negative ΔG/nt = more structured/stable = higher score (invert percentile)
        scores["dg"] = 100.0 - _percentile(mfe / len(cds[:300]), pct["dg"])

    models, acc, wsum = [], 0.0, 0.0
    for k, (label, _w, refc) in _MODELS.items():
        avail = k in scores
        w = weights.get(k, 0.0)
        models.append({"name": k, "label": label, "score": round(scores.get(k, 0.0), 1),
                       "weight": w, "ref": refc, "available": avail,
                       "kind": "anchored"})
        if avail:
            acc += w * scores[k]
            wsum += w
    epi = acc / wsum if wsum else 50.0

    # hard hazard gates from the diagnostic scan
    if scan_result is None:
        from .scanner import scan
        scan_result = scan(transcript, profile=profile)
    hazards: List[Dict] = []
    seen = set()

    def gate(module, factor, reason, cond):
        if module in seen:
            return
        for f in scan_result.flags:
            if f.module == module and f.severity >= Severity.MEDIUM and cond(f):
                hazards.append({"factor": factor, "reason": reason})
                seen.add(module)
                return

    gate("polya", 0.35, "premature poly(A) signal — risk of a truncated transcript",
         lambda f: f.region in ("cds", "5utr") and f.detail.get("context_ok"))
    gate("composition", 0.70, "GC trough — heterochromatin/silencing risk",
         lambda f: f.detail.get("hazard") == "low_gc")
    gate("uorf", 0.85, "upstream AUG competes with the main start", lambda f: True)
    if not cds.upper().startswith("ATG"):
        hazards.append({"factor": 0.15, "reason": "CDS does not begin with ATG"})

    # protein-fate / PTM layer (complementary; drives the cell view routing)
    from .ptm import predict_protein_fate
    fate = predict_protein_fate(transcript)
    if fate.modifier != 1.0:
        hazards.append({"factor": fate.modifier, "reason": fate.modifier_reason})

    for h in hazards:
        epi *= h["factor"]
    epi = max(0.0, min(100.0, epi))

    # opt-in ML model readouts (shown transparently; NOT blended into the index)
    ml_readouts: List[Dict] = []
    if run_ml:
        try:
            from . import plugins
            for p in plugins.REGISTRY:
                if p.NAME == "lineardesign":
                    continue  # its ΔG is already the structural-stability model
                if p.available():
                    for r in p.analyze(transcript):
                        ml_readouts.append({"plugin": r.plugin, "label": r.label,
                                            "value": r.value, "text": r.text,
                                            "citation": r.citation})
        except Exception:
            pass

    if epi >= 65:
        band = "sequence features consistent with STRONG expression"
    elif epi >= 40:
        band = "MODERATE — some sequence-level drag on expression"
    else:
        band = "LIKELY POOR — redesign recommended"

    landscape = [table.weight.get(cds[k:k + 3], 0.0)
                 for k in range(0, len(cds) - len(cds) % 3, 3)]
    nmods = sum(1 for m in models if m["available"]) + len(ml_readouts)
    conf = ("Ensemble of {} model(s), anchored to {} well-expressed native Cr genes. "
            "A MODEL, not a measured yield: uncalibrated pending experimental data. "
            "Directionally validated (Cr-optimised > native, native > foreign)."
            ).format(nmods, ref["n"] if ref else "?")

    return ExpressionResult(epi=epi, band=band, models=models, ml_readouts=ml_readouts,
                            hazards=hazards, landscape=landscape, confidence=conf,
                            fate=fate.to_dict())
