"""Calibration leg -- test the expression index against REAL measured expression data.

Every other validation leg in this project checks direction or specificity against
qualitative literature expectations. This leg is different: it checks the index
against actual published numbers (see data/calibration_anchors.yaml).

What it can and cannot establish, stated plainly:
  CAN  -- catch ORDERING failures. If a construct that really expresses better is
          scored worse, that is a genuine defect regardless of scale.
  CANNOT -- verify magnitude. The index is bounded 0-100 and is not a fold-change
          predictor, so it will never reproduce a 16x ratio. Reporting a magnitude
          mismatch as a failure would be misleading.

This leg rests on TWO independent anchor sets (Frontiers 2025 / NanoLuc; Baier 2018
NAR / mVenus) that corroborate each other on direction and magnitude class. That is
stronger than one -- but both probe the SAME axis, introns. Until a second axis is
anchored (codon adaptation, promoter strength, UTR effects), say "directionally
validated, two corroborating anchors on one axis", never "calibrated".
"""
from __future__ import annotations

import os
import re
from typing import Dict, List, Optional

from ..profiles import PACKAGE_DIR, load_yaml

ANCHORS_PATH = os.path.join(PACKAGE_DIR, "data", "calibration_anchors.yaml")
CORPUS_FASTA = os.path.join(PACKAGE_DIR, "data", "construct_grammar",
                            "moclo_corpus", "corpus.fasta")


def load_anchors(path: str = ANCHORS_PATH) -> Dict:
    with open(path, "r", encoding="utf-8") as fh:
        return load_yaml(fh.read()) or {}


def _read_fasta(path: str) -> Dict[str, str]:
    out, cur = {}, None
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                cur = line[1:].split()[0]
                out[cur] = []
            elif cur:
                out[cur].append(line)
    return {k: "".join(v) for k, v in out.items()}


def _excised_intron() -> Optional[str]:
    """The real bTUB2i1 intron, excised from its real MoClo plasmid at BsaI geometry."""
    if not os.path.isfile(CORPUS_FASTA):
        return None
    import importlib.util
    bc_path = os.path.join(PACKAGE_DIR, "data", "construct_grammar", "build_corpus.py")
    spec = importlib.util.spec_from_file_location("_bc", bc_path)
    bc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bc)
    seqs = _read_fasta(CORPUS_FASTA)
    keys = sorted(k for k in seqs if "Intron" in k)
    if not keys:
        return None
    ex = bc.excise_bsai_part(seqs[keys[0]])
    return ex[0] if ex else None


def _base_cds() -> str:
    """The shipped Cr-codon-optimised GFP literature case."""
    cases = open(os.path.join(PACKAGE_DIR, "data", "cases.yaml"), encoding="utf-8").read()
    m = re.search(r"id: gfp_cr_codon_optimized.*?sequence: \"([ACGT]+)\"", cases, re.S)
    return m.group(1) if m else ""


def run() -> Dict:
    from ..seqio import Transcript
    from ..scanner import scan
    from ..expression import predict_expression

    anchors = load_anchors()
    intron = _excised_intron()
    cds = _base_cds()
    if not intron or not cds:
        return {"available": False,
                "reason": "corpus not fetched locally; run "
                          "synatvis/data/construct_grammar/build_corpus.py first"}

    def epi(seq: str) -> float:
        tx = Transcript(cds=seq, name="calib")
        res = scan(tx, profile="cr_nuclear")
        return predict_expression(tx, res.profile, scan_result=res).epi

    scores = {
        0: epi(cds),
        1: epi(cds[:60] + intron + cds[60:]),
        2: epi(cds[:60] + intron + cds[60:400] + intron + cds[400:]),
    }
    measured = {0: 1.0, 1: 6.0, 2: 16.0}  # lower bounds from the anchor

    violations: List[str] = []
    if scores[1] < scores[0]:
        violations.append("score(1 intron) < score(0 introns), but 1 intron measures 6-9x MORE")
    if scores[2] < scores[1]:
        violations.append("score(2 introns) < score(1 intron), but 2 introns measure >16x MORE")

    anchor_list = anchors.get("anchors", [])
    return {"available": True, "scores": scores, "measured_fold": measured,
            "violations": violations, "intron_bp": len(intron),
            "n_anchors": len(anchor_list),
            "citations": [a.get("citation", "") for a in anchor_list],
            "anchor_id": "intron_mediated_enhancement (2 corroborating sets)"}


def main(profile: str = "cr_nuclear") -> int:
    r = run()
    print("Calibration leg -- expression index vs REAL measured expression")
    if not r["available"]:
        print(f"  SKIPPED: {r['reason']}")
        return 0
    print(f"  anchors: {r['n_anchors']} independent set(s), same axis (introns)")
    for c in r["citations"]:
        print(f"    - {c[:104]}")
    print(f"  real excised intron: {r['intron_bp']} bp")
    print()
    print("  introns   index    real measured")
    for n in (0, 1, 2):
        note = "baseline" if n == 0 else f">={r['measured_fold'][n]:.0f}x more protein"
        print(f"    {n}       {r['scores'][n]:5.1f}    {note}")
    print()
    if r["violations"]:
        print(f"  ORDERING VIOLATIONS ({len(r['violations'])}):")
        for v in r["violations"]:
            print(f"    - {v}")
        print("\n  => the index contradicts real measured data. See "
              "data/calibration_anchors.yaml -> known_defect_exposed.")
    else:
        print("  => ordering consistent with the measured anchor.")
    print("\n  NOTE: both anchors probe the SAME axis (introns). Corroborating, but NOT "
          "sufficient to call the index calibrated; magnitude is deliberately not scored.")
    return 1 if r["violations"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
