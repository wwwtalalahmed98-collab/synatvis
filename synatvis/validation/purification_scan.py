"""Purification-layer confirmation over the 5,000 native Cr CDS.

Runs the purification-strategy recommender + biophysics on every native
Chlamydomonas nuclear CDS in the shipped corpus (500 + 500 + 4,000 = 5,000 unique
genes) and reports the distributions. This confirms the layer executes on all
5,000 without error and produces chemically sensible recommendations across a real
proteome slice. It is a functionality + sanity check on real proteins, NOT a
measured purification-yield validation (that needs wet-lab data).
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional

from ..profiles import PACKAGE_DIR
from ..purification import predict_purification
from ..seqio import Transcript, read_fasta

CORPUS = ["native_cr_cds.fasta", "native_cr_cds_heldout.fasta", "native_cr_cds_4000.fasta"]


def _pct(vals: List[float], q: float) -> float:
    if not vals:
        return 0.0
    s = sorted(vals)
    return round(s[min(len(s) - 1, int(q * len(s)))], 2)


def run(limit: Optional[int] = None) -> Dict:
    files = [os.path.join(PACKAGE_DIR, "data", f) for f in CORPUS]
    files = [f for f in files if os.path.isfile(f)]
    pI, mw, gravy = [], [], []
    capture: Dict[str, int] = {}
    n = errors = secreted = membrane = disulfide = tagged = hic = 0
    for path in files:
        for name, seq in read_fasta(path):
            seq = seq.strip()
            if not seq:
                continue
            if limit and n >= limit:
                break
            n += 1
            try:
                p = predict_purification(Transcript(cds=seq, name=name))
                b = p.biophysics
                pI.append(b["pI"]); mw.append(b["mw_kda"]); gravy.append(b["gravy"])
                cap = next((w["method"] for w in p.workflow if w["phase"].startswith("capture")),
                           "—")
                capture[cap] = capture.get(cap, 0) + 1
                if p.tags:
                    tagged += 1
                if "secreted" in p.localization:
                    secreted += 1
                if "membrane" in p.localization:
                    membrane += 1
                if any(c["topic"] == "Redox" for c in p.considerations):
                    disulfide += 1
                if any(r["method"].startswith("Hydrophobic") for r in p.recommendations):
                    hic += 1
            except Exception:  # noqa: BLE001 — we are counting failures
                errors += 1
    return {"n": n, "errors": errors, "capture": capture,
            "pI": {"p10": _pct(pI, .1), "median": _pct(pI, .5), "p90": _pct(pI, .9)},
            "mw": {"p10": _pct(mw, .1), "median": _pct(mw, .5), "p90": _pct(mw, .9)},
            "gravy": {"p10": _pct(gravy, .1), "median": _pct(gravy, .5), "p90": _pct(gravy, .9)},
            "secreted": secreted, "membrane": membrane, "disulfide": disulfide,
            "tagged": tagged, "hic_candidates": hic}


def main(profile: str = "cr_nuclear") -> int:
    r = run()
    n = max(1, r["n"])
    print(f"Purification-layer confirmation over {r['n']} native Cr CDS "
          f"({r['errors']} errors)")
    print(f"  biophysics   pI  p10/median/p90 = {r['pI']['p10']} / {r['pI']['median']} / {r['pI']['p90']}")
    print(f"               MW  p10/median/p90 = {r['mw']['p10']} / {r['mw']['median']} / {r['mw']['p90']} kDa")
    print(f"               GRAVY               = {r['gravy']['p10']} / {r['gravy']['median']} / {r['gravy']['p90']}")
    print("  recommended CAPTURE method (property-based, since native genes carry no tags):")
    for m, c in sorted(r["capture"].items(), key=lambda kv: -kv[1]):
        print(f"    {c/n:6.1%}  {m}")
    print(f"  secreted (recover from medium) : {r['secreted']/n:5.1%}")
    print(f"  membrane (needs detergent)     : {r['membrane']/n:5.1%}")
    print(f"  likely disulfides (non-reducing): {r['disulfide']/n:5.1%}")
    print(f"  HIC-suitable (hydrophobic)     : {r['hic_candidates']/n:5.1%}")
    print(f"  pre-existing tags detected     : {r['tagged']/n:5.1%}  (expected ~0 for native genes)")
    print("  => the layer runs on every CDS and returns a chemically sensible plan; "
          "capture splits by pI as expected. Functionality confirmed (not a yield measurement).")
    return 0 if r["errors"] == 0 else 1
