"""Leg 1 — specificity from native Cr transcripts (CLAUDE.md §7).

Native, highly-expressed Cr nuclear CDSs are clean negatives BY DEFINITION. Run
the scanner over a few hundred; the per-module medium/high flag rate *is* the
false-positive rate. This is how the TGTAA poly(A) threshold and the GC-trough
threshold get set. The user supplies a Phytozome Cr FASTA; a tiny stub ships so
the harness runs out of the box.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional

from ..flags import Severity
from ..profiles import PACKAGE_DIR, load_profile
from ..scanner import scan
from ..seqio import Transcript, read_fasta

STUB = os.path.join(PACKAGE_DIR, "data", "native_cr_cds.fasta")


def run_specificity(fasta_path: str, profile_name: str = "cr_nuclear") -> Dict:
    profile = load_profile(profile_name)
    records = read_fasta(fasta_path)
    per_module_hits: Dict[str, int] = {}
    per_module_seqs: Dict[str, int] = {}
    n = 0
    for name, seq in records:
        seq = seq.strip()
        if not seq:
            continue
        n += 1
        tx = Transcript(cds=seq, name=name)
        result = scan(tx, profile=profile)
        flagged_modules = set()
        for f in result.flags:
            if f.severity >= Severity.MEDIUM:
                flagged_modules.add(f.module)
        for m in result.module_meta:
            per_module_seqs[m] = per_module_seqs.get(m, 0) + 1
            if m in flagged_modules:
                per_module_hits[m] = per_module_hits.get(m, 0) + 1
    rates = {
        m: (per_module_hits.get(m, 0) / per_module_seqs[m] if per_module_seqs.get(m) else 0.0)
        for m in per_module_seqs
    }
    return {"n_sequences": n, "fp_rate_medium_or_high": rates,
            "hits": per_module_hits, "totals": per_module_seqs}


def main(fasta: Optional[str] = None, profile: str = "cr_nuclear") -> int:
    path = fasta or STUB
    if not os.path.isfile(path):
        print(f"Leg 1: FASTA not found: {path}")
        return 2
    res = run_specificity(path, profile)
    print(f"Leg 1 — specificity over {res['n_sequences']} native Cr CDS "
          f"(profile: {profile})")
    print("  per-module medium/high FALSE-POSITIVE rate (want low):")
    # Deterministic exact-match detectors: a "hit" on a native gene is a REAL
    # motif (e.g. an internal Type IIS site that genuinely needs domestication),
    # so their rate is site prevalence, not a false-positive rate.
    exact_match = {"cloning"}
    for m, r in sorted(res["fp_rate_medium_or_high"].items(),
                       key=lambda kv: -kv[1]):
        bar = "#" * int(round(r * 20))
        kind = "prevalence " if m in exact_match else "FP-rate    "
        print(f"    {m:<14} {r:6.1%}  {kind}{bar}")
    print("  (heuristic modules -> false-positive rate, want low; exact-match")
    print("   modules like 'cloning' -> real site prevalence, expected to be high")
    print("   in native GC-rich genes and NOT a specificity problem.)")
    if res["n_sequences"] < 50:
        print("  NOTE: corpus is tiny — supply more real Cr CDS to set operating points.")
    return 0
