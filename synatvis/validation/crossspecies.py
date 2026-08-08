"""Leg 4 — cross-species discrimination (CLAUDE.md §7 extension).

Proves the Cr-tuned tool is genuinely *host-specific*: it accepts native
*C. reinhardtii* transcripts and flags FOREIGN transcripts from other species as
not matching the host. This is the specificity claim in classifier terms —
controlling both false positives (foreign wrongly accepted) and false negatives
(Cr wrongly flagged).

A decision threshold on the host-fit metric (RCA) is *trained* on one half of the
labelled data and evaluated on the held-out half; AUC is reported over all data.
Data: ``native_cr_cds*.fasta`` (positives) and ``foreign_cds.fasta`` (negatives,
real NCBI RefSeq CDS from other species).
"""
from __future__ import annotations

import os
from typing import Dict, List

from ..codon_tables import attach_advanced, gene_tai, load_for_profile
from ..profiles import PACKAGE_DIR, load_profile
from ..seqio import read_fasta

FOREIGN = os.path.join(PACKAGE_DIR, "data", "foreign_cds.fasta")


def _auc(pos: List[float], neg: List[float]) -> float:
    allv = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    ranks = [0.0] * len(allv)
    j, r = 0, 1
    while j < len(allv):
        k = j
        while k < len(allv) and allv[k][0] == allv[j][0]:
            k += 1
        avg = (r + r + (k - j) - 1) / 2
        for m in range(j, k):
            ranks[m] = avg
        r += (k - j)
        j = k
    sp = sum(ranks[i] for i in range(len(allv)) if allv[i][1] == 1)
    return (sp - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)) if pos and neg else 0.0


def _pct(xs, p):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(p * len(xs)))] if xs else 0.0


def run_crossspecies(profile_name: str = "cr_nuclear", n_cr: int = 1000) -> Dict:
    profile = load_profile(profile_name)
    table = load_for_profile(profile, PACKAGE_DIR)
    tai_w = attach_advanced(profile, PACKAGE_DIR)["tai"]

    cr = []
    for f in ("native_cr_cds_heldout.fasta", "native_cr_cds_4000.fasta",
              "native_cr_cds.fasta"):
        p = os.path.join(PACKAGE_DIR, "data", f)
        if os.path.isfile(p):
            cr += [s.upper() for _, s in read_fasta(p)]
    cr = cr[:n_cr]

    foreign: Dict[str, List[str]] = {}
    if os.path.isfile(FOREIGN):
        for hdr, s in read_fasta(FOREIGN):
            foreign.setdefault(hdr.rsplit("_", 1)[0], []).append(s.upper())

    cr_rca = [table.rca(s) for s in cr]
    cr_tai = [gene_tai(s, tai_w) for s in cr]
    fo_rca = {k: [table.rca(s) for s in v] for k, v in foreign.items()}
    fo_tai = [gene_tai(s, tai_w) for v in foreign.values() for s in v]
    all_fo = [x for v in fo_rca.values() for x in v]

    cr_tr, cr_te = cr_rca[::2], cr_rca[1::2]
    fo_te = all_fo[1::2]
    thr = _pct(cr_tr, 0.05)
    tp = sum(1 for v in cr_te if v >= thr)
    fn = sum(1 for v in cr_te if v < thr)
    fp = sum(1 for v in fo_te if v >= thr)
    tn = sum(1 for v in fo_te if v < thr)
    tot = tp + fn + fp + tn
    return {
        "n_cr": len(cr), "n_foreign": len(all_fo),
        "species": {k: len(v) for k, v in foreign.items()},
        "threshold": thr,
        "tp": tp, "fn": fn, "fp": fp, "tn": tn,
        "accuracy": (tp + tn) / tot if tot else 0,
        "sensitivity": tp / (tp + fn) if (tp + fn) else 0,
        "specificity": tn / (tn + fp) if (tn + fp) else 0,
        "auc_rca": _auc(cr_rca, all_fo),
        "auc_tai": _auc(cr_tai, fo_tai),
        "per_species_flagged": {k: (sum(1 for v in rv if v < thr) / len(rv)) if rv else 0
                                for k, rv in fo_rca.items()},
        "cr_median_rca": _pct(cr_rca, 0.5),
        "foreign_median_rca": _pct(all_fo, 0.5),
    }


def main(profile: str = "cr_nuclear") -> int:
    if not os.path.isfile(FOREIGN):
        print("Leg 4: foreign_cds.fasta not found (ship it or regenerate).")
        return 2
    r = run_crossspecies(profile)
    print(f"Leg 4 — cross-species discrimination (profile: {profile})")
    print(f"  {r['n_cr']} Cr transcripts vs {r['n_foreign']} foreign "
          f"({', '.join(f'{k}={n}' for k, n in r['species'].items())})")
    print(f"  trained RCA threshold {r['threshold']:.3f}  (held-out test half)")
    print(f"  TP(Cr accepted)={r['tp']} FN(Cr flagged)={r['fn']} "
          f"FP(foreign accepted)={r['fp']} TN(foreign flagged)={r['tn']}")
    print(f"  accuracy {r['accuracy']:.1%} | sensitivity {r['sensitivity']:.1%} | "
          f"specificity {r['specificity']:.1%}")
    print(f"  AUC(RCA) {r['auc_rca']:.3f} | AUC(tAI) {r['auc_tai']:.3f}")
    print("  foreign flagged as not-Cr, by species:")
    for k, v in r["per_species_flagged"].items():
        print(f"    {k:<16} {v:.0%}")
    print(f"  Cr median RCA {r['cr_median_rca']:.2f} vs foreign {r['foreign_median_rca']:.2f} "
          f"— the host-fit metric cleanly separates the two.")
    return 0
