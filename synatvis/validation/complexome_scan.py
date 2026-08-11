"""Complex-membership matcher — specificity + true-positive confirmation.

Two checks, mirroring the project's other validation legs:

1. SPECIFICITY: run the name matcher over all 5,000 native Cr CDS in the shipped
   corpus. Their headers are NCBI accession IDs (e.g. ``XM_001700409.2``), not
   gene symbols, so a correct matcher should find ~0 matches here. A high count
   would mean the patterns are too loose and firing on noise.
2. TRUE POSITIVES: run the same matcher over ``data/named_complex_genes.fasta``
   -- 13 real genes fetched directly from the Chlamydomonas reinhardtii
   chloroplast genome (NCBI RefSeq NC_005353.1), covering Rubisco, Photosystem
   I, Photosystem II, and ATP synthase. A correct matcher should find a real
   match for every one of these 13.

This is a functionality + specificity check on real data, not a claim that the
6-complex pattern list is exhaustive or that every Cr gene with a matching name
was itself imaged by cryo-ET (see complexome.py's docstring).
"""
from __future__ import annotations

import os
from typing import Dict, List

from ..complexome import identify_complexes, load_complexes
from ..profiles import PACKAGE_DIR
from ..seqio import read_fasta

NATIVE_CORPUS = ["native_cr_cds.fasta", "native_cr_cds_heldout.fasta", "native_cr_cds_4000.fasta"]
NAMED_GENES = "named_complex_genes.fasta"


def run() -> Dict:
    complexes = load_complexes()

    # 1. specificity over the 5,000-gene native corpus
    n_native = fp = 0
    fp_examples: List[str] = []
    for f in NATIVE_CORPUS:
        path = os.path.join(PACKAGE_DIR, "data", f)
        if not os.path.isfile(path):
            continue
        for header, _seq in read_fasta(path):
            name = header.split()[0]
            n_native += 1
            matches = identify_complexes(name, complexes)
            if matches:
                fp += 1
                if len(fp_examples) < 10:
                    fp_examples.append(f"{name} -> {[m.complex_name for m in matches]}")

    # 2. true positives over the 13 real named chloroplast genes
    path = os.path.join(PACKAGE_DIR, "data", NAMED_GENES)
    n_named = tp = 0
    misses: List[str] = []
    hits: Dict[str, str] = {}
    if os.path.isfile(path):
        for header, _seq in read_fasta(path):
            name = header.split()[0]
            n_named += 1
            matches = identify_complexes(name, complexes)
            if matches:
                tp += 1
                hits[name] = matches[0].complex_name
            else:
                misses.append(name)

    return {
        "n_native": n_native, "false_positives": fp, "fp_examples": fp_examples,
        "n_named": n_named, "true_positives": tp, "misses": misses, "hits": hits,
    }


def main() -> int:
    r = run()
    print(f"Complex-membership matcher — specificity + true-positive check")
    print(f"  Specificity: {r['n_native']} native Cr genes (accession-style names) "
          f"-> {r['false_positives']} matched (expected ~0)")
    if r["fp_examples"]:
        print("    examples:", "; ".join(r["fp_examples"]))
    print(f"  True positives: {r['true_positives']}/{r['n_named']} real named chloroplast "
          f"genes correctly matched to their complex")
    if r["misses"]:
        print("    missed:", ", ".join(r["misses"]))
    ok = r["false_positives"] == 0 and r["true_positives"] == r["n_named"]
    print("  => PASS" if ok else "  => CHECK NEEDED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
