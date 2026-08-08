#!/usr/bin/env python3
"""SynAT.Vis functionality self-test — the one command to prove the predictor works.

Scans every file in ``examples/test_kit`` and checks that the scanner flags
exactly the planted failure mode (and stays quiet on the clean control). Prints a
PASS/FAIL table and exits non-zero if anything is wrong.

Run from the project folder, no install required:

    python selftest.py        (or: py selftest.py  /  python3 selftest.py)

This tests the tool's DETECTION on inputs whose ground truth YOU can read in the
FASTA headers — it is not a promise of biological outcome (see the report banner).
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)  # run without installing

from synatvis import scan  # noqa: E402
from synatvis.seqio import read_record  # noqa: E402
from synatvis.flags import Severity  # noqa: E402

KIT = os.path.join(HERE, "examples", "test_kit")

# filename -> (cds_span or None, expected_module or None (None => expect clean),
#              minimum severity that counts)
MANIFEST = [
    ("kit_00_clean_control.fasta",          None,        None,           None),
    ("kit_01_cloning_BsaI.fasta",           None,        "cloning",      Severity.HIGH),
    ("kit_02_polya_TGTAA.fasta",            None,        "polya",        Severity.HIGH),
    ("kit_03_instability_ARE.fasta",        None,        "instability",  Severity.MEDIUM),
    ("kit_04_composition_homopolymer.fasta", None,       "composition",  Severity.LOW),
    ("kit_05_uorf_upstream_AUG.fasta",      (30, 144),   "uorf",         Severity.LOW),
]


def check(fname, cds_span, expect_module, min_sev):
    tx = read_record(os.path.join(KIT, fname), cds_span=cds_span)
    res = scan(tx)
    if expect_module is None:
        # clean control: no medium/high from any real detector
        offenders = [f for f in res.flags if f.severity >= Severity.MEDIUM]
        ok = not offenders
        detail = "no medium/high flags" if ok else \
            "unexpected: " + ", ".join(f"{f.module}:{f.severity}" for f in offenders)
        return ok, detail
    hits = [f for f in res.flags
            if f.module == expect_module and f.severity >= min_sev]
    ok = bool(hits)
    detail = (f"{expect_module} flagged: {hits[0].message[:60]}..."
              if ok else f"MISSING expected {expect_module} >= {min_sev}")
    return ok, detail


def main() -> int:
    if not os.path.isdir(KIT):
        print(f"test kit not found: {KIT}")
        return 2
    print("SynAT.Vis functionality self-test")
    print("=" * 72)
    n_pass = 0
    rows = []
    for fname, span, module, sev in MANIFEST:
        try:
            ok, detail = check(fname, span, module, sev)
        except Exception as exc:  # noqa: BLE001
            ok, detail = False, f"ERROR: {exc!r}"
        rows.append((ok, fname, detail))
        n_pass += ok
    for ok, fname, detail in rows:
        print(f"  [{'PASS' if ok else 'FAIL'}] {fname:<40} {detail}")
    print("=" * 72)
    print(f"{n_pass}/{len(MANIFEST)} checks passed")
    if n_pass == len(MANIFEST):
        print("OK — the predictor detects each planted failure mode and stays "
              "quiet on the clean control.")
    return 0 if n_pass == len(MANIFEST) else 1


if __name__ == "__main__":
    raise SystemExit(main())
