"""Leg 3 — literature cases, gated by COMPARTMENT (CLAUDE.md §7).

Each case in ``cases.yaml`` is tagged ``{organism, compartment, lineage, module,
expect_flag, mechanism, citation, sequence}``. GATING RULE: a nuclear-specific
detector is validated on nuclear cases only; a chloroplast case is gated OUT of
nuclear-module scoring. Sensitivity/specificity are reported STRATIFIED by
compartment so the gating is visible to a reviewer.

The seed corpus is a SCAFFOLD — the user populates it, preferring single-variable
RESCUE PAIRS (one element changed, expression restored) over orphan failures.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional

from ..flags import Severity
from ..profiles import PACKAGE_DIR, load_profile, load_yaml
from ..scanner import scan
from ..seqio import Transcript

SEED = os.path.join(PACKAGE_DIR, "data", "cases.yaml")


def _load_cases(path: str) -> List[Dict]:
    with open(path, "r", encoding="utf-8") as fh:
        data = load_yaml(fh.read())
    if isinstance(data, dict) and "cases" in data:
        data = data["cases"]
    return [c for c in (data or []) if isinstance(c, dict)]


def run_cases(cases_path: str, profile_name: str = "cr_nuclear") -> Dict:
    profile = load_profile(profile_name)
    active_compartment = profile.get("meta", {}).get("compartment")

    # stratify: compartment -> {tp, fn, tn, fp, gated}
    strata: Dict[str, Dict[str, int]] = {}

    for case in _load_cases(cases_path):
        comp = case.get("compartment", "unknown")
        strata.setdefault(comp, {"tp": 0, "fn": 0, "tn": 0, "fp": 0, "gated": 0})
        seq = (case.get("sequence") or "").strip()
        if not seq:
            continue
        # GATING: only score cases whose compartment matches the active profile
        if comp != active_compartment:
            strata[comp]["gated"] += 1
            continue

        tx = Transcript(cds=seq, name=case.get("id", "case"))
        result = scan(tx, profile=profile)
        module = case.get("module")
        flagged = any(f.module == module and f.severity >= Severity.LOW
                      for f in result.flags)
        expect = bool(case.get("expect_flag"))

        if expect and flagged:
            strata[comp]["tp"] += 1
        elif expect and not flagged:
            strata[comp]["fn"] += 1
        elif not expect and not flagged:
            strata[comp]["tn"] += 1
        else:
            strata[comp]["fp"] += 1

    return {"active_compartment": active_compartment, "strata": strata}


def _rate(num: int, den: int) -> str:
    return f"{(num / den):.0%}" if den else "n/a"


def main(cases: Optional[str] = None, profile: str = "cr_nuclear") -> int:
    path = cases or SEED
    if not os.path.isfile(path):
        print(f"Leg 3: cases file not found: {path}")
        return 2
    res = run_cases(path, profile)
    print(f"Leg 3 — literature cases, gated by compartment "
          f"(active profile compartment: {res['active_compartment']})")
    for comp, s in res["strata"].items():
        gate_note = " [GATED OUT of scoring]" if comp != res["active_compartment"] else ""
        sens = _rate(s["tp"], s["tp"] + s["fn"])
        spec = _rate(s["tn"], s["tn"] + s["fp"])
        print(f"  compartment={comp}{gate_note}")
        print(f"    scored: tp={s['tp']} fn={s['fn']} tn={s['tn']} fp={s['fp']} "
              f"gated={s['gated']}")
        if comp == res["active_compartment"]:
            print(f"    sensitivity={sens}  specificity={spec}")
    print("  Corpus: real sequences (NCBI accessions in cases.yaml). Expand with "
          "more single-variable rescue pairs to tighten the sens/spec estimate.")
    return 0
