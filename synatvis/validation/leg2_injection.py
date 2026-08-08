"""Leg 2 — sensitivity by synthetic injection (CLAUDE.md §7).

Insert a known Cr motif (TGTAA + G-rich context; a Type IIS site; an ARE) into a
clean parent, confirm the right module flags it AND stays silent on the parent.
This tests DETECTION, not biological consequence — labelled as such.

The 'clean parent' is a synthetic GC-rich in-frame ORF built to trip nothing.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, List

from ..flags import Severity
from ..profiles import load_profile
from ..scanner import scan
from ..seqio import Transcript

# A GC-rich, in-frame parent CDS using preferred Cr codons; designed to be clean.
# (Met + a run of Ala/Gly/Leu/Glu preferred codons + stop.)
CLEAN_PARENT = "ATG" + ("GCCGGCCTGGAGAAGGCC" * 8) + "TAA"


@dataclass
class InjectionCase:
    name: str
    expect_module: str
    build: Callable[[str], str]   # parent -> mutated CDS


def _inject(parent: str, motif: str, at_codon: int) -> str:
    """Splice *motif* into *parent* at a codon boundary (keeps frame length-wise)."""
    pos = at_codon * 3
    return parent[:pos] + motif + parent[pos:]


CASES: List[InjectionCase] = [
    InjectionCase(
        "polya_TGTAA_in_Grich",
        "polya",
        # TGTAA followed by a G/C-rich context so it reads as a functional Cr NUE
        lambda p: _inject(p, "TGTAA" + "GGCCGGCGGCCGCGGCCGGC", 4),
    ),
    InjectionCase(
        "cloning_BsaI",
        "cloning",
        lambda p: _inject(p, "GGTCTC", 6),
    ),
    InjectionCase(
        "instability_ARE_cluster",
        "instability",
        lambda p: _inject(p, "ATTTA" + "CG" + "ATTTA", 5),
    ),
]


def run_injection(profile_name: str = "cr_nuclear") -> Dict:
    profile = load_profile(profile_name)
    parent_tx = Transcript(cds=CLEAN_PARENT, name="clean_parent")
    parent_result = scan(parent_tx, profile=profile)
    parent_modules_flagged = {f.module for f in parent_result.flags
                              if f.severity >= Severity.MEDIUM}

    results = []
    for case in CASES:
        mutated = case.build(CLEAN_PARENT)
        tx = Transcript(cds=mutated, name=case.name)
        res = scan(tx, profile=profile)
        detected = any(f.module == case.expect_module and f.severity >= Severity.LOW
                       for f in res.flags)
        parent_silent = case.expect_module not in parent_modules_flagged
        results.append({
            "case": case.name,
            "expect_module": case.expect_module,
            "detected": detected,
            "parent_silent": parent_silent,
            "pass": detected and parent_silent,
        })
    return {"parent_modules_flagged": sorted(parent_modules_flagged),
            "cases": results}


def main(profile: str = "cr_nuclear") -> int:
    res = run_injection(profile)
    print(f"Leg 2 — sensitivity by injection (profile: {profile}) "
          f"[tests DETECTION, not biological consequence]")
    print(f"  clean parent medium/high modules: "
          f"{res['parent_modules_flagged'] or 'none'}")
    all_pass = True
    for c in res["cases"]:
        status = "PASS" if c["pass"] else "FAIL"
        all_pass &= c["pass"]
        print(f"    [{status}] {c['case']:<26} module={c['expect_module']:<12} "
              f"detected={c['detected']} parent_silent={c['parent_silent']}")
    return 0 if all_pass else 1
