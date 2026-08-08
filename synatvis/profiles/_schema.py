"""Profile schema validation (CLAUDE.md §4).

A profile must supply every key the modules read. Validating up front means a
module never has to defend against a missing constant at scan time, and adding a
new compartment/lineage profile fails loudly if it forgets a key.
"""
from __future__ import annotations

from typing import Any, Dict, List, Tuple

# section -> required keys
REQUIRED: Dict[str, List[str]] = {
    "meta": ["host", "compartment", "gating_axis", "validated"],
    "codon": ["table", "rare_weight_threshold", "cluster_window_codons",
              "cluster_min_rare", "low_rca_warn"],
    "composition": ["target_gc", "genome_gc", "window", "step", "gc_trough_low",
                    "gc_trough_warn", "min_trough_len", "homopolymer_min"],
    "polya": ["nue_motifs", "upstream_min", "upstream_max", "gc_context_window",
              "gc_context_min"],
    "splice": ["donor_consensus", "acceptor_consensus", "min_intron",
               "max_intron", "inserted_intron_name"],
    "structure": ["start_context_window", "utr5_max_len", "utr5_struct_window",
                  "max_pairing_frac"],
    "uorf": ["flag_any_uaug", "strong_context_flag"],
    "cloning": ["enzymes", "scope"],
    "instability": ["are_motifs", "are_cluster_min", "are_window"],
    "silencing": ["validated", "guidance"],
}


class ProfileError(ValueError):
    """Raised when a profile is missing a key a module needs."""


def validate(profile: Dict[str, Any]) -> Tuple[bool, List[str]]:
    """Return ``(ok, problems)``; ``ok`` is True iff ``problems`` is empty."""
    problems: List[str] = []
    for section, keys in REQUIRED.items():
        if section not in profile or not isinstance(profile[section], dict):
            problems.append(f"missing section: {section!r}")
            continue
        for key in keys:
            if key not in profile[section] or profile[section][key] is None:
                # empty string / empty list are allowed (used by the plastid stub)
                if key in profile[section] and profile[section][key] in ("", []):
                    continue
                problems.append(f"missing key: {section}.{key}")
    return (not problems), problems


def require_valid(profile: Dict[str, Any]) -> Dict[str, Any]:
    ok, problems = validate(profile)
    if not ok:
        raise ProfileError("invalid profile:\n  " + "\n  ".join(problems))
    return profile
