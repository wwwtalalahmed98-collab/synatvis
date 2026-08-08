"""Stage 0 — ground-truth corpus admission gate for the construct-grammar recognizer.

This module enforces, in code, the criteria frozen in
``data/construct_grammar/inclusion_criteria.yaml`` (human rationale in the
sibling ``INCLUSION_CRITERIA.md``). A :class:`CandidatePart` is a proposed
addition to the training corpus; :func:`evaluate_candidate` is the single gate
every candidate must pass through before Stage 1 (part segmentation) sees it.

Criterion text and citations live in the YAML, not here, so the checkable logic
and the citable rationale can never silently drift apart — this mirrors how
``THRESHOLDS.md`` documents the provenance of every numeric threshold used
elsewhere in the scanner.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

from .profiles import PACKAGE_DIR, load_yaml

CRITERIA_PATH = os.path.join(PACKAGE_DIR, "data", "construct_grammar", "inclusion_criteria.yaml")

_ACCESSION_TYPES_OK = (
    "addgene_plasmid", "genbank_accession",
    "chlamycollection_catalog_entry", "paper_supplementary_sequence_file",
)
_ACCESSION_TYPES_NOT_OK = ("supplementary_table_only", "figure_only", "described_only")


def load_criteria(path: str = CRITERIA_PATH) -> Dict[str, Any]:
    """Load the frozen criteria document (IC-1..IC-5, EX-1..EX-4 + rationale)."""
    with open(path, "r", encoding="utf-8") as fh:
        data = load_yaml(fh.read())
    return data or {}


class Verdict(str, Enum):
    INCLUDE = "INCLUDE"
    CANDIDATE_TIER_ONLY = "CANDIDATE_TIER_ONLY"
    PENDING = "PENDING"
    EXCLUDE = "EXCLUDE"


@dataclass
class CandidatePart:
    """A proposed ground-truth corpus entry, prior to admission.

    Field names mirror ``candidate_record_schema`` in ``inclusion_criteria.yaml``.
    """

    part_id: str
    so_term: Optional[str] = None
    assembly_sites: List[str] = field(default_factory=list)  # e.g. ["BsaI", "BsmBI"]
    fusion_overhangs: Optional[Tuple[str, str]] = None
    validated_hosts: List[str] = field(default_factory=list)  # e.g. ["chlamydomonas_reinhardtii_nuclear"]
    syntax_compliant_only: bool = False  # True if validated only via Phytobrick syntax, not Cr itself
    sequence: str = ""
    accession: str = ""
    accession_type: str = "described_only"
    citation: str = ""
    functional_evidence: str = ""
    sequence_conflict: bool = False  # EX-2: supplementary vs. deposited sequence disagree
    assembly_standard: str = "type_iis"  # "type_iis" | "gateway" | "other"
    is_preprint_only: bool = False  # EX-1: no peer-reviewed record and no deposit


@dataclass
class CriterionCheck:
    id: str
    passed: bool
    reason: str


@dataclass
class CriteriaResult:
    part_id: str
    verdict: Verdict
    checks: List[CriterionCheck]

    def failed(self) -> List[CriterionCheck]:
        return [c for c in self.checks if not c.passed]


def _check_ic1(p: CandidatePart) -> CriterionCheck:
    type_iis = {"bsai", "bsmbi", "bpii"}
    has_sites = any(s.lower() in type_iis for s in p.assembly_sites)
    has_overhangs = p.fusion_overhangs is not None
    ok = has_sites or has_overhangs
    reason = "Type IIS site or Phytobrick overhang present" if ok else \
        "no recognized Type IIS site and no Phytobrick fusion overhang recorded"
    return CriterionCheck("IC-1", ok, reason)


def _check_ic2(p: CandidatePart) -> CriterionCheck:
    """Passes only for direct Cr validation. Syntax-only (non-Cr) validation is
    real evidence for Tier A (junction grammar) but must not count as passing
    IC-2 for primary ground truth -- it routes to CANDIDATE_TIER_ONLY instead.
    """
    cr_hosts = {h for h in p.validated_hosts if "chlamydomonas" in h.lower()}
    ok = bool(cr_hosts)
    if ok:
        reason = "validated directly in Chlamydomonas"
    elif p.syntax_compliant_only and p.validated_hosts:
        reason = "validated only in a Phytobrick-compliant non-Cr system (syntax-tier only)"
    else:
        reason = "no recorded functional validation host"
    return CriterionCheck("IC-2", ok, reason)


def _check_ic3(p: CandidatePart) -> CriterionCheck:
    ok = p.accession_type in _ACCESSION_TYPES_OK and bool(p.accession) and bool(p.sequence)
    if ok:
        reason = f"primary source deposited ({p.accession_type}: {p.accession})"
    elif p.accession_type in _ACCESSION_TYPES_NOT_OK:
        reason = f"only a {p.accession_type.replace('_', ' ')} — not an addressable primary record"
    else:
        reason = "missing accession and/or sequence"
    return CriterionCheck("IC-3", ok, reason)


def _check_ic4(p: CandidatePart) -> CriterionCheck:
    ok = bool(p.so_term)
    reason = f"SO term assigned ({p.so_term})" if ok else "no Sequence Ontology term assignable"
    return CriterionCheck("IC-4", ok, reason)


def _check_ic5(p: CandidatePart) -> CriterionCheck:
    ok = bool(p.citation.strip())
    reason = "independently citable" if ok else "no DOI/PMID/Addgene citation recorded"
    return CriterionCheck("IC-5", ok, reason)


def _check_ex1(p: CandidatePart) -> CriterionCheck:
    violated = p.is_preprint_only and p.accession_type not in _ACCESSION_TYPES_OK
    reason = "preprint-only, no deposited sequence" if violated else "not preprint-only"
    return CriterionCheck("EX-1", not violated, reason)


def _check_ex2(p: CandidatePart) -> CriterionCheck:
    reason = "sequence conflict between supplementary and deposited records" if p.sequence_conflict \
        else "no recorded sequence conflict"
    return CriterionCheck("EX-2", not p.sequence_conflict, reason)


def _check_ex3(p: CandidatePart) -> CriterionCheck:
    violated = p.assembly_standard.lower() not in ("type_iis",)
    reason = f"assembly standard '{p.assembly_standard}' is not Type-IIS/Phytobrick-compatible" \
        if violated else "Type-IIS compatible"
    return CriterionCheck("EX-3", not violated, reason)


def _check_ex4(p: CandidatePart) -> CriterionCheck:
    has_function = bool(p.functional_evidence.strip())
    reason = "no functional evidence recorded (candidate-tier only)" if not has_function \
        else "functional evidence recorded"
    return CriterionCheck("EX-4", has_function, reason)


_CHECKS = (_check_ic1, _check_ic2, _check_ic3, _check_ic4, _check_ic5,
           _check_ex1, _check_ex2, _check_ex3, _check_ex4)


def evaluate_candidate(part: CandidatePart) -> CriteriaResult:
    """Evaluate *part* against every frozen IC/EX criterion.

    Verdict logic (checked in this order):
      EXCLUDE              -- fails EX-2 (conflicting sequence) or EX-3 (non-Type-IIS standard):
                               these are unconditional, independent of deposit status.
      PENDING               -- EX-1 case: preprint-only with no deposited sequence yet.
                               This is exactly why it lacks a primary source, so it takes
                               priority over an IC-3 failure rather than compounding into EXCLUDE.
      EXCLUDE               -- otherwise fails a hard gate: IC-1, IC-3, IC-4, or IC-5.
      CANDIDATE_TIER_ONLY   -- passes structural/provenance gates but fails IC-2 (validated
                               outside Cr only) or EX-4 (no functional evidence at all).
      INCLUDE               -- passes every criterion; eligible for primary ground truth.
    """
    checks = [chk(part) for chk in _CHECKS]
    by_id = {c.id: c for c in checks}

    if not by_id["EX-2"].passed or not by_id["EX-3"].passed:
        verdict = Verdict.EXCLUDE
    elif not by_id["EX-1"].passed:
        verdict = Verdict.PENDING
    elif any(not by_id[i].passed for i in ("IC-1", "IC-3", "IC-4", "IC-5")):
        verdict = Verdict.EXCLUDE
    elif not by_id["IC-2"].passed or not by_id["EX-4"].passed:
        verdict = Verdict.CANDIDATE_TIER_ONLY
    else:
        verdict = Verdict.INCLUDE

    return CriteriaResult(part_id=part.part_id, verdict=verdict, checks=checks)
