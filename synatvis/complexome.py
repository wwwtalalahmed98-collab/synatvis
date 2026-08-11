"""Multiprotein complex-membership context (opt-in, name-based, never fabricated).

Chlamydomonas proteins rarely work alone -- most real biological machines are
teams of proteins (a "complex"). A 2025 large-scale cryo-electron tomography
study (Chromatin-Structure-Rhythms-Lab community dataset, EMPIAR-11830)
directly photographed intact Chlamydomonas cells in 3D and confirmed that
~25 real macromolecular complexes physically exist and assemble as imaged --
not just as annotated gene lists.

What this module does NOT do: identify a complex from raw DNA sequence alone.
There is no reliable way to do that without a genome-scale structural/BLAST
search, which is out of scope for this stdlib-only tool. Instead, it matches
the transcript's declared NAME (from a FASTA header or an explicit gene
symbol) against a curated list of well-established gene-name patterns for
known complex subunits (see data/complexes.yaml). A match means "this gene's
name is a known member of complex X's subunit family, and complex X's real
3D assembly has been directly imaged in this organism" -- a citable identity
+ context hint, never a claim that THIS specific molecule was itself imaged.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional

from .profiles import PACKAGE_DIR, load_yaml

COMPLEXES_PATH = os.path.join(PACKAGE_DIR, "data", "complexes.yaml")


@dataclass
class ComplexMatch:
    complex_name: str
    function: str
    matched_pattern: str
    identity_citation: str
    structural_citation: str


def load_complexes(path: str = COMPLEXES_PATH) -> List[Dict]:
    """Load the curated complex/gene-pattern list (see data/complexes.yaml)."""
    with open(path, "r", encoding="utf-8") as fh:
        data = load_yaml(fh.read())
    return (data or {}).get("complexes", [])


def identify_complexes(name: str, complexes: Optional[List[Dict]] = None) -> List[ComplexMatch]:
    """Match a gene/transcript NAME against known complex-subunit gene patterns.

    Name-based, not sequence-based: only fires when the caller supplies a real
    gene symbol (e.g. a FASTA header like ">rbcL chloroplast CDS"). Never
    guesses a gene identity from raw DNA content.
    """
    if not name:
        return []
    complexes = complexes if complexes is not None else load_complexes()
    out: List[ComplexMatch] = []
    for c in complexes:
        for pat in c.get("gene_patterns", []):
            if re.search(pat, name):
                out.append(ComplexMatch(
                    complex_name=c["name"], function=c["function"],
                    matched_pattern=pat,
                    identity_citation=c["identity_citation"],
                    structural_citation=c["structural_citation"]))
                break
    return out
