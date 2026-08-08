"""OpenFold3 plugin — structure & co-folding (NVIDIA BioNeMo Fold-CP).

OpenFold3 / BioNeMo Fold-CP (NVIDIA Digital Biology, 2026) predicts biomolecular
structure and complexes with per-residue confidence (pLDDT) and predicted TM-score
(pTM). It gives SynAT.Vis's post-translation layer a real folding / QC axis and,
crucially, feeds the self-explanatory structure interpreter (structure_confidence.py)
with actual pLDDT instead of the TOP-IDP sequence proxy. Operates on the PROTEIN, so
this adapter translates the CDS first.

GPU model — delegates to the user's own inference command or BioNeMo NIM endpoint
(``OPENFOLD3_CMD``), which reads a protein FASTA on stdin and prints
``{"mean_plddt": <0-100>, "ptm": <0-1>, "n_low_confidence": <int>}`` JSON.
"""
from __future__ import annotations

import os
from typing import List

from .base import Plugin, PluginResult, command_available, run_json_command


class OpenFold3Plugin(Plugin):
    NAME = "openfold3"
    ENV = "OPENFOLD3_CMD"
    DESCRIPTION = "OpenFold3 structure / co-folding, pLDDT + pTM (NVIDIA BioNeMo Fold-CP)."
    CITATION = "OpenFold3 / NVIDIA BioNeMo Fold-CP 2026"

    def available(self) -> bool:
        return command_available(os.environ.get(self.ENV))

    def analyze(self, transcript) -> List[PluginResult]:
        from ..ptm import translate
        prot = translate(transcript.cds)
        if not prot:
            return []
        data = run_json_command(os.environ[self.ENV], f">query\n{prot}\n")
        out: List[PluginResult] = []
        if "mean_plddt" in data:
            p = float(data["mean_plddt"])
            nlc = data.get("n_low_confidence")
            tail = f"; {nlc} low-confidence region(s)" if nlc is not None else ""
            out.append(PluginResult(
                plugin=self.NAME, label="OpenFold3 fold confidence (pLDDT)",
                value=round(p, 1),
                text=f"mean pLDDT {p:.1f}/100{tail} (>70 confident, <50 likely disordered)",
                note="Real structure prediction; feeds the plain-language structure interpreter. "
                     "Opt-in, not part of the validated core.",
                citation=self.CITATION))
        if "ptm" in data:
            out.append(PluginResult(
                plugin=self.NAME, label="OpenFold3 predicted TM-score (pTM)",
                value=round(float(data["ptm"]), 3),
                text=f"pTM {float(data['ptm']):.3f} (global fold reliability)",
                note="Experimental structure-model readout.", citation=self.CITATION))
        return out
