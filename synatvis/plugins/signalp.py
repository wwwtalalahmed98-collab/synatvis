"""SignalP 6.0 plugin — transformer-based signal-peptide prediction.

SignalP 6.0 (Teufel et al. 2022, Nature Biotechnology) is the field-standard
signal-peptide predictor and the validated upgrade to SynAT.Vis's own hydrophobicity
heuristic in ptm.py. It operates on the PROTEIN, so this adapter translates the CDS
first. Licensing note: SignalP is free for academic use but COMMERCIAL use needs a
paid licence from DTU — relevant if SynAT.Vis is commercialised. Delegates to the
user's installed SignalP command (``SIGNALP6_CMD``) that reads a protein FASTA on
stdin and prints ``{"signal_peptide": <bool>, "probability": <float>,
"cleavage_site": <int>}`` as JSON.
"""
from __future__ import annotations

import os
from typing import List

from .base import Plugin, PluginResult, command_available, run_json_command


class SignalP6Plugin(Plugin):
    NAME = "signalp6"
    ENV = "SIGNALP6_CMD"
    DESCRIPTION = "SignalP 6.0 signal-peptide prediction (validated; upgrades ptm.py heuristic)."
    CITATION = "Teufel et al. 2022, Nat. Biotechnol. (SignalP 6.0)"

    def available(self) -> bool:
        return command_available(os.environ.get(self.ENV))

    def analyze(self, transcript) -> List[PluginResult]:
        from ..ptm import translate
        prot = translate(transcript.cds)
        if not prot:
            return []
        data = run_json_command(os.environ[self.ENV], f">query\n{prot}\n")
        sp = bool(data.get("signal_peptide"))
        prob = float(data.get("probability", 0.0))
        cs = data.get("cleavage_site")
        cs_txt = f", cleavage after residue {cs}" if cs else ""
        return [PluginResult(
            plugin=self.NAME, label="SignalP 6.0 signal peptide",
            value=round(prob, 3),
            text=f"{'signal peptide' if sp else 'no signal peptide'} "
                 f"(p={prob:.3f}{cs_txt})",
            note="Validated protein-level predictor; replaces the ptm.py heuristic when "
                 "installed. Academic-free; commercial use needs a DTU licence.",
            citation=self.CITATION)]
