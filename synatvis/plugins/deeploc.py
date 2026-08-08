"""DeepLoc 2.1 plugin — protein subcellular localisation from a protein LM.

DeepLoc 2.1 (Ødum et al. 2024) predicts subcellular localisation from a protein
language model — the validated upgrade to SynAT.Vis's rule-of-thumb localisation in
ptm.py. Operates on the PROTEIN, so this adapter translates the CDS first. Academic
tool (check terms for commercial use). Delegates to the user's installed DeepLoc
command (``DEEPLOC_CMD``) that reads a protein FASTA on stdin and prints
``{"localization": "<name>", "probability": <float>}`` as JSON.
"""
from __future__ import annotations

import os
from typing import List

from .base import Plugin, PluginResult, command_available, run_json_command


class DeepLoc2Plugin(Plugin):
    NAME = "deeploc2"
    ENV = "DEEPLOC_CMD"
    DESCRIPTION = "DeepLoc 2.1 subcellular localisation (validated; upgrades ptm.py)."
    CITATION = "Odum et al. 2024 (DeepLoc 2.1)"

    def available(self) -> bool:
        return command_available(os.environ.get(self.ENV))

    def analyze(self, transcript) -> List[PluginResult]:
        from ..ptm import translate
        prot = translate(transcript.cds)
        if not prot:
            return []
        data = run_json_command(os.environ[self.ENV], f">query\n{prot}\n")
        loc = data.get("localization", "?")
        prob = float(data.get("probability", 0.0))
        return [PluginResult(
            plugin=self.NAME, label="DeepLoc 2.1 localisation",
            value=round(prob, 3),
            text=f"{loc} (p={prob:.3f})",
            note="Validated protein-LM predictor; refines the ptm.py localisation call "
                 "when installed. Note Cr organelle biology differs from the training set.",
            citation=self.CITATION)]
