"""APARENT2 plugin — deep-learning poly(A) site strength / 3' cleavage.

APARENT2 (Linder & Seelig) is a residual CNN trained on massively-parallel poly(A)
reporter assays; it predicts cleavage-and-polyadenylation strength far more
accurately than a motif rule. It is the natural upgrade to SynAT.Vis's TGTAA
poly(A) module: where the rule flags a candidate signal, APARENT2 scores how
strong/usable that site actually is. Trained on human data (transfer readout, not
Cr-calibrated), so it is opt-in and experimental. Delegates to the user's own
inference command (``APARENT2_CMD``) that reads a FASTA on stdin and prints
``{"polya_strength": <float 0-1>}`` as JSON.
"""
from __future__ import annotations

import os
from typing import List

from .base import Plugin, PluginResult, command_available, run_json_command


class APARENT2Plugin(Plugin):
    NAME = "aparent2"
    ENV = "APARENT2_CMD"
    DESCRIPTION = "APARENT2 poly(A) site strength from a CNN on poly(A) MPRA (human-trained)."
    CITATION = "Linder & Seelig (APARENT2)"

    def available(self) -> bool:
        return command_available(os.environ.get(self.ENV))

    def analyze(self, transcript) -> List[PluginResult]:
        fasta = f">query\n{transcript.cds}\n"
        data = run_json_command(os.environ[self.ENV], fasta)
        if "polya_strength" not in data:
            return []
        s = float(data["polya_strength"])
        return [PluginResult(
            plugin=self.NAME, label="APARENT2 poly(A) site strength",
            value=round(s, 3),
            text=f"strongest internal poly(A) site score {s:.3f} (0-1; higher = more "
                 "likely premature cleavage)",
            note="Experimental; complements the TGTAA rule. Human-trained readout, "
                 "not part of the validated core.",
            citation=self.CITATION)]
