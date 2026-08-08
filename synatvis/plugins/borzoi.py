"""Borzoi plugin — sequence -> RNA-seq coverage (expression, splicing, poly(A)).

Borzoi (Linder, Kelley et al. 2024/2025, Calico) predicts base-resolution RNA-seq
coverage from DNA sequence over long (~0.5 Mb) context, jointly capturing
expression level, splicing and 3' cleavage/polyadenylation. It is trained on
human/mouse — a transfer-learning readout, NOT Cr-calibrated — so it is opt-in and
clearly labelled experimental. Because Borzoi needs genomic context and a specific
inference stack, this adapter delegates to the user's own configured inference
command (``BORZOI_CMD``): a wrapper that reads a FASTA on stdin and prints
``{"expression": <float>, ...}`` as JSON. It stays dormant until you wire it.
"""
from __future__ import annotations

import os
from typing import List

from .base import Plugin, PluginResult, command_available, run_json_command


class BorzoiPlugin(Plugin):
    NAME = "borzoi"
    ENV = "BORZOI_CMD"
    DESCRIPTION = "Borzoi seq->RNA-seq: expression / splicing / poly(A) (human-trained)."
    CITATION = "Linder, Kelley et al. 2024/2025 (Borzoi, Calico)"

    def available(self) -> bool:
        return command_available(os.environ.get(self.ENV))

    def analyze(self, transcript) -> List[PluginResult]:
        fasta = f">query\n{transcript.cds}\n"
        data = run_json_command(os.environ[self.ENV], fasta)
        out: List[PluginResult] = []
        if "expression" in data:
            out.append(PluginResult(
                plugin=self.NAME, label="Borzoi predicted expression (RNA-seq)",
                value=round(float(data["expression"]), 3),
                text=f"predicted expression track {float(data['expression']):.3f} "
                     "(human-trained; relative, not Cr-calibrated)",
                note="Experimental; needs genomic context. Transfer-learning readout, "
                     "not part of the validated core.",
                citation=self.CITATION))
        for key, label in (("polya", "Borzoi poly(A)/3' cleavage signal"),
                           ("splice", "Borzoi splice-usage signal")):
            if key in data:
                out.append(PluginResult(
                    plugin=self.NAME, label=label, value=round(float(data[key]), 3),
                    text=f"{label}: {float(data[key]):.3f}",
                    note="Experimental transfer-learning readout.", citation=self.CITATION))
        return out
