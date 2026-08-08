"""Pangolin plugin — deep-learning splice-site prediction (multi-species).

Pangolin (Zeng & Li 2022, Genome Biology) predicts splice-site usage across four
species and generally outperforms SpliceAI on non-human sequence — the validated
upgrade to SynAT.Vis's cryptic-splice heuristic. Splicing matters in Cr because
introns are assets but cryptic sites are hazards. Pangolin needs a model + (for
variant scoring) genomic context, so this adapter delegates to the user's own
inference command (``PANGOLIN_CMD``) that reads a FASTA on stdin and prints
``{"max_splice_score": <float>, "n_cryptic_sites": <int>}`` as JSON. Opt-in,
experimental (trained on human/mouse/etc., not Cr).
"""
from __future__ import annotations

import os
from typing import List

from .base import Plugin, PluginResult, command_available, run_json_command


class PangolinPlugin(Plugin):
    NAME = "pangolin"
    ENV = "PANGOLIN_CMD"
    DESCRIPTION = "Pangolin splice-site prediction (validated; upgrades the splice heuristic)."
    CITATION = "Zeng & Li 2022, Genome Biology (Pangolin)"

    def available(self) -> bool:
        return command_available(os.environ.get(self.ENV))

    def analyze(self, transcript) -> List[PluginResult]:
        seq = transcript.utr5 + transcript.cds + transcript.utr3
        data = run_json_command(os.environ[self.ENV], f">query\n{seq}\n")
        out: List[PluginResult] = []
        if "max_splice_score" in data:
            s = float(data["max_splice_score"])
            n = int(data.get("n_cryptic_sites", 0))
            out.append(PluginResult(
                plugin=self.NAME, label="Pangolin splice-site usage",
                value=round(s, 3),
                text=f"strongest predicted splice site {s:.3f}; {n} candidate cryptic site(s)",
                note="Validated splice predictor; refines the cryptic-splice heuristic when "
                     "installed. Trained on human/mouse/etc., not Cr — read as a proxy.",
                citation=self.CITATION))
        return out
