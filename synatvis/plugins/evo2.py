"""Evo 2 plugin — genomic foundation model & generative design (NVIDIA BioNeMo).

Evo 2 (Arc Institute; StripedHyena2; hosted in NVIDIA BioNeMo) is a DNA foundation
model trained across the tree of life. For SynAT.Vis it reads the *whole cassette*
as DNA: a sequence log-likelihood that captures pre-transcriptional / regulatory
plausibility, variant-effect scoring for proposed edits, and — the longer game — a
generative engine to redesign a Cr-conditioned cassette rather than only flag one.
Broadly trained, not Cr-specific, so it is an opt-in transfer-learning readout until
fine-tuned; never part of the validated core.

GPU model — delegates to the user's own inference command or BioNeMo NIM endpoint
(``EVO2_CMD``), which reads a DNA FASTA on stdin and prints
``{"log_likelihood": <float>, "variant_effect": <float>, "generative": <bool>}`` JSON.
"""
from __future__ import annotations

import os
from typing import List

from .base import Plugin, PluginResult, command_available, run_json_command


class Evo2Plugin(Plugin):
    NAME = "evo2"
    ENV = "EVO2_CMD"
    DESCRIPTION = "Evo 2 genomic foundation model: likelihood, variant effect, design (BioNeMo)."
    CITATION = "Brixi, Nguyen et al. 2025 (Evo 2, Arc Institute); NVIDIA BioNeMo"

    def available(self) -> bool:
        return command_available(os.environ.get(self.ENV))

    def analyze(self, transcript) -> List[PluginResult]:
        seq = transcript.utr5 + transcript.cds + transcript.utr3
        data = run_json_command(os.environ[self.ENV], f">query\n{seq}\n")
        out: List[PluginResult] = []
        if "log_likelihood" in data:
            out.append(PluginResult(
                plugin=self.NAME, label="Evo 2 sequence log-likelihood",
                value=round(float(data["log_likelihood"]), 3),
                text="genome-model plausibility of the whole cassette (higher = more natural)",
                note="Transfer-learning readout, not Cr-calibrated. Not part of the validated core.",
                citation=self.CITATION))
        if "variant_effect" in data:
            out.append(PluginResult(
                plugin=self.NAME, label="Evo 2 variant-effect score",
                value=round(float(data["variant_effect"]), 3),
                text=f"predicted effect of the proposed edit(s): {float(data['variant_effect']):.3f}",
                note="Experimental; use to rank silent-fix candidates.", citation=self.CITATION))
        if data.get("generative"):
            out.append(PluginResult(
                plugin=self.NAME, label="Evo 2 generative redesign available",
                text="the configured endpoint can propose a Cr-conditioned cassette redesign",
                note="Generative design is a roadmap capability; review any generated sequence.",
                citation=self.CITATION))
        return out
