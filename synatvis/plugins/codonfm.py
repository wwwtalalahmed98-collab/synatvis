"""CodonFM plugin — codon-resolution foundation model (NVIDIA BioNeMo).

CodonFM (Darabi et al. 2025; NVIDIA Digital Biology) is a family of codon-resolution
language models trained on ~130M protein-coding sequences from 20,000+ species. It
is the natural successor to SynAT.Vis's hand-crafted codon metrics (RCA / tAI): a
learned model of codon usage, giving a per-sequence "codon fitness" and — via its
Sparse-Autoencoder features — a mechanistic *reason* a codon region looks off, which
fits the tool's explain-itself design. Trained broadly (not Cr-specific), so its
value is realised by fine-tuning on Chlamydomonas (BioNeMo LoRA recipes); until then
it is an opt-in, transfer-learning readout, never part of the validated core.

GPU foundation model — this adapter delegates to the user's own inference command or
BioNeMo NIM endpoint (``CODONFM_CMD``), which reads a CDS FASTA on stdin and prints
``{"codon_fitness": <float>, "log_likelihood": <float>, "sae_note": "<str>"}`` JSON.
"""
from __future__ import annotations

import os
from typing import List

from .base import Plugin, PluginResult, command_available, run_json_command


class CodonFMPlugin(Plugin):
    NAME = "codonfm"
    ENV = "CODONFM_CMD"
    DESCRIPTION = "CodonFM codon foundation model + SAE interpretability (NVIDIA BioNeMo)."
    CITATION = "Darabi et al. 2025 (CodonFM); NVIDIA BioNeMo"

    def available(self) -> bool:
        return command_available(os.environ.get(self.ENV))

    def analyze(self, transcript) -> List[PluginResult]:
        data = run_json_command(os.environ[self.ENV], f">query\n{transcript.cds}\n")
        out: List[PluginResult] = []
        if "codon_fitness" in data:
            v = float(data["codon_fitness"])
            out.append(PluginResult(
                plugin=self.NAME, label="CodonFM codon fitness (foundation model)",
                value=round(v, 3),
                text=f"learned codon-usage fitness {v:.3f} (higher = more natural across "
                     "the model's 20k-species training)",
                note="Codon-resolution foundation model; transfer-learning readout, not "
                     "Cr-calibrated (fine-tune via BioNeMo LoRA). Not part of the validated core.",
                citation=self.CITATION))
        if "log_likelihood" in data:
            out.append(PluginResult(
                plugin=self.NAME, label="CodonFM sequence log-likelihood",
                value=round(float(data["log_likelihood"]), 3),
                text="mean per-codon log-likelihood (naturalness)",
                note="Experimental foundation-model readout.", citation=self.CITATION))
        if data.get("sae_note"):
            out.append(PluginResult(
                plugin=self.NAME, label="CodonFM SAE feature (interpretability)",
                text=str(data["sae_note"]),
                note="Mechanistic explanation from a Sparse Autoencoder over CodonFM features.",
                citation=self.CITATION))
        return out
