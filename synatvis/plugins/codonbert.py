"""CodonBERT plugin — mRNA language-model naturalness (opt-in, experimental).

CodonBERT (Li et al. 2024, Genome Research) is a codon-level language model trained
on >10M mRNAs; its pseudo-log-likelihood scores how 'natural' a coding sequence is,
which correlates with expression/stability. This adapter runs only when PyTorch +
Transformers are installed AND a local CodonBERT model is configured via the
``CODONBERT_MODEL`` environment variable — so it never fabricates a score. The
analyze() body is the real integration; it stays dormant until you point it at a
model checkpoint.
"""
from __future__ import annotations

import os
from typing import List

from .base import Plugin, PluginResult


class CodonBERTPlugin(Plugin):
    NAME = "codonbert"
    REQUIRES = ["torch", "transformers"]
    ENV = "CODONBERT_MODEL"
    DESCRIPTION = "CodonBERT mRNA-LM naturalness (pseudo-log-likelihood)."
    CITATION = "Li et al. 2024, Genome Research (CodonBERT)"

    def available(self) -> bool:
        if not os.environ.get(self.ENV):
            return False
        try:
            import torch  # noqa: F401
            import transformers  # noqa: F401
            return True
        except Exception:
            return False

    def analyze(self, transcript) -> List[PluginResult]:
        import math
        import torch
        from transformers import AutoTokenizer, AutoModelForMaskedLM

        model_path = os.environ[self.ENV]
        tok = AutoTokenizer.from_pretrained(model_path)
        model = AutoModelForMaskedLM.from_pretrained(model_path).eval()

        cds = transcript.cds.upper()
        codons = [cds[i:i + 3] for i in range(0, len(cds) - len(cds) % 3, 3)]
        text = " ".join(codons)
        enc = tok(text, return_tensors="pt")
        with torch.no_grad():
            logits = model(**enc).logits
        # mean per-codon log-probability of the observed tokens (pseudo-likelihood)
        logp = torch.log_softmax(logits, dim=-1)
        ids = enc["input_ids"][0]
        scores = [logp[0, i, ids[i]].item() for i in range(1, len(ids) - 1)]
        mean_ll = sum(scores) / len(scores) if scores else float("nan")
        return [PluginResult(
            plugin=self.NAME,
            label="CodonBERT mean codon log-likelihood",
            value=round(mean_ll, 3),
            text=f"mean per-codon log-likelihood {mean_ll:.3f} (higher = more natural)",
            note="Experimental ML prediction; not part of the validated core.",
            citation=self.CITATION,
            validated=False,
        )]
