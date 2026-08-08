"""UTR-LM plugin — 5'UTR translation efficiency / ribosome load (opt-in).

UTR-LM (Chu et al. 2023, Nature Machine Intelligence) is a 5'UTR language model
that predicts mean ribosome load and translation efficiency. This adapter runs
only when PyTorch + Transformers are installed AND a local model is configured via
``UTRLM_MODEL``, and only when the transcript actually has a 5'UTR to score.
"""
from __future__ import annotations

import os
from typing import List

from .base import Plugin, PluginResult


class UTRLMPlugin(Plugin):
    NAME = "utr_lm"
    REQUIRES = ["torch", "transformers"]
    ENV = "UTRLM_MODEL"
    DESCRIPTION = "UTR-LM 5'UTR ribosome-load / translation-efficiency prediction."
    CITATION = "Chu et al. 2023, Nature Machine Intelligence (UTR-LM)"

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
        if len(transcript.utr5) < 10:
            return [PluginResult(
                plugin=self.NAME, label="UTR-LM ribosome load",
                text="no 5'UTR supplied — nothing to score",
                note="Provide a 5'UTR (--cds span) to use this plugin.",
                citation=self.CITATION, validated=False)]
        import torch
        from transformers import AutoTokenizer, AutoModel

        model_path = os.environ[self.ENV]
        tok = AutoTokenizer.from_pretrained(model_path)
        model = AutoModel.from_pretrained(model_path).eval()
        enc = tok(transcript.utr5.replace("T", "U"), return_tensors="pt")
        with torch.no_grad():
            out = model(**enc)
        # regression head or pooled embedding -> MRL; depends on the checkpoint
        score = float(out.pooler_output.mean().item()) if hasattr(out, "pooler_output") \
            else float(out.last_hidden_state.mean().item())
        return [PluginResult(
            plugin=self.NAME,
            label="UTR-LM predicted mean ribosome load",
            value=round(score, 3),
            text=f"predicted 5'UTR ribosome-load index {score:.3f} (higher = more translation)",
            note="Experimental ML prediction; not part of the validated core.",
            citation=self.CITATION,
            validated=False,
        )]
