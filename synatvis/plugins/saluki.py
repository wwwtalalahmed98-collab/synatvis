"""Saluki plugin — mRNA half-life prediction (opt-in, experimental).

Saluki (Agarwal & Kelley 2022, Genome Biology) predicts mRNA half-life from
sequence (coding frame + splice sites) with a hybrid CNN/RNN. This adapter runs
only when TensorFlow is installed AND a local Saluki model is configured via
``SALUKI_MODEL``; otherwise it is unavailable. analyze() is the integration point.
"""
from __future__ import annotations

import os
from typing import List

from .base import Plugin, PluginResult


class SalukiPlugin(Plugin):
    NAME = "saluki"
    REQUIRES = ["tensorflow"]
    ENV = "SALUKI_MODEL"
    DESCRIPTION = "Saluki mRNA half-life prediction from sequence."
    CITATION = "Agarwal & Kelley 2022, Genome Biology (Saluki)"

    def available(self) -> bool:
        if not os.environ.get(self.ENV):
            return False
        try:
            import tensorflow  # noqa: F401
            return True
        except Exception:
            return False

    def analyze(self, transcript) -> List[PluginResult]:
        import numpy as np
        import tensorflow as tf

        model = tf.keras.models.load_model(os.environ[self.ENV])
        # Saluki encodes the mRNA as a 6-track one-hot (4 nt + coding-frame + 5' cap);
        # this builds the minimal nt + first-codon-position tracks from the CDS.
        seq = transcript.full.upper()
        onehot = np.zeros((len(seq), 6), dtype="float32")
        base = {"A": 0, "C": 1, "G": 2, "T": 3}
        cds0 = transcript.cds_start
        for i, ch in enumerate(seq):
            if ch in base:
                onehot[i, base[ch]] = 1.0
            if i >= cds0 and (i - cds0) % 3 == 0:
                onehot[i, 4] = 1.0  # coding-frame track
        pred = float(model.predict(onehot[None, ...], verbose=0).ravel()[0])
        return [PluginResult(
            plugin=self.NAME,
            label="Saluki predicted mRNA half-life (relative)",
            value=round(pred, 3),
            text=f"predicted stability score {pred:.3f} (higher = longer-lived mRNA)",
            note="Experimental ML prediction; not part of the validated core. Track "
                 "encoding may need aligning to your Saluki checkpoint.",
            citation=self.CITATION,
            validated=False,
        )]
