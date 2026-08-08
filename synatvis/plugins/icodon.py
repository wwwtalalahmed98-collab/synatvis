"""iCodon plugin — codon-composition mRNA stability prediction.

iCodon (Diez et al. 2022, NAR Genomics & Bioinformatics) predicts mRNA stability
from codon composition using species-specific coefficients. IMPORTANT limitation:
iCodon ships coefficients for human, mouse, zebrafish and Xenopus only — NOT
Chlamydomonas — so its number is a cross-species proxy and must be read as such.
It is opt-in and clearly labelled. Delegates to the user's own iCodon wrapper
(``ICODON_CMD``; e.g. a small Rscript) that reads a FASTA on stdin and prints
``{"stability": <float>, "species": "<name>"}`` as JSON.
"""
from __future__ import annotations

import os
from typing import List

from .base import Plugin, PluginResult, command_available, run_json_command


class ICodonPlugin(Plugin):
    NAME = "icodon"
    ENV = "ICODON_CMD"
    DESCRIPTION = "iCodon codon-based mRNA stability (species proxy; no Cr model)."
    CITATION = "Diez et al. 2022, NAR Genom. Bioinform. (iCodon)"

    def available(self) -> bool:
        return command_available(os.environ.get(self.ENV))

    def analyze(self, transcript) -> List[PluginResult]:
        fasta = f">query\n{transcript.cds}\n"
        data = run_json_command(os.environ[self.ENV], fasta)
        if "stability" not in data:
            return []
        val = float(data["stability"])
        sp = data.get("species", "?")
        return [PluginResult(
            plugin=self.NAME, label="iCodon predicted mRNA stability",
            value=round(val, 3),
            text=f"predicted stability {val:.3f} (model species: {sp})",
            note=f"Experimental cross-species proxy — iCodon has no Chlamydomonas "
                 f"model, so this is borrowed from {sp}. Not part of the validated core.",
            citation=self.CITATION)]
