"""DeepTMHMM plugin — transmembrane topology from a protein LM.

DeepTMHMM (Hallgren et al. 2022) is the validated deep-learning successor to TMHMM
for transmembrane-topology prediction — the upgrade to the Kyte-Doolittle TM
heuristic in ptm.py. Operates on the PROTEIN, so this adapter translates the CDS
first. Delegates to the user's installed command (``DEEPTMHMM_CMD``; e.g. a
pybiolib wrapper) that reads a protein FASTA on stdin and prints
``{"tm_count": <int>, "topology": "<string>"}`` as JSON.
"""
from __future__ import annotations

import os
from typing import List

from .base import Plugin, PluginResult, command_available, run_json_command


class DeepTMHMMPlugin(Plugin):
    NAME = "deeptmhmm"
    ENV = "DEEPTMHMM_CMD"
    DESCRIPTION = "DeepTMHMM transmembrane topology (validated; upgrades ptm.py TM heuristic)."
    CITATION = "Hallgren et al. 2022 (DeepTMHMM)"

    def available(self) -> bool:
        return command_available(os.environ.get(self.ENV))

    def analyze(self, transcript) -> List[PluginResult]:
        from ..ptm import translate
        prot = translate(transcript.cds)
        if not prot:
            return []
        data = run_json_command(os.environ[self.ENV], f">query\n{prot}\n")
        n = int(data.get("tm_count", 0))
        topo = data.get("topology", "")
        topo_txt = f"; topology {topo}" if topo else ""
        return [PluginResult(
            plugin=self.NAME, label="DeepTMHMM transmembrane helices",
            value=float(n),
            text=f"{n} predicted TM helix/helices{topo_txt}",
            note="Validated topology predictor; replaces the Kyte-Doolittle TM heuristic "
                 "in ptm.py when installed.",
            citation=self.CITATION)]
