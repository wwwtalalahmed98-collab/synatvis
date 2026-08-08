"""Tier-B plugin contract (opt-in ML models).

The validated core of SynAT.Vis is rule-based and stdlib-only. Tier-B plugins wrap
heavy, separately-installed models (CodonBERT, Saluki, UTR-LM, LinearDesign) that
predict things the rules approximate — mRNA naturalness, half-life, ribosome load,
structural stability. They are **opt-in and EXPERIMENTAL**: never imported by the
core, never part of the validated result, and clearly labelled unvalidated. A
plugin runs only when its dependency (and, for the ML models, a configured local
model) is present; otherwise it reports itself as unavailable with an install hint.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
from dataclasses import dataclass
from typing import Dict, List, Optional


def _first_token(cmd: str) -> str:
    """First token of a command line, honouring a leading quoted path."""
    cmd = cmd.strip()
    if cmd and cmd[0] in "\"'":
        end = cmd.find(cmd[0], 1)
        return cmd[1:end] if end > 0 else cmd[1:]
    return cmd.split()[0] if cmd else ""


def command_available(env_value: Optional[str]) -> bool:
    """True if ``env_value`` names a runnable command (its first token resolves)."""
    if not env_value:
        return False
    first = _first_token(env_value)
    return bool(first and (shutil.which(first) or os.path.exists(first)))


def run_json_command(cmd: str, input_text: str, timeout: int = 300) -> Dict:
    """Run a user-configured inference command, feed FASTA on stdin, parse JSON stdout.

    This is how the external-model adapters stay honest: the user points an env var
    at *their own* installed model's inference wrapper (any language), which must
    read a FASTA on stdin and print a single JSON object to stdout. We never
    fabricate a score — if the command is absent the plugin is simply unavailable.
    """
    proc = subprocess.run(cmd, shell=True, input=input_text, capture_output=True,
                          text=True, timeout=timeout)
    if proc.returncode != 0:
        raise RuntimeError((proc.stderr or "non-zero exit").strip()[:300])
    out = (proc.stdout or "").strip()
    # tolerate a wrapper that prints logs before the JSON line
    start = out.find("{")
    if start > 0:
        out = out[start:]
    return json.loads(out)


@dataclass
class PluginResult:
    """One prediction from a Tier-B plugin. Always advisory, never validated."""

    plugin: str
    label: str
    value: Optional[float] = None
    text: str = ""
    note: str = ""
    citation: str = ""
    validated: bool = False  # Tier-B is opt-in / experimental by definition


class Plugin:
    """Base class. Subclasses set NAME/REQUIRES and implement available()/analyze()."""

    NAME: str = "plugin"
    REQUIRES: List[str] = []          # pip packages the plugin needs
    ENV: Optional[str] = None         # env var pointing at a local model, if any
    DESCRIPTION: str = ""
    CITATION: str = ""

    def available(self) -> bool:
        return False

    def analyze(self, transcript) -> List[PluginResult]:
        return []

    def install_hint(self) -> str:
        parts = []
        if self.REQUIRES:
            parts.append("pip install " + " ".join(self.REQUIRES))
        if self.ENV:
            placeholder = ("<inference-command reading FASTA on stdin, printing JSON>"
                           if self.ENV.endswith("_CMD") else "<path-to-model>")
            parts.append(f"set {self.ENV}={placeholder}")
        return "  and  ".join(parts)
