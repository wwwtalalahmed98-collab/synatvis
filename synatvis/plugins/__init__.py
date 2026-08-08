"""Tier-B opt-in ML plugins (CLAUDE.md §1 — opt-in, never a default).

The validated core never imports this package. Plugins are advisory, experimental,
and run only when their (heavy) dependencies are installed. Use ``run_available``
to collect predictions from whatever is installed, or ``status`` to list them.
"""
from __future__ import annotations

from typing import Dict, List

from .base import Plugin, PluginResult
from .codonbert import CodonBERTPlugin
from .lineardesign import LinearDesignPlugin
from .saluki import SalukiPlugin
from .utr_lm import UTRLMPlugin
# Tier-1 frontier readouts (opt-in, transfer-learning, not Cr-calibrated)
from .borzoi import BorzoiPlugin
from .aparent2 import APARENT2Plugin
from .icodon import ICodonPlugin
# Tier-3 validated protein-fate / splicing tools (upgrade the ptm.py / splice heuristics)
from .signalp import SignalP6Plugin
from .deeploc import DeepLoc2Plugin
from .deeptmhmm import DeepTMHMMPlugin
from .pangolin import PangolinPlugin
# NVIDIA BioNeMo foundation-model seams (GPU / NIM microservices; opt-in, dormant until wired)
from .codonfm import CodonFMPlugin
from .openfold3 import OpenFold3Plugin
from .evo2 import Evo2Plugin

REGISTRY: List[Plugin] = [
    LinearDesignPlugin(),
    CodonBERTPlugin(),
    SalukiPlugin(),
    UTRLMPlugin(),
    # Tier-1 frontier readouts
    BorzoiPlugin(),
    APARENT2Plugin(),
    ICodonPlugin(),
    # Tier-3 validated fate / splicing tools
    SignalP6Plugin(),
    DeepLoc2Plugin(),
    DeepTMHMMPlugin(),
    PangolinPlugin(),
    # NVIDIA BioNeMo foundation models
    CodonFMPlugin(),
    OpenFold3Plugin(),
    Evo2Plugin(),
]


def status() -> List[Dict]:
    """List every plugin with availability and install hint."""
    return [
        {"name": p.NAME, "available": p.available(), "description": p.DESCRIPTION,
         "requires": p.REQUIRES, "env": p.ENV, "hint": p.install_hint(),
         "citation": p.CITATION}
        for p in REGISTRY
    ]


def available_plugins() -> List[Plugin]:
    return [p for p in REGISTRY if p.available()]


def run_available(transcript) -> List[PluginResult]:
    out: List[PluginResult] = []
    for p in REGISTRY:
        if not p.available():
            continue
        try:
            out.extend(p.analyze(transcript))
        except Exception as exc:  # a plugin failure must never break the scan
            out.append(PluginResult(plugin=p.NAME, label="plugin error",
                                    note=f"{type(exc).__name__}: {exc}", validated=False))
    return out


__all__ = ["Plugin", "PluginResult", "REGISTRY", "status", "available_plugins",
           "run_available"]
