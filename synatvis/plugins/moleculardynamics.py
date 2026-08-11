"""Molecular dynamics plugin -- opt-in slot for a REAL, separately-installed MD engine.

Molecular dynamics (MD) simulates how a protein's atoms move over time (e.g.
GROMACS, OpenMM, NAMD, AMBER). It is a fundamentally different kind of
computation from everything else in SynAT.Vis: it needs a full 3D starting
structure, a force field, and heavy compute (real MD runs take hours-to-weeks
per system on GPUs/clusters, even for small proteins). No simulation is ever
"perfect" -- every MD run is an approximation with a real, known error margin,
true for every lab, not a limitation specific to this tool.

This tool does NOT run MD itself -- that would mean either faking a result or
silently shipping a fragile in-process simulator no one asked to trust. Instead,
following the same honest pattern as every other Tier-B/Tier-3 plugin in this
package, this is a thin adapter to *your own* separately installed MD software.
Nothing runs, and no numbers appear, unless MDSIM_CMD points at a real command
of yours that reads a protein FASTA on stdin and prints:
    {"rmsd_nm": <float>, "radius_of_gyration_nm": <float>,
     "sim_time_ns": <float>, "force_field": <str>, "temperature_k": <float>}
"""
from __future__ import annotations

import os
from typing import List

from .base import Plugin, PluginResult, command_available, run_json_command


class MolecularDynamicsPlugin(Plugin):
    NAME = "moleculardynamics"
    ENV = "MDSIM_CMD"
    DESCRIPTION = ("Opt-in slot for a real, separately-installed MD engine "
                   "(GROMACS/OpenMM/NAMD/AMBER, your own). Never simulates internally.")
    CITATION = "user-supplied MD engine; no internal simulation is performed"

    def available(self) -> bool:
        return command_available(os.environ.get(self.ENV))

    def analyze(self, transcript) -> List[PluginResult]:
        from ..ptm import translate
        prot = translate(transcript.cds)
        if not prot:
            return []
        data = run_json_command(os.environ[self.ENV], f">query\n{prot}\n")
        out: List[PluginResult] = []
        if "rmsd_nm" in data:
            ff = data.get("force_field", "unspecified force field")
            t = data.get("sim_time_ns")
            tail = f" over {t} ns" if t is not None else ""
            out.append(PluginResult(
                plugin=self.NAME, label="MD backbone RMSD",
                value=round(float(data["rmsd_nm"]), 3),
                text=f"RMSD {float(data['rmsd_nm']):.3f} nm{tail} ({ff})",
                note="From YOUR configured MD engine, not run by SynAT.Vis. Every MD result is "
                     "an approximation (force-field- and sampling-dependent) -- never treat a "
                     "single run as a final answer.",
                citation=self.CITATION))
        if "radius_of_gyration_nm" in data:
            out.append(PluginResult(
                plugin=self.NAME, label="MD radius of gyration",
                value=round(float(data["radius_of_gyration_nm"]), 3),
                text=f"Rg {float(data['radius_of_gyration_nm']):.3f} nm (compactness over the run)",
                note="From YOUR configured MD engine.", citation=self.CITATION))
        return out
