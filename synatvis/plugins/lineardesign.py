"""LinearDesign objective plugin — mRNA structural stability (ΔG).

LinearDesign (Zhang et al. 2020, Nature) jointly optimises codon usage and mRNA
secondary-structure stability; more structure lengthens mRNA half-life. Running
the full lattice design needs the LinearDesign package, but its central *objective*
— the folding free energy of the mRNA — is computable directly with ViennaRNA,
which is a light dependency. This plugin reports that objective for the submitted
sequence, so a user can see how structurally stable their current design is.
"""
from __future__ import annotations

from typing import List

from .base import Plugin, PluginResult


class LinearDesignPlugin(Plugin):
    NAME = "lineardesign"
    REQUIRES = ["ViennaRNA"]
    DESCRIPTION = "mRNA structural-stability objective (ΔG) that LinearDesign optimises."
    CITATION = "Zhang et al. 2020, Nature (LinearDesign)"

    def available(self) -> bool:
        try:
            import RNA  # noqa: F401
            return True
        except Exception:
            return False

    def analyze(self, transcript) -> List[PluginResult]:
        import RNA
        seq = transcript.cds[:1000].replace("T", "U")
        if len(seq) < 30:
            return []
        _struct, mfe = RNA.fold(seq)
        per_nt = mfe / len(seq)
        return [PluginResult(
            plugin=self.NAME,
            label="mRNA folding ΔG (5' up to 1000 nt)",
            value=round(per_nt, 3),
            text=f"{mfe:.1f} kcal/mol total, {per_nt:.3f} kcal/mol/nt",
            note="More negative = more structured / more stable mRNA. This is the "
                 "LinearDesign optimisation objective for the current sequence; the full "
                 "package can redesign toward a lower ΔG.",
            citation=self.CITATION,
            validated=False,
        )]
