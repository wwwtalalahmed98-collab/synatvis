"""cloning — Type IIS domestication scanner for algal MoClo (CLAUDE.md §2, §3.5).

Golden Gate / MoClo is established for Cr (Crozet 2018). An internal Type IIS
recognition site (BsaI/BsmBI/BbsI/SapI) on either strand breaks domestication and
must be removed. This validated, low-noise module flags each site and, for
CDS-internal sites, proposes a minimal synonymous fix (remediation §6).
"""
from __future__ import annotations

from typing import Dict, List

from ..flags import Flag, Severity
from ..remediation import synonymous_fix
from ..util import find_all, reverse_complement

NAME = "cloning"
VALIDATED = True
DEFAULT_ON = True
SUMMARY = "Internal Type IIS sites (BsaI/BsmBI/BbsI/SapI) that break MoClo assembly"

EVIDENCE = "Crozet 2018 (algal MoClo); Schroda 2019 — Type IIS domestication is standard for Cr"


def run(tx, profile: Dict) -> List[Flag]:
    cfg = profile["cloning"]
    enzymes: Dict[str, str] = cfg["enzymes"]
    scope = set(cfg.get("scope", ["5utr", "cds", "3utr"]))
    table = profile.get("_codon_table")

    regions = {
        "5utr": (tx.utr5, 0),
        "cds": (tx.cds, tx.cds_start),
        "3utr": (tx.utr3, tx.cds_end),
    }

    flags: List[Flag] = []
    for enzyme, site in enzymes.items():
        site = site.upper()
        rc = reverse_complement(site)
        variants = {"+": site}
        if rc != site:
            variants["-"] = rc
        for region in scope:
            seq, offset = regions[region]
            for strand, pat in variants.items():
                for pos in find_all(seq, pat):
                    abs_start = offset + pos
                    abs_end = abs_start + len(pat)
                    suggestion = None
                    if region == "cds" and table is not None:
                        edit = synonymous_fix(
                            tx.cds, table,
                            cds_start=pos, cds_end=pos + len(pat),
                            avoid=[site, rc],
                        )
                        suggestion = edit.describe()
                    elif region != "cds":
                        suggestion = ("edit the UTR to remove the site "
                                      "(no coding constraint here)")
                    flags.append(Flag(
                        module=NAME,
                        severity=Severity.HIGH,
                        start=abs_start,
                        end=abs_end,
                        region=region,
                        message=(f"{enzyme} site ({pat}, {strand} strand) in {region}. "
                                 f"Domesticate before MoClo assembly."),
                        evidence=EVIDENCE,
                        suggested_edit=suggestion,
                        detail={"enzyme": enzyme, "strand": strand, "site": pat},
                    ))
    return flags
