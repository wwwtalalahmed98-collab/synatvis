"""Scanner orchestrator: profile + transcript -> ranked flags (CLAUDE.md §4).

Loads the active profile's codon table once, attaches it to the profile under a
private key, runs every enabled module, and returns a :class:`ScanResult`. It
emits NO composite score — only a ranked flag list and the metadata the report
needs (CLAUDE.md §1, §9).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from . import modules  # noqa: F401  (import triggers module registration)
from .codon_tables import CodonTable, load_for_profile
from .flags import Flag, registered, sort_flags
from .profiles import load_profile
from .seqio import Transcript


@dataclass
class ScanResult:
    transcript: Transcript
    profile: Dict[str, Any]
    flags: List[Flag]
    module_meta: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    codon_source: str = ""

    def by_module(self) -> Dict[str, List[Flag]]:
        out: Dict[str, List[Flag]] = {}
        for f in self.flags:
            out.setdefault(f.module, []).append(f)
        return out

    def counts(self) -> Dict[str, int]:
        c = {"info": 0, "low": 0, "medium": 0, "high": 0}
        for f in self.flags:
            c[str(f.severity)] += 1
        return c


def _ensure_codon_table(profile: Dict[str, Any]) -> CodonTable:
    if "_codon_table" in profile and isinstance(profile["_codon_table"], CodonTable):
        return profile["_codon_table"]
    base = profile.get("_base_dir", ".")
    table = load_for_profile(profile, base)
    profile["_codon_table"] = table
    return table


def scan(
    transcript: Transcript,
    profile: "str | Dict[str, Any]" = "cr_nuclear",
    only: Optional[List[str]] = None,
    exclude: Optional[List[str]] = None,
) -> ScanResult:
    """Scan *transcript* under *profile* (name or loaded dict).

    ``only`` / ``exclude`` select modules by name; otherwise every ``default_on``
    module runs.
    """
    if isinstance(profile, str):
        profile = load_profile(profile)
    table = _ensure_codon_table(profile)

    only_set = set(only) if only else None
    exclude_set = set(exclude) if exclude else set()

    all_flags: List[Flag] = []
    module_meta: Dict[str, Dict[str, Any]] = {}
    for spec in registered():
        module_meta[spec.name] = {"validated": spec.validated,
                                  "summary": spec.summary,
                                  "default_on": spec.default_on}
        if only_set is not None:
            if spec.name not in only_set:
                continue
        elif not spec.default_on:
            continue
        if spec.name in exclude_set:
            continue
        all_flags.extend(spec.run(transcript, profile))

    return ScanResult(
        transcript=transcript,
        profile=profile,
        flags=sort_flags(all_flags),
        module_meta=module_meta,
        codon_source=getattr(table, "source", ""),
    )
