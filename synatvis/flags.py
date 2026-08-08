"""Flag dataclass, severity ordering, and the module registry (CLAUDE.md §5).

A :class:`Flag` is the atomic unit of the report. Every module returns a list of
them. Coordinates are 0-based, half-open ``[start, end)``, expressed in the
concatenated transcript (5'UTR + CDS + 3'UTR) so they are directly comparable
across modules; ``region`` records which part the flag falls in.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from enum import IntEnum
from typing import Any, Callable, Dict, List, Optional


class Severity(IntEnum):
    """Ranked severity. Higher sorts first in the report."""

    INFO = 0
    LOW = 1
    MEDIUM = 2
    HIGH = 3

    @classmethod
    def parse(cls, value: "str | Severity") -> "Severity":
        if isinstance(value, Severity):
            return value
        return cls[str(value).strip().upper()]

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.name.lower()


@dataclass
class Flag:
    """A single transcript-level red flag.

    Attributes mirror the CLAUDE.md §5 report contract.
    """

    module: str
    severity: Severity
    start: int
    end: int
    region: str  # "5utr" | "cds" | "3utr" | "transcript"
    message: str
    evidence: str = ""  # one-line literature basis
    suggested_edit: Optional[str] = None
    detail: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.severity = Severity.parse(self.severity)

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["severity"] = str(self.severity)
        return d


@dataclass
class ModuleSpec:
    """Registry entry describing one detection module (CLAUDE.md §4)."""

    name: str
    run: Callable[[Any, Dict[str, Any]], List[Flag]]
    validated: bool = True
    default_on: bool = True
    summary: str = ""


_REGISTRY: "Dict[str, ModuleSpec]" = {}
_ORDER: List[str] = []


def register(spec: ModuleSpec) -> ModuleSpec:
    """Register a module. Later scans iterate registration order."""
    if spec.name in _REGISTRY:
        raise ValueError(f"module {spec.name!r} already registered")
    _REGISTRY[spec.name] = spec
    _ORDER.append(spec.name)
    return spec


def registered() -> List[ModuleSpec]:
    return [_REGISTRY[name] for name in _ORDER]


def get(name: str) -> ModuleSpec:
    return _REGISTRY[name]


def sort_flags(flags: List[Flag]) -> List[Flag]:
    """Rank flags: severity desc, then transcript position, then module."""
    return sorted(flags, key=lambda f: (-int(f.severity), f.start, f.module))
