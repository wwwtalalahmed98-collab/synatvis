"""Profile loading (CLAUDE.md §4).

``load_profile('cr_nuclear')`` reads the bundled YAML, validates it against the
schema, and returns a plain dict. YAML is parsed by PyYAML when installed, else
by the bundled minimal parser (``synatvis._yaml``).
"""
from __future__ import annotations

import os
from typing import Any, Dict

from . import _schema

_HERE = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(_HERE)  # the synatvis/ dir, base for data paths


def load_yaml(text: str) -> Any:
    try:
        import yaml  # type: ignore

        return yaml.safe_load(text)
    except Exception:
        from .. import _yaml

        return _yaml.safe_load(text)


def _resolve(name_or_path: str) -> str:
    if os.path.isfile(name_or_path):
        return name_or_path
    candidate = os.path.join(_HERE, name_or_path)
    if os.path.isfile(candidate):
        return candidate
    candidate = os.path.join(_HERE, name_or_path + ".yaml")
    if os.path.isfile(candidate):
        return candidate
    raise FileNotFoundError(f"no profile named {name_or_path!r} in {_HERE}")


def load_profile(name_or_path: str = "cr_nuclear", validate: bool = True) -> Dict[str, Any]:
    """Load and (by default) validate a host profile.

    The returned dict carries a private ``_base_dir`` key so downstream loaders
    (e.g. the codon table) can resolve relative paths.
    """
    path = _resolve(name_or_path)
    with open(path, "r", encoding="utf-8") as fh:
        profile = load_yaml(fh.read())
    if not isinstance(profile, dict):
        raise ValueError(f"profile {path!r} did not parse to a mapping")
    profile["_base_dir"] = PACKAGE_DIR
    profile["_path"] = path
    if validate:
        _schema.require_valid(profile)
    return profile


def available() -> Dict[str, str]:
    """Map profile name -> path for every bundled ``*.yaml`` profile."""
    out = {}
    for fn in os.listdir(_HERE):
        if fn.endswith(".yaml"):
            out[fn[:-5]] = os.path.join(_HERE, fn)
    return out
