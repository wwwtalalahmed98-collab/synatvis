"""Minimal YAML-subset parser — a dependency-free fallback for PyYAML.

SynAT.Vis profiles and ``cases.yaml`` only use a small, well-behaved subset of
YAML: block mappings, block sequences (including sequences of mappings), inline
flow lists ``[a, b, c]``, comments, and scalars (str/int/float/bool/null). This
module parses exactly that subset so the tool runs on a bare standard library.

If PyYAML is installed, :func:`synatvis.profiles.load_yaml` uses it instead;
this is only the fallback. Keep the shipped ``*.yaml`` files inside this subset.
"""
from __future__ import annotations

from typing import Any, List, Tuple


class MiniYAMLError(ValueError):
    """Raised when the input uses YAML features this parser does not support."""


def _strip_comment(line: str) -> str:
    """Remove an unquoted ``#`` comment from *line*, honouring quotes."""
    out = []
    quote = None
    i = 0
    while i < len(line):
        ch = line[i]
        if quote:
            out.append(ch)
            if ch == quote:
                quote = None
        else:
            if ch in ("'", '"'):
                quote = ch
                out.append(ch)
            elif ch == "#" and (i == 0 or line[i - 1] in " \t"):
                break
            else:
                out.append(ch)
        i += 1
    return "".join(out).rstrip()


def _scalar(text: str) -> Any:
    """Coerce a scalar token to a Python value."""
    text = text.strip()
    if text == "" or text in ("~", "null", "None"):
        return None
    if len(text) >= 2 and text[0] in "'\"" and text[-1] == text[0]:
        return text[1:-1]
    low = text.lower()
    if low in ("true", "yes"):
        return True
    if low in ("false", "no"):
        return False
    if text.startswith("[") and text.endswith("]"):
        inner = text[1:-1].strip()
        if not inner:
            return []
        return [_scalar(part) for part in _split_flow(inner)]
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        pass
    return text


def _split_flow(inner: str) -> List[str]:
    """Split a flow-list body on commas that are not inside quotes."""
    parts, buf, quote = [], [], None
    for ch in inner:
        if quote:
            buf.append(ch)
            if ch == quote:
                quote = None
        elif ch in ("'", '"'):
            quote = ch
            buf.append(ch)
        elif ch == ",":
            parts.append("".join(buf))
            buf = []
        else:
            buf.append(ch)
    if buf:
        parts.append("".join(buf))
    return parts


def _looks_structured(content: str) -> bool:
    """True if a single line opens a mapping or sequence (not a bare scalar)."""
    if content.startswith(("'", '"')):
        return False
    if content.startswith("- "):
        return True
    if content.endswith(":"):
        return True
    # a "key: value" mapping line (colon followed by space), but not a flow list
    return ": " in content and not content.startswith("[")


def _tokenize(text: str) -> List[Tuple[int, str]]:
    lines = []
    for raw in text.replace("\t", "    ").splitlines():
        content = _strip_comment(raw)
        if not content.strip():
            continue
        indent = len(content) - len(content.lstrip(" "))
        lines.append((indent, content.strip()))
    return lines


def _parse_block(lines: List[Tuple[int, str]], i: int, indent: int) -> Tuple[Any, int]:
    if i >= len(lines):
        return None, i
    if lines[i][1].startswith("-"):
        return _parse_seq(lines, i, indent)
    return _parse_map(lines, i, indent)


def _parse_seq(lines, i, indent):
    result = []
    while i < len(lines) and lines[i][0] == indent and lines[i][1].startswith("-"):
        content = lines[i][1]
        first = content[1:].strip()  # text after the dash
        sub: List[Tuple[int, str]] = []
        if first:
            sub.append((indent + 2, first))
        i += 1
        while i < len(lines) and lines[i][0] > indent:
            sub.append(lines[i])
            i += 1
        if not sub:
            result.append(None)
        elif len(sub) == 1 and not _looks_structured(sub[0][1]):
            result.append(_scalar(sub[0][1]))
        else:
            base = min(s[0] for s in sub)
            val, _ = _parse_block(sub, 0, base)
            result.append(val)
    return result, i


def _parse_map(lines, i, indent):
    result = {}
    while i < len(lines) and lines[i][0] == indent and not lines[i][1].startswith("-"):
        content = lines[i][1]
        key, sep, rest = content.partition(":")
        if not sep:
            raise MiniYAMLError(f"expected 'key: value', got: {content!r}")
        key = key.strip()
        rest = rest.strip()
        if rest == "":
            sub = []
            i += 1
            while i < len(lines) and lines[i][0] > indent:
                sub.append(lines[i])
                i += 1
            if sub:
                base = min(s[0] for s in sub)
                val, _ = _parse_block(sub, 0, base)
            else:
                val = None
            result[key] = val
        else:
            result[key] = _scalar(rest)
            i += 1
    return result, i


def safe_load(text: str) -> Any:
    """Parse a YAML-subset *text* into Python data structures."""
    lines = _tokenize(text)
    if not lines:
        return None
    value, _ = _parse_block(lines, 0, lines[0][0])
    return value
