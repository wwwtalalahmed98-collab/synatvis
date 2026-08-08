"""Small sequence helpers shared across modules."""
from __future__ import annotations

import re
from typing import Iterator, List, Tuple

_COMP = str.maketrans("ACGTNRYSWKMBDHVacgtnryswkmbdhv",
                      "TGCANYRSWMKVHDBtgcanyrswmkvhdb")

_IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def iupac_to_regex(motif: str) -> str:
    """Turn an IUPAC motif into a regex (U treated as T)."""
    return "".join(_IUPAC.get(b, re.escape(b)) for b in motif.upper().replace("U", "T"))


def gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for b in seq.upper() if b in "GC")
    return gc / len(seq)


def find_all(seq: str, motif: str, iupac: bool = False) -> List[int]:
    """Return 0-based start indices of *motif* in *seq* (overlapping)."""
    seq = seq.upper()
    starts: List[int] = []
    if iupac:
        pat = re.compile(iupac_to_regex(motif))
        i = 0
        while True:
            m = pat.search(seq, i)
            if not m:
                break
            starts.append(m.start())
            i = m.start() + 1
    else:
        motif = motif.upper()
        i = seq.find(motif)
        while i != -1:
            starts.append(i)
            i = seq.find(motif, i + 1)
    return starts


def windows(seq: str, size: int, step: int) -> Iterator[Tuple[int, str]]:
    """Yield ``(start, subseq)`` sliding windows; a short tail window is skipped."""
    if size <= 0 or len(seq) < size:
        return
    for start in range(0, len(seq) - size + 1, step):
        yield start, seq[start:start + size]


def longest_homopolymer(seq: str) -> Tuple[str, int, int]:
    """Return ``(base, start, length)`` of the longest single-base run."""
    best_base, best_start, best_len = "", 0, 0
    i = 0
    seq = seq.upper()
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        if j - i > best_len:
            best_base, best_start, best_len = seq[i], i, j - i
        i = j
    return best_base, best_start, best_len
