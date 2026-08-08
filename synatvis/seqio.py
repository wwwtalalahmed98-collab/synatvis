"""Sequence I/O and the :class:`Transcript` model (CLAUDE.md §4).

A cassette is modelled as three regions: ``utr5``, ``cds``, ``utr3``. FASTA is
parsed with the standard library; GenBank / SnapGene is parsed via Biopython if
available, otherwise a minimal GenBank ``ORIGIN`` fallback covers common cases.

Region assignment when the input is a single unannotated sequence: if the record
contains a single ORF flanked by UTRs, callers may split explicitly. The
convenience :func:`read_record` returns a :class:`Transcript`; when no UTR
annotation is present the whole sequence is treated as CDS with empty UTRs and a
note is attached, because CDS-only scanning is still meaningful.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, Tuple

_DNA = re.compile(r"[^ACGTUN]")


def normalise(seq: str) -> str:
    """Uppercase, strip whitespace, map U->T. Non-ACGTN raises."""
    seq = "".join(seq.split()).upper().replace("U", "T")
    offenders = set(_DNA.findall(seq))
    if offenders:
        raise ValueError(f"non-nucleotide characters in sequence: {sorted(offenders)}")
    return seq


@dataclass
class Transcript:
    """A cassette split into 5'UTR, CDS and 3'UTR (all 5'->3', coding strand)."""

    cds: str
    utr5: str = ""
    utr3: str = ""
    name: str = "transcript"
    notes: List[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.cds = normalise(self.cds)
        self.utr5 = normalise(self.utr5)
        self.utr3 = normalise(self.utr3)

    # -- geometry -------------------------------------------------------
    @property
    def full(self) -> str:
        return self.utr5 + self.cds + self.utr3

    @property
    def cds_start(self) -> int:
        return len(self.utr5)

    @property
    def cds_end(self) -> int:
        return len(self.utr5) + len(self.cds)

    def abs_from_cds(self, pos: int) -> int:
        """Map a CDS-relative index to a transcript-absolute index."""
        return self.cds_start + pos

    def region_at(self, abs_pos: int) -> str:
        if abs_pos < self.cds_start:
            return "5utr"
        if abs_pos < self.cds_end:
            return "cds"
        return "3utr"

    def codons(self) -> Iterator[Tuple[int, str]]:
        """Yield (codon_index, codon) over the CDS, ignoring a ragged tail."""
        cds = self.cds
        for k in range(0, len(cds) - len(cds) % 3, 3):
            yield k // 3, cds[k:k + 3]


# ---------------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------------
def read_fasta(path: str) -> List[Tuple[str, str]]:
    """Return ``[(header, sequence), ...]`` from a FASTA file."""
    records: List[Tuple[str, str]] = []
    header: Optional[str] = None
    chunks: List[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:].strip()
                chunks = []
            elif line.strip():
                chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def parse_fasta_str(text: str) -> List[Tuple[str, str]]:
    records, header, chunks = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:].strip(), []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


# ---------------------------------------------------------------------------
# GenBank / SnapGene (optional Biopython, minimal fallback)
# ---------------------------------------------------------------------------
def _genbank_features_biopython(path: str) -> Optional[Transcript]:
    try:
        from Bio import SeqIO  # type: ignore
    except Exception:
        return None
    rec = next(SeqIO.parse(path, "genbank"))
    seq = str(rec.seq)
    cds_feats = [f for f in rec.features if f.type == "CDS"]
    if not cds_feats:
        return Transcript(cds=seq, name=rec.id or "transcript",
                          notes=["no CDS feature; whole record treated as CDS"])
    f = cds_feats[0]
    start = int(f.location.start)
    end = int(f.location.end)
    return Transcript(
        cds=seq[start:end],
        utr5=seq[:start],
        utr3=seq[end:],
        name=rec.id or "transcript",
    )


def _genbank_origin_fallback(text: str) -> str:
    """Extract the ORIGIN sequence block from raw GenBank text."""
    out: List[str] = []
    in_origin = False
    for line in text.splitlines():
        if line.startswith("ORIGIN"):
            in_origin = True
            continue
        if in_origin:
            if line.startswith("//"):
                break
            out.append(re.sub(r"[^A-Za-z]", "", line))
    return "".join(out)


def read_record(path: str, cds_span: Optional[Tuple[int, int]] = None) -> Transcript:
    """Read one cassette from *path* (FASTA or GenBank) as a Transcript.

    ``cds_span`` optionally gives ``(start, end)`` 0-based half-open CDS bounds in
    the full sequence, letting a plain FASTA carry UTR annotation.
    """
    lower = path.lower()
    if lower.endswith((".gb", ".gbk", ".genbank", ".dna")):
        tx = _genbank_features_biopython(path)
        if tx is not None:
            return tx
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            seq = _genbank_origin_fallback(fh.read())
        if not seq:
            raise ValueError(f"could not extract sequence from {path}")
        return _split(seq, "record", cds_span)

    records = read_fasta(path)
    if not records:
        raise ValueError(f"no FASTA records in {path}")
    header, seq = records[0]
    return _split(seq, header.split()[0] if header else "transcript", cds_span)


def _split(seq: str, name: str, cds_span: Optional[Tuple[int, int]]) -> Transcript:
    if cds_span is not None:
        s, e = cds_span
        return Transcript(cds=seq[s:e], utr5=seq[:s], utr3=seq[e:], name=name)
    return Transcript(
        cds=seq,
        name=name,
        notes=["no UTR annotation supplied; whole sequence treated as CDS "
               "(pass cds_span or split explicitly for UTR-aware checks)"],
    )
