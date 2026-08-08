# SynAT.Vis test kit — expected results

Six tiny FASTA cassettes whose **ground truth is written in each header**. A
working scanner must flag the planted problem in each and stay quiet on the
control. Run them individually with the commands below, or verify all at once
with the bundled checker from the project root:

```
python selftest.py
```

Expected: `6/6 checks passed`.

| File | Command | Expected flag |
|---|---|---|
| `kit_00_clean_control.fasta` | `synatvis scan kit_00_clean_control.fasta` | **nothing** high/medium (1 INFO = silencing guidance) |
| `kit_01_cloning_BsaI.fasta` | `synatvis scan kit_01_cloning_BsaI.fasta` | `cloning` **HIGH** — internal BsaI site `GGTCTC`, with a synonymous fix |
| `kit_02_polya_TGTAA.fasta` | `synatvis scan kit_02_polya_TGTAA.fasta` | `polya` **HIGH** — premature `TGTAA` poly(A) signal in a G/C-rich context |
| `kit_03_instability_ARE.fasta` | `synatvis scan kit_03_instability_ARE.fasta` | `instability` **MEDIUM** — clustered AU-rich elements (`ATTTA`) |
| `kit_04_composition_homopolymer.fasta` | `synatvis scan kit_04_composition_homopolymer.fasta` | `composition` — 12×`A` homopolymer run |
| `kit_05_uorf_upstream_AUG.fasta` | `synatvis scan kit_05_uorf_upstream_AUG.fasta --cds 30:144` | `uorf` — upstream AUG in the 5'UTR |

`--cds START:END` (0-based, half-open) tells the scanner where the coding region
sits so the UTRs are checked; only `kit_05` needs it because it has a 5'UTR.

## How to read a hit

Each flag prints: severity, region, position, a plain-language message, the
literature `evidence`, and a `fix` (a minimal synonymous edit that removes the
motif while preserving the protein — or an honest "no clean synonymous fix").

## What this proves and what it does not

- **Proves:** the detector fires on known motifs and does not invent flags on a
  clean sequence (detection sensitivity + specificity on controlled inputs).
- **Does not prove:** any biological outcome. SynAT.Vis is transcript-level and
  diagnostic only — it never scores expression, and it is silent on protein
  stability. A clean report is not a promise of expression.
