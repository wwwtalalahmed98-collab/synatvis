# SynAT.Vis — Methods & Validation (preprint draft)

> **Draft for you to verify, adapt, and own.** The numbers are from the validation
> runs in this build; re-run them and confirm before submission. Citations are
> given author–year; format them to your target journal. See the note at the end
> about authorship.

## Implementation

SynAT.Vis checks a recombinant cassette designed for nuclear expression in
*Chlamydomonas reinhardtii* and reports why its transcript might fail. It reads a
cassette (FASTA or GenBank) as three regions — 5′UTR, CDS, 3′UTR — and returns a
severity-ranked list of flags. Each flag carries a location, a one-line literature
basis, and, where one exists, a minimal synonymous edit that removes the flagged
motif without changing the protein. The tool reports no composite expression
score by design; it is diagnostic, and inventing an uncalibrated number is the
failure mode it exists to avoid.

Host biology lives in a single YAML profile. The detection modules read every
constant from that profile and hard-code none, so a chloroplast or a second-lineage
profile drops in without changing module code. The core depends only on the Python
standard library (≥ 3.8). Biopython (GenBank input) and ViennaRNA (5′UTR folding)
are optional; when absent, the affected checks fall back to heuristics that label
themselves as heuristic in the report. The build ships with 15 unit tests.

## Host model and codon usage

The nuclear profile encodes the constants that set *C. reinhardtii* apart from
land-plant expectations: coding GC near 66% with a strong GC3 bias, the UGUAA
(TGTAA) near-upstream polyadenylation element in place of AAUAAA, and introns
treated as expression aids rather than cryptic-splice hazards (Merchant et al.
2007; Barahimipour et al. 2015; Shen et al. 2008; Baier et al. 2018). Codon usage
is taken from the Kazusa Codon Usage Database (taxid 3055, GenBank plant division;
846 CDS, 420,455 codons; coding GC 66.3%) and stored as each codon's fraction
within its amino-acid family. Relative adaptiveness is that fraction normalised to
the most-used codon in the family; a sequence's Relative Codon Adaptation (RCA) is
the geometric mean of those values across its sense codons.

To test whether this table reflects the organism rather than one source's sampling,
we recomputed codon usage from independent *C. reinhardtii* nuclear CDS drawn from
NCBI RefSeq. Against 1,000 CDS (435,901 codons) the most-used codon agreed in 20 of
21 amino-acid families (mean absolute per-codon difference 0.02); against a larger
6,930-CDS count (2.9M codons) it agreed in 18 of 21, the differences confined to
near-ties pushed further toward the GC-rich codon by the GC-enriched sample. The
larger RefSeq-derived table is provided as an alternative
(`cr_nuclear_codon_refseq.tsv`); the Kazusa table remains the default because it is
citable and genome-representative. Both describe the same GC3 bias.

## Detection modules

Nine modules run by default, each reading its thresholds from the active profile:

- **codon** — clusters of low-adaptiveness codons; whole-CDS RCA. In *C.
  reinhardtii* codon choice affects both elongation and mRNA stability
  (Barahimipour et al. 2015), so this is one input, never a verdict.
- **composition** — local GC *troughs* and homopolymer runs. The GC rule is
  inverted relative to land plants: low GC is the hazard, because unfavourable GC
  drives heterochromatinisation and transgene silencing (Barahimipour et al.
  2015). High GC is normal and is never flagged.
- **polya** — premature TGTAA polyadenylation signals in a G/C-rich context
  (Shen et al. 2008; Zhao et al. 2014).
- **cloning** — internal Type IIS sites (BsaI, BsmBI, BbsI, SapI) that break
  Golden Gate / MoClo assembly (Crozet et al. 2018).
- **structure, uorf, instability, splice, silencing** — start-codon context and
  5′UTR structure; upstream AUGs; AU-rich elements; a re-scoped Cr splice check
  with no dicot cryptic-intron logic; and a quarantined, unvalidated silencing
  heuristic that reports only sequence-visible risk and points to strain choice
  (UVM4/UVM11) and intron-containing promoters (Baier et al. 2018; Schroda 2019).

For any pattern-based flag inside the CDS, a remediation step searches for the
smallest set of synonymous codon changes that removes the motif and does not
recreate it in the local window, ranked by adaptiveness on the profile's codon
table. When no clean change exists, the tool says so rather than proposing an edit.

## Validation

The harness has three legs, built alongside the modules rather than after them.

**Specificity (Leg 1).** A well-adapted endogenous gene should raise few flags, so
native highly-expressed CDS serve as negatives. We assembled 500 *C. reinhardtii*
nuclear CDS from NCBI RefSeq — full-length, in-frame, GC-rich (> 55%) — and ranked
them by RCA as a proxy for high expression, then measured the fraction raising a
medium- or high-severity flag per module. At the initial thresholds three modules
over-fired: poly(A) 46.0%, codon 44.8%, composition 13.4%, because a single short
motif recurs by chance in a GC-rich genome. We re-set the operating points against
this measurement — poly(A) context stringency 0.55 → 0.85; rare-codon weight
0.15 → 0.10 with cluster size 4 → 5; GC-trough width 40 → 60 nt at a stricter
depth — bringing the rates to 1.6%, 3.2% and 1.4% while the injected-signal
sensitivity held (Leg 2). The cloning module reads 72%, which is not a
false-positive rate: it is a deterministic match for Type IIS sites, and most
native GC-rich genes genuinely carry one that requires domestication before
assembly.

Thresholds fitted to a dataset can memorise it. To test generalisation we re-ran
Leg 1 on 4,500 genes the thresholds were never tuned on: a 500-CDS set and a
4,000-CDS set, each sharing no sequence with the tuning set or each other
(pairwise-disjoint; 5,000 unique CDS in total). The false-positive rates on the
500-CDS set were 1.6% (poly(A)), 1.8% (codon) and 1.0% (composition); on the
4,000-CDS set 2.5%, 4.2% and 2.0%. Both sit within sampling range of the tuning
set. The operating points generalise rather than overfit.

**Sensitivity (Leg 2).** Into a clean, GC-rich parent ORF we inserted three known
motifs — a TGTAA element in a G/C-rich context, a BsaI site, and a pair of AU-rich
elements — and confirmed each was flagged by the intended module while the parent
stayed silent (3/3). This measures detection, not biological consequence.

**Literature cases (Leg 3).** Scoring is stratified by compartment. Two rescue
pairs use real reporter sequences: native GFP (GC 38%) against a
codon-optimised counterpart (GC 62%), and native firefly luciferase (GC 45%)
against an optimised counterpart (GC 65%); in each pair the protein is unchanged
and only synonymous codons differ. The native, AT-rich members raise composition
flags and the optimised members scan clean, reproducing the documented
codon-optimisation rescue (Fuhrmann et al. 1999; Barahimipour et al. 2015). A
native gene carrying a BsaI site is flagged by the cloning module. An AT-rich
plastid gene is gated out of nuclear scoring, since a nuclear detector says nothing
about a chloroplast construct. On the six cases the nuclear detectors gave 3/3
sensitivity and 2/2 specificity. The corpus is small; read this as a consistency
check, not a powered estimate.

### Table 1. Validation summary

| Leg | Data | Result |
|---|---|---|
| 1 — specificity (tuning) | 500 real Cr nuclear CDS (RefSeq) | FP after tuning: poly(A) 1.6%, codon 3.2%, composition 1.4%; cloning 72% (site prevalence) |
| 1 — specificity (held-out) | 4,500 disjoint real CDS (500 + 4,000; 0 shared) | FP: poly(A) 1.6–2.5%, codon 1.8–4.2%, composition 1.0–2.0% — generalises |
| 2 — sensitivity | 3 injected motifs into a clean parent | 3/3 detected, parent silent |
| 3 — literature cases | 6 real cases (2 rescue pairs + BsaI + gated plastid) | nuclear 3/3 sensitivity, 2/2 specificity; plastid gated out |
| codon table | up to 6,930 independent RefSeq CDS (2.9M codons) | preferred codon agrees 18–20/21 families; mean abs. per-codon difference 0.02–0.04 |

## Availability

SynAT.Vis is written in Python (standard library only; ≥ 3.8) with optional
Biopython and ViennaRNA backends. The distribution includes the profile, the codon
table, the 500-CDS specificity sets, the literature cases, the three validation
legs, a self-test, and a plain-language report mode. [Repository / archive DOI to
be added.]

## Limitations

The tool works at the transcript level and is silent on protein folding, stability,
and proteolysis, which are real causes of low yield it cannot see; a clean report is
not a promise of expression. It reports no expression score. The nuclear and
chloroplast compartments are treated as separate, non-transferable profiles. The
silencing module is a sequence-visible heuristic and is not validated, because
silencing in *C. reinhardtii* is driven mainly by chromatin and strain rather than
CDS sequence. "Highly-expressed" in the specificity sets is approximated by codon
adaptation rather than measured RNA abundance; an expression-ranked set is a
planned upgrade. The literature corpus is six cases, and the sensitivity and
specificity figures on it are consistency checks rather than powered estimates.

---

*Note on authorship.* This is a drafting aid, not a finished submission. The
sequences, runs, and numbers are real, but the biological framing, the citation
choices, and the claims are yours to check and defend — a reviewer will expect the
authors to own every sentence. Verify each figure against a fresh run, confirm each
reference, and rewrite anything you would not say in your own words.
