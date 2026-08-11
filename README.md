# SynAT.Vis

**Syn**thetic **A**lgal **T**ranscript **Vis**ualiser / validator — a
transcript-level **red-flag scanner** for recombinant gene cassettes destined for
**nuclear** expression in *Chlamydomonas reinhardtii*.

It reads a designed cassette and reports *why it might fail to produce an intact,
well-translated transcript in the host*, with a suggested synonymous fix per flag.
It is **diagnostic, not predictive**: it emits a ranked flag list and **never a
composite "expression score."** It is **silent on proteolysis / protein
stability**, and a nuclear report says nothing about a chloroplast construct.

The design contract lives in [`CLAUDE.md`](CLAUDE.md) and is the source of truth;
this README is the operator's guide.

## Why this host is not a generic "algae" profile

Every threshold is set by the **host + compartment** pair, and Cr **nuclear**
biology runs *opposite* to land plants on several axes:

| Axis | Cr nuclear (this tool) | Land-plant intuition (wrong here) |
|---|---|---|
| Coding GC | ~68% GC, strongly preferred | ~44% |
| GC hazard | **LOW** GC → silencing | high GC |
| poly(A) signal | **UGUAA (TGTAA)**, G/C-rich context | AAUAAA, AU-rich |
| Introns | **assets** (rbcS2i1 boosts expression) | cryptic-splice hazards |

The chloroplast pulls the other way again (high-AT, Shine-Dalgarno, no
spliceosomal introns) — it ships only as a **stub** profile.

## Install

The **core runs on the Python standard library alone** (Python ≥ 3.8). Optional
backends upgrade specific features and the tool degrades gracefully without them:

```bash
pip install -r requirements.txt   # optional: pyyaml, biopython, ViennaRNA
```

* **PyYAML** — faster/robust profile parsing (else a bundled minimal parser).
* **Biopython** — GenBank / SnapGene input (else FASTA + a minimal GenBank fallback).
* **ViennaRNA** — real 5'UTR start-codon folding (else a labelled GC heuristic).

## Usage

```bash
# scan a CDS-only FASTA
python -m synatvis scan cassette.fasta

# scan a cassette with an annotated CDS span (0-based, half-open) so UTRs are checked
python -m synatvis scan cassette.fasta --cds 51:801

# GenBank/SnapGene input, JSON output, non-zero exit if any HIGH flag
python -m synatvis scan cassette.gb --json --fail-on-high

# add the expression-propensity prediction (opt-in model; always included in --html)
python -m synatvis scan cassette.fasta --expression

# a polished HTML report with the expression gauge + per-position landscape
python -m synatvis scan cassette.fasta --html --out report.html

# an animated 4D translation-time view (ribosomes traversing the landscape)
python -m synatvis scan cassette.fasta --4d --out translation_4d.html

# the animated cell-journey view (transcription -> localisation, inside the cell)
python -m synatvis scan cassette.fasta --cell --out cell_journey.html

# also run any installed Tier-B ML plugins (experimental, opt-in)
python -m synatvis scan cassette.fasta --plugins

# list host profiles / list Tier-B ML plugins and how to enable them
python -m synatvis profiles
python -m synatvis plugins

# run the validation legs (CLAUDE.md §7)
python -m synatvis validate leg1 --fasta synatvis/data/native_cr_cds.fasta
python -m synatvis validate leg2
python -m synatvis validate leg3
```

Library:

```python
from synatvis import scan, Transcript
res = scan(Transcript(cds="ATG...TAA", utr5="...", utr3="..."))
print(res.counts())
for f in res.flags:
    print(f.severity, f.module, f.message, f.suggested_edit)
```

## Modules (CLAUDE.md §3, §8)

| Module | Severity role | Notes |
|---|---|---|
| `codon` | rare-codon clusters, low RCA, **codon optimality** (mRNA stability), **codon-pair bias**, **tAI** | multi-criteria panel (Demissie 2025); optimality→stability (Presnyak 2015); tAI from GtRNAdb tRNA counts (dos Reis 2004) |
| `composition` | GC **troughs**, homopolymers | **inverted**: low GC is the hazard; high GC never flagged |
| `cloning` | Type IIS sites (BsaI/BsmBI/BbsI/SapI) | MoClo domestication; synonymous fix offered |
| `structure` | start codon, 5'UTR length, **mRNA folding ΔG** (start-region + global) | ViennaRNA MFE if present, else a pure-Python Nussinov fallback (labelled heuristic); structured starts impede initiation (Kudla 2009; Zhang 2020) |
| `uorf` | upstream AUGs / uORFs | scanning competition |
| `instability` | AU-rich element clusters | re-tuned to Cr |
| `polya` | premature **TGTAA** in G/C context | **noisiest**; operating point from Leg-1 |
| `splice` | intended-intron sanity + accidental pairs | low priority; **no** dicot cryptic detector, **no** clade gating |
| `silencing` | GC troughs + strain/intron guidance | **QUARANTINED / not validated**; flags marked HEURISTIC |

Every flag carries `module, severity, start, end, region, message, evidence,
suggested_edit, detail`. The report prints a banner (host, compartment, scope,
no-score notice), groups flags by module, ranks by severity, marks the silencing
module HEURISTIC, and emits both text and JSON.

## Remediation (CLAUDE.md §6)

For a CDS-internal, pattern-based flag, the engine computes the **minimal set of
synonymous codon changes** that removes the motif *and* does not recreate it in
the local window, preserving the protein exactly, ranked by Cr adaptiveness. If
no clean fix exists it says so — it never fabricates an edit.

## Validation harness (CLAUDE.md §7)

* **Leg 1 — specificity.** 500 **real** highly-expressed Cr nuclear CDS (NCBI
  RefSeq) are clean negatives; the per-module medium/high flag rate *is* the
  false-positive rate. This **set the operating points** shipped in the profile:
  `polya` 46%→~2%, `codon` 45%→~3%, `composition` 13%→~1% (`cloning`'s ~72% is
  real Type IIS site prevalence, not noise). Regenerate with your own FASTA.
  **Held-out check:** the thresholds (tuned on 500) were re-tested on **4,500
  disjoint genes** they never saw — a 500-CDS set (`data/native_cr_cds_heldout.fasta`,
  FP 1.0–1.8%) and a 4,000-CDS set (`data/native_cr_cds_4000.fasta`, FP 2.0–4.2%).
  All three sets are pairwise-disjoint (5,000 unique real CDS). The rates stay in
  low single digits, so the operating points **generalise, not overfit**. Run:
  `synatvis validate leg1 --fasta synatvis/data/native_cr_cds_4000.fasta`.
* **Leg 2 — sensitivity by injection.** Inject a known motif into a clean parent;
  confirm the right module flags it and stays silent on the parent. Tests
  *detection*, not biological consequence. `validate leg2sweep` extends this into
  **detection curves** — signal strength vs fraction detected — so each tunable
  threshold has a sensitivity side to read against its Leg-1 false-positive rate
  (an operating point, not an FP number alone).
* **Leg 4 — cross-species discrimination.** `validate crossspecies` proves the
  tool is host-*specific*, not flag-everything: **1,000 real Cr transcripts vs
  1,020 foreign** (human, yeast, Arabidopsis; `data/foreign_cds.fasta`). The
  host-fit metric separates them at **AUC 1.000** (accuracy 97%, specificity
  99.8%, foreign flagged 99–100%). See `BACKEND_HOWITWORKS.pdf` for the chart.
* **Threshold provenance.** `THRESHOLDS.md` classifies every threshold as
  *measured* (set from Leg-1 FP), *biology/literature* (fixed, not tunable),
  *heuristic* (unswept — all in modules already ≤0.1% FP), or *policy*, with each
  module's false-positive rate measured across all 5,000 real CDS.
* **Leg 3 — literature cases, gated by compartment.** **Real** rescue pairs (NCBI
  accessions in `cases.yaml`): native AT-rich GFP vs its Cr-codon-optimized twin
  (Fuhrmann 1999; Barahimipour 2015), a real Cr gene with a BsaI site, and a real
  plastid gene **gated out** of nuclear scoring. Currently 100%/100% sens/spec on
  4 cases — expand the corpus to tighten the estimate.

## Data provenance & limitations (CLAUDE.md §2, §10)

* `data/cr_nuclear_codon.tsv` — **authentic** Kazusa *C. reinhardtii* nuclear
  codon usage (taxid 3055, 846 CDS, 420,455 codons, GC 66.3%); the validated,
  citable, genome-representative default. A larger RefSeq-derived table computed
  from **6,930 real Cr CDS (2.9M codons)** ships alongside as
  `data/cr_nuclear_codon_refseq.tsv` — it agrees with Kazusa on the preferred
  codon in 18–20/21 families (mean |Δ| 0.02–0.04; it is GC-enriched by the
  GC>55% nuclear filter). Swap it in and re-run Leg 1 to adopt it.
* `data/native_cr_cds.fasta` — **500 real** Cr nuclear CDS (NCBI RefSeq), ranked
  by codon adaptation as a proxy for highly-expressed. `data/cases.yaml` — **real**
  sequences with NCBI accessions and citations (the optimized-GFP partner is
  produced by the documented codon-optimization method, provenance stated).
* **Tier-A frontier panel** (all Cr-specific, in the `codon` module):
  `data/cr_codon_optimality.tsv` (Ikemura optimality from the 5,000 genes vs
  genome — the mRNA-stability axis; Presnyak 2015), `data/cr_codon_pairs.tsv`
  (Coleman codon-pair scores from the 5,000 genes), and `data/cr_tai_weights.tsv`
  + `data/cr_trna_counts.tsv` (tRNA Adaptation Index from **real GtRNAdb Crein5
  tRNA gene counts**, dos Reis 2004). Thresholds set from the 5,000-gene FP
  (optimality 0.9%, codon-pair 0.3%, tAI 1.9%).
* **mRNA folding ΔG axis** (in the `structure` module): a structured initiation
  window impedes ribosome loading (Kudla 2009; Zhang 2020 / LinearDesign). Uses
  ViennaRNA MFE when installed, else a dependency-free Nussinov base-pairing
  fallback (`structure_energy.py`). The start-region paired-fraction cutoff (0.70,
  ViennaRNA) gives **2.0% FP** on the 5,000 genes; a global 5' ΔG is reported as
  informational (folding stabilises the mRNA — context, not a hazard).
* Silencing is chromatin/strain-driven; the tool flags only sequence-visible risk
  (GC troughs) and points to strain (UVM4/UVM11) + intron mitigation.
* A clean report is **not** a promise of expression.

## Expression-propensity prediction (opt-in)

`scan --expression` adds a **relative expression-propensity index** (0–100) from a
**transparent ensemble of models**, each shown separately and anchored to 5,000
well-expressed native Cr genes: codon adaptation (RCA; Sharp & Li 1987), tRNA
adaptation (tAI; dos Reis 2004), silencing resistance (GC; Barahimipour 2015), and
structural stability (mRNA folding ΔG — the LinearDesign objective, Zhang 2020,
active when ViennaRNA is installed). Installed Tier-B ML models (CodonBERT, Saluki,
UTR-LM) appear as additional readouts. Hard hazard gates (premature poly(A),
GC-trough silencing, uORF, missing start) then scale the ensemble. The HTML report
draws a gauge, per-model bars, and a **per-position adaptiveness landscape**;
`scan --4d` renders that landscape **animated over translation time** — ribosomes
that dwell at slow codons and queue into jams (position × feature × structure × time).

### Protein fate, PTM, and the live cell journey

A complementary **protein-level layer** translates the CDS and reads sequence-encoded
features that decide *where the protein goes*: an N-terminal signal peptide
(secretory pathway), transmembrane helices (Kyte–Doolittle hydropathy), N-glycosylation
sequons (N-X-S/T), and cysteine/disulfide potential — yielding a predicted localisation
(secreted / membrane / cytosolic). These are heuristics (SignalP / DeepLoc / NetNGlyc
concepts); at most one gentle, clearly-labelled modifier (heavy secretory-glycosylation
load) touches the index. `scan --cell` renders an **animated cell-journey view**: an
original schematic *Chlamydomonas* cell (nucleus + pore, ER + ribosomes, Golgi,
cup-shaped chloroplast + pyrenoid, mitochondria, flagella, eyespot) through which the
construct travels in real time — transcription → poly(A)/processing → nuclear export →
translation (ribosome dwell from the per-codon adaptiveness) → folding → ER/Golgi
N-glycosylation → localisation, routed by the predicted fate. It is embedded directly in
the `--html` report, so one scan produces the whole integrative picture. Original
schematic (not microscopy); route and timing reflect the model's predictions.

**Every checkpoint now shows its parameters and scales, with citations.** As the
construct moves through each stage, a live readout panel lists the values it passes —
GC%, RCA/tAI, uORF, premature-poly(A), structure confidence, MW/pI/GRAVY, N-glyc, ELP
aggregation propensity, intein cleavage control — each with a pass/warn/fail mark, a
mini-scale, and its reference (`synatvis/journey.py` builds this in Python so it is
unit-tested, not hard-coded in JS). Each stage and organelle has its own scientific
glyph (DNA helix, mRNA, spliceosome, poly(A), ribosome, glycan tree, ER, Golgi, ELP
thermometer, intein scissors, product vial).

**Downstream recovery is part of the journey.** The view extends past localisation to a
**purification bench** (`synatvis/purification.py`): it detects tags in the construct
(His-tag, ELP `(VPGXG)ₙ`, common affinity tags) and lays out the predicted route —
**non-chromatographic ELP inverse-transition cycling** (aggregation propensity from
guest-residue hydrophobicity & chain length; Meyer & Chilkoti), **self-cleaving intein
tag removal** (N- vs C-terminal cleavage, triggers, and the mutations that switch
cleavage direction; Wood & Camarero), and **affinity/IMAC** — following the ELP-intein
self-cleaving concept of Banki, Feng & Wood 2005. Values are published relationships and
rules, clearly labelled where they need wet-lab calibration; nothing is a fitted model.

**Self-explanatory structure confidence.** `synatvis/structure_confidence.py` adds the
piece a DeepMind-style predictor lacks: it translates per-residue confidence into plain
language ("this region folds into a solid, well-defined shape"; "this stretch is likely
floppy — a good place for a linker/tag"). It reads real pLDDT when an AlphaFold3/Boltz
plugin is wired, otherwise a clearly-labelled sequence order/disorder proxy (TOP-IDP).

It is a **model, not a measured yield** — uncalibrated pending wet-lab data (the
grant's Aim 2), and it stays a separate, clearly-labelled module (the validated
diagnostic core still emits no score). It is *directionally validated*: it ranks
Cr-codon-optimised rescue partners far above their native versions (GFP 0→77,
luciferase 0→81) and native genes above foreign (median 79 vs 0; 99% vs 1% above 50).

## Construct-grammar recognizer (in progress — Stage 0)

A future layer aims to recognize the **standardized architectural grammar** of
synthetic/recombinant gene constructs (MoClo/Golden Gate "common syntax," Patron
et al. 2015) as its own thing, separate from native-gene analysis. Before any
recognizer code is trained, Stage 0 builds a ground-truth corpus under the same
no-fabrication discipline as the rest of the tool:

* **Step 1 — inclusion/exclusion criteria, frozen before collection.**
  `synatvis/data/construct_grammar/inclusion_criteria.yaml` (machine-checkable)
  and its sibling `INCLUSION_CRITERIA.md` (human rationale + citations) define
  five inclusion gates (Type-IIS/Phytobricks compatibility, functional
  validation venue, primary-source deposition, Sequence-Ontology categorizability,
  independent citability) and four exclusion gates (preprint-only, conflicting
  sequence versions, non-Type-IIS standards, no functional evidence). Enforced by
  `synatvis/construct_grammar.py`'s `evaluate_candidate()`, which every candidate
  part must pass before entering Stage 1 (segmentation). See
  `INCLUSION_CRITERIA.md` for why each gate exists and what tool-level failure
  mode it prevents (e.g. IC-3 blocks the same class of silent transcription
  errors pLannotate documents in secondary-source plasmid maps).
* **Step 2 — primary-source sequence retrieval, in progress, partly blocked.**
  `synatvis/data/construct_grammar/step2_catalog_ledger.yaml` logs 5 real,
  independently-verified catalog entries from the Chlamydomonas Resource
  Center (chlamycollection.org — the real distributor, correcting an earlier
  assumption that this was on Addgene). Actual DNA sequence text is not yet in
  hand: the founding toolkit paper (Crozet et al. 2018) is paywalled with no
  free copy; outreach to the real corresponding authors is underway (see the
  ledger's `outreach` section for the confirmed contacts and status).
* **Step 3 — Sequence Ontology vocabulary, done.**
  `synatvis/data/construct_grammar/so_vocabulary.yaml` adopts 13 real SO terms
  (promoter, five_prime_UTR, CDS, three_prime_UTR, terminator, intron,
  signal_peptide, selection_marker, engineered_tag, engineered_region, gene,
  engineered_gene, operator), each independently verified against the real
  Sequence Ontology via the EBI Ontology Lookup Service — not guessed. IC-4 in
  `evaluate_candidate()` now checks actual vocabulary membership, so a
  plausible-looking but non-adopted `SO:XXXXXXX` string correctly fails rather
  than passing on the strength of just looking like an ID.
* Steps 4–11 (provenance ledger schema, cross-toolkit deduplication,
  domestication QC, two-tier syntax/identity split, homology-aware
  generalization split, hard-negative mining from the 5,000-gene native
  corpus, frozen held-out eval set, licensing/MTA audit) are planned but not
  yet executed.

## Multiprotein complex membership (name-based, cryo-ET-anchored)

Most proteins don't work alone — they're part of a physical "team" (a complex).
`synatvis/complexome.py` cross-references a transcript's declared gene name/symbol
(from a FASTA header) against a curated list of well-established gene-name patterns
for known complex subunits (`data/complexes.yaml`): Rubisco (`rbcL`/`RBCS`),
Photosystem II (`psbA–D`), Photosystem I (`psaA`/`psaB`), ATP synthase
(`atpA/B/E/F/H/I`, `ATPC/D/G`), the cytosolic ribosome (`RPL#`/`RPS#`), and
microtubules (`TUA#`/`TUB#`). A match adds a "Multiprotein complex membership"
checkpoint citing both *why* the gene is assigned that identity (classic gene
nomenclature) and the *structural* evidence that the complex physically exists
in this organism — a 2025 large-scale cryo-electron tomography community dataset
that directly imaged ~25 real macromolecular complexes inside intact
*Chlamydomonas* cells (Chromatin-Structure-Rhythms-Lab, EMPIAR-11830). This is a
**name-based identity + context hint, not a structural measurement of the
specific scanned molecule** — the cryo-ET study confirms the complex *family*
was imaged, not this particular sequence. Deliberately conservative: only
patterns with unambiguous, well-established nomenclature are included (v1, 6
complexes); expanding to more (nucleosome, clathrin, proteasome, mitoribosome)
needs separately verified Cr-specific gene symbols first, not a guess. Two of
the six patterns (ribosome `RPL#`/`RPS#`, microtubule `TUA#`/`TUB#`) are
marked `confidence: UNVERIFIED for Chlamydomonas specifically` in the data
file — the naming convention is standard across eukaryotes generally, but a
real matching Cr NCBI record was not found when checked (2026-08-08); treat
matches on those two as lower confidence until one is found.

**Confirmed at scale** (`synatvis validate complexome`): the matcher was run
against the real 5,000-gene native Cr corpus (accession-style headers, not
gene symbols) and produced **0 false positives**, and against 13 real named
genes fetched directly from the Chlamydomonas chloroplast genome (NCBI RefSeq
NC_005353.1 — `data/named_complex_genes.fasta`), correctly matching **13/13**.
What "at scale" means here, honestly: most of the 5,000-gene corpus has no
classic gene symbol at all (NCBI's own Cr nuclear annotation mostly assigns
systematic locus IDs like `CHLRE_17g713450v5`, confirmed by direct lookup, not
assumed) — so this is a specificity check (proving no spurious matches), not
5,000 real complex identifications. There is no "model" being trained here;
this is a deterministic name-pattern lookup, not a statistical/ML component.

## Molecular dynamics — opt-in slot, never simulated internally

Molecular dynamics (MD) predicts how a protein's atoms move over time and is a
fundamentally different computation from everything else here: it needs a full
3D starting structure, a force field, and heavy compute (real runs take hours
to weeks per system on GPUs/clusters). No MD run is ever "perfect" — every one
is a real, known approximation, true for every lab, not a gap specific to this
tool. SynAT.Vis does not run MD itself. `plugins/moleculardynamics.py` follows
the same honest external-command contract as every other Tier-B plugin: set
`MDSIM_CMD` to *your own* installed MD engine's wrapper (GROMACS/OpenMM/NAMD/
AMBER), which reads a protein FASTA on stdin and prints JSON
(`rmsd_nm`, `radius_of_gyration_nm`, `sim_time_ns`, `force_field`). Nothing
runs and no numbers appear unless a real command is configured.

## Construct-grammar recognizer (in progress — Stage 0)

* **HTML report** (`scan --html`) — a polished, self-contained, theme-aware page:
  a colour-coded verdict, plain-language findings, suggested fixes, and cited
  evidence. The drag-and-drop launchers (`SCAN-a-gene`) open it automatically.
  Text (`--text`, default), plain-language (`--plain`), and JSON (`--json`) forms
  are also available.
* **Opt-in ML plugins** (`--plugins`, `synatvis plugins`) — wrappers for trained
  models, all **experimental and never part of the validated result** (the stdlib
  core never imports them):
  * *mRNA / expression:* **LinearDesign** (runs on ViennaRNA), **CodonBERT**,
    **Saluki**, **UTR-LM**, and the frontier readouts **Borzoi**, **APARENT2**,
    **iCodon** (transfer-learning, not Cr-calibrated).
  * *validated protein-fate / splice tools* that supersede the `ptm.py` and splice
    heuristics when installed: **SignalP 6.0**, **DeepLoc 2.1**, **DeepTMHMM**,
    **Pangolin**.
  * *NVIDIA BioNeMo foundation-model seams* (GPU / NIM microservices): **CodonFM**
    (codon foundation model + SAE interpretability), **OpenFold3** (folding pLDDT/pTM),
    **Evo 2** (genomic likelihood, variant effect, generative redesign).
  * *Molecular dynamics* (opt-in slot, never simulated internally): wraps *your own*
    installed MD engine (GROMACS/OpenMM/NAMD/AMBER) — see "Molecular dynamics" above.
  * Each external-model adapter runs the user's *own* installed inference command
    (gated on an env var; reads FASTA on stdin, prints JSON) — so it runs a real
    model and never fabricates a score. See **`SynAT.Vis_AI_Integration_Roadmap.pdf`**
    (frontier-model & commercialisation plan) and
    **`SynAT.Vis_BioNeMo_Integration_Roadmap.pdf`** (the NVIDIA BioNeMo mapping).

## Tests

```bash
python tests/test_synatvis.py     # standalone runner (no pytest needed)
pytest -q                         # or under pytest
```

## Extending to another compartment / lineage

Add a profile YAML under `synatvis/profiles/` supplying every schema key
(`profiles/_schema.py`). **No `modules/*.py` edits are required** — biology lives
in data. The v1 gate is `compartment`; add a `lineage` gate if cases are later
pooled across algal models.
