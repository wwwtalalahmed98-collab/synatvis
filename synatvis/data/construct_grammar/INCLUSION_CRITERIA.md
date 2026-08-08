# Construct-Grammar Corpus — Inclusion/Exclusion Criteria (Stage 0, Step 1)

**Status:** frozen 2026-08-07, before any part collection began.
**Machine-readable source of truth:** [`inclusion_criteria.yaml`](inclusion_criteria.yaml). This document is the human-readable rationale for the same IDs — the two must never drift; if a criterion changes, bump `version` in the YAML and update this file in the same commit.

## Why criteria are frozen before collection, not derived from it

A corpus assembled by "collect what's easy to find, then describe what we collected" always overfits to whichever labs are best at self-promotion or happen to use the most searchable repository. Locking IC-1..IC-5 and EX-1..EX-4 *before* touching the literature means every candidate part is measured against a fixed yardstick, not a yardstick quietly redrawn around whatever was already gathered. This mirrors how `THRESHOLDS.md` values were fixed by measurement policy before being applied to the 5,000-CDS scan, rather than picked after seeing which genes passed.

## The criteria

### Inclusion

| ID | Name | What it checks |
|----|------|-----------------|
| IC-1 | assembly_standard_compatible | Flanked by BsaI/BsmBI/BpiI Type IIS sites, and/or overhangs match the Phytobricks 12-position common syntax (Patron et al. 2015, *New Phytologist*; Weber et al. 2011, *PLoS ONE*). |
| IC-2 | functionally_validated | Function demonstrated in *C. reinhardtii* (nuclear or chloroplast), or in another organism via a Phytobrick-compliant system (Crozet et al. 2018, *ACS Synth Biol*; Lauersen et al. 2015, *Appl Microbiol Biotechnol*). |
| IC-3 | primary_source_deposited | Retrievable from Addgene (plasmid + GenBank map) or an NCBI accession — not just a supplementary table or figure (McGuffie & Barrick 2021, *NAR*, on pLannotate's documented plasmid-annotation errors). |
| IC-4 | categorizable | Assignable an unambiguous Sequence Ontology term (Eilbeck et al. 2005, *Genome Biology*; SBOL3 requires SO-backed roles, McLaughlin et al. 2020). |
| IC-5 | independently_citable | Has a DOI/PMID/Addgene ID a third party can check without emailing the authors. |

### Exclusion

| ID | Name | What it checks |
|----|------|-----------------|
| EX-1 | preprint_without_deposit | Preprint-only, no deposited sequence → tagged `pending`, not permanently rejected. |
| EX-2 | conflicting_sequence | Supplementary sequence disagrees with the deposited one → excluded until manually reconciled; never auto-resolved. |
| EX-3 | non_type_iis_standard | Gateway-only or other non-Type-IIS systems → out of scope for *this* corpus, routed to a separate future recombinant-tag-grammar corpus (Luo et al. 2025 pNX; Somogyi et al. 2019). |
| EX-4 | no_functional_evidence | Deposited sequence, no functional data → allowed only into a secondary "candidate" tier, never primary ground truth. |

## How each criterion changes tool performance (not just corpus hygiene)

The corpus these criteria produce becomes training/validation data for the construct-grammar recognizer (Stages 1–2) and eventually the host-context predictor (Stage 3). Every laxity here becomes a systematic error mode later, not a random one — because the recognizer will learn to reproduce whatever the corpus contains, including its mistakes, consistently.

- **IC-1 (assembly compatibility).** This is the literal definition of the signal being detected: a fixed Type IIS junction grammar. If non-compliant sequences (e.g. traditional restriction-site clones with incidental "BsaI-like" motifs) leak in as positives, the recognizer learns a wider, wrong boundary — and in production it will call spurious junctions inside ordinary native CDS, inflating **false-positive fusion-site calls** on real user-submitted constructs. This risk is exactly why the codebase already treats false positives as a first-class, measured cost (`THRESHOLDS.md`'s per-parameter FP-on-5000-genes column) — a construct-grammar false positive is the same category of error, just at the part-boundary level instead of the codon-usage level.
- **IC-2 (functional validation, Cr vs. non-Cr split).** This is what lets the corpus separate "the fusion-site grammar" (host-agnostic, learnable from any Phytobrick-compliant system) from "which parts actually work in *Chlamydomonas*" (host-specific, must stay Cr-anchored). Collapsing this distinction would let the tool report a *Physcomitrella*-validated promoter as "known-good in Cr" — a false reassurance a researcher could act on directly. This is the same class of hazard SynAT.Vis's core design already guards against by keeping the host profile (`cr_nuclear`) explicit rather than assuming any expression system's data is fungible with any other's.
- **IC-3 (primary-source deposition).** pLannotate's own authors document that plasmid maps compiled from secondary sources (tables, figures) frequently carry transcription errors that persist for years because no one re-derives them from the primary record. If SynAT.Vis's ground truth inherited such an error, every future construct that reuses that "part" would be silently mis-annotated — a compounding error that gets *harder* to detect the more the tool is used, not easier. Requiring a re-fetchable accession is what makes Stage 4 auditing possible at all: a wrong call can always be traced back to its source and corrected.
- **IC-4 (SO categorizability).** Without a controlled vocabulary, the model's output classes drift toward whatever ad hoc labels felt convenient during annotation — which breaks both human interpretability of a flag ("what does 'part_type: misc_2' mean?") and any future interoperability with SBOL-based design tools, which require SO-backed roles by spec. This directly affects whether SynAT.Vis's future output can be handed to a wet-lab collaborator's SBOL-aware assembly software without manual relabeling.
- **IC-5 (independent citability).** Enforces the project's standing rule (no assertion the user can't independently verify) at the corpus level, not just the report level. A part-identity claim baked into training data is *harder to audit* than a live flag in a report, so the citability bar here is intentionally at least as strict as the report-level one.
- **EX-1/EX-2 (preprint-only / conflicting sequence).** Both prevent *unresolved uncertainty from being encoded as certainty*. A model trained on a "probably this sequence" data point has no way to represent that uncertainty at inference time — it will report the same false confidence on real user data that it was trained to project. Tagging as `pending` (EX-1) or blocking until reconciled (EX-2) keeps the ground-truth tier free of hidden guesses.
- **EX-3 (non-Type-IIS out of scope).** Prevents grammar contamination between two structurally different systems (Golden Gate junction grammar vs. Gateway attB/attP recombination, or N-/C-terminal recombinant tag grammar). Mixing them would blur the junction-detection boundary IC-1 is trying to keep sharp. Routing excluded parts to a *named* future corpus (rather than silently dropping them) keeps the tool's scope gap visible and plannable instead of an invisible blind spot discovered later by a confused user.
- **EX-4 (candidate tier for unvalidated deposits).** Keeps "this sequence exists" and "this sequence has been shown to work" as separable claims all the way through training, so the eventual output can distinguish "structurally a valid part" from "known functional part" — a distinction a researcher deciding whether to build on a given part actually needs.

**Net effect on the tool:** these eight criteria are what convert "we found some MoClo parts online" into a corpus whose every entry has a traceable, independently checkable provenance chain and an explicit confidence tier. That is the same property that already makes SynAT.Vis's red-flag scanner trustworthy (measured thresholds, cited sources, no fabricated composite score) — Stage 0 is applying that same discipline one level up, to the training data itself, before any model logic is written.

## Machine enforcement

Encoded as `synatvis.construct_grammar.evaluate_candidate()`, which loads criterion text/citations from `inclusion_criteria.yaml` (so code and citable rationale can never silently diverge) and checks a `CandidatePart` record against every IC/EX rule, returning `INCLUDE` / `CANDIDATE_TIER_ONLY` / `PENDING` / `EXCLUDE` with a per-criterion reason. No part enters Stage 1 (segmentation) without passing through this function.
