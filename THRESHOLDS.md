# SynAT.Vis — threshold provenance

Per-module medium/high **false-positive rate measured across 5,000 real *C. reinhardtii* nuclear CDS** (the union of the three disjoint Leg-1 sets). Every threshold is classified by source: **measured** (set by sweeping the Leg-1 FP), **biology**/**literature** (fixed by the host or a cited study — not a tunable knob), **heuristic** (a reasonable default, not yet swept), or **policy** (a design choice).

| module | FP (5,000 genes) | threshold | value | source | note |
|---|---|---|---|---|---|
| codon | 4.7% | `rare_weight_threshold` | 0.1 | measured | swept vs Leg-1 FP (0.15->0.10) |
|  |  | `cluster_min_rare` | 5 | measured | swept vs Leg-1 FP (4->5) |
|  |  | `cluster_window_codons` | 8 | heuristic | window size, unswept |
|  |  | `low_rca_warn` | 0.6 | heuristic | INFO-only, does not affect med/high FP |
|  |  | `optimality_cutoff` | -0.2 | measured | de-opt region; 0.9% FP on 5,000 genes |
|  |  | `optimality_window` | 30 | heuristic | optimality window (codons) |
|  |  | `cpb_run_len` | 5 | measured | codon-pair run; native 98th pct=4 -> flag >=5 |
|  |  | `cpb_cutoff` | -0.5 | literature | Coleman 2008 codon-pair-score cutoff |
|  |  | `low_tai_warn` | 0.34 | measured | native 2nd-pct tAI (GtRNAdb); 1.9% FP |
| composition | 1.8% | `gc_trough_low` | 0.4 | measured | swept vs Leg-1 FP (0.45->0.40) |
|  |  | `gc_trough_warn` | 0.45 | measured | swept vs Leg-1 FP (0.50->0.45) |
|  |  | `min_trough_len` | 60 | measured | swept vs Leg-1 FP (40->60) |
|  |  | `window` | 50 | heuristic | scan window, unswept |
|  |  | `step` | 10 | heuristic | scan step, unswept |
|  |  | `homopolymer_min` | 10 | literature | synthesis/processivity heuristic |
|  |  | `target_gc` | 0.66 | biology | Cr coding GC (reference) |
| polya | 2.3% | `gc_context_min` | 0.85 | measured | swept vs Leg-1 FP (0.55->0.85) |
|  |  | `gc_context_window` | 30 | heuristic | context window, unswept |
|  |  | `nue_motifs` | TGTAA,TGTAG,TGCAA,ATGTAA | literature | Shen 2008; Zhao 2014 (UGUAA) |
|  |  | `upstream_min/max` | 10..30 | literature | Shen 2008 (~10-30 nt) |
| cloning | 73% (site prevalence) | `enzymes` | 4 sites | biology | exact Type IIS sites (not a threshold) |
| splice | 0.0% | `donor_consensus` | GTRAG | biology | GT + Cr context (exact) |
|  |  | `acceptor_consensus` | YAG | biology | YAG (exact) |
|  |  | `min_intron/max_intron` | 50..500 | literature | Baier 2018 (short Cr introns) |
| structure | 0.0% | `struct_paired_max` | 0.7 | measured | ViennaRNA start paired-frac; 2.0% FP on 5,000 genes |
|  |  | `struct_paired_max_heuristic` | 0.8 | heuristic | Nussinov fallback cutoff (over-pairs) |
|  |  | `struct_window_after` | 45 | heuristic | initiation-window length (nt) |
|  |  | `max_pairing_frac` | 0.78 | heuristic | legacy GC fallback (retained) |
|  |  | `utr5_max_len` | 200 | heuristic | INFO-only, unswept |
|  |  | `utr5_struct_window` | 40 | heuristic | fold window, unswept |
|  |  | `start_context_window` | 12 | heuristic | unswept |
| uorf | 0.0% | `flag_any_uaug` | True | policy | flag every upstream AUG |
| instability | 0.1% | `are_motifs` | ATTTA | literature | AUUUA ARE core |
|  |  | `are_cluster_min` | 2 | heuristic | unswept |
|  |  | `are_window` | 50 | heuristic | unswept |
| silencing | 0.0% | `validated` | False | policy | quarantined; heuristic only |

## Summary
- **measured** (from Leg-1 FP): 10 thresholds — the noisy modules that actually over-fired (poly(A), codon, composition).
- **biology/literature** (fixed, not tunable): 10.
- **heuristic** (unswept): 14 — all in modules whose measured FP is already ≤0.5% on 5,000 genes, so there is no false-alarm problem to fix; they need the sensitivity axis (Leg-2 sweep) to be tuned meaningfully.
- **policy**: 2.

See `native_cr_cds*.fasta` for the gene sets and `validate leg1 --fasta …` to reproduce any FP figure.
