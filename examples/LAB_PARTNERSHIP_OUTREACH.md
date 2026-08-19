# Lab partnership outreach — targets and draft approach

## Why this exists

The expression-propensity index is *directionally validated* with ordering confirmed
on the codon axis, and it has one documented open defect on the intron axis. What it
does **not** have is a single measured expression value produced for the purpose of
calibrating it. No amount of computation fixes that — it requires a wet lab.

This is the one gap on the roadmap that cannot be closed by more work at the keyboard.

## What to actually ask for

Be concrete and small. The failure mode of research outreach is asking for a
collaboration in the abstract; the successful version asks for one specific, cheap thing.

**The minimum useful ask:** measured expression for **two constructs that differ only in
synonymous codons** — the same protein, one Cr-optimised and one not. Even a single
matched pair, with an error estimate, converts one axis of the index from directional
to genuinely anchored.

**The better ask, if they are willing:** three or four constructs spanning a range of
codon adaptation, ideally including one that holds GC constant while varying RCA (the
design Barahimipour et al. 2015 used). That would let the index be checked against a
real gradient rather than a single contrast.

**What you can offer in return, honestly:**
- Free, citable pre-experiment screening of their constructs
- A named acknowledgement, or co-authorship if they contribute experimental design
- The tool is source-available, so nothing about the collaboration locks them in
- You will publish the calibration result whether or not it flatters the tool

That last point matters more than it looks. A lab is being asked to spend bench time
validating someone else's model; the credible version of that offer says the result
gets published either way.

## Realistic target types, in order of likelihood

1. **Groups already publishing Cr expression measurements.** They have the assay
   running. The marginal cost of two more constructs is low. The authors of the papers
   already used as calibration anchors are the obvious first approach — they have
   demonstrably done this exact measurement.
2. **Algal synthetic-biology toolkit developers.** The Chlamydomonas MoClo community
   is small and has already shown willingness to respond — the toolkit-sequence request
   was answered.
3. **Bioinformatics groups with a wet-lab partner.** They may value the tool itself and
   have an easier internal route to a few measurements.
4. **Local university labs.** Geography helps more than it should for a first
   collaboration. A group within travel distance can be visited.

## Draft email

> **Subject:** Two constructs' expression data — calibrating an open Chlamydomonas design tool
>
> Dear Professor [NAME],
>
> I am [Talal Ahmed], co-founder of SynAT.Vis, an open, source-available diagnostic
> that reads a designed gene cassette and reports why it may fail to express in
> *Chlamydomonas reinhardtii* nuclear transformants. It is rule-based and cited
> rather than a trained model, and it is validated on 5,000 real Cr genes.
>
> I am writing with a small, specific request rather than a general collaboration
> proposal.
>
> The tool's expression-propensity module is currently *directionally* validated —
> it correctly ranks codon-optimised constructs above their native counterparts —
> but it has never been checked against a measurement made for that purpose. A
> single matched pair would change that: two constructs encoding the same protein
> and differing only in synonymous codon usage, with measured expression and an
> error estimate.
>
> If your group already runs this assay, the marginal cost of including such a pair
> may be small. I would be glad to:
>
> - screen any constructs you are designing, free and without obligation;
> - acknowledge or co-author, as you prefer;
> - publish the calibration result regardless of whether it supports the tool.
>
> I should be straightforward that this is an early-stage project from a small team,
> and that the tool has known open defects, which are documented publicly in the
> repository rather than hidden.
>
> Repository: https://github.com/wwwtalalahmed98-collab/synatvis
>
> If this is not something your group has capacity for, I would be grateful for any
> pointer to a lab that might.
>
> With thanks for your time,
>
> [Talal Ahmed]
> [email] · [affiliation, if any]

## Notes on sending it

- **One lab at a time**, not a mass mailing. Researchers can tell.
- **Read one of their recent papers first** and reference it in a single specific
  sentence. Generic praise is transparent.
- **Expect no reply from most.** The toolkit-sequence request that succeeded got a
  two-sentence answer — which was all it needed to be.
- **A short, specific ask outperforms a long, impressive one.** Resist the urge to
  explain the whole project.
- **Do not overstate.** Saying "validated on 5,000 genes" is true. Saying "calibrated"
  is not, and a reviewer in this field will notice immediately.

## Status

Not yet sent. Sending is a human action — the drafting is done here; the relationship
is not something this repository can build.
