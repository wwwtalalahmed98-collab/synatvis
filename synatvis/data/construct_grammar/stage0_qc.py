"""Stage 0, Steps 4-10 -- corpus QC, deduplication, and a frozen evaluation split.

Runs over the corpus built by build_corpus.py and produces the artefacts the later
recognizer stages need before any model is trained. Every step reports what it
actually found; nothing is silently repaired.

  Step 4  provenance ledger      -- every record traced to a source, checksum and tier
  Step 5  deduplication          -- exact and near-duplicate detection within/across sources
  Step 6  domestication QC       -- INTERNAL Type IIS sites that a real MoClo part must not carry
  Step 8  homology-aware split   -- clusters kept whole, so near-duplicates cannot straddle
                                    the train/test boundary and inflate apparent accuracy
  Step 9  hard negatives         -- native Cr genes that look most construct-like
  Step 10 frozen evaluation set  -- deterministic, written once, hash-stamped

Step 7 (two-tier syntax/identity split) is already enforced by build_corpus.py's tiers.
Step 11 (licensing/MTA audit) is a human/legal task and is deliberately not automated.

Run:
    python synatvis/data/construct_grammar/stage0_qc.py
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CORPUS_DIR = os.path.join(HERE, "moclo_corpus")
OUT_DIR = os.path.join(HERE, "stage0")

TYPE_IIS = {"BsaI": "GGTCTC", "BpiI": "GAAGAC", "BsmBI": "CGTCTC", "SapI": "GCTCTTC"}
_TRANS = str.maketrans("ACGTN", "TGCAN")
K = 12          # k-mer size for near-duplicate detection
NEAR_DUP = 0.90  # containment at/above this counts as a near-duplicate


def _rc(s):
    return s[::-1].translate(_TRANS)


def read_corpus():
    fa = os.path.join(CORPUS_DIR, "corpus.fasta")
    mf = os.path.join(CORPUS_DIR, "corpus_manifest.json")
    if not (os.path.isfile(fa) and os.path.isfile(mf)):
        return None, None
    seqs, cur = {}, None
    with open(fa, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
            elif cur:
                seqs[cur].append(line)
    seqs = {k: "".join(v) for k, v in seqs.items()}
    with open(mf, encoding="utf-8") as fh:
        man = json.load(fh)
    return seqs, man


# ---------------------------------------------------------------- step 5
def comparable_sequence(seqs):
    """The sequence that should actually be compared: the EXCISED PART, not the plasmid.

    Whole-plasmid similarity is dominated by the shared MoClo acceptor backbone, not by
    the biology. Measured on two genuinely unrelated Cr promoter parts (pCM0-001 PSAD vs
    pCM0-002 AR): whole-plasmid k-mer containment 0.314, excised-part containment 0.000.
    Clustering on plasmids therefore groups records by which vector they were cloned into,
    which is exactly the wrong thing -- it made 111 of 115 cr_primary records collapse into
    a single cluster and wrecked the train/test split. Records that cannot be cleanly
    excised fall back to the full sequence and are flagged.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location("_bc", os.path.join(HERE, "build_corpus.py"))
    bc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bc)
    out, fell_back = {}, []
    for rid, s_ in seqs.items():
        ex = bc.excise_bsai_part(s_)
        if ex and len(ex[0]) >= 60:
            out[rid] = ex[0]
        else:
            out[rid] = s_
            fell_back.append(rid)
    return out, fell_back


def deduplicate(seqs):
    """Exact duplicates by hash; near-duplicates by k-mer containment on excised parts."""
    cmp_seq, fell_back = comparable_sequence(seqs)
    by_hash = {}
    for rid, s in cmp_seq.items():
        by_hash.setdefault(hashlib.sha256(s.encode()).hexdigest(), []).append(rid)
    exact = {h: ids for h, ids in by_hash.items() if len(ids) > 1}

    kmers = {rid: {s[i:i + K] for i in range(0, len(s) - K + 1, 4)}
             for rid, s in cmp_seq.items()}
    ids = sorted(seqs)
    near = []
    for i, a in enumerate(ids):
        ka = kmers[a]
        if not ka:
            continue
        for b in ids[i + 1:]:
            kb = kmers[b]
            if not kb:
                continue
            small, big = (ka, kb) if len(ka) <= len(kb) else (kb, ka)
            c = len(small & big) / len(small)
            if c >= NEAR_DUP:
                near.append({"a": a, "b": b, "containment": round(c, 3)})
    # NOTE: the full pair list is returned for clustering. Truncating here and then
    # clustering on the truncated list would silently produce an INCOMPLETE homology
    # split -- exactly the failure this step exists to prevent. Truncation happens only
    # when writing the JSON report.
    return {"exact_groups": exact, "n_exact_groups": len(exact),
            "near_duplicate_pairs": near, "n_near_duplicate_pairs": len(near),
            "compared_on": "excised part where possible",
            "n_fell_back_to_full_sequence": len(fell_back)}


# ---------------------------------------------------------------- step 6
def domestication_qc(seqs, man):
    """Do the EXCISED PARTS carry Type IIS sites that assembly would cut?

    An earlier version ran this on the whole plasmid and flagged 115/115 cr_primary
    records -- which was meaningless, because a plasmid legitimately carries Type IIS
    sites in its backbone and at the part flanks. The only meaningful question is
    whether the PART ITSELF, once excised, still contains a site. That is a genuine
    domestication failure. Records whose part cannot be cleanly excised are reported
    separately as "not assessable" rather than counted either way.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location("_bc", os.path.join(HERE, "build_corpus.py"))
    bc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bc)
    tier = {r["id"]: r["tier"] for r in man["records"]}

    failed, not_assessable, clean = [], [], 0
    for rid, s_ in seqs.items():
        if tier.get(rid) == "synthetic":
            continue                      # assembled in silico; not a deposited part
        ex = bc.excise_bsai_part(s_)
        if not ex:
            not_assessable.append(rid)
            continue
        insert = ex[0]
        core = insert[4:-4] if len(insert) > 8 else insert   # drop the fusion overhangs
        internal = {}
        for enz, site in TYPE_IIS.items():
            n = len(re.findall(site, core)) + len(re.findall(_rc(site), core))
            if n:
                internal[enz] = n
        if internal:
            failed.append({"id": rid, "tier": tier.get(rid, "?"),
                           "internal_sites": internal, "part_bp": len(core)})
        else:
            clean += 1
    by_tier = {}
    for r in failed:
        by_tier[r["tier"]] = by_tier.get(r["tier"], 0) + 1
    return {"assessed_clean": clean, "n_failed": len(failed), "failed_by_tier": by_tier,
            "n_not_assessable": len(not_assessable),
            "not_assessable_note": "part could not be cleanly excised (not exactly one "
                                   "forward and one reverse BsaI site); neither pass nor fail",
            "examples": failed[:25]}


# ---------------------------------------------------------------- step 8
def homology_clusters(seqs, near_pairs):
    """Union-find over near-duplicate pairs, so related records stay together."""
    parent = {rid: rid for rid in seqs}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for p in near_pairs:
        if p["a"] in parent and p["b"] in parent:
            union(p["a"], p["b"])
    clusters = {}
    for rid in seqs:
        clusters.setdefault(find(rid), []).append(rid)
    return list(clusters.values())


# ---------------------------------------------------------------- step 10
def frozen_split(clusters, man, test_fraction=0.2, seed=20260813):
    """Frozen evaluation split -- and a hard constraint discovered while building it.

    THE SYNTHETIC TIER CANNOT BE EVALUATED AGAINST. Synthetic records are assembled
    FROM the cr_primary parts, so they literally contain cr_primary sequence. Homology
    clustering correctly refuses to separate them: one cluster spans 599 synthetic and
    99 cr_primary records. Any split that put synthetic in test and cr_primary in train
    (or the reverse) would leak, and the leak would look like accuracy.

    So the evaluation set is drawn ONLY from real deposited records (cr_primary and
    syntax_only). The synthetic tier is training augmentation only -- which is what it
    was always for: its value is exact junction labels for segmentation, not held-out
    evidence. This is recorded as a constraint, not worked around.
    """
    import random
    rng = random.Random(seed)
    tier = {r["id"]: r["tier"] for r in man["records"]}

    evaluable = {"cr_primary", "syntax_only"}
    train, test, aug = [], [], []
    by_tier_clusters = {}
    for c in clusters:
        # a cluster is evaluable only if every member is a real deposited record
        ts = {tier.get(x, "?") for x in c}
        if ts <= evaluable:
            key = sorted(ts)[0]
            by_tier_clusters.setdefault(key, []).append(c)
        else:
            aug.extend(c)          # mixed or synthetic-containing -> augmentation only

    for t, cls in sorted(by_tier_clusters.items()):
        cls = sorted(cls, key=lambda c: (-len(c), c[0]))
        rng.shuffle(cls)
        total = sum(len(c) for c in cls)
        target = int(round(test_fraction * total))
        n = 0
        for c in cls:
            if n < target:
                test.extend(c)
                n += len(c)
            else:
                train.extend(c)

    return {"seed": seed, "test_fraction_requested": test_fraction,
            "n_train": len(train), "n_test": len(test),
            "n_augmentation_only": len(aug),
            "test_fraction_actual": round(len(test) / max(1, len(train) + len(test)), 3),
            "train_by_tier": _count(train, tier), "test_by_tier": _count(test, tier),
            "augmentation_by_tier": _count(aug, tier),
            "constraint": "synthetic records contain cr_primary sequence, so they are "
                          "training augmentation only and never appear in the evaluation set",
            "train": sorted(train), "test": sorted(test)}


def _count(ids, tier):
    out = {}
    for i in ids:
        out[tier.get(i, "?")] = out.get(tier.get(i, "?"), 0) + 1
    return out


# ---------------------------------------------------------------- step 9
def hard_negatives(limit=200):
    """Native Cr genes that most resemble a construct: they carry Type IIS sites.

    These are the negatives a junction recognizer will find hardest, because they
    contain the very motif that defines a MoClo junction while being ordinary genes.
    """
    sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..")))
    from synatvis.seqio import read_fasta
    from synatvis.profiles import PACKAGE_DIR
    out = []
    for fname in ["native_cr_cds.fasta", "native_cr_cds_heldout.fasta"]:
        p = os.path.join(PACKAGE_DIR, "data", fname)
        if not os.path.isfile(p):
            continue
        for name, seq in read_fasta(p):
            s = seq.strip().upper()
            if len(s) < 300:
                continue
            n = sum(s.count(v) + s.count(_rc(v)) for v in TYPE_IIS.values())
            if n >= 2:
                out.append({"id": name.split()[0], "type_iis_sites": n, "length": len(s)})
            if len(out) >= limit:
                return out
    return out


def main() -> int:
    seqs, man = read_corpus()
    if not seqs:
        print("corpus not built. Run build_corpus.py first.")
        return 1
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"corpus: {len(seqs)} records, tiers {man['by_tier']}\n")

    print("Step 5 -- deduplication")
    dd = deduplicate(seqs)
    print(f"   exact duplicate groups : {dd['n_exact_groups']}")
    print(f"   near-duplicate pairs   : {dd['n_near_duplicate_pairs']} "
          f"(k={K}, containment >= {NEAR_DUP})")
    print(f"   compared on EXCISED PARTS; {dd['n_fell_back_to_full_sequence']} records "
          f"fell back to full sequence")

    print("\nStep 6 -- domestication QC")
    dq = domestication_qc(seqs, man)
    print(f"   excised parts assessed clean        : {dq['assessed_clean']}")
    print(f"   parts with an INTERNAL Type IIS site: {dq['n_failed']}")
    for t, n in sorted(dq["failed_by_tier"].items()):
        print(f"      {t:14s} {n}")
    print(f"   not assessable (no clean excision)  : {dq['n_not_assessable']}")

    print("\nStep 8 -- homology-aware clustering")
    cl = homology_clusters(seqs, dd["near_duplicate_pairs"])
    multi = [c for c in cl if len(c) > 1]
    print(f"   clusters: {len(cl)}  (of which multi-record: {len(multi)}, "
          f"largest: {max((len(c) for c in cl), default=0)})")

    print("\nStep 10 -- frozen evaluation split (clusters kept whole)")
    sp = frozen_split(cl, man)
    print(f"   train {sp['n_train']}  test {sp['n_test']}  "
          f"(actual test fraction {sp['test_fraction_actual']})")
    print(f"   augmentation-only (never evaluated): {sp['n_augmentation_only']}")
    print(f"   train tiers {sp['train_by_tier']}")
    print(f"   test  tiers {sp['test_by_tier']}")

    print("\nStep 9 -- hard negatives from the native corpus")
    hn = hard_negatives()
    print(f"   native genes carrying >=2 Type IIS sites: {len(hn)}")

    ledger = {
        "step4_provenance": {"sources": man["sources"], "by_tier": man["by_tier"],
                             "n_total": man["n_total"],
                             "n_quarantined": man.get("n_quarantined", 0)},
        "step5_deduplication": {**{k: v for k, v in dd.items() if k != "near_duplicate_pairs"},
                                "near_duplicate_pairs_sample": dd["near_duplicate_pairs"][:200]},
        "step6_domestication_qc": dq,
        "step8_clusters": {"n_clusters": len(cl), "n_multi_record": len(multi)},
        "step9_hard_negatives": {"n": len(hn), "examples": hn[:25]},
        "step10_frozen_split": sp,
        "step11_licensing": "NOT AUTOMATED -- human/legal review required. See "
                            "step2_catalog_ledger.yaml redistribution_decision.",
    }
    path = os.path.join(OUT_DIR, "stage0_ledger.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(ledger, fh, indent=1)
    digest = hashlib.sha256(json.dumps(sp["test"], sort_keys=True).encode()).hexdigest()
    print(f"\nwrote {path}")
    print(f"frozen test-set digest: {digest[:32]}...")
    print("Re-running with the same corpus reproduces this split exactly.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
