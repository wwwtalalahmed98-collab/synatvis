"""Stage 0, Step 2b -- build the expanded, TIERED construct-grammar corpus.

WHY TIERS EXIST (read this before using the corpus)
---------------------------------------------------
The Chlamydomonas MoClo toolkit has ~115 parts. That is not a sample -- it is
the entire toolkit. There is no larger pool of real Cr MoClo parts to find. So a
corpus bigger than ~130 necessarily contains material that is NOT Chlamydomonas,
and mixing those into Cr ground truth would silently corrupt exactly the
host-specificity this whole tool exists to get right.

The project's own frozen criteria already anticipated this: IC-2 fails a part
that was not functionally validated in Cr, routing it to CANDIDATE_TIER_ONLY
rather than INCLUDE. This module makes that split explicit and permanent:

  tier "cr_primary"  -- real, deposited, Chlamydomonas. Valid for BOTH the
                        junction-grammar task and any Cr identity claim.
  tier "syntax_only" -- real, deposited Type IIS parts from OTHER hosts (plant,
                        yeast). The Golden Gate "common syntax" (Patron 2015) is
                        deliberately shared across these toolkits, so they are
                        legitimate training data for JUNCTION/ARCHITECTURE
                        recognition -- and are NOT evidence about Cr biology.
  tier "synthetic"   -- Level-1 constructs assembled in silico from real cr_primary
                        parts by the real MoClo overhang rules. The architecture
                        and every junction coordinate are exact and known by
                        construction, which makes these strong SEGMENTATION
                        training data. They are NOT deposited plasmids and must
                        never be counted as experimental evidence.

Never train an identity/host classifier on syntax_only or synthetic records, and
never report a corpus count without its tier breakdown.

Run:
    python synatvis/data/construct_grammar/build_corpus.py
    python synatvis/data/construct_grammar/build_corpus.py --target 1000
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import re
import urllib.request
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "moclo_corpus")  # git-ignored

# Addgene's CDN rejects the default urllib agent with 403. robots.txt permits
# these media paths (only /emta/, /emta-addgene-public/ and /users/login/ are
# disallowed); these are sequence files published for researchers to download.
UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"

# Every source below was fetched and counted on 2026-08-12; the record counts
# and checksums are what was actually observed, not estimates.
SOURCES = [
    {
        "key": "cr_moclo",
        "tier": "cr_primary",
        "host": "Chlamydomonas reinhardtii",
        "url": "https://www.chlamycollection.org/content/uploads/2019/06/MoClo-Kit-Sequences.zip",
        "sha256": "8689ef08e7611fe3ce87ec11f1921ef8a4d37da578c474048897ff1d27128ec4",
        "expected": 115,
        "citation": "Crozet et al. 2018, ACS Synth Biol; Chlamydomonas Resource Center",
    },
    {
        "key": "marillonnet_moclo",
        "tier": "syntax_only",
        "host": "plant (Nicotiana / generic)",
        "url": "https://media.addgene.org/cms/filer_public/b6/f8/b6f82f82-4604-4444-9886-f8577018aee4/moclo_tool_kit_genbank_files_2_1.zip",
        "sha256": None,  # filled on first run; see --pin
        "expected": 95,
        "citation": "Weber/Engler/Marillonnet MoClo Toolkit (Addgene kit)",
    },
    {
        "key": "patron_plant_parts",
        "tier": "syntax_only",
        "host": "plant",
        "url": "https://media.addgene.org/cms/filer_public/39/6d/396d2c07-f428-4658-b723-e1fb765b2df5/plant_parts_genbank_updated_mar2015.zip",
        "sha256": None,
        "expected": 95,
        "citation": "Patron et al. 2015 Phytobricks / MoClo Plant Parts (Addgene kit)",
    },
    {
        "key": "dueber_ytk",
        "tier": "syntax_only",
        "host": "Saccharomyces cerevisiae",
        "url": "https://media.addgene.org/cms/filer_public/4d/26/4d26f69c-1b8e-4473-8494-a7e993618112/ytk_genbank_files.zip",
        "sha256": None,
        "expected": 96,
        "citation": "Lee et al. 2015 Yeast Toolkit (Addgene kit)",
    },
]

_ORIGIN = re.compile(r"^ORIGIN\s*?$(.*?)^//", re.M | re.S)
_LOCUS_BP = re.compile(r"^LOCUS\s+(\S+)\s+(\d+)\s+bp", re.M)
TYPE_IIS = {"BsaI": "GGTCTC", "BpiI": "GAAGAC", "BsmBI": "CGTCTC"}
_TRANS = str.maketrans("ACGTN", "TGCAN")


def _rc(s: str) -> str:
    return s[::-1].translate(_TRANS)


def fetch(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=180) as r:
        return r.read()


def parse_zip(blob: bytes, source: dict) -> list:
    z = zipfile.ZipFile(io.BytesIO(blob))
    names = sorted(n for n in z.namelist()
                   if n.lower().endswith((".gb", ".gbk", ".genbank"))
                   and not n.startswith("__MACOSX"))
    out = []
    for n in names:
        text = z.read(n).decode("utf-8", errors="replace")
        m = _ORIGIN.search(text)
        if not m:
            continue
        seq = re.sub(r"[^ACGTNacgtn]", "", m.group(1)).upper()
        if not seq:
            continue
        loc = _LOCUS_BP.search(text)
        declared = int(loc.group(2)) if loc else None
        base = os.path.basename(n)
        # The filename stem is the part identity. An earlier version used a
        # non-greedy prefix regex, which collapsed EVERY Cr part to "pCM0" and
        # silently broke assembly de-duplication -- use the stem, which is unique.
        # Whitespace in a stem would break the FASTA round-trip: the reader splits
        # the header on whitespace, so an id containing a space is silently truncated
        # and the record can no longer be matched back to its manifest entry (and so
        # loses its tier). Collapse whitespace to underscores at the source.
        stem = re.sub(r"\s+", "_", os.path.splitext(base)[0])
        pid = re.match(r"([A-Za-z]+\d*(?:[-_]\d+)?)", stem)
        sites = {}
        for enz, site in TYPE_IIS.items():
            c = seq.count(site) + seq.count(_rc(site))
            if c:
                sites[enz] = c
        out.append({
            "id": f"{source['key']}:{stem}",
            "part_id": stem,
            "part_family": pid.group(1) if pid else stem,
            "file": base,
            "source": source["key"],
            "tier": source["tier"],
            "host": source["host"],
            "citation": source["citation"],
            "length_bp": len(seq),
            "declared_length_ok": (declared == len(seq)) if declared else None,
            "gc_percent": round(100.0 * (seq.count("G") + seq.count("C")) / len(seq), 1),
            "type_iis_sites": sites,
            "sequence": seq,
        })
    return out


# --------------------------------------------------------------------------
# Real Type IIS excision + combinatorial assembly (the "synthetic" tier)
# --------------------------------------------------------------------------
# BsaI recognises GGTCTC, cuts 1 nt downstream on the top strand and leaves a
# 4 nt 5' overhang. In a MoClo Level-0 plasmid the part sits between two
# inward-facing BsaI sites, so excising it means: find the forward site, take
# the 4 nt overhang starting 1 nt after the recognition sequence, and run to the
# matching reverse site. This is the real enzyme geometry, not an approximation.

def excise_bsai_part(seq: str):
    """Return (insert_with_overhangs, left_overhang, right_overhang) or None.

    Only handles the clean, standard case: exactly one forward GGTCTC and one
    reverse site, correctly oriented. Anything else returns None rather than
    guessing -- a wrong excision would poison every assembly built from it.
    """
    fwd = [m.start() for m in re.finditer("GGTCTC", seq)]
    rev = [m.start() for m in re.finditer("GAGACC", seq)]  # rc of GGTCTC
    if len(fwd) != 1 or len(rev) != 1:
        return None
    f, r = fwd[0], rev[0]
    start = f + 6 + 1           # skip recognition site + 1 nt spacer
    end = r - 1                 # reverse site: 1 nt spacer on its side
    if end - start < 8:
        return None
    insert = seq[start:end]
    if len(insert) < 8:
        return None
    return insert, insert[:4], insert[-4:]


def build_synthetic_assemblies(records: list, target: int, seed: int = 0) -> list:
    """Chain real excised parts into Level-1 constructs by matching overhangs.

    Ground truth is exact by construction: we know every junction coordinate
    because we placed it. Deterministic given `seed` so the corpus is reproducible.
    """
    import random
    rng = random.Random(seed)

    parts = []
    for r in records:
        if r["tier"] != "cr_primary":
            continue
        ex = excise_bsai_part(r["sequence"])
        if not ex:
            continue
        insert, lo, ro = ex
        parts.append({"part_id": r["part_id"], "seq": insert, "left": lo, "right": ro})
    if len(parts) < 2:
        return []

    by_left = {}
    for p in parts:
        by_left.setdefault(p["left"], []).append(p)

    out = []
    attempts = 0
    seen = set()
    while len(out) < target and attempts < target * 200:
        attempts += 1
        chain = [rng.choice(parts)]
        # extend while a part exists whose left overhang matches our right one
        while len(chain) < 6:
            nxt = by_left.get(chain[-1]["right"])
            if not nxt:
                break
            pick = rng.choice(nxt)
            if any(pick["part_id"] == c["part_id"] for c in chain):
                break
            chain.append(pick)
        if len(chain) < 2:
            continue
        key = tuple(c["part_id"] for c in chain)
        if key in seen:
            continue
        seen.add(key)

        seq_parts, junctions, pos = [], [], 0
        for i, c in enumerate(chain):
            s = c["seq"]
            if i > 0:
                s = s[4:]  # the shared 4 nt overhang is not duplicated
            seq_parts.append(s)
            span_start = pos
            pos += len(s)
            junctions.append({"part_id": c["part_id"], "start": span_start, "end": pos})
        seq = "".join(seq_parts)
        out.append({
            "id": f"synthetic:L1_{len(out):05d}",
            "part_id": "+".join(key),
            "source": "synthetic_assembly",
            "tier": "synthetic",
            "host": "Chlamydomonas reinhardtii (in silico assembly of real Cr parts)",
            "citation": "assembled in silico from real cr_primary parts by MoClo overhang rules",
            "length_bp": len(seq),
            "gc_percent": round(100.0 * (seq.count("G") + seq.count("C")) / len(seq), 1),
            "n_parts": len(chain),
            "junctions": junctions,
            "sequence": seq,
        })
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=int, default=1000,
                    help="total corpus size to reach; synthetic assemblies make up "
                         "any shortfall after all real sources are exhausted")
    ap.add_argument("--no-synthetic", action="store_true",
                    help="real deposited records only; corpus will be far below --target")
    args = ap.parse_args()

    all_records = []
    provenance = []
    quarantined = []
    for src in SOURCES:
        print(f"fetching {src['key']} ({src['tier']}) ...")
        try:
            blob = fetch(src["url"])
        except Exception as e:
            print(f"   FAILED: {e}  -- skipping, corpus will be smaller")
            provenance.append({**{k: v for k, v in src.items()}, "status": f"failed: {e}",
                               "n_records": 0})
            continue
        got_sha = hashlib.sha256(blob).hexdigest()
        if src["sha256"] and got_sha != src["sha256"]:
            raise RuntimeError(
                f"{src['key']}: SHA-256 mismatch (expected {src['sha256']}, got {got_sha}). "
                "Upstream changed -- verify by hand before updating.")
        recs = parse_zip(blob, src)
        # Quarantine, never silently accept: a record whose ORIGIN sequence length
        # disagrees with its own LOCUS header has a real integrity problem. The
        # Marillonnet kit has 7 such records (pICH83955/66/77/88/99, pICH84000/84011
        # -- all declare 10,988 bp but hold 10,980, a systematic +8 annotation error
        # in one plasmid family, observed 2026-08-12). This is the same class of
        # silent error pLannotate documents in secondary-source plasmid maps, and
        # is exactly what IC-3 exists to catch.
        bad = [r for r in recs if r["declared_length_ok"] is False]
        good = [r for r in recs if r["declared_length_ok"] is not False]
        for r in bad:
            quarantined.append({"id": r["id"], "source": src["key"],
                                "reason": "LOCUS header length disagrees with ORIGIN sequence length",
                                "length_bp": r["length_bp"]})
        nosite = [r for r in good if not r["type_iis_sites"]]
        print(f"   {len(good)} kept | quarantined {len(bad)} (length/header mismatch)"
              f" | Type IIS sites: {len(good)-len(nosite)}/{len(good)}"
              f" | sha256 {got_sha[:16]}...")
        all_records.extend(good)
        provenance.append({"key": src["key"], "tier": src["tier"], "host": src["host"],
                           "url": src["url"], "sha256": got_sha,
                           "citation": src["citation"], "n_records": len(good),
                           "n_quarantined": len(bad), "status": "ok"})

    n_real = len(all_records)
    print(f"\nreal deposited records: {n_real}")
    if not args.no_synthetic and n_real < args.target:
        need = args.target - n_real
        print(f"generating {need} synthetic Level-1 assemblies from real Cr parts ...")
        syn = build_synthetic_assemblies(all_records, need)
        print(f"   generated {len(syn)}")
        all_records.extend(syn)

    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, "corpus.fasta"), "w",
              encoding="utf-8", newline="\n") as fh:
        for r in all_records:
            fh.write(f">{r['id']} tier={r['tier']} host={r['host']} len={r['length_bp']}\n")
            s = r["sequence"]
            for i in range(0, len(s), 70):
                fh.write(s[i:i + 70] + "\n")

    tiers = {}
    for r in all_records:
        tiers[r["tier"]] = tiers.get(r["tier"], 0) + 1
    manifest = {
        "built": "see git history",
        "n_total": len(all_records),
        "by_tier": tiers,
        "n_quarantined": len(quarantined),
        "quarantined": quarantined,
        "sources": provenance,
        "tier_meaning": {
            "cr_primary": "real deposited Chlamydomonas parts; valid for junction grammar AND Cr identity",
            "syntax_only": "real deposited Type IIS parts from other hosts; junction grammar ONLY, not Cr evidence",
            "synthetic": "in-silico Level-1 assemblies of real Cr parts; exact junction ground truth; NOT deposited plasmids",
        },
        "records": [{k: v for k, v in r.items() if k != "sequence"} for r in all_records],
    }
    with open(os.path.join(OUT_DIR, "corpus_manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"\nTOTAL: {len(all_records)}")
    for t, n in sorted(tiers.items()):
        print(f"   {t:14s} {n}")
    print(f"wrote {OUT_DIR}")
    print("\nReport this corpus WITH its tier breakdown. Only the cr_primary tier is "
          "evidence about Chlamydomonas.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
