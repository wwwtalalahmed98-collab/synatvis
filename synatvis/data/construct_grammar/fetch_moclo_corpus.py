"""Stage 0, Step 2 -- fetch the real MoClo toolkit sequences from the primary source.

The Chlamydomonas Resource Center publishes the complete Cr MoClo toolkit as a
single downloadable archive of GenBank files. This was confirmed directly by the
toolkit's own authors, who replied to an emailed request with: "You can download
all sequence files directly from the Chlamydomonas Resource Center (CRC) catalog."

An earlier pass through this project concluded no sequence download existed. That
was WRONG, and worth recording as a lesson: that pass only checked the individual
per-part product pages (/product/pcm0-001/ etc.), which genuinely have no download
link. The bulk archive is linked from the toolkit LANDING page instead.

This script downloads nothing into version control. It fetches the archive from
the primary source, verifies its SHA-256 against the value recorded when the
corpus was first retrieved, parses every GenBank record, and writes the corpus to
a local, git-ignored directory. Anyone can reproduce the exact corpus from the
primary source; the third-party sequence files are never redistributed here.

Run:
    python -m synatvis.data.construct_grammar.fetch_moclo_corpus
or:
    python synatvis/data/construct_grammar/fetch_moclo_corpus.py
"""
from __future__ import annotations

import hashlib
import io
import json
import os
import re
import sys
import urllib.request
import zipfile

ARCHIVE_URL = "https://www.chlamycollection.org/content/uploads/2019/06/MoClo-Kit-Sequences.zip"
# SHA-256 of the archive as retrieved 2026-08-12. A mismatch means the upstream
# file changed -- do NOT silently accept it: re-verify the contents, then update
# this constant and the ledger together, so the corpus stays reproducible.
EXPECTED_SHA256 = "8689ef08e7611fe3ce87ec11f1921ef8a4d37da578c474048897ff1d27128ec4"
EXPECTED_RECORDS = 115

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "moclo_corpus")  # git-ignored; see .gitignore

_ORIGIN = re.compile(r"^ORIGIN\s*?$(.*?)^//", re.M | re.S)
_LOCUS_BP = re.compile(r"^LOCUS\s+(\S+)\s+(\d+)\s+bp", re.M)
TYPE_IIS = {"BsaI": "GGTCTC", "BpiI": "GAAGAC", "BsmBI": "CGTCTC"}


def _revcomp(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def download(url: str = ARCHIVE_URL) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as resp:
        return resp.read()


def verify(blob: bytes, expected: str = EXPECTED_SHA256) -> str:
    got = hashlib.sha256(blob).hexdigest()
    if got != expected:
        raise RuntimeError(
            f"SHA-256 mismatch.\n  expected {expected}\n  got      {got}\n"
            "The upstream archive changed. Re-verify its contents by hand before "
            "updating EXPECTED_SHA256 and step2_catalog_ledger.yaml together."
        )
    return got


def parse_archive(blob: bytes) -> list:
    """Parse every .gb record in the archive into a plain dict."""
    z = zipfile.ZipFile(io.BytesIO(blob))
    names = [n for n in z.namelist()
             if n.lower().endswith(".gb") and not n.startswith("__MACOSX")]
    out = []
    for n in sorted(names):
        text = z.read(n).decode("utf-8", errors="replace")
        m = _ORIGIN.search(text)
        if not m:
            raise RuntimeError(f"no ORIGIN block in {n} -- archive format changed")
        seq = re.sub(r"[^ACGTNacgtn]", "", m.group(1)).upper()
        if not seq:
            raise RuntimeError(f"empty sequence in {n}")
        base = os.path.basename(n)
        pid = re.match(r"(pCM0-\d+)", base)
        loc = _LOCUS_BP.search(text)
        declared = int(loc.group(2)) if loc else None
        if declared is not None and declared != len(seq):
            raise RuntimeError(
                f"{base}: LOCUS declares {declared} bp but ORIGIN holds {len(seq)}"
            )
        sites = {}
        for enz, site in TYPE_IIS.items():
            c = seq.count(site) + seq.count(_revcomp(site))
            if c:
                sites[enz] = c
        out.append({
            "part_id": pid.group(1) if pid else None,
            "file": base,
            "length_bp": len(seq),
            "gc_percent": round(100.0 * (seq.count("G") + seq.count("C")) / len(seq), 1),
            "type_iis_sites": sites,
            "sequence": seq,
        })
    return out


def write_corpus(records: list, out_dir: str = OUT_DIR) -> None:
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "moclo_parts.fasta"), "w",
              encoding="utf-8", newline="\n") as fh:
        for r in records:
            fh.write(f">{r['part_id']} {r['file']} len={r['length_bp']}\n")
            s = r["sequence"]
            for i in range(0, len(s), 70):
                fh.write(s[i:i + 70] + "\n")
    manifest = [{k: v for k, v in r.items() if k != "sequence"} for r in records]
    with open(os.path.join(out_dir, "manifest.json"), "w", encoding="utf-8") as fh:
        json.dump({"source_url": ARCHIVE_URL, "sha256": EXPECTED_SHA256,
                   "n_records": len(records), "records": manifest}, fh, indent=2)


def main() -> int:
    print(f"fetching {ARCHIVE_URL}")
    blob = download()
    print(f"  {len(blob)} bytes")
    verify(blob)
    print(f"  SHA-256 verified: {EXPECTED_SHA256}")
    records = parse_archive(blob)
    print(f"  parsed {len(records)} GenBank records")
    if len(records) != EXPECTED_RECORDS:
        print(f"  WARNING: expected {EXPECTED_RECORDS} records, got {len(records)}")
    no_sites = [r["part_id"] for r in records if not r["type_iis_sites"]]
    print(f"  records carrying >=1 Type IIS site: "
          f"{len(records) - len(no_sites)}/{len(records)}")
    if no_sites:
        print(f"  without Type IIS sites (would fail IC-1): {no_sites}")
    write_corpus(records)
    print(f"  wrote corpus to {OUT_DIR}")
    print("\nNote: these are whole PLASMIDS, not excised parts. Cutting each part "
          "out from between its Type IIS sites is Stage 1 (segmentation), not Step 2.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
