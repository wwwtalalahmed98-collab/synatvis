"""Batch-fetch real AlphaFold structures for the genes this tool actually scans.

The shipped native corpus identifies every gene by its RefSeq protein accession
(XP_...), so structures can be tied to the real scanned genes rather than to a
generic demonstration protein:

    corpus FASTA header  ->  XP_ accession  ->  UniProt  ->  AlphaFold structure

Nothing is predicted here and nothing is invented. Genes with no UniProt entry,
or no AlphaFold model, are recorded as misses -- they are not substituted.

Run:
    python -m synatvis.fetch_structures_batch --limit 100
    python -m synatvis.fetch_structures_batch --limit 100 --out structures
"""
from __future__ import annotations

import argparse
import json
import os
import re
import time
import urllib.parse
import urllib.request

from .profiles import PACKAGE_DIR

# Both UniProt and the AlphaFold CDN reject a non-browser agent with 403.
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
CORPORA = ["native_cr_cds.fasta", "native_cr_cds_heldout.fasta", "native_cr_cds_4000.fasta"]
ORGANISM = 3055


def _get(url: str, timeout: int = 60) -> bytes:
    return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout).read()


def corpus_accessions(limit: int = 0) -> list:
    """RefSeq protein accessions from the shipped native corpora, in corpus order."""
    seen, out = set(), []
    for fname in CORPORA:
        path = os.path.join(PACKAGE_DIR, "data", fname)
        if not os.path.isfile(path):
            continue
        with open(path, "r", encoding="utf-8") as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                m = re.search(r"(XP_\d+\.\d+)", line)
                if m and m.group(1) not in seen:
                    seen.add(m.group(1))
                    out.append(m.group(1))
                    if limit and len(out) >= limit:
                        return out
    return out


def refseq_to_uniprot(xp: str):
    """Resolve one RefSeq protein accession to a UniProt entry, verifying the xref."""
    base = xp.split(".")[0]
    q = f"xref:refseq-{base} AND organism_id:{ORGANISM}"
    url = ("https://rest.uniprot.org/uniprotkb/search?query=" + urllib.parse.quote(q)
           + "&fields=accession,protein_name,length,xref_refseq&format=json&size=3")
    try:
        data = json.loads(_get(url))
    except Exception:
        return None
    for r in data.get("results", []):
        xrefs = [x.get("id", "") for x in r.get("uniProtKBCrossReferences", [])
                 if x.get("database") == "RefSeq"]
        # only accept when the RefSeq cross-reference really is this accession
        if not any(x.split(".")[0] == base for x in xrefs):
            continue
        name = ((r.get("proteinDescription", {}).get("recommendedName", {}) or {})
                .get("fullName", {}).get("value", ""))
        return {"accession": r["primaryAccession"], "name": name,
                "length": r.get("sequence", {}).get("length")}
    return None


def alphafold_pdb(accession: str):
    try:
        data = json.loads(_get(f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"))
    except Exception:
        return None
    return data[0] if data else None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit", type=int, default=100, help="how many structures to collect")
    ap.add_argument("--scan", type=int, default=0,
                    help="how many corpus genes to try (default: 3x limit)")
    ap.add_argument("--out", default="structures", help="download directory")
    ap.add_argument("--pause", type=float, default=0.25, help="seconds between requests")
    args = ap.parse_args(argv)

    scan_n = args.scan or args.limit * 3
    accs = corpus_accessions(scan_n)
    print(f"corpus accessions to try: {len(accs)}   target structures: {args.limit}")
    os.makedirs(args.out, exist_ok=True)

    got, manifest = 0, []
    no_uniprot = no_model = 0
    for i, xp in enumerate(accs, 1):
        if got >= args.limit:
            break
        up = refseq_to_uniprot(xp)
        time.sleep(args.pause)
        if not up:
            no_uniprot += 1
            continue
        entry = alphafold_pdb(up["accession"])
        time.sleep(args.pause)
        if not entry or not entry.get("pdbUrl"):
            no_model += 1
            continue
        fname = os.path.basename(entry["pdbUrl"])
        path = os.path.join(args.out, fname)
        try:
            blob = _get(entry["pdbUrl"], timeout=120)
        except Exception:
            no_model += 1
            continue
        with open(path, "wb") as fh:
            fh.write(blob)
        got += 1
        manifest.append({"refseq": xp, "uniprot": up["accession"], "protein": up["name"],
                         "length_aa": up["length"], "file": fname, "bytes": len(blob),
                         "model_version": entry.get("latestVersion")})
        if got % 10 == 0 or got == 1:
            print(f"  [{got:3d}/{args.limit}] {xp} -> {up['accession']}  {up['name'][:44]}")
        time.sleep(args.pause)

    with open(os.path.join(args.out, "structures_manifest.json"), "w", encoding="utf-8") as fh:
        json.dump({"source": "AlphaFold DB via UniProt from SynAT.Vis native corpus",
                   "organism_id": ORGANISM, "n": len(manifest),
                   "tried": i if accs else 0,
                   "no_uniprot_entry": no_uniprot, "no_alphafold_model": no_model,
                   "structures": manifest}, fh, indent=1)

    print(f"\ndownloaded {got} real AlphaFold structures into {args.out}")
    print(f"  tried {i} corpus genes | {no_uniprot} had no UniProt entry | "
          f"{no_model} had no AlphaFold model")
    print("  misses are recorded, not substituted -- see structures_manifest.json")
    return 0 if got else 1


if __name__ == "__main__":
    raise SystemExit(main())
