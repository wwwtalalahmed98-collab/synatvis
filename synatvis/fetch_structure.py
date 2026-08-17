"""Fetch a real 3D structure for a gene, so molecular dynamics can be about YOUR gene.

Molecular dynamics needs a 3D starting structure. This helper finds one, automatically,
for genes that already have one in a public database:

    FASTA header gene name  ->  UniProt accession  ->  AlphaFold predicted structure

It never invents a structure. If no public structure exists for the gene -- which is
the normal case for a NOVEL designed protein -- it says so plainly and points at the
ColabFold route, which has to be run by a human in a browser.

Usage:
    python -m synatvis.fetch_structure RBCS2
    python -m synatvis.fetch_structure --fasta my_gene.fasta
    python -m synatvis.fetch_structure RBCS2 --out ./structures
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import urllib.parse
import urllib.request

UA = {"User-Agent": "SynAT.Vis structure fetcher"}
ORGANISM_ID = 3055  # Chlamydomonas reinhardtii


def _get(url: str, timeout: int = 60) -> bytes:
    return urllib.request.urlopen(
        urllib.request.Request(url, headers=UA), timeout=timeout).read()


def find_uniprot(gene: str, organism_id: int = ORGANISM_ID):
    """Look up a gene symbol in UniProt, restricted to the host organism."""
    q = f"gene:{gene} AND organism_id:{organism_id}"
    url = ("https://rest.uniprot.org/uniprotkb/search?query="
           + urllib.parse.quote(q)
           + "&fields=accession,id,protein_name,gene_names,length&format=json&size=5")
    try:
        data = json.loads(_get(url))
    except Exception as exc:
        print(f"  UniProt lookup failed: {exc}")
        return []
    out = []
    for r in data.get("results", []):
        name = (r.get("proteinDescription", {}).get("recommendedName", {})
                 .get("fullName", {}).get("value", "?"))
        out.append({"accession": r.get("primaryAccession"),
                    "name": name,
                    "length": r.get("sequence", {}).get("length")})
    return out


def alphafold_entry(accession: str):
    """Return the AlphaFold DB record for a UniProt accession, or None."""
    try:
        data = json.loads(_get(f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"))
    except Exception:
        return None
    return data[0] if data else None


def download_structure(entry: dict, out_dir: str) -> str:
    os.makedirs(out_dir, exist_ok=True)
    url = entry.get("pdbUrl")
    if not url:
        raise RuntimeError("AlphaFold entry has no pdbUrl")
    path = os.path.join(out_dir, os.path.basename(url))
    with open(path, "wb") as fh:
        fh.write(_get(url, timeout=120))
    return path


def gene_from_fasta(path: str) -> str:
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                return line[1:].split()[0]
    return ""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gene", nargs="?", help="gene symbol, e.g. RBCS2")
    ap.add_argument("--fasta", help="read the gene name from this FASTA header instead")
    ap.add_argument("--organism", type=int, default=ORGANISM_ID,
                    help=f"NCBI taxonomy id (default {ORGANISM_ID} = C. reinhardtii)")
    ap.add_argument("--out", default="structures", help="download directory")
    args = ap.parse_args(argv)

    gene = args.gene or (gene_from_fasta(args.fasta) if args.fasta else "")
    if not gene:
        ap.error("give a gene symbol, or --fasta whose header starts with one")

    print(f"gene: {gene}   (organism {args.organism})")
    hits = find_uniprot(gene, args.organism)
    if not hits:
        print("\nNo UniProt entry found for that gene symbol in this organism.")
        print("That is EXPECTED for a novel designed protein -- public databases only")
        print("hold structures for proteins that already exist.")
        print("\nTo get a structure for a novel sequence you must run a prediction")
        print("yourself. The free route is ColabFold, which runs in a browser and")
        print("needs a Google account, so it cannot be automated from here:")
        print("  https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb")
        return 1

    print(f"\nUniProt matches ({len(hits)}):")
    for h in hits:
        print(f"  {h['accession']}  {h['name']}  ({h['length']} aa)")

    for h in hits:
        acc = h["accession"]
        entry = alphafold_entry(acc)
        if not entry:
            print(f"\n{acc}: no AlphaFold structure available")
            continue
        path = download_structure(entry, args.out)
        size = os.path.getsize(path)
        print(f"\nDownloaded real AlphaFold structure for {acc}")
        print(f"  model version : {entry.get('latestVersion', '?')}")
        print(f"  created       : {entry.get('modelCreatedDate', '?')}")
        print(f"  residues      : {len(entry.get('uniprotSequence', ''))}")
        print(f"  file          : {path}  ({size:,} bytes)")
        print("\nPoint the MD wrapper's PDB_FILE at that path, then run:")
        print("  python -m synatvis scan <your.fasta> --plugins")
        return 0

    print("\nNo AlphaFold structure for any match. Use the ColabFold route above.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
