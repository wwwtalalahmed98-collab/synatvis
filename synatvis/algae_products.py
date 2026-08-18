"""Prior-art lookup: has this protein already been made in algae?

Every scan is otherwise treated as if nothing had ever been expressed in an algal
system before. This module carries the record forward: given a gene name and/or a
coding sequence, it reports whether the product is a known algal product, or
resembles one, and brings the prior art with it.

Two independent matching routes, deliberately kept separate so the evidence for a
hit is always legible:

  NAME match       -- the transcript's declared gene symbol equals a catalogued gene
                      symbol (case-insensitive, punctuation-normalised). Exact and
                      cheap; this is the trustworthy route.
  SIMILARITY match -- the translated protein shares k-mers with a catalogued
                      reference protein. Only possible for entries that carry a real
                      reference sequence, and reported with its measured score so a
                      weak hit can never masquerade as an identification.

Similarity here is a k-mer containment score, not an alignment. It is a screening
aid: it answers "worth a look?", never "this is that protein". Anything below
STRONG_SIMILARITY is reported as a resemblance only.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .profiles import PACKAGE_DIR, load_yaml

CATALOGUE_PATH = os.path.join(PACKAGE_DIR, "data", "algae_products.yaml")
REFSEQ_PATH = os.path.join(PACKAGE_DIR, "data", "algae_product_refs.fasta")

K = 5                      # peptide k-mer length for the screening score
STRONG_SIMILARITY = 0.60   # at/above this, call it a strong resemblance
MIN_SIMILARITY = 0.25      # below this, do not report at all


@dataclass
class ProductHit:
    product: str
    gene: str
    host: str
    compartment: str
    origin: str
    product_class: str
    application: str
    confidence: str
    match_type: str                 # "name" | "similarity"
    similarity: Optional[float] = None

    def summary(self) -> str:
        if self.match_type == "name":
            return (f"{self.product} — already produced in {self.host} "
                    f"({self.compartment}, {self.origin}).")
        return (f"resembles {self.product} ({self.host}) — "
                f"{self.similarity:.0%} peptide k-mer containment, screening only.")


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


def load_catalogue(path: str = CATALOGUE_PATH) -> List[Dict]:
    with open(path, "r", encoding="utf-8") as fh:
        data = load_yaml(fh.read())
    return (data or {}).get("products", [])


def load_reference_proteins(path: str = REFSEQ_PATH) -> Dict[str, str]:
    """Reference PROTEIN sequences keyed by gene symbol, if the file exists."""
    if not os.path.isfile(path):
        return {}
    out, cur = {}, None
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                cur = line[1:].split()[0]
                out[cur] = []
            elif cur:
                out[cur].append(line.strip())
    return {k: "".join(v).upper() for k, v in out.items() if v}


def _kmers(seq: str, k: int = K) -> set:
    return {seq[i:i + k] for i in range(len(seq) - k + 1)} if len(seq) >= k else set()


def similarity(query_prot: str, ref_prot: str, k: int = K) -> float:
    """Containment of the shorter sequence's k-mers in the longer one (0..1).

    Containment rather than Jaccard, so a short peptide genuinely contained in a
    large protein is not penalised for the length difference.
    """
    a, b = _kmers(query_prot, k), _kmers(ref_prot, k)
    if not a or not b:
        return 0.0
    small, big = (a, b) if len(a) <= len(b) else (b, a)
    return len(small & big) / len(small)


def identify(name: str = "", cds: str = "",
             catalogue: Optional[List[Dict]] = None,
             refs: Optional[Dict[str, str]] = None) -> List[ProductHit]:
    """Return catalogued algal products matching this gene by name or resemblance."""
    catalogue = catalogue if catalogue is not None else load_catalogue()
    refs = refs if refs is not None else load_reference_proteins()
    hits: List[ProductHit] = []

    def mk(entry: Dict, match_type: str, sim=None) -> ProductHit:
        return ProductHit(
            product=entry.get("product", "?"), gene=str(entry.get("gene", "")),
            host=entry.get("host", "?"), compartment=entry.get("compartment", "?"),
            origin=entry.get("origin", "?"), product_class=entry.get("product_class", "?"),
            application=entry.get("application", ""), confidence=entry.get("confidence", "reported"),
            match_type=match_type, similarity=sim)

    # --- name route ---
    qn = _norm(name)
    named = set()
    if qn:
        for ent in catalogue:
            for sym in str(ent.get("gene", "")).replace("/", " ").split():
                if _norm(sym) and _norm(sym) == qn:
                    hits.append(mk(ent, "name"))
                    named.add(id(ent))
                    break

    # --- similarity route ---
    if cds and refs:
        from .ptm import translate
        prot = (translate(cds) or "").rstrip("*")
        if len(prot) >= K:
            for ent in catalogue:
                if id(ent) in named:
                    continue
                for sym in str(ent.get("gene", "")).replace("/", " ").split():
                    ref = refs.get(sym)
                    if not ref:
                        continue
                    s = similarity(prot, ref)
                    if s >= MIN_SIMILARITY:
                        hits.append(mk(ent, "similarity", round(s, 3)))
                    break

    hits.sort(key=lambda h: (h.match_type != "name", -(h.similarity or 0)))
    return hits


def catalogue_stats(catalogue: Optional[List[Dict]] = None) -> Dict:
    catalogue = catalogue if catalogue is not None else load_catalogue()
    by_origin, by_class, by_host = {}, {}, {}
    for e in catalogue:
        by_origin[e.get("origin", "?")] = by_origin.get(e.get("origin", "?"), 0) + 1
        by_class[e.get("product_class", "?")] = by_class.get(e.get("product_class", "?"), 0) + 1
        h = e.get("host", "?").split("(")[0].strip()
        by_host[h] = by_host.get(h, 0) + 1
    return {"n": len(catalogue), "by_origin": by_origin,
            "by_class": by_class, "by_host": by_host}
