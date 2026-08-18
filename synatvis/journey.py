"""Assemble the full gene->protein->purified-protein journey as ordered checkpoints.

This is the single place that gathers every parameter the construct passes through
and attaches, to each checkpoint, the values + scales + a pass/warn/fail mark and a
citation. The cell-journey HTML is then a generic renderer over this structure, so
all the science lives here (and is unit-testable), not in JavaScript.

Nothing here is a fitted model: scores come from the (already validated) scan and
the transparent expression ensemble; biophysics/purification/structure come from
published formulas and rules, each carrying its reference.
"""
from __future__ import annotations

from typing import Dict, List, Optional

from .flags import Severity
from .util import gc_fraction


def _p(label, value, scale=None, status="info", ref="", detail=""):
    return {"label": label, "value": value, "scale": scale, "status": status,
            "ref": ref, "detail": detail}


def _module_status(scan_result, module, bad_sev=Severity.HIGH):
    """(status, message) for a diagnostic module from its flags."""
    worst = None
    for f in scan_result.flags:
        if f.module == module:
            if worst is None or f.severity > worst.severity:
                worst = f
    if worst is None:
        return "ok", ""
    if worst.severity >= bad_sev:
        return "bad", worst.message
    if worst.severity >= Severity.MEDIUM:
        return "warn", worst.message
    return "ok", worst.message


def build_journey(transcript, profile: Dict, scan_result, expression,
                  plddt: Optional[List[float]] = None) -> Dict:
    from .ptm import biophysics, translate
    from .purification import predict_purification
    from .structure_confidence import interpret_structure

    prot = translate(transcript.cds)
    bio = biophysics(prot)
    pur = predict_purification(transcript)
    struct = interpret_structure(transcript, plddt=plddt)
    fate = expression.fate or {}
    loc = fate.get("localization", "cytosolic")

    mscore = {m["name"]: m for m in expression.models}

    def sc(name):
        m = mscore.get(name)
        return round(m["score"], 0) if (m and m.get("available")) else None

    gcpct = round(gc_fraction(transcript.cds) * 100, 1)
    polya_status, polya_msg = _module_status(scan_result, "polya")
    # a 3'UTR poly(A) signal is wanted; only CDS/5'UTR premature ones are hazards
    prem_polya = any(f.module == "polya" and f.region in ("cds", "5utr")
                     and f.severity >= Severity.MEDIUM for f in scan_result.flags)
    clone_status, clone_msg = _module_status(scan_result, "cloning", bad_sev=Severity.MEDIUM)
    uorf_status, uorf_msg = _module_status(scan_result, "uorf", bad_sev=Severity.MEDIUM)
    comp_status, comp_msg = _module_status(scan_result, "composition")
    splice_status, splice_msg = _module_status(scan_result, "splice", bad_sev=Severity.MEDIUM)
    inst_status, inst_msg = _module_status(scan_result, "instability", bad_sev=Severity.MEDIUM)
    struct_status, struct_msg = _module_status(scan_result, "structure", bad_sev=Severity.MEDIUM)

    cp: List[Dict] = []

    # 1 — the DNA construct
    cp.append({"key": "gene", "title": "Gene construct (DNA)", "symbol": "dna",
               "organelle": "nucleus", "level": "pretranscript", "params": [
        _p("Length", f"{len(transcript.cds)} nt"),
        _p("GC content", f"{gcpct}%", scale=gcpct,
           status=("ok" if gcpct >= 55 else "bad"),
           ref="Barahimipour 2015 (low GC silences in Cr)",
           detail="Cr nuclear genes are GC-rich (~66%); low GC is the silencing hazard."),
        _p("Type IIS / cloning sites", "clear" if clone_status == "ok" else "present",
           status=clone_status, ref="MoClo / Golden Gate", detail=clone_msg),
    ]})

    # 2 — transcription
    cp.append({"key": "transcription", "title": "Transcription -> pre-mRNA", "symbol": "rna",
               "organelle": "nucleus", "level": "transcript", "params": [
        _p("Silencing / GC-trough risk", "none" if comp_status == "ok" else "flagged",
           status=comp_status, ref="Schroda 2019; Neupert/Bock UVM strains",
           detail=comp_msg or "Chromatin/strain-driven; use a silencing-deficient strain + intron."),
    ]})

    # 3 — splicing (introns are assets in Cr)
    cp.append({"key": "splice", "title": "Splicing", "symbol": "spliceosome",
               "organelle": "nucleus", "level": "posttranscript", "params": [
        _p("Cryptic-splice risk", "low" if splice_status == "ok" else "check",
           status=splice_status, ref="Zeng & Li 2022 (Pangolin)", detail=splice_msg),
        _p("Introns", "asset in Cr (e.g. rbcS2 intron 1 boosts expression)",
           status="info", ref="Baier 2018"),
    ]})

    # 4 — 3' processing / poly(A)  (Cr NUE = TGTAA, not AAUAAA)
    cp.append({"key": "polya", "title": "3' cleavage & poly(A)", "symbol": "polya",
               "organelle": "nucleusEdge", "level": "posttranscript", "params": [
        _p("Premature poly(A) in CDS", "none" if not prem_polya else "PRESENT",
           status=("ok" if not prem_polya else "bad"),
           ref="Cr NUE = TGTAA (not AAUAAA)", detail=polya_msg if prem_polya else
           "No premature TGTAA signal inside the coding sequence."),
        _p("mRNA-destabilising elements (AREs)", "none" if inst_status == "ok" else "flagged",
           status=inst_status, ref="Barreau 2005", detail=inst_msg),
    ]})

    # 5 — nuclear export
    cp.append({"key": "export", "title": "Nuclear export", "symbol": "pore",
               "organelle": "pore", "level": "posttranscript", "params": [
        _p("Mature mRNA", "capped, spliced, poly-adenylated", status="ok"),
    ]})

    # 6 — initiation / 5'UTR (pre-translation)
    has_atg = transcript.cds.upper().startswith("ATG")
    cp.append({"key": "initiation", "title": "Translation initiation (5'UTR)",
               "symbol": "pore", "organelle": "ribo", "level": "pretranslation", "params": [
        _p("Start codon (ATG)", "present" if has_atg else "MISSING",
           status=("ok" if has_atg else "bad"), ref="Kozak",
           detail=None if has_atg else "CDS does not begin with ATG."),
        _p("Upstream AUG / uORF", "none" if uorf_status == "ok" else "present",
           status=uorf_status, ref="Kozak; Calvo 2009", detail=uorf_msg,
           scale=None),
        _p("5' mRNA start-region structure", "clear" if struct_status == "ok" else "structured",
           status=struct_status, ref="Kudla 2009; LinearDesign (Zhang 2020)",
           detail=struct_msg or "Strong 5' folding can throttle ribosome loading."),
    ]})

    # 7 — translation (the adaptiveness landscape lives here)
    cp.append({"key": "translation", "title": "Translation (elongation)", "symbol": "ribosome",
               "organelle": "ribo", "level": "translation", "params": [
        _p("Codon adaptation (RCA)", sc("rca"), scale=sc("rca"),
           status=("ok" if (sc("rca") or 0) >= 50 else "warn"),
           ref="Sharp & Li 1987; Presnyak 2015",
           detail="Percentile vs 5,000 well-expressed native Cr genes."),
        _p("tRNA adaptation (tAI)", sc("tai"), scale=sc("tai"),
           status=("ok" if (sc("tai") or 0) >= 50 else "warn"), ref="dos Reis 2004"),
    ], "landscape": True})

    # 8 — folding (self-explanatory structure confidence)
    cp.append({"key": "folding", "title": "Folding", "symbol": "fold",
               "organelle": "fold", "level": "posttranslation", "params": [
        _p("Structure confidence", f"{struct.mean_confidence:.0f}/100 ({struct.source})",
           scale=struct.mean_confidence,
           status=("ok" if struct.mean_confidence >= 70 else
                   "warn" if struct.mean_confidence >= 55 else "bad"),
           ref="TOP-IDP (Campen 2008) / pLDDT", detail=struct.plain_summary),
        _p("Size", f"{bio['mw_kda']} kDa ({bio['length']} aa)", status="info"),
        _p("Isoelectric point (pI)", bio["pI"], status="info",
           ref="ProtParam", detail="Set purification buffer pH away from pI to keep it soluble."),
        _p("Hydropathy (GRAVY)", bio["gravy"], status="info",
           ref="Kyte-Doolittle", detail="Positive = hydrophobic (aggregation/solubility watch)."),
    ]})

    # 8a — prior art: has this protein already been made in algae?
    from .algae_products import identify as _identify_products
    _prods = _identify_products(name=transcript.name, cds=transcript.cds)
    if _prods:
        cp.append({"key": "algae_prior_art", "title": "Known algal product",
                   "symbol": "vial", "organelle": "bench", "level": "posttranslation",
                   "params": [
            _p(h.product,
               (f"{h.host} · {h.compartment} · {h.origin}" if h.match_type == "name"
                else f"resembles — {h.similarity:.0%} peptide k-mer containment"),
               status="info",
               ref=f"catalogue confidence: {h.confidence}",
               detail=(h.application + (
                   "" if h.match_type == "name" else
                   " — SIMILARITY is a screening aid, not an identification; confirm by alignment.")))
            for h in _prods[:4]
        ]})

    # 8b — multiprotein complex membership (name-based, cryo-ET-anchored context)
    from .complexome import identify_complexes
    complex_matches = identify_complexes(transcript.name)
    if complex_matches:
        cp.append({"key": "complex", "title": "Multiprotein complex membership",
                   "symbol": "complex", "organelle": "fold", "level": "posttranslation", "params": [
            _p(m.complex_name, m.function, status="info",
               ref=m.structural_citation,
               detail=f"Gene-name match ({m.matched_pattern}); identity basis: "
                      f"{m.identity_citation}. This complex was directly imaged intact "
                      f"in Cr cells by cryo-electron tomography — this specific molecule "
                      f"was not itself imaged, only its complex family.")
            for m in complex_matches
        ]})

    # 8/9 — PTM routing by fate
    if "secreted" in loc or "membrane" in loc:
        cp.append({"key": "er", "title": "ER entry (signal peptide)", "symbol": "er",
                   "organelle": "er", "level": "posttranslation", "params": [
            _p("Signal peptide", "yes" if fate.get("signal_peptide") else "no",
               status="ok" if fate.get("signal_peptide") else "info",
               ref="SignalP 6.0 concept"),
        ]})
        cp.append({"key": "golgi", "title": "Golgi — N-glycosylation", "symbol": "glycan",
                   "organelle": "golgi", "level": "posttranslation", "params": [
            _p("N-glycosylation sequons", fate.get("n_glyc_sites", 0),
               status="info", ref="Gavel & von Heijne 1990",
               detail="Cr N-glycosylation differs from mammalian."),
            _p("Disulfide potential", "yes" if fate.get("disulfide_potential") else "no",
               status="info", ref="oxidising ER lumen"),
        ]})
        dest = "secreted" if "secreted" in loc else "membrane"
        cp.append({"key": "localise", "title": f"Localisation — {loc}", "symbol": "vesicle",
                   "organelle": dest, "level": "posttranslation", "params": [
            _p("Predicted localisation", loc, status="ok", ref="DeepLoc 2.1 concept")]})
    else:
        cp.append({"key": "localise", "title": "Released into cytosol", "symbol": "protein",
                   "organelle": "cyto", "level": "posttranslation", "params": [
            _p("Predicted localisation", loc, status="ok", ref="DeepLoc 2.1 concept")]})

    # ---- purification strategy (recommended, chromatographic + non-chromatographic) ----
    _plabel = {"source": "Source", "capture": "Capture", "intermediate": "Intermediate",
               "polish": "Polish"}
    strat_params = []
    for w in pur.to_dict()["workflow"]:
        lab = _plabel.get(w["phase"].split()[0], w["phase"].title())
        tag = "info" if w["mode"] == "prep" else "ok"
        strat_params.append(_p(lab, w["method"], status=tag, detail=w.get("basis", ""),
                               ref=w.get("mode", "")))
    for cc in pur.considerations:
        if cc["topic"] == "Source":
            continue  # already shown as the workflow's source step
        strat_params.append(_p(cc["topic"], cc.get("tag", ""), status="info", detail=cc["note"]))
    b = pur.biophysics
    if b:
        strat_params.append(_p("Protein for purification",
                               f"{b.get('mw_kda')} kDa · pI {b.get('pI')} · GRAVY {b.get('gravy')}",
                               status="info", ref="ProtParam",
                               detail="pI sets the ion-exchange type; GRAVY sets HIC; MW sets SEC."))
    cp.append({"key": "purification_plan", "title": "Purification strategy (recommended)",
               "symbol": "bead", "organelle": "bench", "level": "posttranslation",
               "params": strat_params})

    # ---- downstream: recovery / purification (post-translation) ----
    if pur.tags:
        cp.append({"key": "lysis", "title": "Cell lysis & clarification", "symbol": "lysis",
                   "organelle": "bench", "level": "posttranslation", "params": [
            _p("Detected tags", ", ".join(pur.tags), status="ok"),
            _p("Strategy", pur.strategy, status="info", ref="Banki 2005 / IMAC")]})
        if pur.elp:
            e = pur.elp
            cp.append({"key": "elp", "title": "ELP inverse transition cycling (ITC)",
                       "symbol": "thermo", "organelle": "bench", "level": "posttranslation", "params": [
                _p("ELP aggregation propensity", f"{e['aggregation_index']:.0f}/100",
                   scale=e["aggregation_index"], status="ok", ref=e["ref"], detail=e["band"]),
                _p("ELP length", f"{e['pentamers']} x VPGXG", status="info",
                   detail=e["levers"]),
                _p("Trigger", "warm above Tt and/or add (NH4)2SO4", status="info",
                   ref="Meyer & Chilkoti 1999")]})
        if pur.affinity:
            a = pur.affinity
            cp.append({"key": "affinity", "title": "Affinity capture", "symbol": "bead",
                       "organelle": "bench", "level": "posttranslation", "params": [
                _p("Resin", a.get("resin", a.get("tag", "?")), status="ok",
                   ref=a.get("ref", "")),
                _p("Elution", a.get("elute", "tag-specific"), status="info")]})
        if pur.intein:
            io = pur.intein
            cp.append({"key": "intein", "title": "Intein self-cleavage (tag removal)",
                       "symbol": "scissors", "organelle": "bench", "level": "posttranslation", "params": [
                _p("Cleavage control", "N- vs C-terminal, trigger-controlled", status="ok",
                   ref="Wood & Camarero 2014",
                   detail="Options: " + "; ".join(
                       f"{k}: {v['cleavage']} via {v['trigger']}"
                       for k, v in io["options"].items())),
                _p("Mutations that switch cleavage",
                   f"{len(io['mutations'])} documented", status="info",
                   detail="; ".join(f"{m['mutation']} -> {m['use']}" for m in io["mutations"]))]})
        cp.append({"key": "pure", "title": "Purified, tagless protein", "symbol": "vial",
                   "organelle": "bench", "level": "posttranslation", "params": [
            _p("Product", "tagless target recovered", status="ok",
               ref="Banki 2005", detail="Non-chromatographic where ELP-intein is used.")]})
    else:
        cp.append({"key": "notag", "title": "No purification tag", "symbol": "bench",
                   "organelle": "bench", "level": "posttranslation", "params": [
            _p("Recovery", "add an ELP-intein (non-chromatographic) or His/affinity tag",
               status="warn", ref="Banki 2005")]})

    return {
        "name": transcript.name, "epi": round(expression.epi, 1), "band": expression.band,
        "loc": loc, "landscape": expression.landscape,
        "biophys": bio, "structure": struct.to_dict(),
        "purification": pur.to_dict(), "checkpoints": cp,
    }
