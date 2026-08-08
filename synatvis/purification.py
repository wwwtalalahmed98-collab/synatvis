"""Downstream purification prediction — non-chromatographic & self-cleaving tags.

Extends the gene->protein journey past localisation into how the protein is
*recovered*, using published relationships (not a fitted model): elastin-like
polypeptide (ELP) inverse-transition cycling, self-cleaving intein tags (with the
mutations that switch N- vs C-terminal cleavage), and affinity capture.

It detects tags in the construct's protein and reports the parameters and protocol
each implies. Every number traces to a formula or a citation; where a value needs
wet-lab calibration it is labelled qualitative. Key references:
  * Banki, Feng & Wood 2005, Nat. Methods  — ELP-intein self-cleaving, non-chromatographic
  * Meyer & Chilkoti 1999/2004               — ELP LCST (Tt) vs guest residue, length, [salt]
  * Urry 1997                                 — guest-residue hydrophobicity / Tt ordering
  * Chong et al. 1997 (IMPACT); Perler InBase — controllable intein cleavage
  * Wood & Camarero 2014, JBC review          — engineered inteins for purification
  * Shah & Muir 2014; Mootz                   — split inteins (Npu DnaE) & ligation
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .ptm import _KD, translate, biophysics, predict_protein_fate

# ELP repeat: (VPGXG)n, guest residue X = anything but Pro (Urry)
_ELP = re.compile(r"(?:VPG[^P]G){2,}")
_HIS = re.compile(r"H{5,}")            # His5/His6+ affinity tag
# a few common purification tag signatures (protein level)
_TAGS = {
    "Strep-II": re.compile(r"WSHPQFEK"),
    "FLAG": re.compile(r"DYKDDDDK"),
    "HA": re.compile(r"YPYDVPDYA"),
    "c-Myc": re.compile(r"EQKLISEEDL"),
    "MBP-linker": re.compile(r"NSSSNNNNNNNNNNLGIEGR"),
}

# engineered inteins used for self-cleaving purification (curated from the refs above)
INTEIN_REF = {
    "Mxe GyrA": {
        "cleavage": "C-terminal", "trigger": "thiol (DTT / MESNA / β-ME)",
        "note": "target on the N-side; thiol induces C-terminal cleavage, leaving a "
                "C-terminal thioester (basis of expressed-protein ligation).",
        "ref": "Chong 1997 (IMPACT); Muir EPL"},
    "Ssp DnaB (mini)": {
        "cleavage": "N-terminal", "trigger": "pH shift + temperature",
        "note": "target on the C-side; low pH / raised temperature triggers N-terminal "
                "cleavage to release a tagless target.",
        "ref": "Mathys 1999; Wood 2000"},
    "Sce VMA": {
        "cleavage": "C-terminal", "trigger": "thiol (DTT)",
        "note": "classic yeast intein; DTT-induced C-terminal cleavage.",
        "ref": "Chong 1997"},
    "Npu DnaE (split)": {
        "cleavage": "trans-splicing / traceless ligation", "trigger": "spontaneous on N+C assembly",
        "note": "fastest naturally split intein; used for traceless protein ligation and "
                "C-terminal labelling rather than simple tag removal.",
        "ref": "Zettler 2009; Shah & Muir 2014"},
}

# mutations that redirect intein chemistry (N- vs C-terminal cleavage control)
INTEIN_MUTATIONS = [
    {"mutation": "Block1 Cys1/Ser1 → Ala",
     "effect": "abolishes N-terminal cleavage (no N-S/N-O acyl shift, no thioester)",
     "use": "force C-terminal-only cleavage", "ref": "Chong 1998; Wood 2014"},
    {"mutation": "C-terminal Asn → Ala (block Asn cyclisation)",
     "effect": "abolishes C-terminal cleavage / splicing",
     "use": "force N-terminal-only cleavage", "ref": "Xu & Perler 1996"},
    {"mutation": "+1 nucleophile (Cys/Ser/Thr of C-extein) → Ala",
     "effect": "blocks branch resolution → no C-terminal cleavage",
     "use": "controllable N-terminal cleavage tags", "ref": "Southworth 1999"},
    {"mutation": "penultimate His → Ala",
     "effect": "slows C-terminal Asn cyclisation (reduces premature cleavage)",
     "use": "tighten trigger control", "ref": "Wood 2014"},
]


@dataclass
class Purification:
    tags: List[str] = field(default_factory=list)
    strategy: str = "standard chromatography"
    elp: Optional[Dict] = None
    intein: Optional[Dict] = None
    affinity: Optional[Dict] = None
    steps: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)
    # --- purification-strategy layer (chromatographic + non-chromatographic) ---
    recommendations: List[Dict] = field(default_factory=list)  # every applicable method
    workflow: List[Dict] = field(default_factory=list)         # ordered Capture->Intermediate->Polish
    considerations: List[Dict] = field(default_factory=list)   # source / redox / detergent flags
    biophysics: Dict = field(default_factory=dict)
    localization: str = ""

    def to_dict(self) -> Dict:
        return {"tags": self.tags, "strategy": self.strategy, "elp": self.elp,
                "intein": self.intein, "affinity": self.affinity,
                "steps": self.steps, "notes": self.notes,
                "recommendations": self.recommendations, "workflow": self.workflow,
                "considerations": self.considerations, "biophysics": self.biophysics,
                "localization": self.localization}


def _iex_by_pi(pI: float) -> Dict:
    """Ion-exchange choice + working pH from the theoretical pI (classic CIPP rule)."""
    if pI <= 6.0:
        return {"method": "Anion exchange (Q)", "mode": "chromatographic", "phase": "capture",
                "basis": f"pI ≈ {pI}: net-negative near pH 7.5–8.0 → binds a Q (quaternary-amine) resin",
                "buffer": "bind at pH 7.5–8.0 in low salt; elute on a NaCl gradient", "ref": "IEX theory"}
    if pI >= 8.0:
        return {"method": "Cation exchange (S)", "mode": "chromatographic", "phase": "capture",
                "basis": f"pI ≈ {pI}: net-positive near pH 5.5–6.5 → binds an S (sulfopropyl) resin",
                "buffer": "bind at pH 5.5–6.5 in low salt; elute on a NaCl gradient", "ref": "IEX theory"}
    return {"method": "Ion exchange (weak — near-neutral pI)", "mode": "chromatographic", "phase": "capture",
            "basis": f"pI ≈ {pI} is near neutral — IEX selectivity is weak; prefer HIC, an affinity tag, or SEC",
            "buffer": "try Q at pH 8.5 or S at pH 5.0; expect modest resolution", "ref": "IEX theory"}


def recommend_purification(bio: Dict, fate, his, other, elp) -> "tuple[List[Dict], List[Dict]]":
    """Full chromatographic + non-chromatographic plan from intrinsic biophysics + tags."""
    pI = bio.get("pI", 7.0); mw = bio.get("mw_kda", 0.0); gravy = bio.get("gravy", 0.0)
    loc = fate.localization
    recs: List[Dict] = []
    cons: List[Dict] = []

    # source / feedstock
    if "secreted" in loc:
        cons.append({"topic": "Source", "tag": "from medium",
                     "note": "Secreted (signal peptide, no TM) — recover from clarified culture medium; "
                             "no cell lysis and far fewer host contaminants."})
    elif fate.tm_count > 0:
        cons.append({"topic": "Source", "tag": "detergent",
                     "note": f"Predicted membrane protein ({fate.tm_count} TM helix/helices) — solubilise "
                             "in a mild detergent (e.g. DDM) before any column."})
    else:
        cons.append({"topic": "Source", "tag": "lyse",
                     "note": "Cytosolic — lyse cells and clarify the lysate before capture."})

    # capture: tag-based first, else property-based
    if his:
        recs.append({"method": "IMAC (Ni-NTA / Co-TALON)", "mode": "chromatographic", "phase": "capture",
                     "basis": f"His{len(his.group(0))} tag detected — immobilised-metal affinity capture",
                     "buffer": "bind low-imidazole; elute on an imidazole gradient (or low pH)",
                     "ref": "Porath 1975; Hochuli 1987"})
    elif other:
        recs.append({"method": f"{other[0]} affinity", "mode": "chromatographic", "phase": "capture",
                     "basis": f"{other[0]} tag detected — tag-specific affinity resin",
                     "buffer": "tag-specific elution", "ref": "tag supplier"})
    elif elp:
        recs.append({"method": "ELP inverse transition cycling (ITC)", "mode": "non-chromatographic",
                     "phase": "capture", "basis": f"ELP tag ({elp['pentamers']}× VPGXG) — salt/temperature "
                     "phase separation; no resin", "buffer": elp["band"], "ref": "Meyer & Chilkoti 1999; Banki 2005"})
    else:
        recs.append(_iex_by_pi(pI))

    # tagless alternative capture (if a tag drove the primary capture)
    if his or other or elp:
        alt = _iex_by_pi(pI); alt["phase"] = "capture (tagless alternative)"
        recs.append(alt)

    # non-chromatographic concentration / partial purification
    recs.append({"method": "Ammonium-sulfate precipitation (salting-out)", "mode": "non-chromatographic",
                 "phase": "capture / concentrate", "basis": f"cheap, scalable concentration by hydrophobicity "
                 f"(GRAVY {gravy}); a good front step before HIC", "buffer": "raise (NH4)2SO4 stepwise; spin; "
                 "resolubilise the pellet", "ref": "Wingfield 1998"})

    # intermediate — HIC when relatively hydrophobic
    if gravy > -0.2:
        recs.append({"method": "Hydrophobic interaction (HIC)", "mode": "chromatographic", "phase": "intermediate",
                     "basis": f"relatively hydrophobic (GRAVY {gravy}) — binds in high (NH4)2SO4, elutes at low "
                     "salt; pairs naturally after salting-out or high-salt IEX", "buffer": "bind ~1 M (NH4)2SO4; "
                     "elute on a descending-salt gradient", "ref": "Queiroz 2001"})

    # polish — SEC always
    recs.append({"method": "Size-exclusion (SEC)", "mode": "chromatographic", "phase": "polish",
                 "basis": f"MW ≈ {mw} kDa — final polishing, buffer exchange, and oligomeric-state check; pick a "
                 f"matrix whose fractionation range spans {mw} kDa", "buffer": "isocratic; removes aggregates",
                 "ref": "SEC theory"})

    # redox / cofactor considerations
    if fate.disulfide_potential:
        cons.append({"topic": "Redox", "tag": "non-reducing",
                     "note": "Even cysteine count → likely disulfide bonds; keep buffers non-reducing (no DTT/"
                             "β-ME). Note this conflicts with thiol-inducible intein cleavage — choose a "
                             "pH/temperature-controlled intein instead."})
    return recs, cons


def _elp_metrics(prot: str) -> Optional[Dict]:
    m = _ELP.search(prot)
    if not m:
        return None
    block = m.group(0)
    n_pent = len(block) // 5
    guests = [block[i + 3] for i in range(0, len(block) - 4, 5)]  # X in VPGXG
    # relative Tt driver: more hydrophobic guest + longer chain => LOWER Tt (aggregates sooner)
    hydro = sum(_KD.get(g, 0.0) for g in guests) / max(1, len(guests))
    # transparent 0-100 "aggregation propensity" index (higher = lower Tt); NOT a temperature
    idx = max(0.0, min(100.0, 50 + 8 * hydro + 0.4 * (n_pent - 20)))
    if idx >= 66:
        band = "low transition temperature — aggregates near/below room temp (easy warm ITC)"
    elif idx >= 40:
        band = "moderate Tt — trigger with mild warming or kosmotropic salt ((NH4)2SO4)"
    else:
        band = "high Tt — needs more salt / higher temperature to phase-separate"
    return {"pentamers": n_pent, "guest_hydropathy": round(hydro, 2),
            "aggregation_index": round(idx, 0), "band": band,
            "levers": "Tt falls with: longer ELP, more hydrophobic guest residue, higher "
                      "protein concentration, and higher kosmotropic salt (Meyer & Chilkoti).",
            "ref": "Meyer & Chilkoti 1999/2004; Urry 1997"}


def _detect_intein(prot: str, notes: List[str]) -> Optional[Dict]:
    # sequence-level intein identification is unreliable; report the design reference so the
    # user selects the intein, rather than guessing one from sequence.
    return {"identified": False,
            "options": INTEIN_REF,
            "mutations": INTEIN_MUTATIONS,
            "note": "Specify the intein used (identity is not reliably callable from "
                    "sequence). The table lists cleavage direction, trigger and the "
                    "mutations that switch N- vs C-terminal cleavage."}


def predict_purification(transcript) -> Purification:
    prot = translate(transcript.cds)
    p = Purification()
    if not prot:
        p.notes.append("empty translation")
        return p

    elp = _elp_metrics(prot)
    his = _HIS.search(prot)
    other = [name for name, rx in _TAGS.items() if rx.search(prot)]

    if his:
        p.tags.append(f"His{len(his.group(0))}")
    if elp:
        p.tags.append(f"ELP ({elp['pentamers']}x VPGXG)")
    p.tags.extend(other)

    # choose a strategy, richest-tag first
    if elp:
        p.elp = elp
        p.strategy = "non-chromatographic — ELP inverse transition cycling (ITC)"
        p.intein = _detect_intein(prot, p.notes)
        p.steps = [
            "Express ELP–(intein)–target fusion; clarify lysate.",
            "Trigger phase transition above Tt: warm and/or add kosmotropic salt "
            "((NH4)2SO4) → ELP-fusion coacervates.",
            "Warm centrifugation → pellet the ELP-fusion; discard soluble contaminants.",
            "Resolubilise the pellet in cold buffer; cold spin removes irreversible aggregates.",
            "Trigger intein self-cleavage (thiol or pH/T per intein) → release tagless target.",
            "Second ITC → aggregate & pellet the ELP–intein; recover pure target in the "
            "supernatant. No columns, no resin.",
        ]
        p.notes.append("Self-cleaving ELP-intein tag: non-chromatographic capture + traceless "
                       "tag removal in one system (Banki, Feng & Wood 2005).")
    elif his:
        p.strategy = "affinity — immobilised-metal (IMAC, Ni-NTA)"
        p.affinity = {"resin": "Ni-NTA / Co-TALON", "bind": "native or denaturing",
                      "elute": "imidazole gradient (or low pH)",
                      "ref": "Porath 1975; Hochuli 1987"}
        p.intein = _detect_intein(prot, p.notes)
        p.steps = [
            "Bind His-tagged protein to Ni-NTA; wash low-imidazole.",
            "Elute with imidazole gradient (or drop pH).",
            "Optional: on-column intein self-cleavage releases tagless target directly, "
            "leaving the tag+intein bound (self-cleaving IMAC).",
        ]
    elif other:
        p.strategy = f"affinity — {other[0]} tag"
        p.affinity = {"tag": other[0], "ref": "tag-specific resin"}
    else:
        p.notes.append("No purification tag detected — the plan below is built from the protein's own "
                       "properties; adding an ELP-intein (non-chromatographic) or His tag would simplify capture.")
        p.intein = None

    # --- purification-strategy layer (works for every protein, tag or not) ---
    bio = biophysics(prot)
    fate = predict_protein_fate(transcript)
    recs, cons = recommend_purification(bio, fate, his, other, elp)
    p.biophysics = bio
    p.localization = fate.localization
    p.recommendations = recs
    p.considerations = cons
    if not p.tags:
        p.strategy = recs[0]["method"] + (" (property-based capture)"
                                          if recs[0]["mode"] == "chromatographic" else "")
    # ordered Capture -> Intermediate -> Polish workflow (one method per phase)
    def _pick(phase_key):
        for r in recs:
            if r["phase"].startswith(phase_key):
                return {"phase": r["phase"], "method": r["method"], "mode": r["mode"], "basis": r["basis"]}
        return None
    src = cons[0] if cons else None
    p.workflow = [w for w in (
        {"phase": "source", "method": src["note"].split(" — ")[0] if src else "clarify feedstock",
         "mode": "prep", "basis": src["note"] if src else ""} if src else None,
        _pick("capture"), _pick("intermediate"), _pick("polish"),
    ) if w]
    return p
