"""Report rendering — text and JSON, NO composite score (CLAUDE.md §1, §5, §9).

The banner states host, compartment, the transcript-level-only scope, the
proteolysis/protein-stability blind spot, and the no-score notice. Flags are
grouped by module and ranked by severity; any module with ``validated=False``
(silencing) is marked HEURISTIC in-line. Both a human-readable and a JSON form
are produced. No probability, no grade.
"""
from __future__ import annotations

import json
from typing import Any, Dict, List

from .flags import Severity, sort_flags
from .scanner import ScanResult

_BANNER_LIMITATIONS = [
    "Transcript-level only: splicing/processing, mRNA stability, translation "
    "initiation, silencing-risk features.",
    "SILENT on proteolysis and protein stability — a real cause of low yield it "
    "cannot see. A clean report is NOT a promise of expression.",
    "Diagnostic, not predictive: NO composite expression score is emitted.",
    "This nuclear report says nothing about a chloroplast construct.",
    "poly(A) is the noisiest module; its operating point comes from Leg-1, not "
    "assumption. The shipped codon table is a starting point (see provenance).",
]


def _banner(result: ScanResult) -> List[str]:
    meta = result.profile.get("meta", {})
    lines = [
        "=" * 74,
        "  SynAT.Vis — Synthetic Algal Transcript red-flag report",
        "=" * 74,
        f"  Host        : {meta.get('host', '?')}",
        f"  Compartment : {meta.get('compartment', '?')}  "
        f"(lineage: {meta.get('lineage', '?')}; gate: {meta.get('gating_axis', '?')})",
        f"  Transcript  : {result.transcript.name}  "
        f"(5'UTR {len(result.transcript.utr5)} / CDS {len(result.transcript.cds)} / "
        f"3'UTR {len(result.transcript.utr3)} nt)",
        f"  Codon table : {result.codon_source}",
        "-" * 74,
        "  SCOPE & LIMITATIONS:",
    ]
    for lim in _BANNER_LIMITATIONS:
        wrapped = _wrap(lim, 68)
        lines.append(f"   - {wrapped[0]}")
        lines.extend(f"     {w}" for w in wrapped[1:])
    for note in result.transcript.notes:
        lines.append(f"   ! input note: {note}")
    lines.append("=" * 74)
    return lines


def _wrap(text: str, width: int) -> List[str]:
    words, line, out = text.split(), "", []
    for w in words:
        if line and len(line) + 1 + len(w) > width:
            out.append(line)
            line = w
        else:
            line = f"{line} {w}".strip()
    if line:
        out.append(line)
    return out or [""]


def render_text(result: ScanResult) -> str:
    lines = _banner(result)
    counts = result.counts()
    lines.append(
        f"  SUMMARY: {len(result.flags)} flag(s) — "
        f"high {counts['high']}, medium {counts['medium']}, "
        f"low {counts['low']}, info {counts['info']}."
    )
    lines.append("")

    by_mod = result.by_module()
    # module display order follows registration order in module_meta
    for name, meta in result.module_meta.items():
        mod_flags = sort_flags(by_mod.get(name, []))
        if not mod_flags:
            continue
        tag = "  [HEURISTIC — NOT VALIDATED]" if not meta["validated"] else ""
        lines.append(f"### {name}{tag}")
        lines.append(f"    {meta['summary']}")
        for f in mod_flags:
            sev = str(f.severity).upper()
            span = f"{f.start}-{f.end}"
            lines.append(f"  [{sev:<6}] {f.region:<10} {span:<12} {f.message}")
            if f.evidence:
                lines.append(f"           evidence: {f.evidence}")
            if f.suggested_edit:
                lines.append(f"           fix     : {f.suggested_edit}")
        lines.append("")

    if not result.flags:
        lines.append("  No flags raised. (Reminder: a clean report is not a promise "
                     "of expression — see limitations above.)")
    return "\n".join(lines)


def render_plugins_text(results) -> str:
    """Format Tier-B plugin predictions as a labelled, clearly-experimental block."""
    if not results:
        return ("\n  (Tier-B ML plugins: none active. Run `synatvis plugins` to see "
                "what to install to enable CodonBERT / Saluki / UTR-LM / LinearDesign.)")
    lines = ["", "### EXPERIMENTAL — Tier-B ML plugins (opt-in, NOT part of the "
             "validated result)"]
    for r in results:
        val = f"  = {r.value}" if r.value is not None else ""
        lines.append(f"  [{r.plugin}] {r.label}{val}")
        if r.text:
            lines.append(f"        {r.text}")
        if r.note:
            lines.append(f"        note: {r.note}")
        if r.citation:
            lines.append(f"        ref : {r.citation}")
    return "\n".join(lines)


def render_expression_text(expr) -> str:
    """Format the expression-propensity ensemble as a labelled block."""
    L = ["", "### PREDICTED EXPRESSION PROPENSITY  (ensemble of models; not a measured yield)"]
    L.append(f"  Index: {expr.epi:.0f}/100  —  {expr.band}")
    for m in expr.models:
        if m["available"]:
            bar = "#" * int(round(m["score"] / 5))
            L.append(f"    {m['label']:<46} {m['score']:>3.0f}/100  {bar}")
        else:
            L.append(f"    {m['label']:<46} (inactive — install ViennaRNA)")
    for r in expr.ml_readouts:
        v = f" = {r['value']}" if r.get("value") is not None else ""
        L.append(f"    [ML: {r['plugin']}] {r['label']}{v}")
    for h in expr.hazards:
        L.append(f"    hazard gate x{h['factor']}: {h['reason']}")
    L.append(f"  {expr.confidence}")
    if expr.fate:
        f = expr.fate
        L.append("")
        L.append("  Protein fate & PTM (protein-level layer; sequence heuristics):")
        L.append(f"    localisation : {f.get('localization')}"
                 f"  |  signal peptide: {'yes' if f.get('signal_peptide') else 'no'}"
                 f"  |  TM: {f.get('tm_count', 0)}"
                 f"  |  N-glyc sequons: {f.get('n_glyc_sites', 0)}"
                 f"  |  Cys: {f.get('cys_count', 0)}")
        for n in f.get("notes", []):
            L.append(f"      - {n}")
    return "\n".join(L)


def _landscape_svg(track, w=680, h=66):
    """Inline SVG of per-codon adaptiveness along the CDS (dips = slow/rare codons)."""
    if not track:
        return ""
    step = max(1, len(track) // 240)
    pts = [sum(track[i:i + step]) / len(track[i:i + step]) for i in range(0, len(track), step)]
    m = len(pts)

    def X(i):
        return (i / (m - 1) * w) if m > 1 else 0.0

    def Y(v):
        return h - 3 - max(0.0, min(1.0, v)) * (h - 9)

    line = "M" + " L".join(f"{X(i):.1f},{Y(v):.1f}" for i, v in enumerate(pts))
    area = f"M0,{h} " + " ".join(f"L{X(i):.1f},{Y(v):.1f}" for i, v in enumerate(pts)) + f" L{w},{h} Z"
    lows = "".join(
        f'<rect x="{X(i):.1f}" y="0" width="{max(2, w / m):.1f}" height="{h}" fill="#d03b3b" opacity="0.14"/>'
        for i, v in enumerate(pts) if v < 0.3)
    return (f'<svg class="land" viewBox="0 0 {w} {h}" preserveAspectRatio="none" '
            f'role="img" aria-label="per-codon adaptiveness along the coding sequence">'
            f'{lows}<path d="{area}" fill="var(--accent)" opacity="0.12"/>'
            f'<path d="{line}" fill="none" stroke="var(--accent)" stroke-width="1.4"/>'
            f'<line x1="0" y1="{Y(0.3):.1f}" x2="{w}" y2="{Y(0.3):.1f}" stroke="#d03b3b" '
            f'stroke-width="0.6" stroke-dasharray="3 3" opacity="0.5"/></svg>')


def render_expression_html(expr) -> str:
    import html as _h
    color = "#0ca30c" if expr.epi >= 65 else ("#e07b39" if expr.epi >= 40 else "#d03b3b")
    axbars = ""
    for m in expr.models:
        if m["available"]:
            axbars += (
                f'<div class="axrow"><span class="axl">{_h.escape(m["label"])}</span>'
                f'<span class="axtrack"><span class="axfill" style="width:{m["score"]:.0f}%"></span></span>'
                f'<span class="axp">{m["score"]:.0f}</span></div>')
        else:
            axbars += (
                f'<div class="axrow off"><span class="axl">{_h.escape(m["label"])}</span>'
                f'<span class="axtrack"></span><span class="axp">—</span></div>')
    for r in expr.ml_readouts:
        v = f' = {r["value"]}' if r.get("value") is not None else ""
        axbars += (f'<div class="axrow"><span class="axl">{_h.escape(r["label"])} '
                   f'<span class="tag exp">{_h.escape(r["plugin"])}</span></span>'
                   f'<span class="axtrack ml"></span><span class="axp">{_h.escape(str(v))}</span></div>')
    haz = ""
    if expr.hazards:
        haz = '<div class="hazrow">' + "".join(
            f'<span class="hchip">×{h["factor"]} {_h.escape(h["reason"].split(" — ")[0])}</span>'
            for h in expr.hazards) + "</div>"
    return (
        '<h2 class="svh" style="--sc:#256abf">Predicted expression propensity'
        ' <span class="svsub">— a model, not a measured yield</span></h2>'
        '<div class="expcard">'
        f'<div class="gauge"><div class="gnum" style="color:{color}">{expr.epi:.0f}<small>/100</small></div>'
        f'<div class="gband" style="color:{color}">{_h.escape(expr.band)}</div>'
        f'<div class="gbar"><div class="gfill" style="width:{expr.epi:.0f}%;background:{color}"></div></div></div>'
        f'<div class="axes">{axbars}</div>{haz}'
        '<div class="landwrap"><div class="landlbl">Per-codon adaptiveness along the coding '
        'sequence — dips (below the dashed line) mark slow, expression-limiting regions</div>'
        f'{_landscape_svg(expr.landscape)}</div>'
        f'<div class="ev">{_h.escape(expr.confidence)}</div></div>')


def render_fate_html(fate) -> str:
    """Protein-fate & PTM panel (complementary protein-level layer)."""
    import html as _h
    if not fate:
        return ""
    loc = fate.get("localization", "cytosolic")
    chips = [("Localisation", loc)]
    if fate.get("signal_peptide"):
        chips.append(("Signal peptide", "yes → secretory pathway"))
    if fate.get("tm_count"):
        chips.append(("Transmembrane", f'{fate["tm_count"]} helix/helices'))
    if fate.get("n_glyc_sites"):
        chips.append(("N-glyc sequons", f'{fate["n_glyc_sites"]} (N-X-S/T)'))
    if fate.get("cys_count"):
        ds = " (possible disulfides)" if fate.get("disulfide_potential") else ""
        chips.append(("Cysteines", f'{fate["cys_count"]}{ds}'))
    chip_html = "".join(
        f'<span class="fchip"><b>{_h.escape(str(v))}</b><span>{_h.escape(k)}</span></span>'
        for k, v in chips)
    notes = "".join(f'<li>{_h.escape(n)}</li>' for n in fate.get("notes", []))
    notes_html = f'<ul class="fnotes">{notes}</ul>' if notes else ""
    modrow = ""
    if fate.get("modifier", 1.0) != 1.0:
        modrow = (f'<div class="hazrow"><span class="hchip">×{fate["modifier"]} '
                  f'{_h.escape(fate.get("modifier_reason",""))}</span></div>')
    return (
        '<h2 class="svh" style="--sc:#2e8d6b">Protein fate &amp; modifications'
        ' <span class="svsub">— translated protein-level layer, sequence heuristics</span></h2>'
        f'<div class="expcard"><div class="fchips">{chip_html}</div>{modrow}{notes_html}'
        '<div class="ev">Heuristic, sequence-based (SignalP / DeepLoc / NetNGlyc concepts). '
        'These govern where the protein goes, not how much is made; verify with the '
        'dedicated predictors for decisions.</div></div>')


_LV_ICON = {"pass": ("&#10003;", "#1f9d55"), "warn": ("!", "#d98a1f"),
            "fail": ("&#10007;", "#cf4b3a")}


def render_levels_html(levels, verdict) -> str:
    """Six-stage check-mark pipeline: how the protein gets made, level by level."""
    import html as _h
    tiles = []
    for lv in levels:
        icon, col = _LV_ICON[lv["status"]]
        sc = (f'<span class="lvscore">{lv["score"]:.0f}</span>'
              if lv.get("score") is not None else "")
        issues = ""
        if lv["issues"]:
            items = "".join(f'<li>{_h.escape(i["label"])}</li>' for i in lv["issues"][:3])
            issues = f'<ul class="lvissues">{items}</ul>'
        else:
            issues = '<div class="lvok">clear</div>'
        tiles.append(
            f'<div class="lvtile" style="--lc:{col}">'
            f'<div class="lvtop"><span class="lvnum">{lv["order"]}</span>'
            f'<span class="lvic" style="color:{col}">{icon}</span>{sc}</div>'
            f'<div class="lvname">{_h.escape(lv["name"])}</div>'
            f'<div class="lvdesc">{_h.escape(lv["desc"])}</div>{issues}</div>')
    vcol = _LV_ICON[verdict["verdict"]][1]
    head = (f'{verdict["passes"]} of {verdict["n"]} levels clear'
            + (f' &middot; first block at <b>{_h.escape(verdict["blocker"])}</b>'
               if verdict["blocker"] else ' &middot; full path to protein is clear'))
    return (
        '<h2 class="svh" style="--sc:#1a3c5e">How the protein gets made — six checkpoints'
        f' <span class="svsub" style="color:{vcol}">{head}</span></h2>'
        f'<div class="lvpipe">{"".join(tiles)}</div>')


def render_levels_text(levels, verdict) -> str:
    L = ["", "### PROTEIN-MAKING CHECKPOINTS (six molecular-biology levels)"]
    mark = {"pass": "[PASS]", "warn": "[WARN]", "fail": "[FAIL]"}
    for lv in levels:
        sc = f"  score {lv['score']:.0f}/100" if lv.get("score") is not None else ""
        L.append(f"  {mark[lv['status']]}  {lv['order']}. {lv['name']:<20} "
                 f"({lv['n_checks']} checks{sc})")
        for i in lv["issues"]:
            L.append(f"           - {i['label']}: {i['value']}")
    v = verdict
    tail = (f"blocked at {v['blocker']}" if v["blocker"] else "full path clear")
    L.append(f"  => {v['passes']}/{v['n']} levels clear; {tail}.")
    return "\n".join(L)


def render_cell_embed(journey) -> str:
    """Embed the animated gene->protein->purified journey as an isolated iframe."""
    from .viz_cell import render_cell_document
    doc = render_cell_document(journey)
    srcdoc = doc.replace("&", "&amp;").replace('"', "&quot;")
    return (
        '<h2 class="svh" style="--sc:#256abf">Live journey — gene to purified protein'
        ' <span class="svsub">— every checkpoint, parameter &amp; scale, with citations</span></h2>'
        f'<iframe class="celliframe" srcdoc="{srcdoc}" '
        'title="animated gene-to-protein journey" loading="lazy"></iframe>')


def render_json(result: ScanResult) -> str:
    payload: Dict[str, Any] = {
        "tool": "SynAT.Vis",
        "host": result.profile.get("meta", {}).get("host"),
        "compartment": result.profile.get("meta", {}).get("compartment"),
        "transcript": {
            "name": result.transcript.name,
            "utr5_len": len(result.transcript.utr5),
            "cds_len": len(result.transcript.cds),
            "utr3_len": len(result.transcript.utr3),
            "notes": result.transcript.notes,
        },
        "codon_table_source": result.codon_source,
        "no_composite_score": True,
        "counts": result.counts(),
        "modules": {
            name: {"validated": meta["validated"], "summary": meta["summary"]}
            for name, meta in result.module_meta.items()
        },
        "flags": [f.to_dict() for f in result.flags],
        "limitations": _BANNER_LIMITATIONS,
    }
    return json.dumps(payload, indent=2)


# ---------------------------------------------------------------------------
# Plain-language report (for non-programmers; `scan --plain`)
# ---------------------------------------------------------------------------
_PLAIN_TITLE = {
    "cloning": "Hidden lab cut-site inside your gene",
    "polya": 'Premature "stop the message here" signal',
    "composition": "Low-GC region (can get switched off) or long repeat",
    "codon": "Genetic 'spelling' the alga reads slowly",
    "uorf": "Extra start signal before your real start",
    "structure": "Start-of-gene readability",
    "instability": "Element that makes the message break down faster",
    "splice": "Possible splicing site",
    "silencing": "Gene-silencing risk (rough guidance only)",
}

_PLAIN_SEVERITY = {
    Severity.HIGH: ("SERIOUS  -  fix these first", "!!"),
    Severity.MEDIUM: ("WORTH CHECKING", "! "),
    Severity.LOW: ("MINOR", "- "),
    Severity.INFO: ("FOR YOUR INFORMATION", "  "),
}

_PLAIN_REGION = {
    "5utr": "the region before the gene",
    "cds": "the coding sequence (the gene itself)",
    "3utr": "the region after the gene",
    "transcript": "the whole message",
}

# ordered longest/most-specific first so replacements don't collide
_GLOSSARY = [
    ("Relative Codon Adaptation", "codon-fit score"),
    ("AU-rich element", "message-destabilising element"),
    ("AU-rich", "message-destabilising"),
    ("Type IIS site", "lab cloning cut-site"),
    ("Type IIS", "lab cloning"),
    ("MoClo assembly", "Golden Gate assembly"),
    ("MoClo", "Golden Gate"),
    ("Domesticate", "Remove"),
    ("domesticate", "remove"),
    ("no clean synonymous fix", "no clean silent fix"),
    ("synonymous:", "silent DNA change (protein stays identical):"),
    ("synonymous", "silent (protein unchanged)"),
    ("heterochromatinisation", "being switched off"),
    ("poly(A) signal", "end-of-message signal (poly-A)"),
    ("poly(A)", "end-of-message (poly-A)"),
    ("NUE", "end-of-message signal"),
    ("uORF", "extra mini-gene"),
    ("uAUG", "extra start signal"),
    ("RCA=", "codon-fit="),
    ("RCA", "codon-fit score"),
    ("CDS", "coding sequence"),
    ("5'UTR", "region before the gene"),
    ("3'UTR", "region after the gene"),
    ("mRNA", "message"),
    ("transcript", "message"),
]


def _plainify(text: str) -> str:
    if not text:
        return text
    for jargon, plain in _GLOSSARY:
        text = text.replace(jargon, plain)
    return text


def render_plain(result: ScanResult) -> str:
    meta = result.profile.get("meta", {})
    L = []
    L.append("=" * 70)
    L.append("  SynAT.Vis - plain-language report")
    L.append("=" * 70)
    L.append(f"  Your design : {result.transcript.name}")
    L.append(f"  Checked for : {meta.get('host', '?')} "
             f"({meta.get('compartment', '?')})")
    L.append("")
    L.append("  We read your DNA and looked for known problems that stop a gene")
    L.append("  from working in this alga. This checks the 'message' stage only -")
    L.append("  not the protein, and not how much you will get. A clean result")
    L.append("  means 'no obvious red flags', NOT 'guaranteed to work'.")
    L.append("")
    for note in result.transcript.notes:
        L.append(f"  Note: {_plainify(note)}")
    c = result.counts()
    L.append(f"  RESULT: {c['high']} serious, {c['medium']} worth checking, "
             f"{c['low']} minor, {c['info']} for information.")
    L.append("")

    by_sev = {}
    for f in result.flags:
        by_sev.setdefault(f.severity, []).append(f)

    for sev in (Severity.HIGH, Severity.MEDIUM, Severity.LOW, Severity.INFO):
        flags = by_sev.get(sev, [])
        if not flags:
            continue
        header, _mark = _PLAIN_SEVERITY[sev]
        L.append("-" * 70)
        L.append(f"  {header}")
        L.append("-" * 70)
        for f in flags:
            title = _PLAIN_TITLE.get(f.module, f.module)
            # the 3'UTR poly(A) signal is the normal, wanted one - don't alarm
            if f.module == "polya" and f.region == "3utr":
                title = "Normal end-of-message signal (expected here - good)"
            heuristic = not result.module_meta.get(f.module, {}).get("validated", True)
            if heuristic:
                title += "  (rough guess)"
            where = _PLAIN_REGION.get(f.region, f.region)
            L.append(f"  * {title}")
            L.append(f"      Where : {where}, around letter {f.start}-{f.end}")
            L.append(f"      What  : {_plainify(f.message)}")
            if f.suggested_edit:
                L.append(f"      Fix   : {_plainify(f.suggested_edit)}")
            elif f.severity == Severity.INFO:
                L.append("      Action: none needed - this is just information")
            else:
                L.append("      Fix   : no automatic fix - may need a manual redesign")
            L.append("")

    if not result.flags:
        L.append("  No red flags found. (Remember: not a promise of expression.)")
        L.append("")
    L.append("=" * 70)
    L.append("  Reminders: cannot see protein stability * no expression score *")
    L.append("  Chlamydomonas nucleus only * a clean report is not a guarantee.")
    L.append("=" * 70)
    return "\n".join(L)


# ---------------------------------------------------------------------------
# HTML report — a polished, self-contained page (the tool's "response")
# ---------------------------------------------------------------------------
_HTML_SEV = {
    Severity.HIGH:   ("Serious", "fix first", "#d03b3b"),
    Severity.MEDIUM: ("Worth checking", "", "#e07b39"),
    Severity.LOW:    ("Minor", "", "#2a78d6"),
    Severity.INFO:   ("For information", "", "#6b7280"),
}


def render_html(result: ScanResult, plugin_results=None, expression=None) -> str:
    import html as _h
    exp_html = render_expression_html(expression) if expression else ""
    fate_html = render_fate_html(expression.fate) if expression else ""
    levels_html = cell_html = ""
    if expression:
        from .journey import build_journey
        from .levels import summarize_levels, overall_verdict
        journey = build_journey(result.transcript, result.profile, result, expression)
        levels = summarize_levels(journey)
        levels_html = render_levels_html(levels, overall_verdict(levels))
        cell_html = render_cell_embed(journey)
    # inline listener that self-sizes the embedded journey iframe (only when embedded)
    autosize_js = ("<script>\n"
                   "window.addEventListener('message', function(e){\n"
                   "  var d = e.data;\n"
                   "  if(!d || typeof d.synatvisHeight !== 'number') return;\n"
                   "  var f = document.querySelector('.celliframe');\n"
                   "  if(f && e.source === f.contentWindow){\n"
                   "    f.style.height = (Math.max(300, d.synatvisHeight) + 2) + 'px';\n"
                   "  }\n});\n</script>") if cell_html else ""

    meta = result.profile.get("meta", {})
    c = result.counts()
    by_sev = {}
    for f in result.flags:
        by_sev.setdefault(f.severity, []).append(f)

    # headline
    if c["high"]:
        headline = (f"{c['high']} serious issue{'s' if c['high'] != 1 else ''} to fix "
                    f"before you build this.")
        hero_class, hero_icon = "bad", "&#9888;"
    elif c["medium"]:
        headline = f"{c['medium']} thing{'s' if c['medium'] != 1 else ''} worth checking."
        hero_class, hero_icon = "warn", "&#9432;"
    else:
        headline = "No red flags for this host. (Not a promise of expression.)"
        hero_class, hero_icon = "good", "&#10003;"

    def pill(n, label, color):
        return (f'<span class="pill" style="--pc:{color}"><b>{n}</b>'
                f'<span>{label}</span></span>')

    pills = "".join([
        pill(c["high"], "serious", "#d03b3b"),
        pill(c["medium"], "worth checking", "#e07b39"),
        pill(c["low"], "minor", "#2a78d6"),
        pill(c["info"], "info", "#6b7280"),
    ])

    cards = []
    for sev in (Severity.HIGH, Severity.MEDIUM, Severity.LOW, Severity.INFO):
        flags = sort_flags(by_sev.get(sev, []))
        if not flags:
            continue
        label, sub, color = _HTML_SEV[sev]
        subtxt = f' <span class="svsub">— {sub}</span>' if sub else ""
        cards.append(f'<h2 class="svh" style="--sc:{color}">{label}{subtxt} '
                     f'<span class="cnt">{len(flags)}</span></h2>')
        for f in flags:
            title = _PLAIN_TITLE.get(f.module, f.module)
            if f.module == "polya" and f.region == "3utr":
                title = "Normal end-of-message signal (expected here — good)"
            heuristic = not result.module_meta.get(f.module, {}).get("validated", True)
            tag = f'<span class="tag">{_h.escape(f.module)}</span>'
            if heuristic:
                tag += '<span class="tag heur">heuristic</span>'
            where = _PLAIN_REGION.get(f.region, f.region)
            fix = ""
            if f.suggested_edit:
                fix = (f'<div class="fix"><span>Suggested fix</span>'
                       f'<code>{_h.escape(_plainify(f.suggested_edit))}</code></div>')
            elif f.severity != Severity.INFO:
                fix = ('<div class="fix none">No automatic fix — may need a manual '
                       'redesign.</div>')
            ev = (f'<div class="ev">{_h.escape(f.evidence)}</div>' if f.evidence else "")
            cards.append(
                f'<div class="card" style="--cc:{color}">'
                f'<div class="ctop"><span class="ctitle">{_h.escape(title)}</span>{tag}</div>'
                f'<div class="loc">{_h.escape(where)} · position {f.start}–{f.end}</div>'
                f'<p class="msg">{_h.escape(_plainify(f.message))}</p>{fix}{ev}</div>')

    plugin_html = ""
    if plugin_results:
        rows = []
        for r in plugin_results:
            v = f' <b>{r.value}</b>' if r.value is not None else ""
            rows.append(
                f'<div class="pcard"><div class="ctop"><span class="ctitle">'
                f'{_h.escape(r.label)}</span>{v}<span class="tag exp">{_h.escape(r.plugin)}</span></div>'
                f'<p class="msg">{_h.escape(r.text)}</p>'
                f'<div class="ev">{_h.escape(r.note)} · {_h.escape(r.citation)}</div></div>')
        plugin_html = (
            '<h2 class="svh" style="--sc:#8a6d1f">Experimental ML predictions'
            ' <span class="svsub">— opt-in, not part of the validated result</span></h2>'
            + "".join(rows))

    notes = "".join(f'<li>{_h.escape(_plainify(n))}</li>' for n in result.transcript.notes)
    notes_html = f'<ul class="notes">{notes}</ul>' if notes else ""

    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>SynAT.Vis report — {_h.escape(result.transcript.name)}</title>
<style>
:root{{--bg:#f4f5f3;--card:#fff;--ink:#14171a;--ink2:#4b5158;--muted:#8a9099;
--line:#e6e8ea;--ring:rgba(20,23,26,.07);--accent:#1a3c5e;
--good:#0ca30c;--warn:#e07b39;--bad:#d03b3b;}}
@media(prefers-color-scheme:dark){{:root{{--bg:#0e1013;--card:#171a1e;--ink:#f2f4f6;
--ink2:#c3c8ce;--muted:#8a9099;--line:#262b31;--ring:rgba(255,255,255,.06);
--accent:#3a6ea5;}}}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--ink);
font-family:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;line-height:1.55;
-webkit-font-smoothing:antialiased}}
.wrap{{max-width:760px;margin:0 auto;padding:0 18px 64px}}
header{{background:var(--accent);color:#fff;margin:0 -18px 0;padding:22px 26px}}
.brand{{font-weight:800;font-size:15px;letter-spacing:.3px;opacity:.9}}
.gene{{font-size:23px;font-weight:750;margin:2px 0 3px}}
.host{{font-size:13px;opacity:.85}}
.hero{{display:flex;gap:14px;align-items:center;background:var(--card);
border:1px solid var(--ring);border-radius:14px;padding:16px 18px;margin:18px 0 6px;
box-shadow:0 1px 2px var(--ring)}}
.hero .ic{{font-size:26px;line-height:1}}
.hero.good{{--h:var(--good)}} .hero.warn{{--h:var(--warn)}} .hero.bad{{--h:var(--bad)}}
.hero .ic{{color:var(--h)}}
.hero .hl{{font-size:16px;font-weight:650}}
.pills{{display:flex;flex-wrap:wrap;gap:8px;margin:12px 0 4px}}
.pill{{display:inline-flex;align-items:center;gap:6px;font-size:12.5px;color:var(--ink2);
background:var(--card);border:1px solid var(--ring);border-radius:999px;padding:4px 11px}}
.pill b{{color:var(--pc);font-size:14px}}
.help{{font-size:12.5px;color:var(--muted);margin:6px 2px 20px}}
.svh{{font-size:13px;text-transform:uppercase;letter-spacing:.6px;color:var(--sc);
margin:26px 0 10px;font-weight:750;display:flex;align-items:center;gap:8px}}
.svh .svsub{{text-transform:none;letter-spacing:0;color:var(--muted);font-weight:500;font-size:12px}}
.svh .cnt{{margin-left:auto;background:var(--sc);color:#fff;border-radius:999px;
font-size:11px;padding:1px 9px}}
.card,.pcard{{background:var(--card);border:1px solid var(--ring);border-left:4px solid var(--cc,#8a6d1f);
border-radius:12px;padding:13px 16px;margin:9px 0;box-shadow:0 1px 2px var(--ring)}}
.ctop{{display:flex;align-items:center;gap:8px;flex-wrap:wrap}}
.ctitle{{font-weight:680;font-size:15px}}
.tag{{font-size:10.5px;font-family:ui-monospace,Consolas,monospace;color:var(--ink2);
background:var(--bg);border:1px solid var(--line);border-radius:5px;padding:1px 6px}}
.tag.heur{{color:#8a6d1f;border-color:#e3d4a1}}
.tag.exp{{color:#8a6d1f;border-color:#e3d4a1}}
.loc{{font-size:12px;color:var(--muted);margin:3px 0 6px}}
.msg{{margin:0 0 8px;font-size:14px;color:var(--ink2)}}
.fix{{background:color-mix(in srgb,var(--good) 8%,transparent);border:1px solid color-mix(in srgb,var(--good) 25%,transparent);
border-radius:8px;padding:8px 11px;font-size:13px}}
.fix span{{display:block;font-weight:650;color:var(--good);font-size:11px;
text-transform:uppercase;letter-spacing:.5px;margin-bottom:3px}}
.fix code{{font-family:ui-monospace,Consolas,monospace;font-size:12px;color:var(--ink);word-break:break-word}}
.fix.none{{background:var(--bg);border-color:var(--line);color:var(--muted);font-style:italic}}
.ev{{font-size:11.5px;color:var(--muted);margin-top:7px;border-top:1px dashed var(--line);padding-top:6px}}
.notes{{font-size:12.5px;color:var(--muted);margin:8px 0}}
footer{{margin-top:34px;padding-top:14px;border-top:1px solid var(--line);
font-size:12px;color:var(--muted)}}
footer b{{color:var(--ink2)}}
.expcard{{background:var(--card);border:1px solid var(--ring);border-radius:14px;
padding:16px 18px;margin:8px 0 4px;box-shadow:0 1px 2px var(--ring)}}
.gauge{{display:flex;flex-direction:column;gap:6px;margin-bottom:12px}}
.gnum{{font-size:34px;font-weight:800;line-height:1}}
.gnum small{{font-size:15px;color:var(--muted);font-weight:600}}
.gband{{font-size:13px;font-weight:650}}
.gbar{{height:9px;background:var(--bg);border:1px solid var(--ring);border-radius:999px;overflow:hidden}}
.gfill{{height:100%;border-radius:999px}}
.axrow{{display:flex;align-items:center;gap:10px;margin:5px 0;font-size:12.5px}}
.axl{{flex:0 0 46%;color:var(--ink2)}}
.axtrack{{flex:1;height:8px;background:var(--bg);border:1px solid var(--ring);border-radius:999px;overflow:hidden}}
.axfill{{display:block;height:100%;background:var(--accent);border-radius:999px}}
.axp{{flex:0 0 40px;text-align:right;color:var(--muted);font-variant-numeric:tabular-nums;font-size:11.5px}}
.axrow.off .axl{{color:var(--muted)}}
.hazrow{{display:flex;flex-wrap:wrap;gap:6px;margin:10px 0 2px}}
.hchip{{font-size:11px;color:#b23b3b;background:color-mix(in srgb,#d03b3b 8%,transparent);
border:1px solid color-mix(in srgb,#d03b3b 25%,transparent);border-radius:6px;padding:2px 8px}}
.landwrap{{margin-top:12px}}
.landlbl{{font-size:11px;color:var(--muted);margin-bottom:3px}}
.land{{width:100%;height:66px;display:block}}
.lvpipe{{display:grid;grid-template-columns:repeat(6,1fr);gap:8px;margin:6px 0 4px}}
.lvtile{{background:var(--card);border:1px solid var(--ring);border-top:4px solid var(--lc);
border-radius:10px;padding:9px 10px;box-shadow:0 1px 2px var(--ring)}}
.lvtop{{display:flex;align-items:center;gap:6px}}
.lvnum{{font-size:10px;font-weight:700;color:var(--muted);background:var(--bg);border:1px solid var(--line);
border-radius:50%;width:16px;height:16px;display:inline-flex;align-items:center;justify-content:center}}
.lvic{{font-size:15px;font-weight:800;line-height:1}}
.lvscore{{margin-left:auto;font-size:11px;font-weight:700;color:var(--muted);font-variant-numeric:tabular-nums}}
.lvname{{font-size:11.5px;font-weight:700;margin-top:5px;color:var(--ink)}}
.lvdesc{{font-size:9.5px;color:var(--muted);line-height:1.35;margin-top:2px}}
.lvissues{{margin:6px 0 0;padding-left:14px;font-size:9.5px;color:var(--warn)}}
.lvok{{margin-top:6px;font-size:9.5px;color:var(--good);font-weight:600}}
@media(max-width:640px){{.lvpipe{{grid-template-columns:repeat(2,1fr)}}}}
.fchips{{display:flex;flex-wrap:wrap;gap:8px}}
.fchip{{display:inline-flex;flex-direction:column;gap:1px;background:var(--bg);
border:1px solid var(--ring);border-radius:9px;padding:6px 11px;min-width:88px}}
.fchip b{{font-size:13.5px;color:var(--ink)}}
.fchip span{{font-size:10.5px;color:var(--muted);text-transform:uppercase;letter-spacing:.4px}}
.fnotes{{font-size:12px;color:var(--ink2);margin:10px 0 2px;padding-left:18px}}
.celliframe{{width:100%;height:900px;border:1px solid var(--ring);border-radius:14px;
background:var(--card);margin:2px 0;transition:height .15s ease}}
</style></head><body><div class="wrap">
<header><div class="brand">SynAT.Vis</div>
<div class="gene">{_h.escape(result.transcript.name)}</div>
<div class="host">{_h.escape(str(meta.get('host','?')))} · {_h.escape(str(meta.get('compartment','?')))} compartment
· {len(result.transcript.utr5)}/{len(result.transcript.cds)}/{len(result.transcript.utr3)} nt (5'UTR/CDS/3'UTR)</div></header>
<div class="hero {hero_class}"><span class="ic">{hero_icon}</span><span class="hl">{_h.escape(headline)}</span></div>
<div class="pills">{pills}</div>
<div class="help">Read the <b>Serious</b> items first. Each finding shows where it is,
what it means in plain terms, and a suggested fix.</div>
{notes_html}
{levels_html}
{exp_html}
{cell_html}
{fate_html}
{''.join(cards) if cards else '<p class="msg">Nothing flagged.</p>'}
{plugin_html}
<footer>Transcript-level checks only · cannot see protein stability · no expression score ·
<i>{_h.escape(str(meta.get('host','host')))}</i> {_h.escape(str(meta.get('compartment','')))} only.
A clean report is <b>not</b> a promise of expression. Generated by SynAT.Vis.</footer>
{autosize_js}
</div></body></html>"""
