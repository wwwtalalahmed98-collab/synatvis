"""Roll the journey checkpoints up into the six canonical molecular-biology levels.

The protein is made through six stages; this module gives each one a single
pass / warn / fail check-mark, computed as the worst status among the real,
cited sub-checks that belong to it (no new scoring — it reuses the data-anchored
journey checkpoints, whose thresholds are calibrated on 5,000 native Cr genes).

    pre-transcription -> transcription -> post-transcription ->
    pre-translation -> translation -> post-translation
"""
from __future__ import annotations

from typing import Dict, List

# ordered; each: key, short label, one-line meaning
LEVELS = [
    ("pretranscript", "Pre-transcription", "DNA construct & chromatin/silencing context"),
    ("transcript", "Transcription", "primary transcript is made and not silenced"),
    ("posttranscript", "Post-transcription", "splicing, poly(A), stability, export"),
    ("pretranslation", "Pre-translation", "5'UTR, start codon, initiation"),
    ("translation", "Translation", "codon usage & elongation"),
    ("posttranslation", "Post-translation", "folding, PTM, localisation, recovery"),
]

_RANK = {"info": 0, "ok": 1, "warn": 2, "bad": 3}
_STATUS = {0: "pass", 1: "pass", 2: "warn", 3: "fail"}


def summarize_levels(journey: Dict) -> List[Dict]:
    """Return the six levels in order, each with a status and its evidence."""
    buckets: Dict[str, List[Dict]] = {k: [] for k, _, _ in LEVELS}
    for cp in journey.get("checkpoints", []):
        buckets.setdefault(cp.get("level", "posttranslation"), []).append(cp)

    out: List[Dict] = []
    for i, (key, name, desc) in enumerate(LEVELS):
        cps = buckets.get(key, [])
        worst = 0
        checks: List[Dict] = []
        issues: List[Dict] = []
        for cp in cps:
            for p in cp.get("params", []):
                st = p.get("status", "info")
                worst = max(worst, _RANK.get(st, 0))
                checks.append(p)
                if st in ("warn", "bad"):
                    issues.append({"label": p["label"], "value": p["value"],
                                   "status": st, "detail": p.get("detail", ""),
                                   "ref": p.get("ref", "")})
        # a level with only informational checks still passes
        status = _STATUS[worst] if checks else "pass"
        # representative 0-100 score: mean of the level's scaled checks, if any
        scaled = [p["scale"] for p in checks
                  if isinstance(p.get("scale"), (int, float))]
        score = round(sum(scaled) / len(scaled), 0) if scaled else None
        out.append({
            "key": key, "name": name, "desc": desc, "order": i + 1,
            "status": status, "n_checks": len(checks), "score": score,
            "stages": [cp["title"] for cp in cps], "issues": issues,
        })
    return out


def overall_verdict(levels: List[Dict]) -> Dict:
    """Headline: how many levels pass, and the first blocking level if any."""
    fails = [l for l in levels if l["status"] == "fail"]
    warns = [l for l in levels if l["status"] == "warn"]
    passes = sum(1 for l in levels if l["status"] == "pass")
    if fails:
        verdict, blocker = "fail", fails[0]["name"]
    elif warns:
        verdict, blocker = "warn", warns[0]["name"]
    else:
        verdict, blocker = "pass", None
    return {"verdict": verdict, "passes": passes, "warns": len(warns),
            "fails": len(fails), "blocker": blocker, "n": len(levels)}
