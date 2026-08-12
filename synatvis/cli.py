"""Command-line interface (CLAUDE.md §8).

Examples
--------
  python -m synatvis scan cassette.fasta --cds 51:801
  python -m synatvis scan cassette.gb --json
  python -m synatvis profiles
  python -m synatvis validate leg1 --fasta synatvis/data/native_cr_cds.fasta
"""
from __future__ import annotations

import argparse
import sys
from typing import List, Optional, Tuple

from . import __version__
from .profiles import available, load_profile
from .report import render_json, render_plain, render_text
from .scanner import scan
from .seqio import read_record


def _parse_cds(spec: Optional[str]) -> Optional[Tuple[int, int]]:
    if not spec:
        return None
    sep = ":" if ":" in spec else "-"
    a, b = spec.split(sep)
    return int(a), int(b)


def _cmd_scan(args) -> int:
    cds_span = _parse_cds(args.cds)
    tx = read_record(args.input, cds_span=cds_span)
    result = scan(tx, profile=args.profile,
                  only=args.only.split(",") if args.only else None,
                  exclude=args.exclude.split(",") if args.exclude else None)
    plugin_res = None
    if getattr(args, "plugins", False):
        from .plugins import run_available
        plugin_res = run_available(tx)
    expr = None
    if (getattr(args, "expression", False) or getattr(args, "html", False)
            or getattr(args, "four_d", False) or getattr(args, "cell", False)):
        from .expression import predict_expression
        expr = predict_expression(tx, result.profile, scan_result=result)

    if getattr(args, "cell", False):
        import os as _os
        from .journey import build_journey
        from .viz_cell import render_cell_document
        journey = build_journey(tx, result.profile, result, expr)
        out = args.out or (_os.path.splitext(args.input)[0] + "_cell.html")
        with open(out, "w", encoding="utf-8") as fh:
            fh.write(render_cell_document(journey))
        print(f"Cell-journey view written: {out}")
    elif getattr(args, "four_d", False):
        import os as _os
        from .viz4d import render_4d_html
        out = args.out or (_os.path.splitext(args.input)[0] + "_4d.html")
        with open(out, "w", encoding="utf-8") as fh:
            fh.write(render_4d_html(tx.name, expr.landscape, expr.epi, expr.band))
        print(f"4D translation view written: {out}")
    elif getattr(args, "html", False):
        import os as _os
        from .report import render_html
        out = args.out or (_os.path.splitext(args.input)[0] + "_report.html")
        with open(out, "w", encoding="utf-8") as fh:
            fh.write(render_html(result, plugin_res, expression=expr))
        print(f"HTML report written: {out}")
    elif args.json:
        print(render_json(result))
    elif args.plain:
        print(render_plain(result))
    else:
        print(render_text(result))
        if getattr(args, "expression", False) and expr is not None:
            from .journey import build_journey
            from .levels import summarize_levels, overall_verdict
            from .report import render_expression_text, render_levels_text
            _lv = summarize_levels(build_journey(tx, result.profile, result, expr))
            print(render_levels_text(_lv, overall_verdict(_lv)))
            print(render_expression_text(expr))
        if plugin_res is not None:
            from .report import render_plugins_text
            print(render_plugins_text(plugin_res))
    high = result.counts()["high"]
    return 1 if (args.fail_on_high and high) else 0


def _cmd_plugins(args) -> int:
    from .plugins import status
    print("Tier-B ML plugins (opt-in, experimental — never part of the validated core):\n")
    for p in status():
        mark = "AVAILABLE" if p["available"] else "not installed"
        print(f"  [{mark:<13}] {p['name']:<13} {p['description']}")
        if not p["available"] and p["hint"]:
            print(f"                    enable: {p['hint']}")
        if p["citation"]:
            print(f"                    ref   : {p['citation']}")
    return 0


def _cmd_profiles(args) -> int:
    for name, path in available().items():
        try:
            p = load_profile(name)
            meta = p.get("meta", {})
            status = "validated" if meta.get("validated") else "STUB/unvalidated"
            print(f"{name:<16} {meta.get('host','?')} / {meta.get('compartment','?')}"
                  f"  [{status}]")
        except Exception as exc:  # pragma: no cover
            print(f"{name:<16} <error: {exc}>")
    return 0


def _cmd_validate(args) -> int:
    from .validation import (leg1_specificity, leg2_injection, leg2_sweep,
                             leg3_cases, crossspecies, purification_scan, complexome_scan,
                             calibration)

    if args.leg == "leg1":
        return leg1_specificity.main(args.fasta, args.profile)
    if args.leg == "leg2":
        return leg2_injection.main(args.profile)
    if args.leg == "leg2sweep":
        return leg2_sweep.main(args.profile)
    if args.leg == "leg3":
        return leg3_cases.main(args.cases, args.profile)
    if args.leg == "crossspecies":
        return crossspecies.main(args.profile)
    if args.leg == "purification":
        return purification_scan.main(args.profile)
    if args.leg == "complexome":
        return complexome_scan.main()
    if args.leg == "calibration":
        return calibration.main(args.profile)
    print("unknown leg", file=sys.stderr)
    return 2


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="synatvis", description=__doc__)
    p.add_argument("--version", action="version", version=f"SynAT.Vis {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("scan", help="scan a cassette for transcript-level red flags")
    s.add_argument("input", help="FASTA or GenBank/SnapGene file")
    s.add_argument("--profile", default="cr_nuclear")
    s.add_argument("--cds", help="CDS span in the full sequence, e.g. 51:801 (0-based, half-open)")
    s.add_argument("--only", help="comma-separated module names to run")
    s.add_argument("--exclude", help="comma-separated module names to skip")
    s.add_argument("--json", action="store_true", help="emit JSON")
    s.add_argument("--plain", action="store_true",
                   help="plain-language report for non-programmers")
    s.add_argument("--expression", action="store_true",
                   help="add the expression-propensity prediction (opt-in model; always in --html)")
    s.add_argument("--plugins", action="store_true",
                   help="also run any installed Tier-B ML plugins (experimental)")
    s.add_argument("--html", action="store_true",
                   help="write a polished HTML report instead of text")
    s.add_argument("--4d", dest="four_d", action="store_true",
                   help="write an animated 4D translation-time view (HTML)")
    s.add_argument("--cell", action="store_true",
                   help="write the animated cell-journey view (transcription -> localisation, HTML)")
    s.add_argument("--out", help="output path for --html (default: <input>_report.html)")
    s.add_argument("--fail-on-high", action="store_true",
                   help="exit non-zero if any HIGH flag is present")
    s.set_defaults(func=_cmd_scan)

    pr = sub.add_parser("profiles", help="list available host profiles")
    pr.set_defaults(func=_cmd_profiles)

    pl = sub.add_parser("plugins", help="list Tier-B ML plugins and how to enable them")
    pl.set_defaults(func=_cmd_plugins)

    v = sub.add_parser("validate", help="run a validation leg (CLAUDE.md §7)")
    v.add_argument("leg", choices=["leg1", "leg2", "leg2sweep", "leg3", "crossspecies",
                                   "purification", "complexome", "calibration"])
    v.add_argument("--profile", default="cr_nuclear")
    v.add_argument("--fasta", help="native Cr CDS FASTA (leg1)")
    v.add_argument("--cases", help="cases.yaml (leg3)")
    v.set_defaults(func=_cmd_validate)
    return p


def _force_utf8_stdout() -> None:
    """Avoid UnicodeEncodeError on cp1252 Windows consoles (report uses '—', etc.)."""
    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is not None:
            try:
                reconfigure(encoding="utf-8")
            except Exception:
                pass


def main(argv: Optional[List[str]] = None) -> int:
    _force_utf8_stdout()
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
