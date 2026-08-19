"""SynAT.Vis smoke + behaviour tests.

Runs under pytest, or standalone: ``python tests/test_synatvis.py``.
Covers: profile validity, YAML fallback, codon adaptiveness, the inverted GC
rule, TGTAA poly(A) detection, Type IIS detection + remediation, the no-score
contract, and all three validation legs including compartment gating.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from synatvis import scan  # noqa: E402
from synatvis.seqio import Transcript  # noqa: E402
from synatvis.profiles import load_profile  # noqa: E402
from synatvis.profiles._schema import validate  # noqa: E402
from synatvis.codon_tables import load_for_profile  # noqa: E402
from synatvis.remediation import synonymous_fix  # noqa: E402
from synatvis.report import render_json, render_text  # noqa: E402
from synatvis.flags import Severity  # noqa: E402
from synatvis import _yaml  # noqa: E402
from synatvis.construct_grammar import (  # noqa: E402
    CandidatePart, Verdict, evaluate_candidate, load_criteria,
    load_so_vocabulary, valid_so_ids,
)


def _table():
    p = load_profile("cr_nuclear")
    return load_for_profile(p, p["_base_dir"])


def test_yaml_fallback_parses_profile():
    p = load_profile("cr_nuclear")
    assert p["meta"]["host"].startswith("Chlamydomonas")
    assert p["meta"]["compartment"] == "nuclear"
    assert "TGTAA" in p["polya"]["nue_motifs"]
    assert p["cloning"]["enzymes"]["BsaI"] == "GGTCTC"
    assert p["composition"]["target_gc"] > 0.6  # GC-rich host


def test_yaml_parses_list_of_mappings():
    data = _yaml.safe_load(
        "cases:\n  - id: a\n    n: 1\n  - id: b\n    n: 2\n"
    )
    assert data["cases"][0]["id"] == "a"
    assert data["cases"][1]["n"] == 2


def test_both_profiles_schema_valid():
    for name in ("cr_nuclear", "cr_chloroplast"):
        ok, problems = validate(load_profile(name, validate=False))
        assert ok, f"{name}: {problems}"


def test_codon_adaptiveness_gc_biased():
    t = _table()
    # preferred GC-ending codons should have weight 1.0 in their family
    assert t.weight["GCC"] == 1.0            # Ala
    assert t.weight["CTG"] == 1.0            # Leu
    assert t.weight["GGC"] == 1.0            # Gly
    assert t.weight["GCA"] < t.weight["GCC"]  # AT-ending is less adaptive
    rca = t.rca("ATGGCCGGCCTGGAGAAGTAA")
    assert 0.0 < rca <= 1.0


def test_clean_parent_no_high_flags():
    parent = "ATG" + ("GCCGGCCTGGAGAAGGCC" * 6) + "TAA"
    res = scan(Transcript(cds=parent, name="clean"))
    highs = [f for f in res.flags if f.severity == Severity.HIGH]
    assert not highs, [f.message for f in highs]


def test_low_gc_is_the_hazard_not_high():
    # a GC-rich CDS must NOT be flagged by composition (high GC is normal)
    gc_rich = "ATG" + ("GCCGGCCTGGAGAAGGCC" * 6) + "TAA"
    res = scan(Transcript(cds=gc_rich), only=["composition"])
    assert not [f for f in res.flags if f.detail.get("hazard") == "low_gc"]
    # an AT-rich stretch DOES trip a GC trough
    at_rich = "ATG" + ("AAATTTAAATTT" * 8) + "TAA"
    res2 = scan(Transcript(cds=at_rich), only=["composition"])
    assert any(f.detail.get("hazard") == "low_gc" for f in res2.flags)


def test_polya_tgtaa_detected_in_gc_context():
    # strong (>85%) G/C context, matching the Leg-1-derived operating point
    cds = "ATGGCCGGCCTGGAG" + "TGTAA" + "GGCCGGCGGCCGCGGCCGGCGGCCGGCGGC" + "TAA"
    res = scan(Transcript(cds=cds), only=["polya"])
    assert any(f.severity >= Severity.MEDIUM for f in res.flags), \
        [f.message for f in res.flags]


def test_polya_weak_context_not_flagged_high():
    # a TGTAA in weak (non-GC-rich) context must NOT raise a medium/high flag
    cds = "ATGGCCGGCCTGGAG" + "TGTAA" + "GAGAAGATAAAGAATAAGAAT" + "GCCTAA"
    res = scan(Transcript(cds=cds), only=["polya"])
    assert not any(f.severity >= Severity.MEDIUM for f in res.flags)


def test_cloning_bsai_detected_and_remediated():
    cds = "ATGGGTCTCGCCAAGGACCAGTTCATCACCGTGCCCAACTAA"  # GGTCTC at pos 3
    res = scan(Transcript(cds=cds), only=["cloning"])
    bsa = [f for f in res.flags if f.detail.get("enzyme") == "BsaI"]
    assert bsa and bsa[0].severity == Severity.HIGH
    assert bsa[0].suggested_edit and "synonymous" in bsa[0].suggested_edit


def test_remediation_removes_motif_and_preserves_protein():
    t = _table()
    cds = "ATGGGTCTCGCCAAGTAA"  # GGTCTC (BsaI) across codons 1-2
    edit = synonymous_fix(cds, t, cds_start=3, cds_end=9,
                          avoid=["GGTCTC", "GAGACC"])
    assert edit.ok, edit.reason
    new_cds = cds[:edit.cds_start] + edit.new_window + \
        cds[edit.cds_start + len(edit.old_window):]
    assert "GGTCTC" not in new_cds
    assert t.translate(cds) == t.translate(new_cds)  # synonymous


def test_expression_model_directional_and_bounded():
    from synatvis.expression import predict_expression
    from synatvis.profiles import load_profile, load_yaml, PACKAGE_DIR
    import os
    p = load_profile("cr_nuclear")
    cases_path = os.path.join(PACKAGE_DIR, "data", "cases.yaml")
    cases = {c["id"]: c for c in load_yaml(open(cases_path, encoding="utf-8").read())["cases"]}

    def epi(cid):
        return predict_expression(Transcript(cds=cases[cid]["sequence"]), p).epi

    # bounded 0-100
    r = predict_expression(Transcript(cds="ATG" + "GCCGGCCTGGAGAAGGCC" * 6 + "TAA"), p)
    assert 0.0 <= r.epi <= 100.0 and r.landscape and r.confidence
    # directional: Cr-optimised rescue partner scores well above the native version
    assert epi("gfp_cr_codon_optimized") > epi("gfp_native_at_rich") + 20
    assert epi("luc_cr_codon_optimized") > epi("luc_native_at_rich") + 20


def test_render_html_is_self_contained():
    from synatvis.report import render_html
    res = scan(Transcript(cds="ATGGGTCTCGCCAAGGACCAGTTCATCACCGTGCCCAACTAA"))  # has BsaI
    html = render_html(res)
    assert html.lstrip().lower().startswith("<!doctype html>") and "</html>" in html
    # self-contained: no external scripts, stylesheets, images, or fetches
    for forbidden in ("<script", "<link", "src=", "http://", "https://"):
        assert forbidden not in html
    assert "SUGGESTED FIX".lower() in html.lower() or "cloning" in html


def test_plugins_status_and_core_stays_torch_free():
    import sys as _sys
    from synatvis import plugins
    names = {p["name"] for p in plugins.status()}
    assert {"lineardesign", "codonbert", "saluki", "utr_lm",
            "borzoi", "aparent2", "icodon",
            "signalp6", "deeploc2", "deeptmhmm", "pangolin",
            "codonfm", "openfold3", "evo2"} <= names
    # every command-gated frontier/validated tool must be unavailable without its ENV
    cmd_gated = ("borzoi", "aparent2", "icodon", "signalp6", "deeploc2",
                 "deeptmhmm", "pangolin", "codonfm", "openfold3", "evo2")
    assert not any(p["available"] and p["name"] in cmd_gated for p in plugins.status())
    scan(Transcript(cds="ATGGCCGGCTAA"))
    assert "torch" not in _sys.modules  # the core never imports a heavy ML dep


def test_json_command_plugin_runs_real_command(tmp_path, monkeypatch):
    """A command-gated plugin runs the user's command and never fabricates a score."""
    import sys as _sys
    from synatvis import plugins
    from synatvis.plugins.icodon import ICodonPlugin
    script = tmp_path / "fake_icodon.py"
    script.write_text("import sys,json; sys.stdin.read();"
                      "print(json.dumps({'stability':0.42,'species':'zebrafish'}))")
    monkeypatch.setenv("ICODON_CMD", f'"{_sys.executable}" "{script}"')
    p = ICodonPlugin()
    assert p.available() is True
    res = p.analyze(Transcript(cds="ATGGCCGGCGCCTAA"))
    assert res and res[0].value == 0.42 and res[0].validated is False
    monkeypatch.delenv("ICODON_CMD")
    assert p.available() is False           # dormant again without the ENV


def test_report_has_no_composite_score():
    res = scan(Transcript(cds="ATGGCCGGCTAA"))
    js = render_json(res)
    assert '"no_composite_score": true' in js
    # no numeric verdict key of any kind
    for forbidden in ('"score":', '"grade":', '"probability":', '"expression_score":'):
        assert forbidden not in js
    txt = render_text(res)
    assert "no composite expression score" in txt.lower()
    assert "proteolysis" in txt.lower()


def test_silencing_marked_heuristic():
    res = scan(Transcript(cds="ATG" + "AAATTTAAATTT" * 8 + "TAA"))
    assert res.module_meta["silencing"]["validated"] is False
    txt = render_text(res)
    if any(f.module == "silencing" for f in res.flags):
        assert "HEURISTIC" in txt


def test_structure_energy_fold_backend_agnostic():
    from synatvis.structure_energy import fold
    _, paired_hp, _ = fold("GGGGGGGGGGAAAACCCCCCCCCC")   # perfect hairpin
    _, paired_ss, _ = fold("AAAAAAAAAAAAAAAAAAAAAAAA")    # no pairs possible
    assert paired_hp > 0.5 and paired_ss < 0.1


def test_codon_panel_tables_load_and_tai_discriminates():
    from synatvis.codon_tables import attach_advanced, gene_tai
    from synatvis.profiles import PACKAGE_DIR
    p = load_profile("cr_nuclear")
    adv = attach_advanced(p, PACKAGE_DIR)
    assert adv["optimality"] and adv["pairs"] and adv["tai"]  # all three tables present
    assert adv["tai"]["CTG"] > adv["tai"]["ATA"]              # preferred > avoided
    opt_cds = "ATG" + "GCCGGCCTGGAGAAGGCC" * 6 + "TAA"        # Cr-optimal
    at_cds = "ATG" + "GCAGGATTAGAAAAAGCA" * 6 + "TAA"          # AT-ending synonyms
    assert gene_tai(opt_cds, adv["tai"]) > gene_tai(at_cds, adv["tai"])


def test_codon_panel_flags_low_tai_and_stays_quiet_on_optimal():
    at_cds = "ATG" + "GCAGGATTAGAAAAAGCA" * 6 + "TAA"
    res = scan(Transcript(cds=at_cds), only=["codon"])
    assert any("tRNA adaptation" in f.message for f in res.flags)
    opt_cds = "ATG" + "GCCGGCCTGGAGAAGGCC" * 6 + "TAA"
    res2 = scan(Transcript(cds=opt_cds), only=["codon"])
    assert not any("tRNA adaptation" in f.message for f in res2.flags)


def test_leg2_sweep_codon_curve():
    from synatvis.validation.leg2_sweep import sweep_codon
    cmin, rows = sweep_codon(load_profile("cr_nuclear"))
    d = dict(rows)
    # detection turns on exactly at the cluster_min_rare operating point
    assert d[0] == 0.0 and d[cmin] == 1.0 and d[8] == 1.0
    assert all(v == 0.0 for k, v in rows if k < cmin)


def test_leg2_injection_all_pass():
    from synatvis.validation.leg2_injection import run_injection
    res = run_injection("cr_nuclear")
    for c in res["cases"]:
        assert c["pass"], c


def test_leg1_specificity_runs():
    from synatvis.validation.leg1_specificity import run_specificity, STUB
    res = run_specificity(STUB, "cr_nuclear")
    assert res["n_sequences"] >= 1
    assert "polya" in res["fp_rate_medium_or_high"]


def test_crossspecies_discrimination():
    import os
    from synatvis.validation.crossspecies import run_crossspecies, FOREIGN
    if not os.path.isfile(FOREIGN):
        return  # foreign set not shipped in this build
    r = run_crossspecies(n_cr=300)
    # the host-fit metric must separate Cr from foreign strongly
    assert r["auc_rca"] > 0.95
    assert r["specificity"] > 0.90
    assert r["cr_median_rca"] > r["foreign_median_rca"] + 0.2


def test_leg3_compartment_gating():
    from synatvis.validation.leg3_cases import run_cases, SEED
    res = run_cases(SEED, "cr_nuclear")
    assert res["active_compartment"] == "nuclear"
    # the chloroplast case must be gated OUT (not scored under the nuclear profile)
    assert res["strata"]["chloroplast"]["gated"] >= 1
    assert res["strata"]["chloroplast"]["tp"] == 0
    # nuclear rescue pairs should score
    nuc = res["strata"]["nuclear"]
    assert nuc["tp"] >= 1 and nuc["tn"] >= 1


def test_ptm_signal_and_localization():
    from synatvis.ptm import predict_protein_fate, translate
    # strong hydrophobic leader + N-glyc sequons -> secreted
    sec = "ATG" + "AAG" + "CTGCTGCTGCTGGCCCTGCTGGTGCTGGCCTGC" + "GCC" \
          + "AACGCCACC" * 4 + "GCCGGCGAGGCC" * 20 + "TAA"
    f = predict_protein_fate(Transcript(cds=sec))
    assert f.signal_peptide is True
    assert f.localization == "secreted"
    assert len(f.n_glyc_sites) >= 3
    # a globular cytosolic protein (GFP-like N-terminus) must NOT be called secreted
    gfp = ("ATGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGAC"
           "GTAAACGGCCACAAGTTCAGCGTG") + "GCCGGC" * 10 + "TAA"
    g = predict_protein_fate(Transcript(cds=gfp))
    assert g.signal_peptide is False
    assert g.localization == "cytosolic"


def test_ptm_translate_stops_at_stop():
    from synatvis.ptm import translate
    assert translate("ATGGCCTAAGGG") == "MA"


def test_cell_view_self_contained_and_organelles():
    import re
    from synatvis.viz_cell import render_cell_document
    journey = {"name": "demo", "epi": 78.0, "band": "STRONG", "loc": "secreted",
               "landscape": [0.2, 0.9, 0.5, 1.0], "biophys": {}, "structure": {},
               "purification": {"tags": []},
               "checkpoints": [{"key": "gene", "title": "Gene", "symbol": "dna",
                                "organelle": "nucleus",
                                "params": [{"label": "GC", "value": "66%", "scale": 66,
                                            "status": "ok", "ref": "x", "detail": ""}]}]}
    html = render_cell_document(journey)
    assert not re.search(r'https?://(?!www\.w3\.org)', html)   # self-contained
    assert html.count("{") == html.count("}")
    assert "chloroplast" in html and "nucleus" in html and "Golgi" in html


def test_journey_routes_golgi_only_when_secretory():
    from synatvis.scanner import scan as _scan
    from synatvis.expression import predict_expression
    from synatvis.journey import build_journey

    def keys_for(cds, name):
        tx = Transcript(name=name, cds=cds)
        res = _scan(tx, profile="cr_nuclear")
        expr = predict_expression(tx, res.profile, scan_result=res)
        return [c["key"] for c in build_journey(tx, res.profile, res, expr)["checkpoints"]]

    secreted = ("ATGAAGCTGCTGCTGCTGGCCCTGCTGGTGCTGGCCTGCGCC"
                + "AACGCCACC" * 3 + "GCCGGCGAGGCC" * 16 + "TAA")
    cytosolic = "ATG" + "GCCGGCGAGGCC" * 16 + "TAA"
    assert "golgi" in keys_for(secreted, "sec")
    assert "golgi" not in keys_for(cytosolic, "cyt")


def test_biophysics_reasonable():
    from synatvis.ptm import biophysics, translate
    # insulin-ish acidic vs basic control
    acidic = translate("ATG" + "GAAGAAGAAGAT" * 5 + "TAA")
    basic = translate("ATG" + "AAGAAGCGCCGC" * 5 + "TAA")
    assert biophysics(acidic)["pI"] < 6.0
    assert biophysics(basic)["pI"] > 8.0
    assert biophysics(acidic)["mw_kda"] > 0


def test_purification_detects_elp_his_and_intein():
    from synatvis.purification import predict_purification, INTEIN_MUTATIONS
    elp = "GTGCCGGGCGTGGGC" * 8          # (VPGVG)x8
    his = "CATCATCATCATCATCAT"           # His6
    tx = Transcript(cds="ATG" + his + elp + "GCCGGCGCC" * 8 + "TAA")
    pur = predict_purification(tx)
    assert any(t.startswith("His") for t in pur.tags)
    assert any(t.startswith("ELP") for t in pur.tags)
    assert "ITC" in pur.strategy and pur.elp and pur.elp["pentamers"] == 8
    assert pur.intein and len(pur.intein["mutations"]) == len(INTEIN_MUTATIONS)
    # no-tag construct falls back
    p2 = predict_purification(Transcript(cds="ATGGCCGGCGCCTAA"))
    assert not p2.tags


def test_purification_recommends_by_biophysics():
    from synatvis.purification import predict_purification
    acidic = predict_purification(Transcript(cds="ATG" + "GAAGAAGATGAAGAAGATGAA" * 10 + "TAA"))
    basic = predict_purification(Transcript(cds="ATG" + "AAGAAGCGCAAGCGCAAG" * 10 + "TAA"))

    def cap(p):
        return next(w["method"] for w in p.workflow if w["phase"].startswith("capture"))

    assert "Anion" in cap(acidic)          # low pI -> Q
    assert "Cation" in cap(basic)           # high pI -> S
    # a complete plan offers both chromatographic and non-chromatographic options + a polish step
    modes = {r["mode"] for r in acidic.recommendations}
    assert "chromatographic" in modes and "non-chromatographic" in modes
    assert any(w["phase"] == "polish" for w in acidic.workflow)


def test_purification_layer_confirms_over_corpus():
    from synatvis.validation.purification_scan import run
    r = run(limit=150)
    assert r["n"] >= 120 and r["errors"] == 0            # runs on every CDS, no crashes
    assert sum(r["capture"].values()) == r["n"]          # every gene gets a capture method


def test_journey_has_purification_plan_checkpoint():
    from synatvis.scanner import scan as _scan
    from synatvis.expression import predict_expression
    from synatvis.journey import build_journey
    tx = Transcript(name="d", cds="ATG" + "GCCGGCGAGGCC" * 14 + "TAA")
    res = _scan(tx, profile="cr_nuclear")
    expr = predict_expression(tx, res.profile, scan_result=res)
    J = build_journey(tx, res.profile, res, expr)
    plan = [c for c in J["checkpoints"] if c["key"] == "purification_plan"]
    assert plan and any(p["label"] == "Capture" for p in plan[0]["params"])


def test_build_journey_checkpoints_and_params():
    from synatvis.scanner import scan as _scan
    from synatvis.expression import predict_expression
    from synatvis.journey import build_journey
    elp = "GTGCCGGGCGTGGGC" * 8
    his = "CATCATCATCATCATCAT"
    leader = "ATGAAGCTGCTGCTGCTGGCCCTGCTGGTGCTGGCCTGCGCC"
    tx = Transcript(name="elp demo", cds=leader + his + elp + "GCCGGCGAGGCC" * 12 + "TAA")
    res = _scan(tx, profile="cr_nuclear")
    expr = predict_expression(tx, res.profile, scan_result=res)
    J = build_journey(tx, res.profile, res, expr)
    keys = [c["key"] for c in J["checkpoints"]]
    for expected in ("gene", "initiation", "translation", "folding", "elp", "intein", "pure"):
        assert expected in keys, expected
    # every checkpoint carries at least one parameter, a status, and a level
    for c in J["checkpoints"]:
        assert c["params"]
        assert c["level"] in ("pretranscript", "transcript", "posttranscript",
                              "pretranslation", "translation", "posttranslation")
        for p in c["params"]:
            assert p["status"] in ("ok", "warn", "bad", "info")
    assert J["structure"]["plain_summary"]           # self-explanatory text present


def test_six_levels_summary_and_verdict():
    from synatvis.scanner import scan as _scan
    from synatvis.expression import predict_expression
    from synatvis.journey import build_journey
    from synatvis.levels import summarize_levels, overall_verdict, LEVELS
    # premature poly(A) inside the CDS -> must fail at the post-transcription level
    tx = Transcript(name="polya", cds="ATGAAGCTGCTGGCCCTGCTGTGTAAGGCC"
                    + "GGCGAGGCC" * 6 + "TAA")
    res = _scan(tx, profile="cr_nuclear")
    expr = predict_expression(tx, res.profile, scan_result=res)
    levels = summarize_levels(build_journey(tx, res.profile, res, expr))
    assert len(levels) == 6
    assert [l["key"] for l in levels] == [k for k, _, _ in LEVELS]  # ordered
    for l in levels:
        assert l["status"] in ("pass", "warn", "fail")
    post = next(l for l in levels if l["key"] == "posttranscript")
    assert post["status"] == "fail"
    assert any("poly(A)" in i["label"] for i in post["issues"])
    v = overall_verdict(levels)
    assert v["fails"] >= 1 and v["blocker"]


def test_cell_view_renders_from_journey():
    import re
    from synatvis.scanner import scan as _scan
    from synatvis.expression import predict_expression
    from synatvis.journey import build_journey
    from synatvis.viz_cell import render_cell_document
    tx = Transcript(name="d", cds="ATG" + "GCCGGCGAGGCC" * 12 + "TAA")
    res = _scan(tx, profile="cr_nuclear")
    expr = predict_expression(tx, res.profile, scan_result=res)
    html = render_cell_document(build_journey(tx, res.profile, res, expr))
    assert not re.search(r'https?://(?!www\.w3\.org)', html)   # self-contained
    assert html.count("{") == html.count("}")
    assert "PURIFICATION BENCH" in html and "checkpoints" in html


def test_complexome_matches_known_gene_names():
    from synatvis.complexome import identify_complexes
    assert [m.complex_name for m in identify_complexes("rbcL")] == ["Rubisco"]
    assert [m.complex_name for m in identify_complexes("RBCS2")] == ["Rubisco"]
    assert [m.complex_name for m in identify_complexes("psbA")] == ["Photosystem II"]
    assert [m.complex_name for m in identify_complexes("RPL3")] == ["Cytosolic ribosome (80S)"]
    assert identify_complexes("some_unrelated_gene_123") == []
    assert identify_complexes("") == []


def test_complexome_checkpoint_appears_only_on_a_real_match():
    from synatvis.scanner import scan as _scan
    from synatvis.expression import predict_expression
    from synatvis.journey import build_journey
    cds = "ATG" + "GCCGGCGAGGCC" * 12 + "TAA"
    tx_match = Transcript(name="rbcL", cds=cds)
    tx_nomatch = Transcript(name="unrelated_transcript", cds=cds)
    res_m = _scan(tx_match, profile="cr_nuclear")
    expr_m = predict_expression(tx_match, res_m.profile, scan_result=res_m)
    J_m = build_journey(tx_match, res_m.profile, res_m, expr_m)
    keys_m = [c["key"] for c in J_m["checkpoints"]]
    assert "complex" in keys_m
    complex_cp = next(c for c in J_m["checkpoints"] if c["key"] == "complex")
    assert complex_cp["params"][0]["label"] == "Rubisco"

    res_n = _scan(tx_nomatch, profile="cr_nuclear")
    expr_n = predict_expression(tx_nomatch, res_n.profile, scan_result=res_n)
    J_n = build_journey(tx_nomatch, res_n.profile, res_n, expr_n)
    keys_n = [c["key"] for c in J_n["checkpoints"]]
    assert "complex" not in keys_n


def test_moleculardynamics_plugin_runs_users_own_command_and_never_fabricates(tmp_path, monkeypatch):
    import sys as _sys
    from synatvis.plugins.moleculardynamics import MolecularDynamicsPlugin
    script = tmp_path / "fake_md.py"
    script.write_text("import sys,json; sys.stdin.read();"
                      "print(json.dumps({'rmsd_nm':0.21,'radius_of_gyration_nm':1.4,"
                      "'sim_time_ns':50,'force_field':'AMBER99SB-ILDN'}))")
    monkeypatch.setenv("MDSIM_CMD", f'"{_sys.executable}" "{script}"')
    p = MolecularDynamicsPlugin()
    assert p.available() is True
    res = p.analyze(Transcript(cds="ATGGCCGGCGCCTAA"))
    assert any(r.label == "MD backbone RMSD" and r.value == 0.21 and r.validated is False
              for r in res)
    monkeypatch.delenv("MDSIM_CMD")
    assert p.available() is False  # no MD engine configured -> plugin is simply absent, nothing faked


def test_complexome_scan_real_corpus_specificity_and_true_positives():
    from synatvis.validation.complexome_scan import run
    r = run()
    assert r["n_native"] == 5000
    assert r["false_positives"] == 0
    assert r["n_named"] == 13
    assert r["true_positives"] == 13


def test_moclo_fetcher_parses_real_archive_if_present():
    """Offline test: if the corpus has been fetched locally, sanity-check it.

    Skips (rather than fails) when the git-ignored corpus is absent, so the suite
    stays runnable on a clean checkout with no network.
    """
    import json as _json
    import importlib.util as _ilu
    from synatvis.profiles import PACKAGE_DIR
    cg = os.path.join(PACKAGE_DIR, "data", "construct_grammar")
    # the data dir is not a package, so load the fetcher by path
    spec = _ilu.spec_from_file_location("fetch_moclo_corpus",
                                        os.path.join(cg, "fetch_moclo_corpus.py"))
    F = _ilu.module_from_spec(spec)
    spec.loader.exec_module(F)
    manifest = os.path.join(cg, "moclo_corpus", "manifest.json")
    if not os.path.isfile(manifest):
        return  # corpus not fetched on this machine; nothing to check
    with open(manifest, encoding="utf-8") as fh:
        data = _json.load(fh)
    assert data["sha256"] == F.EXPECTED_SHA256
    assert data["n_records"] == F.EXPECTED_RECORDS
    recs = data["records"]
    # every real MoClo part must carry a Type IIS site -- this is IC-1
    assert all(r["type_iis_sites"] for r in recs)
    assert all(r["length_bp"] > 0 for r in recs)


def test_expanded_corpus_is_tiered_and_synthetic_labels_are_exact():
    """Offline: if the expanded corpus was built locally, check its integrity.

    The synthetic tier's junction labels are exact by construction, so they must
    tile each sequence with no gaps or overlaps. A failure here means generated
    training labels are wrong -- worse than having no synthetic data at all.
    """
    import json as _json
    from synatvis.profiles import PACKAGE_DIR
    corp = os.path.join(PACKAGE_DIR, "data", "construct_grammar", "moclo_corpus")
    man_path = os.path.join(corp, "corpus_manifest.json")
    if not os.path.isfile(man_path):
        return  # corpus not built on this machine
    with open(man_path, encoding="utf-8") as fh:
        man = _json.load(fh)

    tiers = man["by_tier"]
    assert "cr_primary" in tiers and tiers["cr_primary"] > 0
    # tiers must stay separate and honestly labelled
    assert set(tiers) <= {"cr_primary", "syntax_only", "synthetic"}
    assert sum(tiers.values()) == man["n_total"]

    ids = [r["id"] for r in man["records"]]
    assert len(ids) == len(set(ids)), "duplicate record ids"

    # only Chlamydomonas records may carry the cr_primary tier
    for r in man["records"]:
        if r["tier"] == "cr_primary":
            assert "hlamydomonas" in r["host"]
        if r["tier"] == "syntax_only":
            assert "hlamydomonas" not in r["host"]

    seqs, cur = {}, None
    with open(os.path.join(corp, "corpus.fasta"), encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
            elif cur:
                seqs[cur].append(line)
    seqs = {k: "".join(v) for k, v in seqs.items()}
    assert len(seqs) == man["n_total"]

    for r in man["records"]:
        if r["tier"] != "synthetic":
            continue
        s = seqs[r["id"]]
        assert len(s) == r["length_bp"]
        js = r["junctions"]
        assert js[0]["start"] == 0 and js[-1]["end"] == len(s)
        for a, b in zip(js, js[1:]):
            assert a["end"] == b["start"], "junction spans must tile exactly"


def test_splice_deliberate_intron_specificity_on_native_genes():
    """Guard the measured specificity win: an absolute GC>=0.60 rule flagged 41.3%
    of real intron-less native Cr genes as containing a 'deliberate intron'. The
    measured rule (length + GC relative to flanks + polypyrimidine tract) cut that
    to 9.1%. This asserts we never drift back toward the old behaviour.
    """
    from synatvis.seqio import Transcript, read_fasta
    from synatvis.scanner import scan as _scan
    from synatvis.flags import Severity
    from synatvis.profiles import PACKAGE_DIR
    path = os.path.join(PACKAGE_DIR, "data", "native_cr_cds.fasta")
    if not os.path.isfile(path):
        return
    n = hit = 0
    for name, seq in read_fasta(path):
        seq = seq.strip()
        if len(seq) < 200:
            continue
        n += 1
        res = _scan(Transcript(cds=seq, name=name), profile="cr_nuclear")
        if any(f.module == "splice" and f.severity == Severity.INFO for f in res.flags):
            hit += 1
        if n >= 200:
            break
    assert n > 0
    rate = hit / n
    # measured 9.1% over 992 genes; 20% leaves headroom for sampling noise while
    # still failing loudly if the old ~41% behaviour ever returns
    assert rate < 0.20, f"deliberate-intron false-positive rate regressed to {rate:.1%}"


def test_splice_does_not_flag_a_clean_cr_optimised_cds():
    """The Cr-codon-optimised GFP literature case contains no intron. It used to
    draw a false 'deliberate intron' call, which is what made the intron-masking
    fix delete real coding sequence."""
    import re as _re
    from synatvis.seqio import Transcript
    from synatvis.scanner import scan as _scan
    from synatvis.flags import Severity
    from synatvis.profiles import PACKAGE_DIR
    cases = open(os.path.join(PACKAGE_DIR, "data", "cases.yaml"), encoding="utf-8").read()
    m = _re.search(r"id: gfp_cr_codon_optimized.*?sequence: \"([ACGT]+)\"", cases, _re.S)
    assert m
    res = _scan(Transcript(cds=m.group(1), name="gfp"), profile="cr_nuclear")
    calls = [f for f in res.flags
             if f.module == "splice" and f.severity == Severity.INFO]
    assert not calls, f"false deliberate-intron call on an intron-less CDS: {calls}"


def test_fetch_structure_reads_gene_from_fasta_header(tmp_path):
    """Offline check of the structure fetcher's parsing. The network paths
    (UniProt / AlphaFold) are deliberately not tested here -- a unit test must not
    depend on two external services being reachable."""
    from synatvis.fetch_structure import gene_from_fasta
    p = tmp_path / "g.fasta"
    p.write_text(">RBCS2 some description here\nATGGCC\n")
    assert gene_from_fasta(str(p)) == "RBCS2"
    p2 = tmp_path / "empty.fasta"
    p2.write_text("ATGGCC\n")
    assert gene_from_fasta(str(p2)) == ""


def test_algae_product_catalogue_loads_and_is_well_formed():
    from synatvis.algae_products import load_catalogue, catalogue_stats
    cat = load_catalogue()
    assert len(cat) >= 25, "catalogue failed to parse (stdlib YAML fallback is strict)"
    for e in cat:
        assert e.get("product") and e.get("gene")
        assert e.get("origin") in ("transgenic", "native")
        assert e.get("confidence") in ("verified_sequence", "reported", "indicative")
    s = catalogue_stats()
    assert s["by_origin"].get("transgenic", 0) > 0
    assert s["by_origin"].get("native", 0) > 0


def test_algae_product_name_and_similarity_matching():
    """Name matching must be exact; similarity must link the SAME protein encoded by
    DIFFERENT DNA, and must stay silent on unrelated sequence."""
    import re as _re
    from synatvis.algae_products import identify
    from synatvis.profiles import PACKAGE_DIR

    named = identify(name="RBCS2")
    assert named and named[0].match_type == "name"

    cases = open(os.path.join(PACKAGE_DIR, "data", "cases.yaml"), encoding="utf-8").read()
    native_gfp = _re.search(r'id: gfp_native_at_rich.*?sequence: "([ACGT]+)"',
                            cases, _re.S).group(1)
    # native AT-rich GFP shares no meaningful DNA identity with the Cr-optimised
    # reference, but encodes the same protein -- similarity must still find it
    sim = [h for h in identify(name="unnamed_query", cds=native_gfp)
           if h.match_type == "similarity"]
    assert sim, "similarity route failed to link the same protein from different DNA"
    assert sim[0].similarity >= 0.6

    rnd = "ATG" + "GCTAGCACCTGACTGATCGATCGTACGATCAGCTAGCATCGATCGATCGTAGCTAGCTA" * 4 + "TAA"
    assert not identify(name="definitely_not_catalogued", cds=rnd)


def test_journey_reports_known_algal_product():
    from synatvis.scanner import scan as _scan
    from synatvis.expression import predict_expression
    from synatvis.journey import build_journey
    tx = Transcript(name="RBCS2", cds="ATG" + "GCCGGCGAGGCC" * 12 + "TAA")
    res = _scan(tx, profile="cr_nuclear")
    expr = predict_expression(tx, res.profile, scan_result=res, run_ml=False)
    J = build_journey(tx, res.profile, res, expr)
    assert "algae_prior_art" in [c["key"] for c in J["checkpoints"]]


def test_calibration_anchors_are_real_and_cited():
    from synatvis.validation.calibration import load_anchors
    a = load_anchors()
    assert len(a["anchors"]) >= 3, "expected at least three anchor sets"
    axes = set()
    for anc in a["anchors"]:
        assert anc["citation"], "an anchor without a citation is not usable"
        assert anc["measurements"], "an anchor must carry real measured values"
        ms = anc["measurements"]
        if any("introns" in m for m in ms):
            # intron-axis schema: fold activity keyed on intron count
            axes.add("intron")
            by_n = {}
            for m in ms:
                v = m.get("relative_protein_activity")
                if v is not None:      # some entries are ordering evidence only
                    by_n.setdefault(m["introns"], []).append(v)
            assert by_n.get(0) and min(by_n[0]) == 1.0, "baseline must be 1.0"
            assert min(by_n[max(by_n)]) > 1.0, "more introns must measure above baseline"
        else:
            # codon-axis schema: ranked variants with measured RCA
            axes.add("codon")
            ranks = [m["rank"] for m in ms]
            rcas = [m["rca_percent"] for m in ms]
            assert len(set(ranks)) == len(ranks), "variant ranks must be distinct"
            # rank order must agree with RCA order -- that is the paper's finding
            paired = sorted(zip(ranks, rcas))
            assert [r for _, r in paired] == sorted(rcas), \
                "measured rank must be monotonic in RCA"
    assert len(axes) >= 2, "anchors must span at least two independent axes"
    assert "known_defect_exposed" in a
    # the honesty guard: the word "calibrated" must never be claimed outright
    assert "NOT enough to call the index" in a["status"] or "never \"calibrated\"" in a["status"]


def test_calibration_leg_runs_and_reports_the_known_ordering_defect():
    """The leg must actually detect the documented defect, not quietly pass.

    If this starts failing because there are no violations, the expression model
    was changed -- update calibration_anchors.yaml's known_defect_exposed block
    rather than deleting this test.
    """
    from synatvis.validation.calibration import run
    r = run()
    if not r["available"]:
        return  # corpus not fetched on this machine
    assert set(r["scores"]) == {0, 1, 2}
    assert r["scores"][1] >= r["scores"][0], "1 intron should not score below baseline"
    # the two-intron ordering defect is real and currently expected
    assert r["violations"], "expected the documented two-intron ordering violation"


def test_so_vocabulary_loads_and_contains_real_verified_terms():
    vocab = load_so_vocabulary()
    ids = valid_so_ids()
    assert len(vocab["terms"]) >= 10
    # spot-check a few terms independently verified against the real Sequence
    # Ontology (EBI OLS) rather than guessed
    for so_id, label in [("SO:0000167", "promoter"), ("SO:0000316", "CDS"),
                         ("SO:0000141", "terminator"), ("SO:0000188", "intron")]:
        assert so_id in ids
        entry = next(t for t in vocab["terms"] if t["id"] == so_id)
        assert entry["label"] == label


def test_construct_grammar_unknown_so_term_fails_ic4_even_if_real_looking():
    """A plausible-looking 'SO:XXXXXXX' string that just isn't in the adopted
    vocabulary yet must fail IC-4 -- IC-4 checks vocabulary membership, not
    merely that some string was supplied."""
    part = CandidatePart(
        part_id="mystery_part",
        so_term="SO:9999999",  # not a real adopted term
        assembly_sites=["BsaI"],
        validated_hosts=["chlamydomonas_reinhardtii_nuclear"],
        accession="Addgene_11111",
        accession_type="addgene_plasmid",
        sequence="ATGGCC" * 20,
        citation="doi:10.1000/xyz",
        functional_evidence="confirmed",
    )
    result = evaluate_candidate(part)
    assert result.verdict is Verdict.EXCLUDE
    assert "IC-4" in {c.id for c in result.failed()}


def test_construct_grammar_criteria_yaml_parses_and_matches_ids():
    data = load_criteria()
    ic_ids = {c["id"] for c in data["inclusion"]}
    ex_ids = {c["id"] for c in data["exclusion"]}
    assert ic_ids == {"IC-1", "IC-2", "IC-3", "IC-4", "IC-5"}
    assert ex_ids == {"EX-1", "EX-2", "EX-3", "EX-4"}
    for c in data["inclusion"] + data["exclusion"]:
        assert c.get("rationale"), f"{c['id']} missing rationale"


def test_construct_grammar_well_formed_moclo_cr_part_includes():
    part = CandidatePart(
        part_id="CDS_mCherry_Cr",
        so_term="SO:0000316",  # CDS
        assembly_sites=["BsaI"],
        validated_hosts=["chlamydomonas_reinhardtii_nuclear"],
        accession="Addgene_112345",
        accession_type="addgene_plasmid",
        sequence="ATGGCC" * 20,
        citation="Crozet et al. 2018, ACS Synth Biol, doi:10.1021/acssynbio.7b00203",
        functional_evidence="fluorescence confirmed in Cr nuclear transformants",
        assembly_standard="type_iis",
    )
    result = evaluate_candidate(part)
    assert result.verdict is Verdict.INCLUDE
    assert all(c.passed for c in result.checks)


def test_construct_grammar_missing_deposited_sequence_fails_ic3():
    part = CandidatePart(
        part_id="Promoter_HSP70A_RBCS2",
        so_term="SO:0000167",  # promoter
        assembly_sites=["BsaI"],
        validated_hosts=["chlamydomonas_reinhardtii_nuclear"],
        accession="",
        accession_type="supplementary_table_only",
        sequence="",
        citation="Some Paper 2020",
        functional_evidence="reported functional in a supplementary table",
    )
    result = evaluate_candidate(part)
    assert result.verdict is Verdict.EXCLUDE
    failed_ids = {c.id for c in result.failed()}
    assert "IC-3" in failed_ids


def test_construct_grammar_gateway_only_part_fails_ex3():
    part = CandidatePart(
        part_id="GatewayVector_attB1",
        so_term="SO:0000316",
        assembly_sites=[],
        validated_hosts=["chlamydomonas_reinhardtii_nuclear"],
        accession="NC_000001",
        accession_type="genbank_accession",
        sequence="ATGGCC" * 20,
        citation="doi:10.1000/xyz",
        functional_evidence="confirmed",
        assembly_standard="gateway",
    )
    result = evaluate_candidate(part)
    assert result.verdict is Verdict.EXCLUDE
    failed_ids = {c.id for c in result.failed()}
    assert "IC-1" in failed_ids  # no Type IIS sites/overhangs recorded
    assert "EX-3" in failed_ids


def test_construct_grammar_non_cr_validation_is_candidate_tier_only():
    part = CandidatePart(
        part_id="Terminator_Nos_Phytobrick",
        so_term="SO:0000141",  # terminator
        assembly_sites=["BsaI"],
        validated_hosts=["nicotiana_benthamiana"],
        syntax_compliant_only=True,
        accession="Addgene_99999",
        accession_type="addgene_plasmid",
        sequence="ATGGCC" * 20,
        citation="Engler et al. 2014",
        functional_evidence="confirmed in N. benthamiana transient assay",
    )
    result = evaluate_candidate(part)
    assert result.verdict is Verdict.CANDIDATE_TIER_ONLY


def test_construct_grammar_preprint_only_is_pending():
    part = CandidatePart(
        part_id="CDS_novel_reporter",
        so_term="SO:0000316",
        assembly_sites=["BsmBI"],
        validated_hosts=["chlamydomonas_reinhardtii_nuclear"],
        accession="",
        accession_type="described_only",
        sequence="",
        citation="bioRxiv 2026.01.01.000000",
        functional_evidence="reported in preprint",
        is_preprint_only=True,
    )
    result = evaluate_candidate(part)
    assert result.verdict is Verdict.PENDING


def test_construct_grammar_conflicting_sequence_excludes():
    part = CandidatePart(
        part_id="CDS_conflicting",
        so_term="SO:0000316",
        assembly_sites=["BsaI"],
        validated_hosts=["chlamydomonas_reinhardtii_nuclear"],
        accession="Addgene_55555",
        accession_type="addgene_plasmid",
        sequence="ATGGCC" * 20,
        citation="doi:10.1000/abc",
        functional_evidence="confirmed",
        sequence_conflict=True,
    )
    result = evaluate_candidate(part)
    assert result.verdict is Verdict.EXCLUDE
    assert "EX-2" in {c.id for c in result.failed()}


def test_construct_grammar_catalog_metadata_without_sequence_still_excludes():
    """A real, verified catalog entry (e.g. pCM0-001) with no DNA sequence text
    retrieved yet must still fail IC-3 -- a citable source is not the same as
    having the actual sequence in hand."""
    part = CandidatePart(
        part_id="pCM0-001",
        so_term="SO:0000167",
        assembly_sites=["BsaI"],
        validated_hosts=["chlamydomonas_reinhardtii_nuclear"],
        accession="pCM0-001",
        accession_type="chlamycollection_catalog_entry",
        sequence="",  # not yet retrieved
        citation="Crozet et al. 2018, ACS Synth Biol",
        functional_evidence="distributed as a validated MoClo toolkit part",
    )
    result = evaluate_candidate(part)
    assert result.verdict is Verdict.EXCLUDE
    assert "IC-3" in {c.id for c in result.failed()}


def _run_standalone() -> int:
    import inspect
    fns = [v for k, v in sorted(globals().items())
           if k.startswith("test_") and callable(v)]
    failed = skipped = 0
    for fn in fns:
        if inspect.signature(fn).parameters:  # needs pytest fixtures (tmp_path, monkeypatch)
            print(f"SKIP {fn.__name__} (run under pytest for fixture-based tests)")
            skipped += 1
            continue
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except Exception as exc:  # noqa: BLE001
            failed += 1
            print(f"FAIL {fn.__name__}: {exc!r}")
    ran = len(fns) - skipped
    print(f"\n{ran - failed}/{ran} passed ({skipped} skipped — need pytest)")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(_run_standalone())
