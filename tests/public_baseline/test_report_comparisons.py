"""Report derivations preserve original evidence and use exact declared arms."""
import hashlib
from pathlib import Path

import pandas as pd
import pytest

from katana_next.reports.comparisons import _rank, _rs_eligibility, _safe, _cross_patterns, _direction_summary, _Views, build_comparison_views


def _contrast():
    return {"contrast_id": "B_vs_A", "numerator": {"sample_ids": ["B1", "B2", "B3"]},
            "denominator": {"sample_ids": ["A1", "A2", "A3"]}}


def _events():
    return pd.DataFrame([
        dict(row_id="r" + str(i), event_id="e" + str(i), event_type="SE", gene_id="g" + str(i), gene_symbol="G" + str(i),
             comparison_id="B_vs_A", fdr=q, dpsi=dpsi, significant=selected,
             stable_significant=True, IJC_SAMPLE_1=reads, SJC_SAMPLE_1="0,0,0",
             IJC_SAMPLE_2="20,20,20", SJC_SAMPLE_2="0,0,0", IncLevel1=".2,.3,.4", IncLevel2=".4,.5,.6")
        for i, q, dpsi, selected, reads in [(0, 1e-250, .3, True, "10,10,10"),
                                           (1, .001, -.4, True, "20,20,9"),
                                           (2, .9, .1, False, "20,20,20")]
    ])


def test_every_sample_count_gate_rejects_low_third_replicate_and_old_flags():
    frame = _events()
    eligible, reason = _rs_eligibility(frame, _contrast())
    assert eligible.tolist() == [True, False, True]
    assert reason == "original_counts_all_declared_samples_ge10"
    assert frame.stable_significant.all()
    for bad in ["10,10", "10,10,nan", "10,10,-1", "10,10,NA"]:
        broken = frame.copy()
        broken.loc[0, "IJC_SAMPLE_1"] = bad
        assert not _rs_eligibility(broken, _contrast())[0].iloc[0]


def test_ranking_keeps_subnormal_FDR_and_uses_absolute_effect_then_ID():
    frame = pd.DataFrame({"row_id": ["b", "a", "c", "d"], "fdr": [.01, .01, 5e-310, None], "dpsi": [.2, -.3, .1, .8]})
    original = frame.copy(deep=True)
    ranked = _rank(frame, "fdr", "dpsi")
    assert ranked.row_id.tolist() == ["c", "a", "b", "d"]
    assert ranked.iloc[0].fdr == 5e-310
    pd.testing.assert_frame_equal(frame, original)


def test_complete_report_views_keep_both_policies_units_and_native_values(tmp_path, monkeypatch):
    from katana_next.reports import comparisons
    artifacts = []
    for analysis, name, frame in [
        ("differential.as", "Differential AS events", _events()),
        ("differential.deg", "Differentially expressed genes", pd.DataFrame({
            "row_id": ["g0", "g1"], "gene_id": ["g0", "g1"], "gene_symbol": ["G0", "G1"],
            "comparison_id": ["B_vs_A"] * 2, "log2fc": [.3, -2.], "deg_fdr": [.5, .001], "significant": [False, True]})),
        ("differential.dtu", "DTU and exon-usage evidence", pd.DataFrame({
            "row_id": ["tx", "ex"], "gene_id": ["g0", "g1"], "gene_symbol": ["G0", "G1"],
            "transcript_id": ["tx1", "E001"], "comparison_id": ["B_vs_A"] * 2,
            "source_kind": ["dtu_result", "dexseq_result"], "dtu_fdr": [.01, .02],
            "gene_screen_fdr": [.01, None], "significant": [True, True]})),
    ]:
        path = tmp_path / (analysis + ".parquet")
        frame.to_parquet(path)
        artifacts.append({"analysis_id": analysis, "name": name, "path": str(path), "parameters": {}})
    original = {a["path"]: hashlib.sha256(Path(a["path"]).read_bytes()).hexdigest() for a in artifacts}
    # Test the complete data path while one separately tested volcano provides
    # actual image rendering. Other plotting calls still execute and validate.
    def capture(self, fig, identity, chapter, title, comparison, caption, policy=None):
        import matplotlib.pyplot as plt
        plt.close(fig)
        self.panels.append({"id": identity, "chapter": chapter, "contrast_id": comparison, "policy": policy})
    monkeypatch.setattr(comparisons._Views, "figure", capture)
    result = build_comparison_views({"contrasts": [_contrast()], "assets": []}, artifacts, tmp_path / "report")
    assert result["new_scientific_fits"] == 0 and result["source_invariance"]["unchanged"]
    assert result["policy_data"]["significant"]["B_vs_A"]["events"] == 2
    assert result["policy_data"]["stable_significant"]["B_vs_A"]["events"] == 1
    tops = {x["analysis"] + "-" + str(x["policy"]): x for x in result["top10"]}
    assert tops["DAS-stable_significant"]["rows"][0]["fdr"] == 1e-250
    assert len(tops["DAS-significant"]["rows"]) == 2
    assert tops["DTU-None"]["unit"] == "transcript"
    assert tops["DEXSeq-None"]["unit"] == "exon_bin"
    assert "dpsi" not in tops["DTU-None"]["columns"]
    as_table = next(t for t in result["tables"] if t["name"] == "DAS-complete-original-and-report-flags")
    saved = pd.read_parquet(tmp_path / "report" / as_table["path"])
    pd.testing.assert_frame_equal(saved[list(_events())], _events())
    assert saved.report_rs_das.tolist() == [True, False, False]
    for entry in result["gene_sets"]:
        foreground = pd.read_csv(tmp_path / "report" / entry["foreground"], sep="\t")
        background = pd.read_csv(tmp_path / "report" / entry["background"], sep="\t")
        assert set(foreground.gene_id).issubset(background.gene_id)
    assert original == {name: hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in original}


def test_volcano_writes_real_png_svg_and_preserves_zero_FDR(tmp_path):
    from katana_next.reports.comparisons import _Views, _volcano
    frame = _events()
    frame.loc[0, "fdr"] = 0.
    before = frame.copy(deep=True)
    view = _Views({"contrasts": [_contrast()]}, [], tmp_path / "report")
    _volcano(view, frame, "B_vs_A", "DAS", "dpsi", "fdr", frame.significant, policy="stable_significant")
    png = tmp_path / "report" / view.panels[0]["path"]
    assert png.read_bytes().startswith(b"\x89PNG")
    assert "<svg" in png.with_suffix(".svg").read_text(encoding="utf-8")
    pd.testing.assert_frame_equal(frame, before)


def test_distinct_comparison_names_cannot_create_duplicate_html_IDs():
    values = ["B:A", "B A", "B-A", "B" * 180 + "A", "B" * 180 + "B"]
    assert len({_safe(value) for value in values}) == len(values)
    assert _safe("B_vs_A") == "B_vs_A"


def test_exact_upset_keeps_every_pattern_member_and_unique_gene_exon_unit(tmp_path):
    frame = pd.DataFrame({"comparison_id": ["A", "B", "A", "B"],
        "exact_feature_id": ["g1::E1", "g1::E1", "g2::E1", "g2::E2"],
        "gene_id": ["g1", "g1", "g2", "g2"], "gene_symbol": ["G1", "G1", "G2", "G2"],
        "dtu_fdr": [1e-250, .01, .02, .03]})
    before = frame.copy(deep=True)
    view = _Views({}, [], tmp_path / "report")
    result = _cross_patterns(view, frame, "exact_feature_id", ["A", "B"], "DEXSeq-features")
    assert result["objects"] == 3 and result["patterns"] == 3
    counts = pd.read_csv(tmp_path / "report" / result["intersection_counts"], sep="\t")
    assert counts.object_count.sum() == 3
    assert set(counts.included_comparisons) == {"A", "B", "A; B"}
    members = pd.read_parquet(tmp_path / "report" / result["intersection_members"])
    assert members.dtu_fdr.tolist() == before.dtu_fdr.tolist()
    assert len(members) == 4
    assert view.panels[0]["chapter"] == 4 and view.panels[0]["contrast_id"] == "all"
    assert all(t.get("chapter") == 4 for t in view.tables)
    browser = (tmp_path / "report" / result["complete_browser"]).read_text(encoding="utf-8")
    assert "g2::E1" in browser and "g2::E2" in browser
    assert '../../../REPORT.html' in browser
    pd.testing.assert_frame_equal(frame, before)


def test_shared_event_direction_distinguishes_opposition_missing_and_zero():
    frame = pd.DataFrame({"comparison_id": ["A", "A", "B", "B", "C"],
        "event_id": ["e1", "zero", "e1", "zero", "unshared"], "dpsi": [.3, 0., -.2, .5, .1]})
    original = frame.copy(deep=True)
    summary = _direction_summary(frame, ["A", "B", "C"]).set_index(["comparison_a", "comparison_b"])
    row = summary.loc[("A", "B")]
    assert row.shared_selected_events == 2 and row.direction_comparable_events == 1
    assert row.opposite_direction_events == 1 and row.same_direction_fraction == 0
    assert pd.isna(summary.loc[("A", "C"), "same_direction_fraction"])
    assert summary.loc[("A", "C"), "shared_selected_events"] == 0
    pd.testing.assert_frame_equal(frame, original)
