"""Display decisions use existing native FDR, never a new correction family."""
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

from katana_next.reports.motif_display import (
    apply_motif_significance_display, motif_display_summary, significant_motif_ranking,
)


def _ranking(values):
    return pd.DataFrame([dict(model_id=f"m{i:02}", q_new_window_family=q,
                              raw_P_new_paired_swaps=q / 3 if pd.notna(q) else np.nan,
                              comparison_id="B_vs_A", display_rank=i + 1)
                         for i, q in enumerate(values)])


def test_only_significant_models_enter_main_top10_without_invented_values():
    native = _ranking([.001] * 12 + [.05, .051, np.nan])
    before = native.copy(deep=True)
    selected = significant_motif_ranking(native)
    assert selected.model_id.tolist() == [f"m{i:02}" for i in range(10)]
    assert selected.q_new_window_family.eq(.001).all()
    summary = motif_display_summary(native, "B_vs_A")
    assert summary["models"] == 15
    assert summary["significant_models"] == 13
    assert summary["displayed_models"] == 10
    assert summary["model_level_fdr_created"] is False
    pd.testing.assert_frame_equal(native, before)


@pytest.mark.parametrize("values,status", [([.6] * 44, "no_significant_motif"),
                                          ([np.nan] * 3, "not_estimable")])
def test_repeated_nonsignificant_or_nonestimable_native_fdr_gives_compact_notice(values, status):
    summary = motif_display_summary(_ranking(values), "B_vs_A")
    assert summary["status"] == status
    assert summary["significant_models"] == summary["displayed_models"] == 0
    assert summary["main_model_ids"] == []
    assert "No " in summary["notice"]
    assert "outputs" in summary["notice"]


@pytest.mark.parametrize("column", ["q_new_window_family", "raw_P_new_paired_swaps"])
@pytest.mark.parametrize("invalid", [-.1, 1.1, np.inf])
def test_invalid_native_statistics_are_not_hidden_by_display_filter(column, invalid):
    ranking = _ranking([.1])
    ranking.loc[0, column] = invalid
    with pytest.raises(ValueError, match="finite"):
        significant_motif_ranking(ranking)


def _report(tmp_path, values):
    ranking = _ranking(values)
    path = tmp_path / "ranking.tsv"
    ranking.to_csv(path, sep="\t", index=False, float_format="%.17g")
    manifest = {"panels": [dict(id=f"card-{i}", contrast_id="B_vs_A", output_only=False,
        metadata=dict(motif_view_kind="top_motif", model_id=f"m{i:02}")) for i in range(len(values))],
        "complete_results": [{"name": "all-motif-ranked-recorded-tests | B_vs_A",
                              "contrast_id": "B_vs_A", "source": "ranking.tsv"}]}
    return manifest, path


def test_existing_report_filter_preserves_native_files_and_all_output_panels(tmp_path):
    manifest, source = _report(tmp_path, [.02, .05, .051, np.nan])
    before, bytes_before = deepcopy(manifest), source.read_bytes()
    filtered, audit = apply_motif_significance_display(manifest, tmp_path)
    assert manifest == before
    assert source.read_bytes() == bytes_before
    assert [p["id"] for p in filtered["panels"] if not p["output_only"]] == ["card-0", "card-1"]
    assert len(filtered["panels"]) == 4
    assert audit["passed"] and audit["main_sequence_cards"] == 2
    assert audit["native_statistics_unchanged"]
    assert filtered["motif_summaries"][0]["ranking_path"] == "ranking.tsv"


def test_existing_report_zero_significance_is_completed_without_main_cards(tmp_path):
    manifest, _ = _report(tmp_path, [.5925147523709168] * 10)
    filtered, audit = apply_motif_significance_display(manifest, tmp_path)
    assert all(p["output_only"] for p in filtered["panels"])
    assert audit["main_sequence_cards"] == 0
    assert filtered["motif_summaries"][0]["status"] == "no_significant_motif"


def test_missing_significant_sequence_card_cannot_be_silently_omitted(tmp_path):
    manifest, _ = _report(tmp_path, [.01, .1])
    manifest["panels"] = manifest["panels"][1:]
    with pytest.raises(ValueError, match="omits or duplicates"):
        apply_motif_significance_display(manifest, tmp_path)


def test_foreign_ranking_scope_rejected(tmp_path):
    manifest, path = _report(tmp_path, [.01])
    path.write_text(path.read_text().replace("B_vs_A", "OTHER"))
    with pytest.raises(ValueError, match="another comparison"):
        apply_motif_significance_display(manifest, tmp_path)
