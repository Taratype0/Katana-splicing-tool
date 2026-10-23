"""Exercise the real frozen 64MiB display limit with generated display records."""
import hashlib
import json

import pandas as pd
import pytest

from katana_next.reports.motif_native_index import render_budgeted_consequences


def test_real_oversize_frozen_display_uses_complete_file_index_without_rewriting(tmp_path):
    names = ["ANALYSIS_PLAN_BEFORE_SCORING.json", "FROZEN_SELECTED_EVENTS.tsv", "MATCHING_ALL_SELECTED_WITH_LOSSES.tsv",
        "ACTUAL_RNA_FORMS.tsv.gz", "ACTUAL_LOCAL_RNA.fa", "ACTUAL_LOCUS_HITS.tsv.gz", "PAIRED_FORM_INSTANCE_CHANGES.tsv.gz",
        "PER_REPLICATE_SEQUENCE_OPPORTUNITY.tsv.gz", "PAIRED_FORM_NORMALIZED_SCORES.tsv.gz", "NEW_CLUSTER_FWER.tsv",
        "MAX_CLUSTER_NULL.npy", "RESULTS.json"]
    for name in names:
        (tmp_path / name).write_bytes(b"Generated display fixture only; no biological claim.\n")
    # This is an intentionally large display field, not invented biological
    # statistics. It trips the actual frozen serializer's byte-size limit.
    tests = pd.DataFrame({"model_id": ["fixture"], "raw_P_new_paired_swaps": [.01], "q_new_window_family": [.02],
                          "generated_display_field": ["X" * (65 * 1024 * 1024)]})
    tests.to_csv(tmp_path / "NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz", sep="\t", index=False)
    before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in tmp_path.iterdir()}
    summary = {"comparison_id": "generated_display_only", "inference_status": "fixture_not_computed"}
    path = tmp_path / "MOTIFTRACE_EVENT_CONSEQUENCES.html"
    render_budgeted_consequences(pd.DataFrame(), pd.DataFrame(), tests, pd.DataFrame(), summary, path)
    assert path.stat().st_size < 20000
    page = path.read_text()
    payload = json.JSONDecoder().raw_decode(page.split("const D=", 1)[1])[0]
    assert payload["native_statistics_modified"] is False
    assert payload["statistical_family_split"] is False
    assert payload["presentation"] == "complete_native_files_index"
    for name, digest in before.items():
        assert name in payload["available_files"]
        assert hashlib.sha256((tmp_path / name).read_bytes()).hexdigest() == digest
        assert payload["file_identities"][name]["sha256"] == digest
    recovered = pd.read_csv(tmp_path / "NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz", sep="\t")
    pd.testing.assert_frame_equal(tests, recovered)


def test_other_native_render_errors_are_not_hidden(tmp_path):
    with pytest.raises(ValueError):
        render_budgeted_consequences(pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(),
                                    {"invalid": float("inf")}, tmp_path / "report.html")
    assert not (tmp_path / "report.html").exists()
