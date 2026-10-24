"""Positive comparison-complete acceptance of real rendered generated reports.

The supplied statistical records are controlled display fixtures, not claims
of native computation; the full generic Nextflow acceptance covers computation.
"""
import copy
import json
from pathlib import Path

import pandas as pd
import pytest

from katana_next.flow.motif_acceptance import required_motif_report_gaps, _empty_selection
from katana_next.flow.model import file_sha256
from katana_next.flow.presentation import build_presentation_assets
from katana_next.reports import build_report
from test_flow_presentation import _fixture_collection
from test_report_multicontrast_motif_acceptance import _controlled_family


@pytest.fixture(scope="module")
def accepted(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("positive-motif-gate")
    ids = ["B_vs_A", "C_vs_A"]
    project = {"project_id": "synthetic-positive-acceptance", "name": "Generated acceptance only",
        "metadata": {"synthetic": True}, "annotation": {"genome_build": "generated"}, "samples": [],
        "contrasts": [{"contrast_id": cid, "numerator": {"label": cid[0], "sample_ids": []},
                       "denominator": {"label": "A", "sample_ids": []}} for cid in ids]}
    tasks, products = [], {}
    for index, cid in enumerate(ids):
        source, _proof = _controlled_family(tmp_path / "sources" / cid, cid, index * 3)
        tasks.append({"id": cid, "workflow_id": "prioritize.motiftrace_consequences", "dependencies": [],
                      "scope": {"contrast_ids": [cid]}})
        products[cid] = list(source.iterdir()), {}
    root = tmp_path / "collection"
    plan, receipt_pairs, transport = _fixture_collection(root, project, tasks, products)
    views = build_presentation_assets(root, plan, receipt_pairs, transport)
    report_dir = tmp_path / "report"
    report = build_report({"project_manifest": str(root / "project.json"), "artifacts": [],
                           "extra_views": views["extra_views"], "dummy": True,
                           "report_require_plot_data": True}, report_dir)
    receipts = {key: value for key, (_path, value) in receipt_pairs.items()}
    def transported(receipt, original):
        record = next(record for record in receipt["outputs"] if record["path"] == original)
        member = next(row for row in transport if row["source"]["path"] == original)
        path = root / member["collection_path"]
        assert file_sha256(path) == record["sha256"]
        return path
    return root, plan, report, report_dir, project, receipts, transported


def check(case, *, plan=None, report=None, receipts=None):
    root, base_plan, base_report, directory, project, base_receipts, transport = case
    return required_motif_report_gaps(root, plan or base_plan, base_report if report is None else report,
        directory, project=project, receipts=base_receipts if receipts is None else receipts, transported=transport)


def cards(report, cid="B_vs_A"):
    return [p for p in report["panels"] if p.get("contrast_id") == cid and
            p.get("metadata", {}).get("motif_view_kind") == "top_motif"]


def test_two_rendered_comparisons_have_full_families_and_top10(accepted):
    assert check(accepted) == []
    for cid in ("B_vs_A", "C_vs_A"):
        assert len(cards(accepted[2], cid)) == 7
    (accepted[3].parent / "positive-gate.receipt.json").write_text(json.dumps({
        "passed": True, "scope": "generated controlled display records, no biological claims",
        "comparisons": 2, "declared_models_per_comparison": 12, "main_motif_cards": 14,
        "native_same_row_P_FDR_effect": True, "complete_family_attachments": True}, indent=2))


def test_missing_whole_comparison_without_unavailable_cannot_pass(accepted):
    report = copy.deepcopy(accepted[2])
    report["panels"] = [p for p in report["panels"] if p.get("contrast_id") != "B_vs_A"]
    report["unavailable"] = []
    assert any(g.get("contrast_id") == "B_vs_A" for g in check(accepted, report=report))


@pytest.mark.parametrize("kind", ["omitted", "duplicate", "pooled"])
def test_each_expected_comparison_requires_one_independent_task(accepted, kind):
    plan = copy.deepcopy(accepted[1])
    if kind == "omitted":
        plan["tasks"] = plan["tasks"][1:]
    elif kind == "duplicate":
        plan["tasks"].append(copy.deepcopy(plan["tasks"][0]))
    else:
        plan["tasks"][0]["scope"]["contrast_ids"] = ["B_vs_A", "C_vs_A"]
    assert check(accepted, plan=plan)


@pytest.mark.parametrize("kind", ["card", "output_only", "model", "rank", "family", "scope", "link_scope", "native_link", "image", "full_family"])
def test_partial_or_substituted_motif_display_is_rejected(accepted, kind):
    report = copy.deepcopy(accepted[2])
    panel = cards(report)[0]
    if kind == "card": report["panels"].remove(panel)
    elif kind == "output_only": panel["output_only"] = True
    elif kind == "model": panel["metadata"]["model_id"] = "different"
    elif kind == "rank": panel["metadata"]["display_rank"] = 12
    elif kind == "family": panel["metadata"]["statistical_family_split"] = True
    elif kind == "scope": panel["metadata"]["das_scope"]["recorded_event_set"] = "significant"
    elif kind == "link_scope": panel["data_links"][0]["contrast_id"] = "C_vs_A"
    elif kind == "native_link": panel["data_links"] = [r for r in panel["data_links"] if r.get("role") != "full_results"]
    elif kind == "image": panel["metadata"]["motif_image_sha256"] = "0" * 64
    elif kind == "full_family": report["complete_results"] = []
    assert check(accepted, report=report)


def test_different_test_P_cannot_be_combined_with_best_FDR(accepted):
    report = copy.deepcopy(accepted[2])
    ranking_link = next(l for l in cards(report)[0]["data_links"] if l["label"] == "All ranked models and native test statistics")
    ranking = pd.read_csv(accepted[3] / ranking_link["path"], sep="\t", float_precision="round_trip")
    ranking.loc[0, "raw_P_new_paired_swaps"] = 1e-9
    changed = accepted[3] / "deliberately-corrupted-ranking.tsv"
    ranking.to_csv(changed, sep="\t", index=False, float_format="%.17g")
    for panel in cards(report):
        for link in panel["data_links"]:
            if link["label"] == ranking_link["label"]:
                link.update(path=changed.name, portable_sha256=file_sha256(changed))
    for table in report["complete_results"]:
        if table.get("name") == "all-motif-ranked-recorded-tests | B_vs_A":
            table["source"] = changed.name
    assert any("same actual selected native test row" in g["reason"] for g in check(accepted, report=report))


def _zero_foreground(tmp_path, *, positive=False):
    from katana_next.science.bulk import READ_SUPPORT_POLICY_VERSION
    source = tmp_path / "events.parquet"
    pd.DataFrame({"comparison_id": ["B_vs_A"], "significant": [positive], "stable_significant": [positive],
                  "measurement_eligible": [True], "read_support_policy": [READ_SUPPORT_POLICY_VERSION],
                  "criterion_min_reads": [10]}).to_parquet(source)
    selection = {"source_artifact_id": "full-as", "source_rows": 1, "selected_rows": 0,
                 "event_set": "stable_significant", "contrast_ids": ["B_vs_A"]}
    audit = tmp_path / "audit.parquet"
    pd.DataFrame([selection]).to_parquet(audit)
    task = {"selection": {"analysis_id": "differential.as", "event_set": "stable_significant", "contrast_ids": ["B_vs_A"]}}
    receipt = {"status": "not_applicable", "selection": selection, "execution_evidence": selection,
        "artifacts": [{"artifact_id": "full-as", "analysis_id": "differential.as", "path": str(source), "contrast_ids": ["B_vs_A"]}],
        "generated_artifacts": [{"entity_type": "eligibility_audit", "path": str(audit)}]}
    owner = {"status": "completed", "outputs": [{"path": str(source)}]}
    return task, receipt, {"as": owner}, lambda receipt, path: Path(path)


def test_real_zero_RS_foreground_is_allowed_only_with_full_source_and_audit(tmp_path):
    values = _zero_foreground(tmp_path)
    _empty_selection(*values, "B_vs_A")
    values[1]["generated_artifacts"] = []
    with pytest.raises(ValueError, match="eligibility audit"):
        _empty_selection(*values, "B_vs_A")


def test_zero_claim_with_actual_positive_RS_foreground_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="contradicts"):
        _empty_selection(*_zero_foreground(tmp_path, positive=True), "B_vs_A")


def test_completed_sequence_with_no_matched_pairs_keeps_models_and_explicit_NA(tmp_path):
    from test_motiftrace_execution_budget import _generated
    from katana_next.science.motiftrace_execution import analyze_consequences
    frame, ids, _models, fasta = _generated(tmp_path)
    frame = frame.loc[frame.stable_significant].iloc[:2].copy()
    ids = set(frame.event_id)
    # A long homopolymer has no hit in generated ACGT RNA; missing matching
    # controls makes the native statistical family genuinely zero-pair.
    models = {"absent_sequence": {"kind": "consensus", "pattern": "A" * 20, "width": 20}}
    source = tmp_path / "native"
    analyze_consequences(frame, ids, models, fasta, source, permutations=99)
    project = {"project_id": "generated-empty-matching", "samples": [], "contrasts": [{"contrast_id": "B_vs_A",
               "numerator": {"label": "B", "sample_ids": []}, "denominator": {"label": "A", "sample_ids": []}}]}
    tasks = [{"id": "motif", "workflow_id": "prioritize.motiftrace_consequences", "scope": {"contrast_ids": ["B_vs_A"]}, "dependencies": []}]
    root = tmp_path / "collection"
    plan, pairs, transport = _fixture_collection(root, project, tasks, {"motif": (list(source.iterdir()), {})})
    views = build_presentation_assets(root, plan, pairs, transport)
    report_dir = tmp_path / "report"
    report = build_report({"project_manifest": str(root / "project.json"), "artifacts": [], "extra_views": views["extra_views"],
                           "dummy": True, "report_require_plot_data": True}, report_dir)
    def transported(receipt, original):
        row = next(r for r in transport if r["source"]["path"] == original)
        return root / row["collection_path"]
    actual = required_motif_report_gaps(root, plan, report, report_dir, project=project,
                                      receipts={key: value for key, (_, value) in pairs.items()}, transported=transported)
    assert actual == []
    assert len(cards(report)) == 0
    assert report["motif_summaries"][0]["status"] == "not_estimable"
    assert report["motif_summaries"][0]["significant_models"] == 0
    assert len(pd.read_csv(source / "NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz", sep="\t")) == 0


def test_no_significant_native_motifs_require_truthful_notice_and_complete_outputs(tmp_path):
    source, _ = _controlled_family(tmp_path/"source", "B_vs_A", 0)
    native = source/"NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz"
    frame = pd.read_csv(native, sep="\t")
    frame["q_new_window_family"] = .7
    frame["raw_P_new_paired_swaps"] = .2
    frame.to_csv(native, sep="\t", index=False)
    complete = json.loads((source/"COMPLETE.json").read_text())
    complete["members"][native.name] = {"sha256":file_sha256(native), "bytes":native.stat().st_size}
    (source/"COMPLETE.json").write_text(json.dumps(complete))
    project = {"project_id":"controlled-nonsignificant", "samples":[], "contrasts":[{"contrast_id":"B_vs_A",
               "numerator":{"label":"B","sample_ids":[]}, "denominator":{"label":"A","sample_ids":[]}}]}
    tasks = [{"id":"motif", "workflow_id":"prioritize.motiftrace_consequences", "dependencies":[],
              "scope":{"contrast_ids":["B_vs_A"]}}]
    root = tmp_path/"collection"
    plan,pairs,transport = _fixture_collection(root, project, tasks, {"motif":(list(source.iterdir()),{})})
    views = build_presentation_assets(root, plan, pairs, transport)
    report_dir = tmp_path/"report"
    report = build_report({"project_manifest":str(root/"project.json"), "artifacts":[],
                           "extra_views":views["extra_views"], "dummy":True, "report_require_plot_data":True}, report_dir)
    def transported(receipt, original):
        return root/next(row["collection_path"] for row in transport if row["source"]["path"] == original)
    def gaps(value):
        return required_motif_report_gaps(root,plan,value,report_dir,project=project,
               receipts={key:value for key,(_,value) in pairs.items()},transported=transported)
    assert not gaps(report)
    assert cards(report) == []
    assert report["motif_summaries"][0]["status"] == "no_significant_motif"
    text = (report_dir/"REPORT.html").read_text("utf-8")
    assert "No motif has a recorded native test with FDR" in text
    assert "top10-motif-" not in text
    missing = copy.deepcopy(report); missing["motif_summaries"] = []
    assert any("truthful significance summary" in g["reason"] for g in gaps(missing))
    forged = copy.deepcopy(report); forged["motif_summaries"][0]["significant_models"] = 1
    assert any("contradicts" in g["reason"] for g in gaps(forged))
    dropped = copy.deepcopy(report); dropped["complete_results"] = [r for r in dropped["complete_results"] if "NEW_WINDOW_AND_FULL_SEQUENCE_TESTS" not in r["name"]]
    assert any("complete scientific family" in g["reason"] for g in gaps(dropped))
