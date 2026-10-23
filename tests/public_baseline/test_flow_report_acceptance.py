"""Report-gate contracts on real rendering with synthetic lifecycle receipts.

These tests do not claim native MAEI scientific execution or full-profile E2E
acceptance. They exercise unchanged native-table transport, actual figure/data
rendering and exact task/applicability proof validation independently of R fits.
"""
from pathlib import Path
import shutil

import pandas as pd
import pytest

from katana_next.flow.acceptance import required_report_gaps, verify_mapp_not_applicable
from katana_next.flow.mapp_views import build_mapp_views
from katana_next.flow.model import file_sha256, read_json, seal_plan, write_json
from katana_next.flow.runtime import canonical_hash, file_record


def _fixture(tmp_path, *, parallel=False, not_applicable=False):
    from katana_next.full_fixture import create_full_fixture
    source = create_full_fixture(tmp_path / "study")
    root = tmp_path / "collection"
    root.mkdir()
    shutil.copyfile(source, root / "project.json")
    project = read_json(root / "project.json")
    cid = project["contrasts"][0]["contrast_id"]
    if parallel:
        project["assets"].append({"asset_id": "declared-but-absent-parallel", "kind": "motiftrace_parallel_results",
                                  "path": str(tmp_path / "missing-parallel.json")})
        write_json(root / "project.json", project)
    task = {"id": "native", "workflow_id": "flow.mapp.maei", "adapter": "mapp.maei", "dependencies": [],
            "parameters": {}, "scope": {"contrast_ids": [cid], "site": "3ss", "window_start": -50, "window_end": 0}}
    plan = seal_plan({"project_manifest": str(source), "project_manifest_sha256": file_sha256(source),
                      "contrast_ids": [cid], "configuration": {"design_factors": ["condition"]}, "tasks": [task]})
    task = plan["tasks"][0]
    write_json(root / "plan.json", plan)
    original = tmp_path / "native" / "AUTHOR_NATIVE_MIXTURE_P_BONFERRONI.tsv"
    original.parent.mkdir()
    pd.DataFrame({"activities_A1": [-.2, .4], "activities_B1": [.3, -.1],
                  "combined.Zscore": [2., 3.], "corr_pval_A1": [.1, .02], "corr_pval_B1": [.04, .2]},
                 index=["window|AAA", "window|CCC"]).to_csv(original, sep="\t")
    products, status, evidence = [original], "completed", {"native_fitting_performed": True}
    if not_applicable:
        original.unlink()
        audit = original.with_name("complete_eligibility.tsv")
        pd.DataFrame([{"event_id": "SE:1", "gene_id": "G1", "status": "no_exact_annotated_included_or_skipped_transcript"}]).to_csv(audit, sep="\t", index=False)
        contrast = project["contrasts"][0]
        proof = {"schema": "katana.flow.mapp_applicability/v1", "status": "not_applicable", "reason_code": "no_exact_mapped_se",
                 "comparison_id": cid, "site": "3ss", "window_start": -50, "window_end": 0,
                 "sample_ids": contrast["numerator"]["sample_ids"] + contrast["denominator"]["sample_ids"],
                 "complete_measured_universe": True, "native_fitting_performed": False, "eligible_events": 0,
                 "measured_events": 1, "eligibility_table": str(audit), "eligibility_table_sha256": file_sha256(audit)}
        proof_path = write_json(original.with_name("applicability.json"), proof)
        products, status = [audit, proof_path], "not_applicable"
        evidence = {"native_fitting_performed": False, "applicability_proof": str(proof_path)}
    receipt = {"schema": "katana.flow.task_receipt/v1", "task_id": task["id"], "plan_id": plan["plan_id"],
               "task_contract_sha256": canonical_hash(task), "source_project_manifest_sha256": plan["project_manifest_sha256"],
               "status": status, "execution_evidence": evidence,
               "outputs": [file_record(path) for path in products], "synthetic_lifecycle_fixture": True}
    receipt["receipt_sha256"] = canonical_hash(receipt)
    receipt_path = write_json(root / "evidence/native/task.receipt.json", receipt)
    transport = []
    for path in products:
        destination = root / "products" / path.name
        destination.parent.mkdir(exist_ok=True)
        shutil.copyfile(path, destination)
        transport.append({"source": file_record(path), "collection_path": destination.relative_to(root).as_posix()})
    write_json(root / "transport-provenance.json", {"products": transport})
    collection = {"project_manifest": "project.json", "task_receipts": [{"task_id": "native",
                  "collection_path": "evidence/native/task.receipt.json", "sha256": file_sha256(receipt_path)}]}
    receipts = {"native": (receipt_path, receipt)}
    views = build_mapp_views(root, plan, receipts, transport)
    from katana_next.reports import build_report
    report_dir = tmp_path / "report"
    actual = build_report({"project_manifest": str(root / "project.json"), "artifacts": [],
                           "extra_views": views["extra_views"], "report_require_plot_data": True}, report_dir)
    return root, plan, receipts, collection, actual, report_dir


def _mapp_gaps(result):
    return [gap for gap in result if gap.get("component") in {"mapp_native", "mapp_parallel"}]


def test_verified_native_activity_render_does_not_require_absent_optional_parallel(tmp_path):
    root, plan, receipts, collection, report, out = _fixture(tmp_path)
    gaps = required_report_gaps(root, plan, report, out, collection_receipt=collection)
    assert not _mapp_gaps(gaps)
    panel = next(panel for panel in report["panels"] if panel.get("metadata", {}).get("component") == "mapp_native")
    assert panel["metadata"]["fit_task_ids"] == ["native"]
    assert {link["role"] for link in panel["data_links"]} == {"plot_data", "full_results"}
    optional = next(row for row in report["unavailable"] if row.get("metadata", {}).get("component") == "mapp_parallel")
    assert optional["metadata"]["requirement"] == "conditional_import"
    assert "no new model is fitted" not in (out / "REPORT.html").read_text(encoding="utf-8")
    assert (out / "REPORT.zip").stat().st_size > 0


def test_native_receipt_verification_accepts_unicode_resource_paths(tmp_path):
    root, plan, receipts, collection, report, out = _fixture(tmp_path / "项目")
    assert not _mapp_gaps(required_report_gaps(root, plan, report, out, collection_receipt=collection))


def test_missing_required_native_fit_fails_even_when_optional_import_is_not_requested(tmp_path):
    root, plan, receipts, collection, report, out = _fixture(tmp_path)
    path, receipt = receipts["native"]
    receipt["execution_evidence"]["native_fitting_performed"] = False
    receipt.pop("receipt_sha256")
    receipt["receipt_sha256"] = canonical_hash(receipt)
    write_json(path, receipt)
    collection["task_receipts"][0]["sha256"] = file_sha256(path)
    gaps = _mapp_gaps(required_report_gaps(root, plan, report, out, collection_receipt=collection))
    assert any(gap["component"] == "mapp_native" and "not actually completed" in gap["reason"] for gap in gaps)


def test_supplied_but_missing_parallel_collection_is_a_required_failure(tmp_path):
    root, plan, receipts, collection, report, out = _fixture(tmp_path, parallel=True)
    gaps = required_report_gaps(root, plan, report, out, collection_receipt=collection)
    assert any(gap["component"] == "mapp_parallel" and "no planned import" in gap["reason"] for gap in gaps)
    assert not any(gap["component"] == "mapp_native" for gap in gaps)


def test_one_native_fit_panel_cannot_cover_another_missing_window(tmp_path):
    root, plan, receipts, collection, report, out = _fixture(tmp_path)
    panel = next(panel for panel in report["panels"] if panel.get("metadata", {}).get("component") == "mapp_native")
    panel["metadata"]["fit_task_ids"].append("not-a-completed-window")
    gaps = _mapp_gaps(required_report_gaps(root, plan, report, out, collection_receipt=collection))
    assert any("every completed window" in gap["reason"] for gap in gaps)


def test_optional_caption_words_cannot_disable_required_typed_clip_gap(tmp_path):
    root, plan, receipts, collection, report, out = _fixture(tmp_path)
    report["unavailable"] = [{"chapter": 7, "metadata": {"component": "clip"}, "reason": "Missing matched binding results"},
                              {"chapter": 7, "title": "Unrelated attachment", "reason": "mapp bootstrap interaction in a file name"}]
    gaps = required_report_gaps(root, plan, report, out, collection_receipt=collection)
    assert [gap["component"] for gap in gaps] == ["clip"]


def test_only_complete_exact_zero_eligibility_proof_is_not_applicable(tmp_path):
    root = tmp_path / "own"
    root.mkdir()
    table = root / "eligibility.tsv"
    pd.DataFrame([{"event_id": "SE:1", "gene_id": "G1", "status": "no_exact_annotated_included_or_skipped_transcript"}]).to_csv(table, sep="\t", index=False)
    proof = {"schema": "katana.flow.mapp_applicability/v1", "status": "not_applicable", "reason_code": "no_exact_mapped_se",
             "comparison_id": "B_vs_A", "site": "3ss", "window_start": -50, "window_end": 0,
             "complete_measured_universe": True, "native_fitting_performed": False, "eligible_events": 0,
             "measured_events": 1, "eligibility_table": str(table), "eligibility_table_sha256": file_sha256(table)}
    path = write_json(root / "applicability.json", proof)
    receipt = {"status": "not_applicable", "execution_evidence": {"native_fitting_performed": False,
               "applicability_proof": str(path)}}
    task = {"scope": {"contrast_ids": ["B_vs_A"], "site": "3ss", "window_start": -50, "window_end": 0}}
    assert verify_mapp_not_applicable(task, receipt, Path)["measured_events"] == 1
    task["scope"]["contrast_ids"] = ["C_vs_A"]
    with pytest.raises(ValueError, match="different comparison"):
        verify_mapp_not_applicable(task, receipt, Path)
    task["scope"]["contrast_ids"] = ["B_vs_A"]
    proof["complete_measured_universe"] = False
    write_json(path, proof)
    with pytest.raises(ValueError, match="complete zero-eligible"):
        verify_mapp_not_applicable(task, receipt, Path)


def test_zero_exact_eligible_targets_are_disclosed_without_fabricating_native_activity(tmp_path):
    root, plan, receipts, collection, report, out = _fixture(tmp_path, not_applicable=True)
    gaps = required_report_gaps(root, plan, report, out, collection_receipt=collection)
    assert not _mapp_gaps(gaps)
    assert not any(panel.get("metadata", {}).get("component") == "mapp_native" for panel in report["panels"])
    note = next(row for row in report["unavailable"] if row.get("task_id") == "native")
    assert note["status"] == "not_applicable" and note["metadata"]["eligible_events"] == 0
    assert not any(Path(row["path"]).name == "AUTHOR_NATIVE_MIXTURE_P_BONFERRONI.tsv"
                   for row in receipts["native"][1]["outputs"])
