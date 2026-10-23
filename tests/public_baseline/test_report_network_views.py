"""Network views preserve populations, matched contrasts and frozen corrections."""
from hashlib import sha256
import json

import numpy as np
import pandas as pd
import pytest

from katana_cosplice.trait_adapter import safe_trait_association
from katana_next.reports.network_views import build_network_assets, matched_module_associations


def _data(tmp_path):
    src = tmp_path / "source"
    src.mkdir()
    sample_ids = ["B1", "B2", "B3", "A1", "A2", "A3", "C1", "C2", "C3"]
    samples = pd.DataFrame({"sample_id": sample_ids, "biological_sample_id": sample_ids,
                            "condition": ["B"] * 3 + ["A"] * 3 + ["C"] * 3})
    scores = pd.DataFrame({"M01": [1.2, 1.1, 1.4, -.2, -.4, -.3, 90, 92, 96],
                           "M02": [.3, .2, .5, .8, .9, 1.2, 70, 72, 74]}, index=sample_ids)
    assignments = pd.DataFrame({"event_key": ["e1", "e2", "e3", "e4", "e5"],
        "module": ["M01", "M01", "M02", "M02", "grey"],
        "network_connectivity": [4., 3., 2., 1., .1],
        "module_membership": [.9, .7, .8, .6, np.nan],
        "stability_score": [.5, .4, .3, .2, .1]})
    annotation = pd.DataFrame({"event_key": assignments.event_key,
        "gene_symbol": ["G1", "G1", "G2", "G3", "G4"],
        "event_type": ["SE", "RI", "A3SS", "MXE", "A5SS"]})
    psi = pd.DataFrame(np.linspace(.1, .9, 45).reshape(5, 9), index=assignments.event_key, columns=sample_ids)
    for stem, frame in (("assignments", assignments), ("annotation", annotation), ("samples", samples)):
        frame.to_csv(src / (stem + ".tsv"), sep="\t", index=False)
    scores.to_csv(src / "scores.tsv", sep="\t")
    psi.to_csv(src / "psi.tsv", sep="\t")
    stages = {key: {"status": "finished"} for key in ("network_observed", "associate_observed", "annotate_modules")}
    (src / "completion.json").write_text(json.dumps({"status": "completed", "completion_scope": "observed_network_baseline_only",
        "baseline_job": "123", "actual_sacct": "123|COMPLETED|0:0|256G\n", "stage_snapshot": stages}))
    (src / "sets.gmt").write_text("term1\tdescription\tG1\tG2\nterm2\tdescription\tG4\tGENOME_ONLY_GENE\n")
    contrasts = [{"contrast_id": "B_vs_A", "numerator": {"sample_ids": sample_ids[:3]}, "denominator": {"sample_ids": sample_ids[3:6]}}]
    return src, samples, scores, contrasts


def test_matched_scores_exclude_other_groups_and_preserve_full_production_family(tmp_path):
    _, samples, scores, contrasts = _data(tmp_path)
    complete, selected, audits = matched_module_associations(scores, samples, contrasts)
    assert len(complete) == 6 and len(selected) == 2
    assert audits[0]["adapter_audit"]["BH_family_size"] == 6
    expected_samples = samples.iloc[:6].copy()
    expected_samples["declared_contrast_arm"] = ["numerator"] * 3 + ["denominator"] * 3
    expected = safe_trait_association(scores.iloc[:6].rename_axis("sample_id"), expected_samples, factor="declared_contrast_arm", unit_column="biological_sample_id").condition_tests
    pd.testing.assert_frame_equal(complete[expected.columns], expected)
    assert selected.set_index("module").loc["M01", "effect"] == pytest.approx(1.5333333333333332)
    assert "C1" not in str(audits)


def test_constant_scores_stay_not_computed_and_no_relabelled_vsrest(tmp_path):
    _, samples, scores, contrasts = _data(tmp_path)
    scores.iloc[:6, 0] = 1
    complete, selected, audits = matched_module_associations(scores, samples, contrasts)
    assert complete.empty and selected.empty
    assert audits[0]["status"] == "not_computed"
    assert "Constant" in audits[0]["reason"]


def test_network_views_keep_all_events_grey_measured_background_and_sources(tmp_path):
    src, _, _, contrasts = _data(tmp_path)
    before = {p.name: sha256(p.read_bytes()).hexdigest() for p in src.iterdir()}
    out = tmp_path / "views"
    result = build_network_assets(assignments_path=src / "assignments.tsv", scores_path=src / "scores.tsv",
        psi_path=src / "psi.tsv", annotation_path=src / "annotation.tsv", samples_path=src / "samples.tsv",
        completion_manifest=src / "completion.json", contrasts=contrasts, output_dir=out,
        libraries={"fixed_mouse": src / "sets.gmt"}, project_id="test-project")
    assert result["receipt"]["events"] == 5
    assert result["receipt"]["unassigned_events"] == 1
    assert result["receipt"]["bootstrap_accepted"] is False
    background = pd.read_csv(out / "tables/enrichment_background.tsv", sep="\t")
    assert set(background.gene_symbol) == {"G1", "G2", "G3", "G4"}
    assert "GENOME_ONLY_GENE" not in set(background.gene_symbol)
    tested = pd.read_csv(out / "tables/module_enrichment.tsv.gz", sep="\t")
    assert len(tested) == 4
    assert tested.background_genes.eq(4).all()
    assert len(pd.read_csv(out / "tables/all_module_members.tsv", sep="\t")) == 5
    assert len(pd.read_csv(out / "tables/all_event_psi.tsv.gz", sep="\t")) == 5
    page = (out / "browser/ALL_MODULES.html").read_text(encoding="utf-8")
    assert "e5" in page and "60 rows per page" in page
    assert "global adjacency sum" in page
    assert len(list((out / "figures").glob("*.png"))) >= 10
    assert before == {p.name: sha256(p.read_bytes()).hexdigest() for p in src.iterdir()}
    assert all(a["project_id"] == "test-project" for a in result["artifacts"])


def test_native_running_receipt_cannot_create_views(tmp_path):
    src, _, _, contrasts = _data(tmp_path)
    proof = json.loads((src / "completion.json").read_text())
    proof["actual_sacct"] = "123|RUNNING|0:0|256G\n"
    (src / "completion.json").write_text(json.dumps(proof))
    out = tmp_path / "views"
    with pytest.raises(ValueError, match="Native baseline"):
        build_network_assets(assignments_path=src / "assignments.tsv", scores_path=src / "scores.tsv",
            psi_path=src / "psi.tsv", annotation_path=src / "annotation.tsv", samples_path=src / "samples.tsv",
            completion_manifest=src / "completion.json", contrasts=contrasts, output_dir=out)
    assert not out.exists()


def _recovered_data(tmp_path):
    """Independent-native collection fixture with an explicitly FAILED source."""
    from katana_next.workflows.module_functions import enrich_modules
    src, samples, scores, contrasts = _data(tmp_path)
    phase = src / "postprocessing"
    phase.mkdir()
    full, selected, audits = matched_module_associations(scores, samples, contrasts)
    assignments = pd.read_csv(src / "assignments.tsv", sep="\t")
    annotation = pd.read_csv(src / "annotation.tsv", sep="\t")
    tested, background = enrich_modules(assignments, annotation, {"fixed_mouse": src / "sets.gmt"})
    full.to_csv(phase / "trait_full_family.tsv", sep="\t", index=False)
    selected.to_csv(phase / "trait_selected.tsv", sep="\t", index=False)
    (phase / "trait_audit.json").write_text(json.dumps(audits))
    tested.to_csv(phase / "module_enrichment.tsv.gz", sep="\t", index=False)
    background.to_csv(phase / "module_background.tsv", sep="\t", index=False)
    assignments.merge(annotation, on="event_key")[["module", "gene_symbol"]].drop_duplicates().to_csv(phase / "module_input_genes.tsv", sep="\t", index=False)

    def identity(path):
        return {"path": path.relative_to(src).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path.read_bytes()).hexdigest()}

    (src / "PLAN.json").write_text(json.dumps({"source_job_id": "123", "design": {"contrasts": contrasts}}))
    (src / "IMPORT_PROVENANCE.json").write_text(json.dumps({"schema": "katana.verified_failed_plot_numerical_import/v1",
        "status": "verified_import", "source_job_id": "123", "fit_network_rerun": False}))
    (phase / "completion.receipt.json").write_text(json.dumps({
        "schema": "katana.network_verified_numeric_postprocessing/v1", "status": "completed",
        "completion_scope": "observed_numerical_import_and_postprocessing_only", "native_job_id": "456",
        "source_import_manifest_sha256": identity(src / "IMPORT_PROVENANCE.json")["sha256"],
        "plan_sha256": identity(src / "PLAN.json")["sha256"], "declared_contrasts": [c["contrast_id"] for c in contrasts],
        "source_failed_status_preserved": True, "frozen_kernels_unchanged": True, "fit_network_rerun": False,
        "files": [identity(path) for path in sorted(phase.iterdir())]}))
    (src / "native_completion.receipt.json").write_text(json.dumps({
        "schema": "katana.network_numeric_recovery_native_acceptance/v1", "passed": True,
        "completion_scope": "observed_numerical_import_and_postprocessing_only", "native_job_id": "456",
        "actual_sacct": {"returncode": 0, "stdout": "456|COMPLETED|0:0|\n"},
        "source_actual_sacct": {"returncode": 0, "stdout": "123|FAILED|1:0|\n"},
        "plan_sha256": identity(src / "PLAN.json")["sha256"], "source_failed_status_preserved": True,
        "fit_network_rerun": False, "module_bootstrap_accepted": False,
        "files": [identity(path) for path in sorted(src.rglob("*")) if path.is_file()]}))
    return src, contrasts, full, selected, tested


def test_verified_failed_plot_import_uses_completed_statistics_without_refitting(tmp_path, monkeypatch):
    import katana_next.reports.network_views as views
    import katana_next.workflows.module_functions as functions
    src, contrasts, full, selected, tested = _recovered_data(tmp_path)
    before = {str(path): sha256(path.read_bytes()).hexdigest() for path in src.rglob("*") if path.is_file()}

    def forbidden(*args, **kwargs):
        raise AssertionError("A completed family must not be fitted a second time")

    monkeypatch.setattr(views, "matched_module_associations", forbidden)
    monkeypatch.setattr(functions, "enrich_modules", forbidden)
    out = tmp_path / "views"
    result = views.build_network_assets(assignments_path=src / "assignments.tsv", scores_path=src / "scores.tsv",
        psi_path=src / "psi.tsv", annotation_path=src / "annotation.tsv", samples_path=src / "samples.tsv",
        completion_manifest=src / "native_completion.receipt.json", contrasts=contrasts, output_dir=out)
    assert (out / "tables/matched_module_tests_full_family.tsv").read_bytes() == (src / "postprocessing/trait_full_family.tsv").read_bytes()
    assert (out / "tables/module_enrichment.tsv.gz").read_bytes() == (src / "postprocessing/module_enrichment.tsv.gz").read_bytes()
    assert result["receipt"]["new_trait_fits"] == 0 and result["receipt"]["new_module_ora"] is False
    assert result["receipt"]["native_baseline_job"] is None
    assert result["receipt"]["native_postprocessing_job"] == "456"
    assert result["receipt"]["original_source_native_job"] == "123"
    assert result["receipt"]["bootstrap_accepted"] is False
    assert any(row["title"] == "Original network dendrogram" for row in result["unavailable"])
    assert all(row["metadata"]["completion_scope"] == "observed_numerical_import_and_postprocessing_only" for row in result["artifacts"])
    assert before == {str(path): sha256(path.read_bytes()).hexdigest() for path in src.rglob("*") if path.is_file()}


@pytest.mark.parametrize("error", ["running", "modified_table", "different_arms", "source_relabelled"])
def test_recovered_network_rejects_incomplete_native_or_wrong_source_scope(tmp_path, error):
    from katana_next.reports.network_views import _completed_network_scope
    src, contrasts, *_ = _recovered_data(tmp_path)
    proof_path = src / "native_completion.receipt.json"
    proof = json.loads(proof_path.read_text())
    if error == "running":
        proof["actual_sacct"] = "456|RUNNING|0:0|\n"
    elif error == "source_relabelled":
        proof["source_actual_sacct"] = "123|COMPLETED|0:0|\n"
    elif error == "modified_table":
        with (src / "postprocessing/trait_selected.tsv").open("a") as stream:
            stream.write("changed")
    else:
        contrasts[0]["numerator"]["sample_ids"] = ["C1", "C2", "C3"]
    proof_path.write_text(json.dumps(proof))
    with pytest.raises(ValueError):
        _completed_network_scope(proof_path, [src / "assignments.tsv"], contrasts)
