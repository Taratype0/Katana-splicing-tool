"""Detection universes are evidence, not a reference-gene or significant list."""
from types import SimpleNamespace

import pandas as pd
import pytest

from katana_next.science.enrichment import detected_gene_background, library_overrepresentation
from katana_next.workflows.interpretation import _candidate_background, _enrichment_groups


def test_reference_annotation_and_significance_alone_cannot_define_background():
    for frame in (pd.DataFrame({"gene_id": ["G1", "G2"], "gene_type": ["protein_coding"] * 2}),
                  pd.DataFrame({"gene_id": ["G1"], "significant": [True]})):
        with pytest.raises(ValueError, match="reference gene list is not a measured background"):
            detected_gene_background(frame)


def test_background_keeps_detected_nonsignificant_genes_and_native_statistics():
    frame = pd.DataFrame({"gene_id": ["G1", "G2", "G3"],
                          "detected": [True, True, False], "tested": [True, True, True],
                          "pvalue": [.0001, .98, .0002], "fdr": [.001, .99, .001]})
    original = frame.copy(deep=True)
    kept, audit = detected_gene_background(frame)
    assert kept.gene_id.tolist() == ["G1", "G2"]
    assert kept.pvalue.tolist() == [.0001, .98]
    assert audit.background_source.eq("detected_union").all()
    assert audit.background_eligible.tolist() == [True, True, False]
    pd.testing.assert_frame_equal(frame, original)


def test_assay_scope_does_not_borrow_detection_from_another_assay():
    frame = pd.DataFrame({"gene_id": ["G1", "G2", "G3"],
                          "as_available": [True, False, False],
                          "deg_available": [False, True, True],
                          "deg_detected": [False, True, False]})
    assert detected_gene_background(frame, "das")[0].gene_id.tolist() == ["G1"]
    assert detected_gene_background(frame, "deg")[0].gene_id.tolist() == ["G2"]
    assert detected_gene_background(frame)[0].gene_id.tolist() == ["G1", "G2"]
    with pytest.raises(ValueError, match="No DTU detection"):
        detected_gene_background(frame, "dtu")


def test_serialized_false_and_missing_flags_are_not_truthy():
    frame = pd.DataFrame({"gene_id": ["G1", "G2", "G3", "G4"],
                          "detected": ["True", "False", None, "0"]})
    assert detected_gene_background(frame)[0].gene_id.tolist() == ["G1"]
    with pytest.raises(ValueError, match="Invalid detected/tested flag"):
        detected_gene_background(pd.DataFrame({"detected": ["unknown"]}))


def test_gene_deduplication_and_library_mapping_do_not_use_genome_genes():
    from katana_next.science.enrichment import canonical_gene_mapping, map_gene_sets
    rows = pd.DataFrame({"gene_id": ["G1", "G1", "G2", "GENOME_ONLY"],
                         "gene_symbol": ["A", "A", "B", "C"],
                         "as_available": [True, True, True, False]})
    kept, _ = detected_gene_background(rows, "das")
    ids, aliases, _ = canonical_gene_mapping(kept)
    sets, _ = map_gene_sets({"GO / term": {"A", "B", "C"}}, aliases)
    result = library_overrepresentation({"G1"}, set(ids), sets)
    assert result.background_genes.eq(2).all()
    assert result.term_genes_in_background.eq(2).all()
    assert result.selected_genes.eq(1).all()


def test_selected_cell_type_does_not_admit_other_cell_type_background(monkeypatch):
    import katana_next.workflows.interpretation as module
    frame = pd.DataFrame({"row_id": ["c::A::G1", "c::A::G2", "c::B::G1"],
                          "comparison_id": ["c"] * 3, "cell_type": ["A", "A", "B"],
                          "gene_id": ["G1", "G2", "G1"], "as_available": [True] * 3})
    selected = frame.iloc[[0]].assign(canonical_gene_id="G1")
    monkeypatch.setattr(module, "annotated_gene_source", lambda *_: frame.copy())
    monkeypatch.setattr(module, "_background_artifact", lambda *args: object())
    monkeypatch.setattr(module, "native_background_source_provenance", lambda *args: {"test_only": True})
    kept, audit, _, _ = _candidate_background(SimpleNamespace(), {"selection_id": "s", "foreground_source": "das"}, object(), selected)
    assert kept.cell_type.eq("A").all()
    assert kept.gene_id.tolist() == ["G1", "G2"]
    assert kept.selected_input.tolist() == [True, False]
    assert audit.background_source.eq("das").all()


def test_same_contrast_cell_types_form_distinct_enrichment_scopes():
    frame = pd.DataFrame({"comparison_id": ["c", "c", "c"],
                          "cell_type": ["A", "B", "A"], "gene_id": ["G1", "G1", "G2"]})
    groups = list(_enrichment_groups(frame))
    assert len(groups) == 2
    assert [len(group) for _, group in groups] == [2, 1]
    assert groups[0][0] != groups[1][0]


def test_legacy_background_function_retains_historical_behavior():
    # Native motif kernels call the frozen compatibility function. New strict
    # gene-background policy must not change those numerical inputs silently.
    from katana_next.science.enrichment import enrichment_background
    frame = pd.DataFrame({"gene_id": ["G1", "G2"]})
    kept, _ = enrichment_background(frame)
    assert kept.equals(frame)


def test_failed_native_models_cannot_be_rescued_by_stale_numeric_pvalues():
    frame = pd.DataFrame({"gene_id": ["G1", "G2", "G3"],
                          "model_status": ["fit", "optimizer_failed", "not_computed"],
                          "pvalue": [.9, .01, .001], "measurement_eligible": [True] * 3})
    assert detected_gene_background(frame)[0].gene_id.tolist() == ["G1"]


def test_candidate_adapter_respects_native_as_measurement_without_changing_effects():
    from katana_next.science.enrichment import attach_native_gene_measurement_flags
    candidates = pd.DataFrame({"comparison_id": ["c"] * 3, "gene_id": ["G1", "G2", "G3"],
                               "as_available": [True] * 3, "as_supported": [False, True, False],
                               "best_as_fdr": [.9, .001, .95]})
    events = pd.DataFrame({"comparison_id": ["c"] * 4, "gene_id": ["G1", "G1", "G2", "G3"],
                           "measurement_eligible": [True, True, False, True], "fdr": [.9, .8, .001, float("nan")]})
    original = candidates.copy(deep=True)
    enriched = attach_native_gene_measurement_flags(candidates, events, pd.DataFrame(), pd.DataFrame())
    assert enriched.as_eligible.tolist() == [True, False, False]
    assert detected_gene_background(enriched, "das")[0].gene_id.tolist() == ["G1"]
    pd.testing.assert_frame_equal(enriched[original.columns], original)
    pd.testing.assert_frame_equal(candidates, original)


@pytest.mark.parametrize("workflow", ["prioritize.enrichment", "prioritize.rbp_support"])
def test_api_fabricated_candidate_flags_cannot_replace_native_sources(tmp_path, workflow):
    from test_supervisor_feedback import _api_catalog_fixture
    from katana_next.domain.models import JobStatus
    from katana_next.workflows.base import ArtifactWriter
    from katana_next.workflows.bulk import CandidateEvidenceWorkflow
    runtime, client, context, project, _, comparisons = _api_catalog_fixture(tmp_path)
    contrast = comparisons[0]
    frame = pd.DataFrame({
        "row_id": [f"{contrast}::GENE{i}" for i in range(1, 5)],
        "comparison_id": [contrast] * 4,
        "gene_id": [f"GENE{i}" for i in range(1, 5)],
        "gene_symbol": [f"GENE{i}" for i in range(1, 5)],
        "as_available": [False, True, True, False],
        "as_supported": [False, True, True, False],
        "deg_available": [True, True, False, False],
        "deg_supported": [True, False, False, False],
        "evidence_dimensions": [1, 1, 1, 0],
    })
    writer = ArtifactWriter(context, "prioritize.evidence", "test-detected-background", {"contrast_ids": [contrast]}, CandidateEvidenceWorkflow.version)
    writer.table("Candidate evidence matrix", "candidate_evidence.parquet", frame, entity_type="gene")
    writer.finish()
    started = client.post(f"/api/projects/{project}/workflows/{workflow}/run", json={"parameters": {
        "contrast_ids": [contrast], "foreground_source": "deg"}})
    assert started.status_code == 200, started.text
    finished = runtime.jobs.wait(started.json()["job_id"], timeout=90)
    assert finished.status == JobStatus.FAILED
    assert "traceable complete native parent" in finished.error


def test_candidate_workflow_propagates_actual_imported_measurement_evidence(tmp_path):
    from katana_next.api.app import PlatformRuntime
    from katana_next.domain.models import JobStatus
    from katana_next.full_fixture import create_full_fixture
    manifest = create_full_fixture(tmp_path / "fixture")
    runtime = PlatformRuntime(application_root=tmp_path / "source")
    context = runtime.projects.open_project(manifest)
    project = context.manifest.project_id
    for workflow in ("differential.as", "differential.deg", "differential.dtu", "prioritize.evidence"):
        submitted = runtime.workflows.run(workflow, project, parameters={})
        finished = runtime.jobs.wait(submitted.job_id, timeout=120)
        assert finished.status == JobStatus.COMPLETE, (workflow, finished.error)
    source = next(a for a in context.catalog.list_artifacts(project, analysis_id="prioritize.evidence", kind="table")
                  if a.name == "Candidate evidence matrix")
    candidates = pd.read_parquet(source.path)
    assert {"as_eligible", "deg_tested", "dtu_tested", "dexseq_tested"} <= set(candidates)
    events = pd.read_parquet(next(a.path for a in context.catalog.list_artifacts(project, analysis_id="differential.as", kind="table")
                                 if a.name == "Differential AS events"))
    expected = events.loc[events.measurement_eligible & events.fdr.between(0, 1), ["comparison_id", "gene_id"]].drop_duplicates()
    actual = candidates.loc[candidates.as_eligible, ["comparison_id", "gene_id"]].drop_duplicates()
    assert set(map(tuple, actual.to_numpy())) == set(map(tuple, expected.to_numpy()))
    # Standard complete native imports still have a working per-assay route.
    for workflow, name in (("prioritize.enrichment", "Gene-set enrichment"),
                           ("prioritize.rbp_support", "RBP declared-set enrichment")):
        submitted = runtime.workflows.run(workflow, project, parameters={"foreground_source": "deg"})
        completed = runtime.jobs.wait(submitted.job_id, timeout=120)
        assert completed.status == JobStatus.COMPLETE, completed.error
        result = next(a for a in context.catalog.list_artifacts(project, analysis_id=workflow, kind="table") if a.name == name)
        data = pd.read_parquet(result.path)
        assert data.background_source.eq("deg").all() and data.selected_genes.eq(9).all()
        if workflow == "prioritize.enrichment":
            assert data.tested_genes_before_annotation.eq(30).all()
        else:
            assert data.background_genes.eq(30).all()
        provenance = result.metadata["native_background_provenance"]
        assert provenance["delivery_completeness"] == "not_independently_verified"
        assert all(item["sha256"] for item in provenance["native_inputs"])
    from katana_next.science.enrichment import native_background_source_provenance
    with pytest.raises(ValueError, match="foregrounds only"):
        native_background_source_provenance(context, source.model_copy(update={"metadata": {**source.metadata, "result_scope": "significant_only"}}))
    pd.read_parquet(source.path).head(2).to_parquet(source.path, index=False)
    with pytest.raises(ValueError, match="content changed"):
        native_background_source_provenance(context, source)


def test_numeric_pvalues_and_model_status_are_not_detection_proof():
    hit_only = pd.DataFrame({"gene_id": ["G1"], "significant": [True], "model_status": ["fit"], "pvalue": [.001]})
    with pytest.raises(ValueError, match="reference gene list is not a measured background"):
        detected_gene_background(hit_only)


def test_complete_scope_declaration_cannot_certify_an_arbitrary_hit_table(tmp_path):
    from test_supervisor_feedback import _api_catalog_fixture
    from katana_next.workflows.base import ArtifactWriter
    from katana_next.science.enrichment import native_background_source_provenance
    _, _, context, _, _, _ = _api_catalog_fixture(tmp_path)
    writer = ArtifactWriter(context, "custom.significant_gene_list", "test-declaration-only", {}, "1.0",
                            common_metadata={"result_scope": "complete", "certified": True})
    artifact = writer.table("Declared complete", "hits.parquet", pd.DataFrame({"row_id": ["G1"], "gene_id": ["G1"], "detected": [True]}), entity_type="gene")
    writer.finish()
    with pytest.raises(ValueError, match="foregrounds only"):
        native_background_source_provenance(context, artifact)
