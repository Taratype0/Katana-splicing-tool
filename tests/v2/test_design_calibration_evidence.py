import hashlib
import json

import numpy as np
import pandas as pd
import pytest

from katana_next.client import KatanaProject
from katana_next.full_fixture import create_full_fixture
from katana_next.domain.models import InputAsset
from katana_next.science.psi_calibration import psi_intervals, length_calibration_grid, audit_rmats_lengths
from katana_next.science.versioned_evidence import protein_feature_support, binding_support


def test_length_corrected_psi_and_missingness():
    result = psi_intervals([20, 0, 0, np.nan], [10, 0, 10, 3], 200, 100)
    assert result["psi"][0] == .5
    assert np.isnan(result["psi"][1]) and np.isnan(result["psi_lower"][1])
    assert np.isnan(result["psi"][3])
    assert result["psi"][2] == 0
    assert np.all(result["psi_lower"][[0, 2]] < result["psi_upper"][[0, 2]])
    for counts in (-1, np.inf, .5):
        with pytest.raises(ValueError):
            psi_intervals(counts, 10, 100, 100)
    with pytest.raises(ValueError):
        psi_intervals(1, 1, 0, 1)


def test_caller_length_parity_and_known_truth():
    frame = pd.DataFrame([dict(IncFormLen=200, SkipFormLen=100, IJC_SAMPLE_1="20,0", SJC_SAMPLE_1="10,0", IncLevel1="0.5,NA",
                               IJC_SAMPLE_2="10,1", SJC_SAMPLE_2="10,0", IncLevel2="0.333,1.0")])
    assert audit_rmats_lengths(frame).caller_parity.all()
    a = length_calibration_grid(events=20, seed=7)
    pd.testing.assert_frame_equal(a, length_calibration_grid(events=20, seed=7))
    assert len(a) == 60
    assert set(a.length_ratio) == {.25, .5, 1, 2, 4}


def test_design_edit_preview_revision_backup_and_cache(tmp_path):
    path = create_full_fixture(tmp_path/"fixture")
    with KatanaProject(path) as project:
        before = project.context.manifest_path.read_bytes()
        original = project.context.snapshot.manifest_fingerprint
        exported = project.editable_design()
        edit = {k: exported[k] for k in ("revision", "samples")}
        edit["reason"] = "Correct declared genotype"
        edit["samples"][0]["attributes"]["genotype"] = "MUT"
        result = project.revise_design(edit)
        assert not result["applied"]
        assert project.context.manifest_path.read_bytes() == before
        assert project.jobs.list() == []
        old_contrasts = project.context.manifest.contrasts
        result = project.revise_design(edit, apply=True)
        assert result["applied"]
        assert project.context.snapshot.manifest_fingerprint != original
        assert project.context.manifest.contrasts == old_contrasts
        from pathlib import Path
        assert Path(result["backup"]).read_bytes() == before
        with pytest.raises(ValueError, match="changed"):
            project.revise_design(edit, apply=True)
        fresh = project.editable_design()
        bad = {"revision": fresh["revision"], "samples": fresh["samples"][:-1], "reason": "bad"}
        with pytest.raises(ValueError, match="exactly once"):
            project.revise_design(bad)


def resource(path, kind, metadata):
    return InputAsset(asset_id=kind, kind=kind, path=path, metadata={"species": "Homo sapiens", "version": "fixture1",
        "source_url": "https://example.org/test-fixture", "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), **metadata})


def test_exact_sequence_protein_evidence_not_invented_domain_loss(tmp_path):
    path = tmp_path/"uniprot.json"
    path.write_text(json.dumps({"primaryAccession": "test", "organism": {"scientificName": "Homo sapiens"}, "sequence": {"value": "MK"},
        "entryAudit": {"sequenceVersion": 1, "entryVersion": 5},
        "features": [{"type": "Domain", "location": {"start": {"value": 1}, "end": {"value": 2}}, "description": "Fixture domain"}]}))
    asset = resource(path, "uniprot_protein_json", {})
    features = pd.DataFrame([{"transcript_id": "T1", "protein_status": "validated_annotated_ORF", "protein_sha256": hashlib.sha256(b"MK").hexdigest()},
                             {"transcript_id": "T2", "protein_status": "validated_annotated_ORF", "protein_sha256": "different"}])
    result, audit = protein_feature_support(features, [asset], species="Homo sapiens")
    assert result.mapping_status.tolist() == ["exact_protein_sequence", "not_evaluable_no_exact_sequence_match"]
    assert result.entry_version.iloc[0] == 5
    with pytest.raises(ValueError, match="species mismatch"):
        protein_feature_support(features, [asset], species="Mus musculus")
    path.write_text("{}")
    with pytest.raises(ValueError, match="SHA256"):
        protein_feature_support(features, [asset], species="Homo sapiens")


def test_binding_uses_exact_build_strand_and_reports_all_selected(tmp_path):
    bed = tmp_path/"peaks.bed"
    bed.write_text("chr1\t99\t110\tp1\t10\t+\nchr1\t105\t120\tp2\t10\t-\n")
    asset = resource(bed, "rbp_binding_peaks", {"genome_build": "GRCh38", "coordinate_system": "0-based-half-open",
        "assay": "eCLIP", "biosample": "K562", "target": "PTBP1", "accession": "test", "peak_calling": "fixture"})
    events = pd.DataFrame([dict(row_id="e1", event_id="e1", chr="chr1", strand="+", event_coordinates="100"),
                           dict(row_id="e2", event_id="e2", chr="chr1", strand="+", event_coordinates="500"),
                           dict(row_id="e3", event_id="e3", chr="chr1", strand="+", event_coordinates=""),
                           dict(row_id="e4", event_id="e4", chr="chr1", strand=".", event_coordinates="100")])
    result, _ = binding_support(events, [asset], species="Homo sapiens", genome_build="GRCh38", window_bp=5)
    assert result.overlapping_peaks.iloc[0] == 1
    assert json.loads(result.overlapping_peak_records.iloc[0])[0]["name"] == "p1"
    assert result.status.tolist() == ["overlap_in_reference_biosample", "no_overlap_in_resource", "not_evaluable", "not_evaluable"]
    assert result.not_evaluable_reason.iloc[3] == "missing_or_unresolved_strand"
    assert pd.isna(result.overlapping_peaks.iloc[3])
    with pytest.raises(ValueError, match="genome build"):
        binding_support(events, [asset], species="Homo sapiens", genome_build="GRCm39")


def test_every_declared_selection_route_has_an_explicit_contract(tmp_path):
    from katana_next.services.selection_routes import SELECTION_TARGETS, compatible_workflows, validate_route
    from katana_next.workflows.registry import WorkflowRegistry
    from katana_next.workflows.base import ArtifactWriter
    from katana_next.domain.models import SelectionSet
    from katana_next.workflows.sources import selection_rows
    with KatanaProject(create_full_fixture(tmp_path/"paths")) as project:
        targets = {w.definition.workflow_id for w in WorkflowRegistry().all() if any(p.key == "selection_id" for p in w.definition.parameters)}
        assert targets == SELECTION_TARGETS
        cases = [
            ("gene", "prioritize.evidence", {"gene_id": "G1", "gene_symbol": "Gene"}, {"prioritize.enrichment", "prioritize.rbp_support", "prioritize.selected_gene_views", "validate.perturbation"}),
            ("splice_event", "differential.as", {"gene_id": "G1", "event_id": "e1", "event_type": "SE", "chr": "chr1", "strand": "+", "comparison_id": "A", "IncLevel1": ".2", "IncLevel2": ".3", "IJC_SAMPLE_1": "1,2", "SJC_SAMPLE_1": "3,4", "IJC_SAMPLE_2": "5,6", "SJC_SAMPLE_2": "7,8", "IncFormLen": 200, "SkipFormLen": 100}, SELECTION_TARGETS-{"single_cell.sashimi"}),
            ("cell_splice_event", "single_cell.quantify", {"event_id": "chr1:100-200:.", "cluster_id": "c1"}, {"single_cell.sashimi"}),
            ("cosplice_module", "programs.cosplice", {"module": "M1"}, set()),
            ("module_member", "programs.diagnostics", {"event_key": "e1", "module": "M1", "gene_id": "G1"}, {"prioritize.enrichment", "prioritize.rbp_support", "prioritize.selected_gene_views", "validate.perturbation"}),
        ]
        for index, (entity, workflow, columns, allowed) in enumerate(cases):
            writer = ArtifactWriter(project.context, workflow, f"routing-{index}", {}, WorkflowRegistry().get(workflow).version)
            frame = pd.DataFrame([{**columns, "row_id": "wanted"}, {**columns, "row_id": "not_selected"}])
            a = writer.table(f"Routing {index}", "rows.parquet", frame, entity_type=entity)
            writer.finish()
            assert set(compatible_workflows(a)) == allowed
            selection = SelectionSet(project_id=project.project_id, name="frozen", entity_type=entity, source_artifact_id=a.artifact_id, entity_ids=("wanted",))
            project.projects.save_selection(selection)
            parameters = {"selection_id": selection.selection_id}
            assert selection_rows(project.context, parameters, a, frame).row_id.tolist() == ["wanted"]
            for target in targets:
                if target in allowed:
                    validate_route(project.context, target, parameters)
                else:
                    with pytest.raises(ValueError, match="cannot be sent"):
                        validate_route(project.context, target, parameters)
            bad = selection.model_copy(update={"entity_ids": ("absent",)})
            with pytest.raises(ValueError, match="missing"):
                project.projects.save_selection(bad)
        assert project.jobs.list() == []


def test_indexed_python_bam_retains_long_spliced_overhang(tmp_path):
    import bamnostic
    from katana_next.sashimi_fixture import write_indexed_bam
    from katana_next.science.single_cell_sashimi import indexed_reads
    bam = tmp_path/"long.bam"
    write_indexed_bam(bam)
    with bamnostic.AlignmentFile(str(bam), "rb") as handle:
        reads = list(indexed_reads(handle, "chr1", 210, 230))
    assert len(reads) == 4
    assert all(r.reference_start == 100 for r in reads)


def test_saved_api_design_refreshes_ui_payload_and_engine_metadata(tmp_path):
    from fastapi.testclient import TestClient
    from katana_next.api.app import PlatformRuntime, create_app
    from katana_next.services.design_editor import sync_cosplice_design
    from types import SimpleNamespace
    runtime = PlatformRuntime(application_root=tmp_path/"app")
    try:
        client = TestClient(create_app(runtime))
        opened = client.post("/api/projects/open", json={"manifest_path": str(create_full_fixture(tmp_path/"fixture"))}).json()
        project_id = opened["snapshot"]["project_id"]
        url = f"/api/projects/{project_id}/design/edit"
        exported = client.get(url).json()
        edit = {k: exported[k] for k in ("revision", "samples")}
        edit["reason"] = "Verify saved metadata reaches the engine and UI"
        edit["samples"][0]["condition"] = "revised"
        edit["samples"][0]["batch"] = "  "
        result = client.post(url+"?apply=true", json=edit)
        assert result.status_code == 200, result.text
        payload = result.json()["project"]
        assert payload["manifest"]["samples"][0]["condition"] == "revised"
        assert payload["snapshot"]["manifest_fingerprint"] != opened["snapshot"]["manifest_fingerprint"]
        assert runtime.jobs.list() == []
        context = runtime.projects.get(project_id)
        metadata = tmp_path/"metadata.tsv"
        pd.DataFrame([{"sample_id": s.sample_id, "condition": "old", "bam": "original.bam"} for s in context.manifest.samples]).to_csv(metadata, sep="\t", index=False)
        config = SimpleNamespace(sample_metadata=metadata)
        sync_cosplice_design(config, context.manifest, tmp_path/"engine")
        aligned = pd.read_csv(config.sample_metadata, sep="\t")
        assert aligned.condition.iloc[0] == "revised"
        assert aligned.bam.eq("original.bam").all()
        assert pd.isna(aligned.batch.iloc[0])
    finally:
        runtime.jobs.shutdown()
