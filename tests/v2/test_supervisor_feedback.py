from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from katana_next.science.cell_embedding import _feature_panel, build_cell_embedding, cell_feature_values
from katana_next.science.binding_enrichment import paired_binding_enrichment
from katana_next.science.coverage_bias import read_gene_body_coverage, summarize_coverage
from katana_next.science.enrichment_views import annotate_categories, enrichment_bubble_plot, read_go_categories
from katana_next.resources import sha256
from katana_next.workflows.registry import WorkflowRegistry


def test_rseqc_wide_profile_has_explicit_sample_mapping(tmp_path):
    path = tmp_path / "coverage.tsv"
    pd.DataFrame([["alignment.bam", *range(1, 101)]], columns=["Percentile", *map(str, range(1, 101))]).to_csv(path, sep="\t", index=False)
    profile = read_gene_body_coverage(path, {"alignment.bam": "s1"})
    curves, metrics = summarize_coverage(profile, pd.DataFrame({"sample_id": ["s1"], "condition": ["WT"]}))
    assert len(curves) == 100
    assert curves.relative_coverage.mean() == pytest.approx(1)
    assert metrics.three_to_five_ratio.iloc[0] > 1
    assert "pass/fail" in metrics.interpretation.iloc[0]


@pytest.mark.parametrize("problem", ["negative", "nonfinite", "duplicate", "truncated", "unknown_sample"])
def test_coverage_rejects_uninterpretable_input(tmp_path, problem):
    frame = pd.DataFrame({"sample_id": "s1", "percentile": np.arange(1, 101), "coverage": 1.})
    if problem == "negative":
        frame.loc[0, "coverage"] = -1
    elif problem == "nonfinite":
        frame.loc[0, "coverage"] = np.inf
    elif problem == "duplicate":
        frame = pd.concat([frame, frame.iloc[:1]])
    elif problem == "truncated":
        frame = frame.iloc[40:60]
    path = tmp_path / "input.tsv"
    frame.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError):
        summarize_coverage(read_gene_body_coverage(path), pd.DataFrame({"sample_id": ["s2" if problem == "unknown_sample" else "s1"]}))


def test_zero_coverage_is_unmeasurable_not_good_qc():
    frame = pd.DataFrame({"sample_id": "s1", "percentile": np.arange(1, 101), "coverage": 0.})
    curves, metrics = summarize_coverage(frame, pd.DataFrame({"sample_id": ["s1"]}))
    assert curves.relative_coverage.isna().all()
    assert metrics.profile_status.iloc[0] == "no_coverage"
    assert np.isnan(metrics.three_to_five_ratio.iloc[0])


def test_go_namespace_annotation_preserves_statistics(tmp_path):
    path = tmp_path / "go.obo"
    path.write_text("[Term]\nid: GO:0000001\nnamespace: biological_process\nalt_id: GO:0000002\n\n[Term]\nid: GO:0000003\nnamespace: molecular_function\nis_obsolete: true\n", encoding="ascii")
    assert read_go_categories(path) == {"GO:0000001": "GO BP", "GO:0000002": "GO BP"}
    frame = pd.DataFrame({"term": ["go / GO:0000002 test", "kegg / pathway"], "overlap_genes": [3, 2], "selected_genes": [10, 10], "pvalue": [.01, .2], "fdr": [.03, .3]})
    manifest = SimpleNamespace(assets=[SimpleNamespace(kind="go_ontology_obo", path=path, metadata={"sha256": sha256(path)}), SimpleNamespace(kind="gene_sets_gmt", asset_id="kegg", metadata={"ontology": "KEGG"})])
    result = annotate_categories(frame, manifest)
    assert result.ontology_category.tolist() == ["GO BP", "KEGG"]
    assert result.gene_ratio.tolist() == [.3, .2]
    pd.testing.assert_frame_equal(frame, result[frame.columns])
    manifest.assets[0].metadata["sha256"] = "wrong"
    with pytest.raises(ValueError, match="SHA256"):
        annotate_categories(frame, manifest)


def test_bubble_category_filter_is_display_only():
    frame = pd.DataFrame({"term": ["go / bp", "go / mf"], "comparison_id": ["A", "A"], "ontology_category": ["GO BP", "GO MF"], "passes_display_overlap": [True, True], "fdr": [.01, .4], "gene_ratio": [.2, .1], "minus_log10_fdr": [2., .398], "overlap_genes": [4, 2]})
    original = frame.copy(deep=True)
    chart = enrichment_bubble_plot(frame, category="GO BP")
    assert chart["data"]["values"][0]["term"] == "go / bp"
    assert len(chart["data"]["values"]) == 1
    pd.testing.assert_frame_equal(frame, original)


def test_unmeasured_event_psi_is_missing_not_zero():
    counts = pd.DataFrame({"cell_id": ["a", "b"], "event_id": ["e", "e"], "inclusion": [0, 5], "exclusion": [0, 5]})
    values = cell_feature_values(counts)
    assert np.isnan(values.value.iloc[0])
    assert values.value.iloc[1] == .5


@pytest.mark.parametrize("value", [-1, .5, np.inf])
def test_cell_counts_must_be_valid_integer_evidence(value):
    with pytest.raises(ValueError, match="integer evidence"):
        cell_feature_values(pd.DataFrame({"cell_id": ["a"], "event_id": ["e"], "inclusion": [value], "exclusion": [4]}))


def test_sparse_expression_variance_includes_implicit_zeros():
    values = cell_feature_values(pd.DataFrame({"cell_id": ["a", "a", "b", "b"], "gene_id": ["g1", "g2", "g2", "g3"], "count": [3, 3, 3, 3]}), mode="expression")
    features, panel = _feature_panel(values, {"a", "b"}, min_fraction=.5, max_features=5, expression=True)
    assert "g1" in panel and "g3" in panel
    assert "g2" not in panel
    assert features.set_index("feature_id").loc["g1", "measured_fraction"] == .5


def test_feature_selection_uses_training_cells_not_heldout_donor():
    values = pd.DataFrame({"cell_id": ["a", "b", "z"], "feature_id": ["e", "e", "e"], "value": [.5, .5, 1.]})
    _, panel = _feature_panel(values, {"a", "b"}, min_fraction=.2, max_features=5)
    assert panel == []


def test_cell_embedding_retains_donors_and_technical_baseline():
    metadata = pd.DataFrame([{"cell_id": f"c{i}", "donor_id": f"d{i // 8}", "cell_type": "A" if i % 2 else "B", "condition": "HC"} for i in range(24)])
    rng = np.random.default_rng(17)
    counts = pd.DataFrame([{"cell_id": row.cell_id, "event_id": f"e{j}", "inclusion": int(rng.integers(2, 19)), "exclusion": 20} for row in metadata.itertuples() for j in range(8)])
    result = build_cell_embedding(counts, metadata, n_neighbors=3)
    assert len(result.coordinates) == 24
    assert result.donor_validation.status.eq("evaluated").all()
    assert result.donor_validation.depth_missingness_baseline.notna().all()
    assert result.coordinates.embedding_source.eq("splicing").all()
    permuted = metadata.copy()
    permuted["cell_type"] = "other_label"
    repeated = build_cell_embedding(counts, permuted, n_neighbors=3, evaluate_cell_types=False)
    np.testing.assert_allclose(result.coordinates[["umap_1", "umap_2"]], repeated.coordinates[["umap_1", "umap_2"]])


def test_registry_integrates_coverage_and_cell_embeddings():
    workflows = {d["workflow_id"]: d for d in WorkflowRegistry().definitions()}
    assert workflows["quality.read_coverage_bias"]["workspace_id"] == "quality"
    assert workflows["single_cell.cell_embedding"]["workspace_id"] == "single-cell"
    assert "ontology_category" in {p["key"] for p in workflows["prioritize.enrichment"]["parameters"]}
    assert "prioritize.binding_enrichment" in workflows


def binding_fixture():
    pairs = pd.DataFrame({"foreground_event": ["f1", "f2", "f3"], "background_event": ["b1", "b2", "b3"], "status": "matched"})
    evidence = pd.DataFrame({"resource_id": "clip1", "comparison_id": "mut_vs_wt", "RBP": "MBNL1", "biosample": "reference_tissue", "accession": "experiment1", "assay": "eCLIP", "event_id": ["f1", "b1", "f2", "b2", "f3", "b3"], "status": ["overlap_in_reference_biosample", "no_overlap_in_resource"] * 3})
    return evidence, pairs


def test_clip_enrichment_counts_independent_matched_pairs():
    evidence, pairs = binding_fixture()
    result = paired_binding_enrichment(evidence, pairs).iloc[0]
    assert result.evaluable_gene_pairs == 3
    assert result.foreground_hit_genes == 3
    assert result.background_hit_genes == 0
    assert result.enrichment_pvalue == .125
    assert result.enrichment_fdr == .125


def test_unevaluable_binding_is_not_a_negative_hit():
    evidence, pairs = binding_fixture()
    evidence.loc[evidence.event_id.isin(["b2", "b3"]), "status"] = "not_evaluable"
    result = paired_binding_enrichment(evidence, pairs).iloc[0]
    assert result.evaluable_gene_pairs == 1
    assert result.unevaluable_gene_pairs == 2
    assert pd.isna(result.enrichment_pvalue) and pd.isna(result.enrichment_fdr)


def test_clip_enrichment_rejects_pooled_comparisons():
    evidence, pairs = binding_fixture()
    evidence.loc[0, "comparison_id"] = "other"
    with pytest.raises(ValueError, match="separately"):
        paired_binding_enrichment(evidence, pairs)


def test_coverage_workflow_runs_through_real_api(tmp_path):
    import json
    from fastapi.testclient import TestClient
    from katana_next.api.app import PlatformRuntime, create_app
    from katana_next.full_fixture import create_full_fixture
    from katana_next.domain.models import JobStatus
    manifest = create_full_fixture(tmp_path / "fixture")
    payload = json.loads(manifest.read_text())
    coverage = tmp_path / "coverage.tsv"
    pd.DataFrame([{"sample_id": sample["sample_id"], "percentile": p, "coverage": 1.} for sample in payload["samples"] for p in range(1, 101)]).to_csv(coverage, sep="\t", index=False)
    payload["assets"].append({"asset_id": "coverage", "kind": "gene_body_coverage", "path": str(coverage)})
    manifest.write_text(json.dumps(payload))
    runtime = PlatformRuntime(application_root=tmp_path / "app")
    client = TestClient(create_app(runtime))
    opened = client.post("/api/projects/open", json={"manifest_path": str(manifest)})
    assert opened.status_code == 200, opened.text
    project = opened.json()["snapshot"]["project_id"]
    started = client.post(f"/api/projects/{project}/workflows/quality.read_coverage_bias/run", json={"parameters": {}, "force": False})
    assert started.status_code == 200, started.text
    finished = runtime.jobs.wait(started.json()["job_id"], timeout=60)
    assert finished.status == JobStatus.COMPLETE, finished
    artifacts = client.get(f"/api/projects/{project}/artifacts").json()
    table = next(a for a in artifacts if a["name"] == "Read coverage bias metrics")
    metrics = pd.read_parquet(table["path"])
    assert metrics.profile_status.eq("measurable").all()
    assert metrics.three_to_five_ratio.eq(1).all()


def test_comparison_selection_does_not_reuse_genes_across_conditions():
    from katana_next.workflows.interpretation import _comparison_selection_mask
    frame = pd.DataFrame({"comparison_id": ["A", "B", "A", "B"], "canonical_gene_id": ["G1", "G1", "G2", "G2"]})
    selected = frame.iloc[[0, 3]]
    assert _comparison_selection_mask(frame, selected).tolist() == [True, False, False, True]


def test_expression_overlay_keeps_missing_libraries_distinct_from_sparse_zeros():
    from katana_next.workflows.cell_embedding import expression_gene_overlay
    expression = pd.DataFrame({"cell_id": ["a", "a", "b"], "gene_id": ["G1", "G2", "G2"], "count": [5, 5, 10]})
    result = expression_gene_overlay(expression, pd.Series(["a", "b", "no_library"]), "G1")
    assert result.iloc[0] == pytest.approx(np.log1p(5000))
    assert result.iloc[1] == 0
    assert pd.isna(result.iloc[2])
    with pytest.raises(ValueError, match="unambiguous"):
        expression_gene_overlay(expression, pd.Series(["a"]), "unknown")


def test_category_extension_preserves_original_membership(tmp_path):
    import json
    from katana_next.resources import add_ontology_categories
    obo, gmt = tmp_path / "go.obo", tmp_path / "library.gmt"
    obo.write_text("[Term]\nid: GO:0000001\nnamespace: biological_process\n", encoding="ascii")
    gmt.write_text("term\tsource\tG1\n", encoding="ascii")
    bundle, extended = tmp_path / "resources.json", tmp_path / "extended.json"
    payload = {"assets": [{"asset_id": "lib", "kind": "gene_sets_gmt", "path": str(gmt), "metadata": {"sha256": sha256(gmt), "semantics": "KEGG declared test library"}}],
               "sources": [{"path": str(obo), "sha256": sha256(obo), "url": "https://purl.obolibrary.org/obo/go/go-basic.obo"}]}
    bundle.write_text(json.dumps(payload))
    original = sha256(bundle)
    add_ontology_categories(bundle, extended)
    result = json.loads(extended.read_text())
    assert sha256(bundle) == original
    assert result["assets"][0]["metadata"]["sha256"] == sha256(gmt)
    assert result["assets"][0]["metadata"]["ontology"] == "KEGG"
    assert result["category_extension"]["membership_changed"] is False
    assert result["assets"][1]["kind"] == "go_ontology_obo"
    with pytest.raises(FileExistsError):
        add_ontology_categories(bundle, extended)


def _api_catalog_fixture(tmp_path, *, clip=False, rs_proof=True):
    import json
    from fastapi.testclient import TestClient
    from katana_next.api.app import PlatformRuntime, create_app
    from katana_next.full_fixture import create_full_fixture
    from katana_next.workflows.base import ArtifactWriter
    manifest = create_full_fixture(tmp_path / "fixture")
    payload = json.loads(manifest.read_text())
    if clip:
        fasta = tmp_path / "genome.fa"
        fasta.write_text(">chr1\n" + "ACGT" * 2500 + "\n", encoding="ascii")
        Path(str(fasta) + ".fai").write_text("chr1\t10000\t6\t10000\t10001\n", encoding="ascii")
        payload["annotation"]["genome_fasta"] = str(fasta)
        payload["metadata"]["species"] = "Mus_musculus"
        peaks = tmp_path / "clip.bed"
        peaks.write_text("chr1\t990\t1010\tpeak1\t1\t+\nchr1\t1490\t1510\tpeak2\t1\t+\n", encoding="ascii")
        payload["assets"].append({"asset_id": "clip", "kind": "rbp_binding_peaks", "path": str(peaks), "metadata": {
            "source_url": "https://example.org/test-fixture-not-biological-evidence", "version": "test1", "sha256": sha256(peaks),
            "species": "Mus_musculus", "genome_build": payload["annotation"]["genome_build"], "coordinate_system": "0-based-half-open",
            "assay": "eCLIP", "biosample": "test_fixture", "target": "RBP1", "accession": "fixture1", "peak_calling": "fixture_only"}})
    manifest.write_text(json.dumps(payload))
    runtime = PlatformRuntime(application_root=tmp_path / "app")
    client = TestClient(create_app(runtime))
    opened = client.post("/api/projects/open", json={"manifest_path": str(manifest)})
    assert opened.status_code == 200, opened.text
    project = opened.json()["snapshot"]["project_id"]
    context = runtime.projects.get(project)
    comparisons = [c.contrast_id for c in context.manifest.contrasts[:2]]
    rows = [{"row_id": f"{comparison}::e{i}", "comparison_id": comparison, "gene_id": f"GENE{i+1}", "gene_symbol": f"GENE{i+1}",
             "event_id": f"e{i}", "event_type": "SE", "chr": "chr1", "strand": "+", "exonStart_0base": 1000 + i*500,
             "exonEnd": 1100 + i*500, "upstreamES": 800 + i*500, "upstreamEE": 900 + i*500,
             "downstreamES": 1200 + i*500, "downstreamEE": 1300 + i*500, "model_status": "fit", "pvalue": .5}
            for comparison in comparisons for i in range(6)]
    frame = pd.DataFrame(rows)
    if clip and rs_proof:
        # This fixture tests an RS-required foreground. Derive its evidence
        # from declared per-arm sample counts with the production QC adapter;
        # do not merely attach an RS label to a generic imported event catalog.
        from katana_next.science.bulk import apply_stable_event_qc, significance_columns
        contrasts = {c.contrast_id: c for c in context.manifest.contrasts}
        for index, row in frame.iterrows():
            contrast = contrasts[row.comparison_id]
            n1, n2 = len(contrast.numerator.sample_ids), len(contrast.denominator.sample_ids)
            frame.loc[index, "numerator_expected_samples"] = n1
            frame.loc[index, "denominator_expected_samples"] = n2
            for suffix, number, psi in (("1", n1, ".6"), ("2", n2, ".4")):
                frame.loc[index, f"IncLevel{suffix}"] = ",".join([psi]*number)
                frame.loc[index, f"IJC_SAMPLE_{suffix}"] = ",".join(["12"]*number)
                frame.loc[index, f"SJC_SAMPLE_{suffix}"] = ",".join(["8"]*number)
            foreground = row.event_id in {"e0", "e1"}
            frame.loc[index, "pvalue"] = .001 if foreground else .5
            frame.loc[index, "fdr"] = .01 if foreground else .5
            frame.loc[index, "dpsi"] = .2 if foreground else .02
            frame.loc[index, "abs_dpsi"] = .2 if foreground else .02
        frame = significance_columns(apply_stable_event_qc(frame, {}), {})
        assert frame.measurement_eligible.all()
        assert frame.stable_significant.eq(frame.event_id.isin({"e0", "e1"})).all()
    writer = ArtifactWriter(context, "test.measured_event_catalog", "frozen", {}, "1.0")
    artifact = writer.table("Complete test event catalog", "events.parquet", frame, entity_type="splicing_event")
    writer.finish()
    chosen = [f"{comparisons[0]}::e0", f"{comparisons[0]}::e1"]
    saved = client.post(f"/api/projects/{project}/selections", json={"name": "Exact comparison A", "entity_type": "splicing_event",
        "source_artifact_id": artifact.artifact_id, "entity_ids": chosen, "filters": {}})
    assert saved.status_code == 200, saved.text
    return runtime, client, context, project, saved.json()["selection_id"], comparisons


def test_unverified_event_catalog_cannot_define_gene_enrichment_background(tmp_path):
    from katana_next.domain.models import JobStatus
    runtime, client, context, project, selection, comparisons = _api_catalog_fixture(tmp_path)
    for workflow, name in [("prioritize.enrichment", "Gene-set enrichment"), ("prioritize.rbp_support", "RBP declared-set enrichment")]:
        started = client.post(f"/api/projects/{project}/workflows/{workflow}/run", json={"parameters": {"selection_id": selection}, "force": False})
        assert started.status_code == 200, started.text
        finished = runtime.jobs.wait(started.json()["job_id"], timeout=90)
        assert finished.status == JobStatus.FAILED
        assert "foregrounds only" in finished.error


def test_clip_enrichment_workflow_runs_exact_matching_and_provenance(tmp_path):
    from katana_next.domain.models import JobStatus
    runtime, client, context, project, selection, comparisons = _api_catalog_fixture(tmp_path, clip=True)
    started = client.post(f"/api/projects/{project}/workflows/prioritize.binding_enrichment/run", json={"parameters": {"selection_id": selection}, "force": False})
    assert started.status_code == 200, started.text
    finished = runtime.jobs.wait(started.json()["job_id"], timeout=90)
    assert finished.status == JobStatus.COMPLETE, finished.error
    artifacts = context.catalog.list_artifacts(project, analysis_id="prioritize.binding_enrichment")
    result = pd.read_parquet(next(a.path for a in artifacts if a.name == "Matched-event CLIP enrichment"))
    assert result.comparison_id.eq(comparisons[0]).all()
    assert result.evaluable_gene_pairs.eq(2).all()
    assert result.foreground_hit_genes.eq(2).all()
    assert result.background_hit_genes.eq(0).all()
    assert result.enrichment_pvalue.eq(.25).all()
    assert next(a for a in artifacts if a.name == "Experimental CLIP provenance").path.is_file()


def test_clip_enrichment_rejects_catalog_without_rs_proof_before_submission(tmp_path):
    runtime, client, context, project, selection, comparisons = _api_catalog_fixture(tmp_path, clip=True, rs_proof=False)
    try:
        source = next(a for a in context.catalog.list_artifacts(project) if a.name == "Complete test event catalog")
        before = source.path.read_bytes()
        started = client.post(f"/api/projects/{project}/workflows/prioritize.binding_enrichment/run",
            json={"parameters": {"selection_id": selection}, "force": False})
        assert started.status_code == 422, started.text
        assert "lacks stable_significant" in started.json()["detail"]
        assert not runtime.jobs.list(project)
        assert source.path.read_bytes() == before
    finally:
        runtime.jobs.shutdown()


def test_cell_embedding_overlay_uses_same_coordinates_and_expression_mode_runs(tmp_path):
    import json
    from fastapi.testclient import TestClient
    from katana_next.api.app import PlatformRuntime, create_app
    from katana_next.domain.models import JobStatus
    from katana_next.full_fixture import create_full_fixture
    manifest = create_full_fixture(tmp_path / "fixture")
    payload = json.loads(manifest.read_text())
    metadata = pd.DataFrame([{"cell_id": f"c{i}", "sample_id": f"d{i // 8}", "donor_id": f"d{i // 8}",
                             "cell_type": "A" if i % 2 else "B", "condition": "test"} for i in range(24)])
    rng = np.random.default_rng(17)
    counts = pd.DataFrame([{"cell_id": c, "event_id": f"e{j}", "inclusion": int(rng.integers(2, 19)), "exclusion": 20}
                           for c in metadata.cell_id for j in range(8)])
    expression = pd.DataFrame([{"cell_id": c, "gene_id": f"G{j}", "count": int(rng.integers(1, 50))}
                               for c in metadata.cell_id for j in range(8)])
    payload["assets"] = [a for a in payload["assets"] if a["kind"] not in {"single_cell_metadata", "single_cell_junction_counts"}]
    for kind, frame in [("single_cell_metadata", metadata), ("single_cell_junction_counts", counts), ("single_cell_expression_counts", expression)]:
        path = tmp_path / f"{kind}.tsv"
        frame.to_csv(path, sep="\t", index=False)
        payload["assets"].append({"asset_id": kind, "kind": kind, "path": str(path)})
    manifest.write_text(json.dumps(payload))
    runtime = PlatformRuntime(application_root=tmp_path / "app")
    client = TestClient(create_app(runtime))
    opened = client.post("/api/projects/open", json={"manifest_path": str(manifest)})
    project = opened.json()["snapshot"]["project_id"]
    context = runtime.projects.get(project)
    results = []
    for parameters in [{}, {"expression_overlay_gene": "G1"}, {"mode": "expression", "expression_overlay_gene": "G1"}]:
        started = client.post(f"/api/projects/{project}/workflows/single_cell.cell_embedding/run", json={"parameters": parameters, "force": False})
        assert started.status_code == 200, started.text
        finished = runtime.jobs.wait(started.json()["job_id"], timeout=90)
        assert finished.status == JobStatus.COMPLETE, finished.error
        artifacts = [context.catalog.get_artifact(a) for a in finished.artifact_ids]
        frame = pd.read_parquet(next(a.path for a in artifacts if a.name == "Cell-level UMAP coordinates"))
        results.append(frame)
        if parameters:
            assert frame.gene_expression_logCP10K.notna().all()
            spec = next(a for a in artifacts if a.name.endswith("gene_expression_logCP10K"))
            assert spec.path.is_file()
    np.testing.assert_array_equal(results[0][["umap_1", "umap_2"]], results[1][["umap_1", "umap_2"]])
    assert results[2].embedding_source.eq("expression").all()


def test_missing_project_id_is_stable_without_modifying_manifest(tmp_path):
    import json
    from katana_next.client import KatanaProject
    from katana_next.domain.models import JobStatus
    from katana_next.full_fixture import create_full_fixture
    manifest = create_full_fixture(tmp_path / "fixture")
    payload = json.loads(manifest.read_text())
    del payload["project_id"]
    manifest.write_text(json.dumps(payload))
    original = sha256(manifest)
    with KatanaProject(manifest) as project:
        identity = project.project_id
        result = project.run("differential.deg")
        assert result.status == JobStatus.COMPLETE
        count = len(project.artifacts())
    with KatanaProject(manifest) as project:
        assert project.project_id == identity
        assert len(project.artifacts()) == count
        assert project.run("differential.deg").status == JobStatus.CACHE_HIT
    assert sha256(manifest) == original
