import numpy as np
import pandas as pd
import pytest

from katana_next.science.blind_clustering import cluster_latent, reveal_partition
from katana_next.science.coverage_bias import read_gene_body_coverage, summarize_coverage
from katana_next.workflows.registry import WorkflowRegistry


def latent_fixture():
    rng = np.random.default_rng(17)
    x = np.concatenate([rng.normal(center, .2, (20, 4)) for center in (-5, 0, 5)])
    frame = pd.DataFrame(x, columns=[f"pc_{i + 1}" for i in range(4)])
    frame.insert(0, "cell_id", [f"c{i:03d}" for i in range(60)])
    return frame


def test_clustering_is_label_blind_order_invariant_and_seeded():
    frame = latent_fixture()
    original, metrics = cluster_latent(frame)
    frame["cell_type"] = "ignored"
    frame["donor_id"] = np.random.default_rng(5).integers(0, 5, len(frame))
    permuted, repeated = cluster_latent(frame.sample(frac=1, random_state=3))
    pd.testing.assert_frame_equal(original, permuted)
    assert metrics == repeated
    assert metrics["clusters"] == 3 and not metrics["labels_used"]


@pytest.mark.parametrize("problem", ["duplicate", "missing", "infinite", "constant", "2D_umap"])
def test_bad_latent_inputs_rejected(problem):
    frame = latent_fixture()
    if problem == "duplicate":
        frame.loc[1, "cell_id"] = frame.loc[0, "cell_id"]
    elif problem == "missing":
        frame.loc[1, "cell_id"] = None
    elif problem == "infinite":
        frame.loc[0, "pc_1"] = np.inf
    elif problem == "constant":
        frame.loc[:, [f"pc_{i + 1}" for i in range(4)]] = 0
    else:
        frame = frame.rename(columns={f"pc_{i + 1}": f"umap_{i + 1}" for i in range(4)})
    with pytest.raises(ValueError):
        cluster_latent(frame)


@pytest.mark.parametrize("resolution", [0, -1, np.nan, np.inf])
def test_resolution_is_validated(resolution):
    with pytest.raises(ValueError):
        cluster_latent(latent_fixture(), resolution=resolution)


def test_reveal_is_separate_and_requires_exact_frozen_cells():
    frame = latent_fixture()
    assignments, _ = cluster_latent(frame)
    annotations = pd.DataFrame({"cell_id": frame.cell_id, "cell_type": np.repeat(["A", "B", "C"], 20)})
    _, agreement = reveal_partition(assignments, annotations)
    assert agreement["adjusted_rand_index"] == 1
    with pytest.raises(ValueError, match="exactly"):
        reveal_partition(assignments, annotations.iloc[:-1])


def test_resolution_and_umap_spacing_are_distinct_registered_controls():
    definition = next(d for d in WorkflowRegistry().definitions() if d["workflow_id"] == "single_cell.cell_embedding")
    keys = {p["key"] for p in definition["parameters"]}
    assert {"min_dist", "pca_components", "cluster_cells", "clustering_resolution", "clustering_neighbors"} <= keys


def test_qualimap_requires_declared_identity_and_preserves_metric_semantics(tmp_path):
    path = tmp_path / "coverage.txt"
    pd.DataFrame({"#Transcript position": np.arange(101), "Transcript coverage profile": np.arange(1, 102)}).to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="sample_id"):
        read_gene_body_coverage(path)
    frame = read_gene_body_coverage(path, sample_id="s1")
    curves, metrics = summarize_coverage(frame, pd.DataFrame({"sample_id": ["s1"]}))
    assert len(curves) == 101
    assert np.isnan(metrics.three_to_five_ratio.iloc[0])
    assert metrics.position_end_to_start_ratio.iloc[0] > 1
    assert "not native Qualimap" in metrics.end_ratio_semantics.iloc[0]


def test_native_terminal_summary_preserves_reported_ratio(tmp_path):
    from katana_next.science.coverage_bias import read_qualimap_summary
    path = tmp_path / "rnaseq_qc_results.txt"
    path.write_text("5' bias = 0.86\n3' bias = 0.71\n5'-3' bias = 1.07\n5'-3' bias region size = 100\n5'-3' bias number of top transcripts = 1000\nprotocol = strand-specific-reverse\ngff file = reference.gtf\n")
    record = read_qualimap_summary(path, sample_id="s1")
    assert record["qualimap_five_to_three_bias"] == 1.07
    assert record["qualimap_five_to_three_bias"] != record["qualimap_five_prime_bias"] / record["qualimap_three_prime_bias"]
    assert record["bias_window_bp"] == 100
    with pytest.raises(ValueError, match="sample_id"):
        read_qualimap_summary(path, sample_id=None)
    path.write_text("5' bias = 0.86")
    with pytest.raises(ValueError, match="Missing"):
        read_qualimap_summary(path, sample_id="s1")


@pytest.mark.parametrize("value", [-.1, 1.1, np.nan])
def test_umap_spacing_is_validated(value):
    from katana_next.science.cell_embedding import build_cell_embedding
    metadata = pd.DataFrame({"cell_id": ["c1"], "donor_id": ["d1"], "cell_type": ["A"], "condition": ["HC"]})
    with pytest.raises(ValueError, match="min_dist"):
        build_cell_embedding(pd.DataFrame(), metadata, min_dist=value)


def test_clustering_api_outputs_and_unchanged_umap(tmp_path):
    import json
    from fastapi.testclient import TestClient
    from katana_next.api.app import PlatformRuntime, create_app
    from katana_next.domain.models import JobStatus
    from katana_next.full_fixture import create_full_fixture
    manifest = create_full_fixture(tmp_path / "fixture")
    payload = json.loads(manifest.read_text())
    metadata = pd.DataFrame([{"cell_id": f"c{i}", "sample_id": f"d{i // 8}", "donor_id": f"d{i // 8}",
        "cell_type": "A" if i % 2 else "B", "condition": "HC"} for i in range(24)])
    rng = np.random.default_rng(17)
    counts = pd.DataFrame([{"cell_id": c, "event_id": f"e{j}", "inclusion": int(rng.integers(2, 19)), "exclusion": 20}
        for c in metadata.cell_id for j in range(8)])
    payload["assets"] = [a for a in payload["assets"] if a["kind"] not in {"single_cell_metadata", "single_cell_junction_counts"}]
    for kind, frame in (("single_cell_metadata", metadata), ("single_cell_junction_counts", counts)):
        path = tmp_path / f"{kind}.tsv"
        frame.to_csv(path, sep="\t", index=False)
        payload["assets"].append({"asset_id": kind, "kind": kind, "path": str(path)})
    manifest.write_text(json.dumps(payload))
    runtime = PlatformRuntime(application_root=tmp_path / "app")
    client = TestClient(create_app(runtime))
    project = client.post("/api/projects/open", json={"manifest_path": str(manifest)}).json()["snapshot"]["project_id"]
    frames = []
    for enabled in (False, True):
        response = client.post(f"/api/projects/{project}/workflows/single_cell.cell_embedding/run", json={"parameters": {"cluster_cells": enabled, "evaluate_cell_types": False}, "force": False})
        assert response.status_code == 200, response.text
        result = runtime.jobs.wait(response.json()["job_id"], timeout=90)
        assert result.status == JobStatus.COMPLETE, result.error
        context = runtime.projects.get(project)
        artifacts = [context.catalog.get_artifact(a) for a in result.artifact_ids]
        frames.append(pd.read_parquet(next(a.path for a in artifacts if a.name == "Cell-level UMAP coordinates")))
        if enabled:
            assert frames[-1].cluster_id.notna().all()
            diagnostics = pd.read_parquet(next(a.path for a in artifacts if a.name == "Label-free Leiden diagnostics"))
            assert diagnostics.labels_used.eq(False).all()
            assert not any(a.name == "Post-clustering annotation agreement" for a in artifacts)
    np.testing.assert_array_equal(frames[0][["umap_1", "umap_2"]], frames[1][["umap_1", "umap_2"]])
