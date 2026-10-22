from pathlib import Path
from types import SimpleNamespace
import json
import numpy as np
import pandas as pd
import pytest

from katana_next.resources import go_gene_sets
from katana_next.science.enrichment import declared_gene_sets
from katana_next.domain.models import InputAsset, AnnotationReference
from katana_next.workflows.intersection_plots import event_membership_sets
from katana_next.science.isoform import load_quant_support, isoform_contrast_summary, dominant_isoform_switches
from katana_next.full_fixture import create_full_fixture
from katana_next.api.app import PlatformRuntime


def test_go_not_annotations_and_parent_propagation(tmp_path):
    obo = tmp_path / "go.obo"
    obo.write_text('[Term]\nid: GO:0003723\nname: RNA binding\n\n[Term]\nid: GO:0000001\nname: child\nis_a: GO:0003723 ! RNA binding\n')
    gaf = tmp_path / "mouse.gaf"
    rows = []
    for symbol, qualifier, taxon in (("Good", "enables", "10090"), ("Negative", "NOT|enables", "10090"), ("Human", "enables", "9606"), ("Modern", "enables", "10090")):
        fields = ["MGI", "MGI:1", symbol, qualifier, "GO:0000001", "PMID:1", "IDA", "", "F", "", "", "protein", f"taxon:{taxon}", "20260909", "MGI", "", ""]
        rows.append("\t".join(fields))
    rows[-1] = rows[-1].replace("taxon:10090", "NCBITaxon:10090")
    gaf.write_text("\n".join(rows))
    sets, audit = go_gene_sets(obo, gaf, "10090")
    assert sets["GO:0003723 RNA binding"] == {"Good", "Modern"}
    assert audit["excluded_annotations"] == 2


def test_all_libraries_loaded_and_species_checked(tmp_path):
    files = [tmp_path / "one.gmt", tmp_path / "two.gmt"]
    for file in files:
        file.write_text("same_term\tcustom\tGene1\n")
    assets = [InputAsset(asset_id=str(i), kind="gene_sets_gmt", path=p, metadata={"species": "Mus_musculus"}) for i, p in enumerate(files)]
    manifest = SimpleNamespace(assets=assets, annotation=AnnotationReference(genome_build="GRCm39"))
    sets, audit = declared_gene_sets(manifest, "gene_sets_gmt")
    assert len(sets) == len(audit) == 2
    manifest.annotation = AnnotationReference(genome_build="GRCh38")
    with pytest.raises(ValueError, match="species"):
        declared_gene_sets(manifest, "gene_sets_gmt")


def test_opposite_direction_not_same_direction_member():
    events = pd.DataFrame({"comparison_id": ["A", "B", "A", "B"], "event_id": ["x", "x", "z", "z"],
                           "dpsi": [.2, -.2, 0, 0], "stable_significant": [True] * 4})
    exact, _ = event_membership_sets(events)
    directional, detail = event_membership_sets(events, "event_and_direction")
    assert len(exact["A"] & exact["B"]) == 2
    assert not directional["A"] & directional["B"]
    assert len(detail) == 2
    assert set(detail.event_id) == {"x"}


@pytest.mark.parametrize("problem", ["missing", "nan", "negative", "schema", "duplicate", "catalog"])
def test_quant_invalid_never_zero_or_silently_skipped(tmp_path, problem):
    manifest = create_full_fixture(tmp_path / "fixture")
    runtime = PlatformRuntime(application_root=tmp_path / "app")
    context = runtime.projects.open_project(manifest)
    sample = context.manifest.samples[-1]
    path = sample.assets["quant_sf"]
    frame = pd.read_csv(path, sep="\t")
    if problem == "missing":
        path.unlink()
    else:
        if problem == "nan":
            frame.loc[0, "TPM"] = np.nan
        elif problem == "negative":
            frame.loc[0, "NumReads"] = -1
        elif problem == "schema":
            frame = frame.drop(columns="TPM")
        elif problem == "duplicate":
            frame = pd.concat([frame, frame.iloc[:1]])
        elif problem == "catalog":
            frame = frame.iloc[1:]
        frame.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError):
        load_quant_support(context, tmp_path / "map.parquet")


def test_unmeasured_transcript_mean_is_not_zero():
    arm = lambda samples: SimpleNamespace(sample_ids=samples, label=str(samples))
    contrast = SimpleNamespace(contrast_id="A", numerator=arm(["n"]), denominator=arm(["d"]))
    context = SimpleNamespace(manifest=SimpleNamespace(contrasts=[contrast]))
    frame = pd.DataFrame({"sample_id": ["n", "d"], "transcript_id": ["a", "b"], "gene_id": ["g", "g"], "gene_symbol": ["G", "G"], "tpm": [1., 2.]})
    assert isoform_contrast_summary(context, frame).delta_tpm.isna().all()


def test_zero_expression_is_not_isoform_switch():
    arm = lambda samples: SimpleNamespace(sample_ids=samples, label=str(samples))
    context = SimpleNamespace(manifest=SimpleNamespace(contrasts=[SimpleNamespace(contrast_id="A", numerator=arm(["n"]), denominator=arm(["d"]))]))
    frame = pd.DataFrame({"sample_id": ["n", "d"], "transcript_id": ["a", "b"], "gene_id": ["g", "g"], "gene_symbol": ["G", "G"], "tpm": [0., 0.]})
    assert dominant_isoform_switches(context, frame).empty


def test_new_workflows_have_explicit_selection_and_cache_contract():
    from katana_next.workflows.registry import WorkflowRegistry
    registry = WorkflowRegistry()
    for name in ("prioritize.selected_gene_views", "prioritize.event_motifs"):
        definition = registry.get(name).definition
        assert "selection_id" in {p.key for p in definition.parameters}
        assert definition.limitations
    assert "membership_mode" in {p.key for p in registry.get("differential.cross_comparison").definition.parameters}


def test_official_r_tom_preserves_reliability_identity():
    from katana_cosplice.benchmark import _rscript_with_wgcna
    from katana_cosplice.r_wgcna import exact_wgcna_labels
    executable, diagnostic = _rscript_with_wgcna()
    if executable is None:
        pytest.skip(diagnostic)
    rng = np.random.default_rng(94)
    latent = rng.normal(size=(3, 24))
    matrix = pd.DataFrame(np.repeat(latent, 20, axis=0) + rng.normal(0, .08, (60, 24)),
                          index=[f"event-{i}" for i in range(60)])
    options = dict(soft_power=4, min_module_size=10, merge_cut_height=.25,
                   deep_split=2, pam_respects_dendro=False, rscript=executable, return_tom=True)
    labels, plain = exact_wgcna_labels(matrix, reliability_power=0, **options)
    _, ones = exact_wgcna_labels(matrix, event_reliability=pd.Series(1., index=matrix.index),
                                 reliability_power=1, **options)
    reliability = pd.Series(1., index=matrix.index)
    reliability.iloc[:20] = .05
    _, weighted = exact_wgcna_labels(matrix, event_reliability=reliability.iloc[::-1],
                                     reliability_power=1, **options)
    assert np.allclose(plain, ones)
    assert np.isfinite(weighted).all() and np.allclose(weighted, weighted.T)
    assert weighted[:20, :20][~np.eye(20, dtype=bool)].mean() < plain[:20, :20][~np.eye(20, dtype=bool)].mean()
    assert len(set(labels) - {"grey"}) >= 2


def test_resource_ui_api_preserves_old_project_and_starts_no_jobs(tmp_path):
    from fastapi.testclient import TestClient
    from katana_next.api.app import create_app
    from katana_next.resources import sha256
    source = create_full_fixture(tmp_path / "old")
    original = source.read_bytes()
    gmt = tmp_path / "resources.gmt"
    gmt.write_text("custom\tdescription\tGENE1\n")
    bundle = tmp_path / "resources.json"
    bundle.write_text(json.dumps({"species": "user_defined", "assets": [{"asset_id": "custom",
                      "kind": "gene_sets_gmt", "path": str(gmt), "metadata": {"sha256": sha256(gmt)}}]}))
    runtime = PlatformRuntime(application_root=tmp_path / "app")
    client = TestClient(create_app(runtime))
    payload = {"manifest_path": str(source), "resource_bundle_path": str(bundle), "output_root": str(tmp_path / "new")}
    response = client.post("/api/projects/configure-resources", json=payload)
    assert response.status_code == 200, response.text
    assert source.read_bytes() == original
    assert not runtime.jobs.list(response.json()["snapshot"]["project_id"])
    assert client.post("/api/projects/configure-resources", json=payload).status_code == 422
    gmt.write_text("tampered\tdescription\tGENE1\n")
    payload["output_root"] = str(tmp_path / "tampered")
    assert client.post("/api/projects/configure-resources", json=payload).status_code == 422
    assert not (tmp_path / "tampered" / "katana.project.json").exists()


def test_cosplice_uses_all_project_interpretation_resources(tmp_path):
    from katana_cosplice.config import load_config
    from katana_next.workflows.programs import configure_interpretation_resources
    runtime = PlatformRuntime(application_root=tmp_path / "app")
    context = runtime.projects.open_project(create_full_fixture(tmp_path / "fixture"))
    config = load_config(next(a.path for a in context.manifest.assets if a.kind == "cosplice_config"))
    root = tmp_path / "engine"
    root.mkdir()
    configure_interpretation_resources(config, context.manifest, root)
    assert config.interpretation.gene_sets_gmt == root / "declared_gene_sets.gmt"
    assert config.rbp_gene_set == root / "declared_rbp_membership.tsv"
    assert config.interpretation.gene_sets_gmt.stat().st_size > 0
    assert {"gene_symbol", "gene_id", "gene_set"}.issubset(pd.read_csv(config.rbp_gene_set, sep="\t"))


def test_pathway_background_does_not_enrich_annotation_coverage():
    from katana_next.science.enrichment import library_overrepresentation
    tested = {f"g{i}" for i in range(100)}
    annotated = {f"g{i}" for i in range(30)}
    result = library_overrepresentation({"g1", "g2"}, tested, {"GO / root": annotated})
    assert result.loc[0, "background_genes"] == 30
    assert result.loc[0, "pvalue"] == 1.
    assert result.loc[0, "tested_genes_before_annotation"] == 100


def test_superseded_workflow_is_not_a_valid_downstream_source():
    from katana_next.workflows.sources import source_complete
    artifact = SimpleNamespace(analysis_id="prioritize.enrichment", metadata={
        "artifact_contract_version": "3.0-correctness", "workflow_version": "2.1.0"})
    assert not source_complete(None, artifact)
