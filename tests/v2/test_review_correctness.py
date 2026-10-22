from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from katana_next.science.bulk import cross_comparison_summary, build_candidate_evidence
from katana_next.science.enrichment import canonical_gene_mapping, map_gene_sets, overrepresentation
from katana_next.science.experimental_design import observation_units, bootstrap_blocks, donor_covariates
from katana_next.science.single_cell_programs import load_declared_cell_embedding
from katana_next.science.single_cell_benchmark import _shared_program_matrix
from katana_next.api.app import PlatformRuntime
from katana_next.full_fixture import create_full_fixture
from katana_next.domain.models import JobStatus, SelectionSet
from katana_next.workflows.sources import resolve_table
from katana_next.workflows.interpretation import _selection_genes


def test_missing_psi_preserves_replicate_and_read_alignment(tmp_path):
    from katana_next.science.bulk import load_bulk_events, parse_number_list, apply_stable_event_qc
    values = parse_number_list("0.1,NA,0.9")
    assert len(values) == 3 and np.isnan(values[1]) and values[2] == .9
    manifest = create_full_fixture(tmp_path / "fixture")
    runtime = PlatformRuntime(application_root=tmp_path / "app")
    context = runtime.projects.open_project(manifest)
    asset = next(a for a in context.manifest.assets if a.kind == "rmats_result")
    frame = pd.read_csv(asset.path, sep="\t")
    frame.loc[0, "IncLevel1"] = "0.1,NA,0.9,NA,NA"
    frame.loc[0, "IJC_SAMPLE_1"] = "10,200,30,0,0"
    frame.loc[0, "SJC_SAMPLE_1"] = "1,2,3,0,0"
    frame.to_csv(asset.path, sep="\t", index=False)
    data = load_bulk_events(context, {"contrast_ids": [asset.metadata["contrast_id"]]})
    event = data.events.loc[(data.events.source_path == str(asset.path)) & data.events.source_row.eq(0)].iloc[0]
    contrast = next(c for c in context.manifest.contrasts if c.contrast_id == asset.metadata["contrast_id"])
    arm = "denominator" if asset.metadata.get("reverse_direction") else "numerator"
    names = asset.metadata.get(f"{arm}_sample_ids") or getattr(contrast, arm).sample_ids
    assert data.psi_matrix.loc[event.event_id, names[0]] == .1
    assert np.isnan(data.psi_matrix.loc[event.event_id, names[1]])
    assert data.psi_matrix.loc[event.event_id, names[2]] == .9
    assert data.read_matrix.loc[event.event_id, names[2]] == 33
    qc = apply_stable_event_qc(data.events, {"stability_profile": "custom", "min_reads": 15})
    row = qc.loc[qc.row_id.eq(event.row_id)].iloc[0]
    assert row[f"{arm}_supported_replicates"] == 1


def test_replicate_mapping_rejects_extra_values():
    from katana_next.science.bulk import _sample_names
    asset = SimpleNamespace(metadata={})
    contrast = SimpleNamespace(contrast_id="test", numerator=SimpleNamespace(sample_ids=["s1", "s2"]))
    with pytest.raises(ValueError, match="mapped samples"):
        _sample_names(asset, contrast, "numerator", 3)


def test_static_facet_width_is_per_panel():
    from katana_next.workflows.figure_exports import _static_layout
    spec = {"width": "container", "encoding": {"column": {"field": "metric"}},
            "data": {"values": [{"metric": str(i)} for i in range(4)]}}
    _static_layout(spec)
    assert 200 <= spec["width"] <= 300


def test_large_pair_plot_preserves_all_counts_in_bounded_density():
    from katana_next.workflows.support_plots import pairwise_display
    frame = pd.DataFrame({"value_a": np.linspace(-10, 10, 20001), "value_b": np.zeros(20001)})
    plotted, mark, encoding, mode = pairwise_display(frame, "transcript_id")
    assert sum(row["count"] for row in plotted) == len(frame)
    assert len(plotted) <= 3600 and mark == "rect" and "density" in mode


def test_custom_delimited_gene_sets(tmp_path):
    from katana_next.science.enrichment import read_gene_sets
    path = tmp_path / "custom_rbp.tsv"
    pd.DataFrame({"set": ["custom", "custom"], "gene": ["A", "B"]}).to_csv(path, sep="\t", index=False)
    assert read_gene_sets(path) == {"custom": {"A", "B"}}


def test_catalog_lock_returns_actionable_conflict(tmp_path, monkeypatch):
    from fastapi.testclient import TestClient
    from katana_next.api.app import create_app
    from katana_next.infrastructure.catalog import ArtifactCatalog, CatalogBusyError
    import duckdb
    catalog = ArtifactCatalog(tmp_path / "results")
    def locked(_path):
        raise duckdb.IOException("The process cannot access the file because it is being used by another process.")
    monkeypatch.setattr(duckdb, "connect", locked)
    with pytest.raises(CatalogBusyError):
        catalog._connect()
    runtime = PlatformRuntime(application_root=tmp_path / "app")
    monkeypatch.setattr(runtime.projects, "open_project", lambda _: catalog._connect())
    with TestClient(create_app(runtime)) as client:
        response = client.post("/api/projects/open", json={"manifest_path": str(tmp_path / "project.json")})
    assert response.status_code == 409 and "Let that operation finish" in response.json()["detail"]


def test_complete_mut_denominator_and_missing_counterpart():
    rows = [dict(comparison_id="MUT", event_id=e, stable_significant=True, dpsi=.2) for e in ("e1", "e2", "e3")]
    rows += [dict(comparison_id="INH", event_id=e, stable_significant=s, dpsi=.3) for e, s in [("e1", True), ("e2", False)]]
    summary, members = cross_comparison_summary(pd.DataFrame(rows))
    row = summary.iloc[0]
    assert row.significant_a == 3
    assert row.containment_a_in_b == pytest.approx(1 / 3)
    assert row.testable_containment_a_in_b == .5
    assert row.a_significant_missing_in_b == 1
    assert members.set_index("event_id").loc["e3", "counterpart_status"] == "not tested in B"


def test_one_gene_one_member_and_bh_family_before_display_filter():
    frame = pd.DataFrame({"gene_id": ["ENSG0001.1", "ENSG0001", "ENSG0002"], "gene_symbol": ["A", "A", "B"]})
    ids, aliases, _ = canonical_gene_mapping(frame)
    sets, audit = map_gene_sets({"hit": {"A", "ENSG0001.5"}, "zero": {"B"}}, aliases)
    assert len(set(ids)) == 2
    assert sets["hit"] == {"ENSG0001"}
    assert audit.iloc[0].mapped_unique_genes == 1
    table = overrepresentation({"ENSG0001"}, set(ids), sets, minimum_overlap=1).set_index("term")
    assert table.loc["hit", "selected_genes"] == 1
    assert table.loc["hit", "background_genes"] == 2
    assert table.loc["hit", "pvalue"] == .5
    assert table.loc["hit", "fdr"] == 1
    assert table.loc["zero", "pvalue"] == 1
    assert not table.loc["zero", "passes_display_overlap"]
    changed_display = overrepresentation({"ENSG0001"}, set(ids), sets, minimum_overlap=2).set_index("term")
    pd.testing.assert_series_equal(table.fdr, changed_display.fdr)


def test_ambiguous_symbols_are_not_silently_mapped():
    frame = pd.DataFrame({"gene_id": ["ID1", "ID2"], "gene_symbol": ["shared", "shared"]})
    _, aliases, audit = canonical_gene_mapping(frame)
    assert "shared" not in aliases
    assert audit.set_index("alias").loc["shared", "mapping_status"] == "ambiguous"


def test_paired_donor_observations_preserved_and_resampled_as_blocks():
    frame = pd.DataFrame({"donor_id": ["D1", "D1", "D2", "D2"], "condition": ["WT", "MUT"] * 2,
                          "event_id": ["e1"] * 4, "psi_posterior": [.1, .9, .2, .8],
                          "total": [10, 90, 20, 80], "measurement_se": [.1] * 4, "batch": ["b"] * 4})
    units = observation_units(frame)
    assert units.observation_id.nunique() == 4
    matrix, _ = _shared_program_matrix(units, pd.Index(units.observation_id), min_unit_reads=1,
                                      min_unit_fraction=0, min_psi_variance=0, max_events=10)
    assert matrix.shape == (1, 4)
    assert matrix.iloc[0].tolist() == [.1, .9, .2, .8]
    for seed in range(10):
        drawn = units.iloc[bootstrap_blocks(units, np.random.default_rng(seed))]
        assert drawn.groupby("donor_id")["condition"].apply(lambda s: (s == "WT").sum() == (s == "MUT").sum()).all()
    covariates, paired, batch = donor_covariates(units)
    assert paired and not batch and covariates.shape == (4, 1)
    with pytest.raises(ValueError, match="Duplicate"):
        observation_units(pd.concat([frame, frame.iloc[:1]]))


def test_manifold_keeps_all_declared_dimensions():
    frame = pd.DataFrame({"cell_id": ["c1", "c2"], "expression_pca_1": [1., 2.],
                          "expression_pca_2": [0., 0.], "expression_pca_3": [-5., 5.]})
    context = SimpleNamespace(manifest=SimpleNamespace(assets=[]))
    embedding, source = load_declared_cell_embedding(context, frame)
    assert list(embedding) == ["cell_id", "state_1", "state_2", "state_3"]
    assert "dimensions=3" in source


def test_deg_only_candidates_keep_rank_and_missing_as():
    deg = pd.DataFrame([dict(comparison_id="A", gene_id="G1", gene_symbol="gene1", log2fc=2., deg_fdr=.001)])
    result = build_candidate_evidence(pd.DataFrame(), deg, pd.DataFrame(), {})
    assert len(result) == 1
    assert result.iloc[0].gene_symbol == "gene1"
    assert not result.iloc[0].as_available
    assert result.iloc[0].deg_supported
    assert result.iloc[0].ranking_score == 1
    assert result.columns[-1] == "ranking_score"


def test_scope_policy_and_selection_remain_bound(tmp_path: Path):
    manifest = create_full_fixture(tmp_path / "fixture")
    runtime = PlatformRuntime(application_root=tmp_path / "app")
    context = runtime.projects.open_project(manifest)
    project_id = context.manifest.project_id

    def run(name, params=None):
        started = runtime.workflows.run(name, project_id, parameters=params)
        finished = runtime.jobs.wait(started.job_id, timeout=90)
        assert finished.status in {JobStatus.COMPLETE, JobStatus.CACHE_HIT}, finished.error
        return [context.catalog.get_artifact(i) for i in finished.artifact_ids]

    quality = run("quality.stable_as", {"stability_profile": "custom", "min_reads": 7, "min_direction_consistency": .6})
    q = next(a for a in quality if a.name == "Stable-AS event QC")
    das = run("differential.as")
    full = next(a for a in das if a.name == "Differential AS events")
    assert full.metadata["measurement_policy_source"] == q.artifact_id
    full_data = pd.read_parquet(full.path)
    assert full_data.criterion_min_reads.eq(7).all()
    run("prioritize.evidence")
    original = resolve_table(context, "prioritize.evidence", "Candidate evidence matrix")
    candidates = pd.read_parquet(original.path)
    selection = SelectionSet(project_id=project_id, name="original", entity_type="gene", source_artifact_id=original.artifact_id,
                             entity_ids=tuple(candidates.row_id.head(3)))
    context.catalog.save_selection(selection)
    subset = [context.manifest.contrasts[0].contrast_id]
    run("differential.as", {"contrast_ids": subset})
    run("prioritize.evidence", {"contrast_ids": subset})
    assert resolve_table(context, "differential.as", "Differential AS events").artifact_id == full.artifact_id
    genes, source = _selection_genes(context, {"selection_id": selection.selection_id})
    assert source.artifact_id == original.artifact_id and genes
    # A file from an unfinished run cannot enter a report just because it exists.
    from katana_next.workflows.reporting import ReproducibilityReportWorkflow
    from katana_next.workflows.base import ArtifactWriter
    writer = ArtifactWriter(context, "test.incomplete", "incomplete", {}, "1.0")
    partial = writer.table("Partial", "partial.parquet", pd.DataFrame({"a": [1]}), entity_type="test")
    assert partial.artifact_id not in {a.artifact_id for a in ReproducibilityReportWorkflow().parent_artifacts(context, {})}
    nested_writer = ArtifactWriter(context, "test.nested", "nested-run", {}, "1.0")
    nested = nested_writer.directory / "plots" / "event" / "counts.tsv"
    nested.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"counts": [3]}).to_csv(nested, sep="\t", index=False)
    nested_writer.existing_file("Nested counts", nested, entity_type="read_evidence")
    for event_name in ("event-a", "event-b"):
        sidecar = nested_writer.directory / "plots" / event_name / "summary.tsv"
        sidecar.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"counts": [3]}).to_csv(sidecar, sep="\t", index=False)
        primary = sidecar.with_suffix(".txt")
        primary.write_text(event_name, encoding="utf-8")
        nested_writer.existing_file(event_name, primary, entity_type="read_evidence", metadata={"exports": {"tsv": str(sidecar)}})
    nested_writer.finish()
    before = len(runtime.jobs.list(project_id))
    run("differential.cross_comparison")
    assert len(runtime.jobs.list(project_id)) == before + 1
    cross = context.catalog.list_artifacts(project_id, analysis_id="differential.cross_comparison")
    assert any("Venn" in a.name for a in cross)
    plots = [a for a in cross if a.metadata.get("exports", {}).get("png")]
    assert len(plots) == 3
    assert all(Path(a.metadata["exports"]["png"]).stat().st_size > 1000 for a in plots)
    report = run("reports.bundle")
    html = next(a for a in report if a.name == "Katana evidence report (HTML)").path.read_text(encoding="utf-8")
    assert "data:image/png;base64," in html
    import zipfile
    archive = next(a for a in report if a.path.suffix == ".zip")
    with zipfile.ZipFile(archive.path) as handle:
        assert any(name.endswith(".png") for name in handle.namelist())
        assert any(name.endswith("parameters.json") for name in handle.namelist())
        assert "artifacts/test.nested/nested-run/parameters.json" in handle.namelist()
        assert "artifacts/test.nested/nested-run/plots/event-a/summary.tsv" in handle.namelist()
        assert "artifacts/test.nested/nested-run/plots/event-b/summary.tsv" in handle.namelist()
        assert len(handle.namelist()) == len(set(handle.namelist()))
    fingerprint, cached = runtime.workflows.compatible_cache("differential.cross_comparison", project_id)
    assert len(cached) == 7
    assert any(a.name == "Intersection member identities and directions" for a in cached)
    # A partially missing figure is not a complete cache hit.
    Path(plots[0].metadata["exports"]["png"]).unlink()
    assert not runtime.workflows.compatible_cache("differential.cross_comparison", project_id)[1]


def test_full_length_without_barcode_uses_explicit_single_cell_mapping(tmp_path, monkeypatch):
    import sys
    from katana_next.science.single_cell import extract_junction_choices_from_bams
    metadata = tmp_path / "metadata.tsv"
    pd.DataFrame([dict(cell_id="cell1", sample_id="s1", donor_id="d1", cell_type="T", condition="WT", assay="full_length")]).to_csv(metadata, sep="\t", index=False)
    bam = tmp_path / "sample.bam"
    bam.touch()

    class Read:
        reference_name = "chr1"
        reference_start = 0
        mapping_quality = 60
        cigartuples = [(0, 10), (3, 100), (0, 10)]

        def __init__(self, name):
            self.query_name = name

        def get_tag(self, _):
            raise KeyError()

    class Alignment:
        def __init__(self, *_):
            pass

        def fetch(self, **_):
            return iter([Read("r1"), Read("r1"), Read("r2")])

        def close(self):
            pass

    monkeypatch.setitem(sys.modules, "pysam", SimpleNamespace(AlignmentFile=Alignment))
    context = SimpleNamespace(manifest=SimpleNamespace(
        assets=[SimpleNamespace(kind="single_cell_metadata", path=metadata)],
        samples=[SimpleNamespace(sample_id="s1", assets={"single_cell_bam": bam})],
        annotation=SimpleNamespace(gtf=None)))
    counts = extract_junction_choices_from_bams(context, tmp_path / "out" / "counts.tsv")
    assert counts.cell_id.tolist() == ["cell1"]
    assert counts.inclusion.tolist() == [2]
    assert counts.count_mode.tolist() == ["read_name"]
    assert not list((tmp_path / "out").glob(".junction-dedup-*"))


def test_four_set_intersections_are_exact_and_plottable(tmp_path):
    from katana_next.workflows.intersection_plots import exact_intersections, render_intersections
    sets = {"A": {"1", "2", "3"}, "B": {"2", "3"}, "C": {"3"}, "D": {"4"}}
    rows = exact_intersections(sets)
    assert rows.event_count.sum() == 4
    assert rows.loc[rows.A & rows.B & rows.C, "event_count"].iloc[0] == 1
    output = render_intersections(sets, {name: name for name in sets}, tmp_path)
    assert Path(output["png"]).stat().st_size > 1000


@pytest.mark.parametrize("bad", [-1, float("inf"), float("nan")])
def test_numeric_workflow_parameters_reject_invalid_values(bad):
    from katana_next.workflows.bulk import CrossComparisonWorkflow
    with pytest.raises(ValueError):
        CrossComparisonWorkflow().normalize_parameters({"min_reads": bad})
