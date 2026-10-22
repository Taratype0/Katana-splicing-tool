import numpy as np
import pandas as pd
import pytest

from katana_next.science.single_cell_isoforms import analyze_isoform_usage, validate_isoform_inputs
from katana_next.workflows.registry import WorkflowRegistry


def fixture(donors=3):
    meta, counts = [], []
    for condition in ("WT", "MUT"):
        for donor in range(donors):
            for cell in range(6):
                identity = f"{condition}-{donor}-{cell}"
                meta.append(dict(cell_id=identity, donor_id=f"{condition}-{donor}", condition=condition, cell_type="neuron"))
                first = (3 if condition == "WT" else 16) + donor
                for transcript, count in (("T1", first), ("T2", 25-first)):
                    counts.append(dict(cell_id=identity, gene_id="G1", transcript_id=transcript, count=count))
    return pd.DataFrame(counts), pd.DataFrame(meta)


def run(counts, metadata, **kwargs):
    return analyze_isoform_usage(counts, metadata, numerator="MUT", denominator="WT", rarefaction_draws=5, **kwargs)


def test_usage_direction_and_donor_unit():
    counts, metadata = fixture()
    result = run(counts, metadata)
    assert len(result.diversity) == 6
    assert result.differential.iloc[0].numerator_donors == 3
    assert result.differential.iloc[0].model_status == "fit"
    assert result.effects.set_index("transcript_id").loc["T1", "delta_fraction"] > .4
    assert np.isclose(result.effects.delta_fraction.sum(), 0)
    assert result.usage.groupby(["donor_id", "gene_id"]).fraction.sum().eq(1).all()
    assert result.differential.fdr.between(0, 1).all()
    assert result.effects.lodo_evaluable.eq(6).all()
    assert result.effects.lodo_direction_consistency.eq(1).all()


def test_two_donors_remain_descriptive():
    result = run(*fixture(2))
    assert result.differential.model_status.eq("insufficient_donors").all()
    assert result.differential.fdr.isna().all()
    assert len(result.effects) == 2
    assert result.effects.lodo_status.eq("not_estimable").all()


def test_absent_gene_is_not_zero_usage():
    counts, meta = fixture()
    counts.loc[counts.cell_id.str.startswith("MUT-0"), "count"] = 0
    result = run(counts, meta)
    assert result.usage.loc[result.usage.donor_id.eq("MUT-0"), "fraction"].isna().all()
    assert result.differential.iloc[0].numerator_donors == 2


def test_no_measured_gene_keeps_exclusions():
    result = run(*fixture(), min_gene_count=100000)
    assert result.effects.empty
    assert result.diagnostics.measured.eq(False).all()
    assert result.diversity.renyi2_bits.isna().all()


def test_entropy_of_mean_cells_is_not_pooled_entropy():
    counts, metadata = fixture()
    cell = counts.cell_id.eq("MUT-0-0")
    counts.loc[cell & counts.transcript_id.eq("T1"), "count"] = 1000
    result = run(counts, metadata)
    row = result.diversity.set_index("donor_id").loc["MUT-0"]
    assert row.renyi2_bits < row.mean_cell_fraction_renyi2_bits


def test_technical_cells_do_not_inflate_donor_count():
    counts, metadata = fixture(2)
    extra_meta = metadata.assign(cell_id=metadata.cell_id + "-extra")
    extra_counts = counts.assign(cell_id=counts.cell_id + "-extra")
    result = run(pd.concat([counts, extra_counts]), pd.concat([metadata, extra_meta]))
    assert result.differential.numerator_donors.eq(2).all()
    assert result.differential.fdr.isna().all()


def test_many_transcripts_pool_tail_without_losing_counts(monkeypatch):
    counts, metadata = fixture()
    extra = [dict(cell_id=cell, gene_id="G1", transcript_id=f"tail{i:02}", count=1)
             for cell in metadata.cell_id for i in range(30)]
    counts = pd.concat([counts, pd.DataFrame(extra)], ignore_index=True)
    def checked_fit(matrix, groups, covariates=None):
        assert matrix.shape == (6, 30)
        assert np.all(matrix.sum(axis=1) == 330)
        return dict(pvalue=.5, status="fit", likelihood_ratio=1)
    monkeypatch.setattr("katana_next.science.single_cell_isoforms.dirichlet_multinomial_two_group", checked_fit)
    result = run(counts, metadata)
    assert result.differential.tested_transcript_choices.eq(30).all()
    assert len(result.differential.iloc[0].pooled_transcripts.split(";")) == 3
    assert result.effects.transcript_id.nunique() == 32


def test_rarefaction_reproducible_and_no_upsampling():
    a = run(*fixture())
    b = run(*fixture())
    pd.testing.assert_frame_equal(a.diversity, b.diversity)
    assert a.diversity.rarefied_renyi2_mean.between(0, 1).all()
    assert run(*fixture(), rarefaction_depth=10000).diversity.rarefaction_status.eq("not_measurable").all()


@pytest.mark.parametrize("bad", [np.nan, -1, .25, np.inf])
def test_invalid_counts_rejected(bad):
    counts, meta = fixture()
    counts["count"] = counts["count"].astype(float)
    counts.loc[0, "count"] = bad
    with pytest.raises(ValueError):
        validate_isoform_inputs(counts, meta)


def test_unknown_cells_and_duplicates_rejected():
    counts, meta = fixture()
    with pytest.raises(ValueError, match="Duplicate"):
        validate_isoform_inputs(pd.concat([counts, counts.iloc[:1]]), meta)
    with pytest.raises(ValueError, match="absent"):
        validate_isoform_inputs(counts, meta.iloc[1:])


def test_single_cell_registration_independent_from_cosplice():
    workflow = WorkflowRegistry().get("single_cell.isoform_usage")
    assert workflow.definition.workspace_id == "single-cell"
    assert workflow.parent_artifacts(None, {}) == []
    assert "CoSplice" not in workflow.definition.method


def test_workflow_exports_and_cache_are_isolated(tmp_path):
    import json
    from pathlib import Path
    from fastapi.testclient import TestClient
    from katana_next.api.app import PlatformRuntime, create_app
    from katana_next.domain.models import JobStatus
    from katana_next.full_fixture import create_full_fixture

    manifest = create_full_fixture(tmp_path / "fixture")
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    metadata_path = Path(next(a["path"] for a in payload["assets"] if a["kind"] == "single_cell_metadata"))
    counts, meta = fixture()
    meta["sample_id"] = meta.donor_id
    meta.to_csv(metadata_path, sep="\t", index=False)
    counts_path = tmp_path / "isoform_counts.tsv"
    counts.to_csv(counts_path, sep="\t", index=False)
    runtime = PlatformRuntime(application_root=tmp_path / "application")
    client = TestClient(create_app(runtime))
    opened = client.post("/api/projects/open", json={"manifest_path": str(manifest)})
    assert opened.status_code == 200, opened.text
    project_id = opened.json()["snapshot"]["project_id"]
    endpoint = f"/api/projects/{project_id}/workflows/single_cell.isoform_usage"
    params = dict(counts_path=str(counts_path), numerator="MUT", denominator="WT", rarefaction_draws=3)
    started = client.post(endpoint + "/run", json={"parameters": params})
    assert started.status_code == 200, started.text
    done = runtime.jobs.wait(started.json()["job_id"], timeout=120)
    assert done.status == JobStatus.COMPLETE, done
    artifacts = client.get(f"/api/projects/{project_id}/artifacts").json()
    assert len(artifacts) == 7
    assert {a["analysis_id"] for a in artifacts} == {"single_cell.isoform_usage"}
    for artifact in artifacts:
        exports = artifact["metadata"]["exports"]
        extension = "png" if artifact["kind"] == "plot_spec" else "tsv"
        assert Path(exports[extension]).stat().st_size > 0
    cache = client.post(endpoint + "/cache", json={"parameters": params})
    assert cache.json()["cache_hit"] is True
    changed = client.post(endpoint + "/cache", json={"parameters": params | {"rarefaction_depth": 10}})
    assert changed.json()["cache_hit"] is False
    empty = client.post(endpoint + "/run", json={"parameters": params | {"min_gene_count": 100000}})
    empty_done = runtime.jobs.wait(empty.json()["job_id"], timeout=120)
    assert empty_done.status == JobStatus.COMPLETE, empty_done
