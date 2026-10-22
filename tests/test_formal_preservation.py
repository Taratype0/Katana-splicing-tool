from __future__ import annotations

import pandas as pd
import numpy as np
import pytest

from katana_cosplice.formal_preservation import (
    prepare_formal_preservation_inputs,
)


def test_native_constant_null_is_not_a_negative_preservation_result():
    from katana_cosplice.preservation import _zscore, _category, _empirical_p
    from katana_cosplice.config import PreservationConfig
    assert np.isnan(_zscore(1., np.ones(20)))
    assert _category(float("nan"), PreservationConfig()) == "not_estimable"
    assert np.isnan(_empirical_p(1., np.ones(20)))
    assert np.isnan(_empirical_p(float("nan"), np.arange(20)))


def _native_fixture(tmp_path):
    from katana_cosplice.config import ExternalCohortConfig, PreservationConfig
    rng = np.random.default_rng(71)
    reference = pd.DataFrame(rng.uniform(size=(18, 12)), index=[f"E{i:02}" for i in range(18)])
    target = pd.DataFrame(rng.uniform(size=(18, 10)), index=reference.index,
                          columns=[f"T{i}" for i in range(10)])
    target_path, samples = tmp_path / "target.tsv", tmp_path / "samples.tsv"
    target.to_csv(target_path, sep="\t")
    pd.DataFrame({"sample_id": target.columns}).to_csv(samples, sep="\t", index=False)
    assignments = pd.DataFrame({"event_key": reference.index[:6], "consensus_module": "M01"})
    return reference, assignments, ExternalCohortConfig("target", target_path, samples), PreservationConfig(
        permutations=30, min_shared_events=3)


@pytest.mark.parametrize("explicit_grey", [True, False])
def test_native_background_includes_unassigned_events(tmp_path, explicit_grey):
    from katana_cosplice.preservation import evaluate_external_module_preservation
    reference, assignments, cohort, settings = _native_fixture(tmp_path)
    if explicit_grey:
        assignments = pd.concat([assignments, pd.DataFrame({"event_key": reference.index[6:],
                                                            "consensus_module": "grey"})])
    result = evaluate_external_module_preservation(reference, assignments, [cohort], settings)
    row = result.module_preservation.iloc[0]
    assert row.background_events == 18 and row.background_only_events == 12
    assert row.native_null_status == "available" and np.isfinite(row.Zsummary_wgcna_like)
    assert set(result.module_preservation.consensus_module) == {"M01"}
    assert result.event_mapping.background_only.sum() == 12


def test_native_background_requires_quality_in_both_cohorts(tmp_path):
    from katana_cosplice.preservation import evaluate_external_module_preservation
    reference, assignments, cohort, settings = _native_fixture(tmp_path)
    reference.loc["E06"] = np.nan
    reference.loc["E07"] = 1.
    reference.loc["E08"] = np.inf
    target = pd.read_csv(cohort.psi_matrix, sep="\t", index_col=0)
    target.loc["E09"] = np.nan
    target.loc["E10"] = .4
    target.to_csv(cohort.psi_matrix, sep="\t")
    result = evaluate_external_module_preservation(reference, assignments, [cohort], settings)
    assert result.module_preservation.iloc[0].background_events == 13
    assert not set(["E06", "E07", "E08", "E09", "E10"]) & set(result.event_mapping.event_key)


def test_native_draws_are_paired_and_independent_of_traversal_order(tmp_path):
    from katana_cosplice.preservation import evaluate_external_module_preservation
    from katana_cosplice.config import ExternalCohortConfig
    reference, assignments, cohort, settings = _native_fixture(tmp_path)
    other = ExternalCohortConfig("other", cohort.psi_matrix, cohort.sample_metadata)
    first = evaluate_external_module_preservation(reference, assignments, [cohort, other], settings)
    second = evaluate_external_module_preservation(reference.iloc[::-1], assignments.iloc[::-1],
                                                   [other, cohort], settings)
    pd.testing.assert_frame_equal(first.module_preservation.sort_values("cohort").reset_index(drop=True),
                                  second.module_preservation.sort_values("cohort").reset_index(drop=True))
    changed = assignments.copy()
    changed.event_key = reference.index[6:12]
    third = evaluate_external_module_preservation(reference, changed, [cohort], settings)
    pd.testing.assert_frame_equal(first.null_summary.iloc[:1].reset_index(drop=True), third.null_summary)


def test_native_matrix_cache_does_not_return_old_path_contents(tmp_path):
    from katana_cosplice.preservation import _read_matrix
    _, _, cohort, _ = _native_fixture(tmp_path)
    old = _read_matrix(cohort.psi_matrix).copy()
    (old + .123456789).to_csv(cohort.psi_matrix, sep="\t")
    np.testing.assert_allclose(_read_matrix(cohort.psi_matrix), old + .123456789)


@pytest.mark.parametrize("problem", ["duplicate_reference", "duplicate_assignment", "unknown_event", "missing_label"])
def test_native_rejects_ambiguous_event_identity(tmp_path, problem):
    from katana_cosplice.preservation import evaluate_external_module_preservation
    reference, assignments, cohort, settings = _native_fixture(tmp_path)
    if problem == "duplicate_reference":
        reference = pd.concat([reference, reference.iloc[:1]])
    elif problem == "duplicate_assignment":
        assignments = pd.concat([assignments, assignments.iloc[:1]])
    elif problem == "unknown_event":
        assignments.loc[0, "event_key"] = "not-in-reference"
    else:
        assignments.loc[0, "consensus_module"] = None
    with pytest.raises(ValueError):
        evaluate_external_module_preservation(reference, assignments, [cohort], settings)


def test_native_stage_cache_tracks_method_and_never_runs_upstream(tmp_path, monkeypatch):
    import json
    import katana_cosplice.pipeline as engine
    from katana_cosplice.config import CoSpliceConfig
    reference, assignments, cohort, settings = _native_fixture(tmp_path)
    output = tmp_path / "engine"
    reference_path = output / "04_network/observed/selected_event_matrix.tsv"
    consensus_path = output / "04c_consensus_modules/consensus_event_module_assignments.tsv"
    reference_path.parent.mkdir(parents=True)
    consensus_path.parent.mkdir(parents=True)
    reference.to_csv(reference_path, sep="\t")
    assignments.to_csv(consensus_path, sep="\t", index=False)
    pipeline = engine.CoSplicePipeline(CoSpliceConfig(
        project_name="cache-test", rmats_root=None, sample_metadata=cohort.sample_metadata,
        output_dir=output, external_cohorts=[cohort], preservation=settings))
    result = pipeline.preserve_modules()
    first = pipeline.status()["stages"]["preserve_modules"]["fingerprint"]
    parameters = json.loads(result["parameters"].read_text())
    assert parameters["native_method_version"] == engine.PRESERVATION_METHOD_VERSION
    pipeline.preserve_modules()
    assert pipeline.status()["stages"]["preserve_modules"]["status"] == "cache_hit"
    monkeypatch.setattr(engine, "PRESERVATION_METHOD_VERSION", "test-next-revision")
    pipeline.preserve_modules()
    state = pipeline.status()["stages"]
    assert set(state) == {"preserve_modules"}
    assert state["preserve_modules"]["status"] == "finished"
    assert state["preserve_modules"]["fingerprint"] != first


def test_native_all_not_estimable_plot_is_explicit(tmp_path):
    from katana_cosplice.preservation import evaluate_external_module_preservation, render_preservation_plots
    reference, assignments, cohort, settings = _native_fixture(tmp_path)
    assignments = pd.DataFrame({"event_key": reference.index, "consensus_module": "M01"})
    result = evaluate_external_module_preservation(reference, assignments, [cohort], settings)
    paths = render_preservation_plots(result, tmp_path / "plots")
    assert len(paths) == 2 and all(path.stat().st_size > 1000 for path in paths)


def test_native_single_module_cannot_use_its_own_members_as_null_universe(tmp_path):
    from katana_cosplice.preservation import evaluate_external_module_preservation
    from katana_cosplice.config import ExternalCohortConfig, PreservationConfig
    matrix = pd.DataFrame(np.random.default_rng(19).uniform(size=(6, 5)),
                          index=[f"E{i}" for i in range(6)], columns=[f"S{i}" for i in range(5)])
    target, samples = tmp_path/"target.tsv", tmp_path/"samples.tsv"
    matrix.to_csv(target, sep="\t")
    pd.DataFrame({"sample_id": matrix.columns}).to_csv(samples, sep="\t", index=False)
    assignments = pd.DataFrame({"event_key": matrix.index, "consensus_module": "M01"})
    result = evaluate_external_module_preservation(matrix, assignments,
        [ExternalCohortConfig("same", target, samples)], PreservationConfig(permutations=8, min_shared_events=3))
    row = result.module_preservation.iloc[0]
    assert row.native_null_status == "module_exhausts_permutation_universe"
    assert row.preservation_class == "not_estimable"
    assert np.isnan(row.Zsummary_wgcna_like) and np.isnan(row.density_empirical_pvalue)


def test_formal_single_module_exhausted_null_is_rejected_before_r():
    matrix = pd.DataFrame(np.random.default_rng(9).uniform(size=(8, 5)), index=[f'e{i}' for i in range(8)])
    assignments = pd.DataFrame({'event_key':matrix.index,'consensus_module':'M1'})
    with pytest.raises(ValueError,match='permutation universe'):
        prepare_formal_preservation_inputs(matrix,matrix,assignments,min_module_size=3)
    assignments = assignments.iloc[:3]
    ref,_,labels = prepare_formal_preservation_inputs(matrix,matrix,assignments,min_module_size=3)
    assert len(ref) == 8 and labels.wgcna_module.eq('grey').sum() == 5


def test_formal_preservation_keeps_common_background() -> None:
    reference = pd.DataFrame(
        {
            "s1": [0.1, 0.2, 0.3, 0.4],
            "s2": [0.2, 0.3, 0.5, 0.7],
            "s3": [0.3, 0.5, 0.7, 0.9],
        },
        index=["A", "B", "C", "D"],
    )
    target = pd.DataFrame(
        {
            "t1": [0.2, 0.1, 0.4, 0.8],
            "t2": [0.3, 0.2, 0.6, 0.6],
            "t3": [0.5, 0.4, 0.8, 0.3],
        },
        index=reference.index,
    )
    assignments = pd.DataFrame(
        {
            "event_key": ["A", "B"],
            "consensus_module": ["M01", "M01"],
        }
    )

    aligned_reference, aligned_target, labels = (
        prepare_formal_preservation_inputs(
            reference, target, assignments, min_module_size=2
        )
    )

    assert aligned_reference.index.tolist() == ["A", "B", "C", "D"]
    assert aligned_target.index.tolist() == ["A", "B", "C", "D"]
    assert labels["wgcna_module"].tolist() == ["M01", "M01", "grey", "grey"]
    assert labels["netrep_module"].tolist() == ["M01", "M01", "0", "0"]
