"""Fixed informative-read support over every sample in the selected contrast."""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from katana_next.science.bulk import (
    DEFAULT_MIN_READS,
    READ_SUPPORT_POLICY_VERSION,
    _declared_sample_names,
    _sample_names,
    apply_stable_event_qc,
    load_bulk_events,
    significance_columns,
)


def event(**changes):
    row = dict(
        IncLevel1="0.8,0.7,0.9", IncLevel2="0.2,0.3,0.1",
        IJC_SAMPLE_1="5,5,5", SJC_SAMPLE_1="5,5,5",
        IJC_SAMPLE_2="4,4,4", SJC_SAMPLE_2="6,6,6",
        numerator_expected_samples=3, denominator_expected_samples=3,
        numerator_replicates=3, denominator_replicates=3,
        reverse_direction=False, FDR=0.031, PValue=0.007,
        IncLevelDifference=0.6, fdr=0.031, dpsi=0.6, abs_dpsi=0.6,
        replicate_pairwise_direction_consistency=0.1,
        leave_one_out_direction_consistency=0.0, max_measurement_se=0.8,
    )
    row.update(changes)
    return row


def evaluate(row, parameters=None):
    return apply_stable_event_qc(pd.DataFrame([row]), parameters or {}).iloc[0]


def test_ten_counts_in_every_sample_passes_without_direction_or_se_gates():
    row = evaluate(event())
    assert DEFAULT_MIN_READS == 10
    assert row.read_support_policy == READ_SUPPORT_POLICY_VERSION
    assert row.measurement_eligible and row.stable_event
    assert row.passes_read_counts and row.passes_data_usability
    assert not row.measurement_se_gate_enabled
    assert "differential_repeat_supported" not in row
    assert "criterion_max_measurement_se" not in row


def test_one_low_sample_out_of_three_fails_even_if_two_samples_pass():
    row = evaluate(event(IJC_SAMPLE_1="5,5,4"))
    assert not row.measurement_eligible
    assert row.numerator_supported_replicates == 2
    assert row.numerator_low_read_samples == 1
    assert "numerator:sample_3:counts_below_10" in row.read_support_failure_reason


@pytest.mark.parametrize("column,value", [
    ("IJC_SAMPLE_1", "5,NA,5"),
    ("SJC_SAMPLE_1", "5,inf,5"),
    ("IJC_SAMPLE_1", "5,-1,5"),
    ("IJC_SAMPLE_1", "5"),
    ("SJC_SAMPLE_1", "5,5,5,5"),
    ("SJC_SAMPLE_1", ""),
    ("IJC_SAMPLE_2", None),
    ("SJC_SAMPLE_2", "6,bad,6"),
])
def test_missing_invalid_or_misaligned_count_vectors_cannot_pass(column, value):
    row = evaluate(event(**{column: value}))
    assert not row.measurement_eligible
    assert row.read_support_failure_reason != "passed"


@pytest.mark.parametrize("psi", ["0.8,NA,0.9", "0.8,1.01,0.9", "0.8,-0.1,0.9", "0.8,inf,0.9"])
def test_invalid_psi_is_unassessable_despite_sufficient_counts(psi):
    row = evaluate(event(IncLevel1=psi))
    assert row.passes_read_counts
    assert not row.passes_data_usability
    assert not row.passes_read_support
    assert row.numerator_invalid_samples == 1


def test_missing_psi_vector_keeps_count_support_separate_from_data_usability():
    row = evaluate(event(IncLevel1=""))
    assert row.passes_read_counts
    assert not row.passes_data_usability
    assert not row.measurement_eligible
    assert row.numerator_read_count_supported_samples == 3
    assert row.numerator_supported_replicates == 0


def test_zero_count_and_exact_ten_count_boundaries():
    assert evaluate(event(IJC_SAMPLE_1="0,0,0", SJC_SAMPLE_1="10,10,10")).measurement_eligible
    row = evaluate(event(IJC_SAMPLE_1="0,0,0", SJC_SAMPLE_1="0,0,0"))
    assert not row.measurement_eligible and row.numerator_low_read_samples == 3


def test_one_vs_one_passes_measurement_qc_without_claiming_model_eligibility():
    row = event(IncLevel1="0.4", IncLevel2="0.6", IJC_SAMPLE_1="5", SJC_SAMPLE_1="5",
                IJC_SAMPLE_2="10", SJC_SAMPLE_2="0", numerator_expected_samples=1, denominator_expected_samples=1)
    result = evaluate(row)
    assert result.measurement_eligible
    assert result.numerator_supported_replicates == 1


def test_reverse_mapping_keeps_declared_arm_and_failure_identity():
    result = evaluate(event(reverse_direction=True, IJC_SAMPLE_1="5,5,4"))
    assert not result.measurement_eligible
    assert result.denominator_supported_replicates == 2
    assert result.numerator_supported_replicates == 3
    assert "denominator:sample_3:counts_below_10" in result.read_support_failure_reason


def test_declared_three_samples_cannot_be_reduced_to_two_vector_entries():
    result = evaluate(event(IncLevel1="0.8,0.7", IJC_SAMPLE_1="5,5", SJC_SAMPLE_1="5,5"))
    assert not result.measurement_eligible
    assert result.numerator_expected_samples == 3
    assert "PSI_vector_length_2_expected_3" in result.read_support_failure_reason


def test_helper_without_declared_counts_infers_psi_length_but_never_an_empty_arm():
    row = event()
    row.pop("numerator_expected_samples")
    row.pop("denominator_expected_samples")
    assert evaluate(row).measurement_eligible
    row.update(IncLevel1="", IJC_SAMPLE_1="", SJC_SAMPLE_1="")
    result = evaluate(row)
    assert not result.measurement_eligible
    assert "numerator:empty_arm" in result.read_support_failure_reason


@pytest.mark.parametrize("count", [0, -1, None, float("nan"), 2.5, True])
def test_invalid_declared_sample_count_is_not_inferred_away(count):
    result = evaluate(event(numerator_expected_samples=count))
    assert not result.measurement_eligible
    assert "invalid_declared_sample_count" in result.read_support_failure_reason


@pytest.mark.parametrize("parameters", [
    {"min_reads": 20}, {"min_reads": 0}, {"stability_profile": "stringent"},
    {"stability_profile": "custom"}, {"min_supported_replicates": 3},
    {"min_direction_consistency": 0.8}, {"apply_measurement_se_gate": True},
    {"max_measurement_se": 0.1}, {"read_support_policy": "unknown"},
])
def test_nondefault_legacy_options_are_rejected_instead_of_silently_reactivated(parameters):
    with pytest.raises(ValueError):
        evaluate(event(), parameters)


def test_inert_legacy_defaults_are_recorded_but_do_not_enforce_two_replicates():
    result = evaluate(event(IncLevel1="0.4", IncLevel2="0.6", IJC_SAMPLE_1="5", SJC_SAMPLE_1="5",
                            IJC_SAMPLE_2="5", SJC_SAMPLE_2="5", numerator_expected_samples=1, denominator_expected_samples=1),
                      {"stability_profile": "balanced", "min_supported_replicates": 2,
                       "min_direction_consistency": .75, "max_measurement_se": .12,
                       "apply_measurement_se_gate": False})
    assert result.measurement_eligible
    assert json.loads(result.ignored_legacy_qc_options)["min_supported_replicates"] == 2


def test_qc_preserves_all_native_statistics_counts_psi_and_input_rows():
    original = pd.DataFrame([event(), event(IJC_SAMPLE_1="5,4,5")])
    before = original.copy(deep=True)
    result = apply_stable_event_qc(original, {})
    pd.testing.assert_frame_equal(original, before)
    pd.testing.assert_frame_equal(result[before.columns], before)
    selected = significance_columns(result, {"fdr": .05, "min_abs_dpsi": .1})
    assert selected.significant.tolist() == [True, True]
    assert selected.stable_significant.tolist() == [True, False]
    assert selected.FDR.tolist() == original.FDR.tolist()


def test_empty_table_never_creates_a_passed_event():
    result = apply_stable_event_qc(pd.DataFrame(), {})
    assert result.empty
    assert result.measurement_eligible.dtype == bool


@pytest.mark.parametrize("mapping", [["a1"], ["a1", "a1"], ["a1", "other"], []])
def test_asset_mapping_cannot_shrink_duplicate_or_replace_selected_contrast_samples(mapping):
    contrast = SimpleNamespace(contrast_id="selected", numerator=SimpleNamespace(sample_ids=["a1", "a2"]))
    asset = SimpleNamespace(metadata={"numerator_sample_ids": mapping})
    with pytest.raises(ValueError, match="every declared"):
        _declared_sample_names(asset, contrast, "numerator")


def test_asset_mapping_allows_a_complete_reordering():
    contrast = SimpleNamespace(contrast_id="selected", numerator=SimpleNamespace(sample_ids=["a1", "a2"]))
    asset = SimpleNamespace(metadata={"numerator_sample_ids": ["a2", "a1"]})
    assert _sample_names(asset, contrast, "numerator", 2) == ["a2", "a1"]


def test_loader_records_manifest_scope_rebuilds_old_cache_and_preserves_psi(tmp_path):
    from katana_next.domain.models import AnnotationReference, Contrast, ContrastArm, InputAsset, ProjectManifest, Sample
    from katana_next.infrastructure.manifest_store import ManifestStore
    from katana_next.services.project_service import ProjectService

    inputs = tmp_path / "inputs"
    inputs.mkdir()
    source = inputs / "SE.MATS.JC.txt"
    native = {key: value for key, value in event().items() if key in {
        "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2",
        "SJC_SAMPLE_1", "SJC_SAMPLE_2", "FDR", "PValue", "IncLevelDifference",
    }}
    native.update(ID="event-a", GeneID="GENE_A", geneSymbol="A")
    pd.DataFrame([native]).to_csv(source, sep="\t", index=False)
    samples = tuple(Sample(sample_id=f"{arm}{number}", condition=arm, biological_sample_id=f"bio-{arm}{number}") for arm in ("a", "b") for number in range(1, 4))
    contrast = Contrast(contrast_id="selected", label="A vs B", numerator=ContrastArm(label="A", sample_ids=["a1", "a2", "a3"]), denominator=ContrastArm(label="B", sample_ids=["b1", "b2", "b3"]))
    asset = InputAsset(asset_id="events", kind="rmats_result", path=source, metadata={"contrast_id": "selected", "event_type": "SE", "numerator_sample_ids": ["a3", "a2", "a1"]})
    manifest = ProjectManifest(project_id="synthetic", name="Synthetic", input_root=inputs, output_root=tmp_path / "results", annotation=AnnotationReference(genome_build="synthetic"), samples=samples, contrasts=(contrast,), assets=(asset,))
    path = ManifestStore().save(manifest, tmp_path / "project.json")
    context = ProjectService(application_root=tmp_path / "application").open_project(path)
    loaded = load_bulk_events(context, {})
    row = loaded.events.iloc[0]
    assert row.numerator_expected_samples == row.denominator_expected_samples == 3
    assert loaded.psi_matrix.loc[row.event_id, "a3"] == .8
    assert loaded.psi_matrix.loc[row.event_id, "a1"] == .9
    assert evaluate(row.to_dict()).measurement_eligible
    cache_metadata = next((manifest.output_root / ".katana/cache/bulk_events").glob("*/cache.json"))
    old = json.loads(cache_metadata.read_text())
    old["loader_version"] = "3.0-na-preserving"
    cache_metadata.write_text(json.dumps(old))
    rebuilt = load_bulk_events(context, {})
    assert rebuilt.source_cache_status == "rebuild"
    pd.testing.assert_frame_equal(rebuilt.psi_matrix, loaded.psi_matrix)


def test_missing_sjc_is_not_padded_to_zero_to_create_passed_support():
    row = evaluate(event(IJC_SAMPLE_1="10,10,10", SJC_SAMPLE_1=""))
    assert not row.passes_read_counts
    assert not row.measurement_eligible
    assert row.numerator_read_count_supported_samples == 0
