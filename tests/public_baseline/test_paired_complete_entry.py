"""Paired biological units survive raw-input validation and downstream inference."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import ttest_rel

from katana_next.domain.models import ProjectManifest
from katana_next.flow.network_inputs import network_scope_identity
from katana_next.reports.network_views import matched_module_associations

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("paired_validator", ROOT / "pipeline/bin/katana_validate_inputs.py")
validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)


def samples():
    return [dict(sample=f"{arm}{i}", condition=arm, subject=f"mouse{i}",
                 single_end="false", strandedness="reverse", read_length="150")
            for arm, order in (("treated", [3, 1, 2]), ("control", [2, 3, 1])) for i in order]


def contrast(**overrides):
    return dict(contrast="treated_vs_control", treatment="treated", control="control", paired="true", **overrides)


def test_subject_order_not_csv_order_and_explicit_model():
    rows, audit = validator.validate_contrasts([contrast()], samples(), design="~ condition", analysis_purpose="analysis")
    assert rows[0]["design"] == "~ subject + condition"
    assert audit[0]["rmats_group1_samples"] == ["treated1", "treated2", "treated3"]
    assert audit[0]["rmats_group2_samples"] == ["control1", "control2", "control3"]
    assert audit[0]["pair_ids"] == ["mouse1", "mouse2", "mouse3"]


@pytest.mark.parametrize("problem", ["missing", "duplicate", "mismatch", "unpaired", "conflicting_formula"])
def test_invalid_pairing_cannot_silently_fall_back(problem):
    records, cmp = samples(), contrast()
    if problem == "missing": records[0]["subject"] = ""
    if problem == "duplicate": records[0]["subject"] = records[1]["subject"]
    if problem == "mismatch": records[0]["subject"] = "other_mouse"
    if problem == "unpaired": cmp["paired"] = "false"
    if problem == "conflicting_formula": cmp["design"] = "~ condition"
    with pytest.raises(validator.InputError):
        validator.validate_contrasts([cmp], records, design="~ condition", analysis_purpose="analysis")


def paired_project():
    rows = samples()
    return ProjectManifest.model_validate(dict(name="pairs", input_root=".", output_root=".",
        annotation={"genome_build": "test"},
        samples=[dict(sample_id=r["sample"], biological_sample_id=r["subject"], condition=r["condition"]) for r in rows],
        contrasts=[dict(contrast_id="paired", label="paired", numerator=dict(label="treated", sample_ids=[r["sample"] for r in rows if r["condition"]=="treated"]),
                        denominator=dict(label="control", sample_ids=[r["sample"] for r in rows if r["condition"]=="control"]))],
        metadata={"comparison_designs": {"paired": {"paired": True, "formula": "~ subject + condition"}}}))


def test_network_counts_pairs_separately_from_observations():
    manifest = paired_project()
    identity = network_scope_identity(manifest, {"network_scope": "contrast", "contrast_ids": ["paired"]})
    assert identity["n_samples"] == 6 and identity["n_independent_units"] == 3
    assert identity["paired"] is True
    unmarked = manifest.model_copy(update={"metadata": {}})
    with pytest.raises(ValueError, match="independent"):
        network_scope_identity(unmarked, {"network_scope": "contrast", "contrast_ids": ["paired"]})


def test_module_blocked_pvalue_and_direction_match_paired_test():
    manifest = paired_project()
    metadata = pd.DataFrame([s.model_dump() for s in manifest.samples])
    values = {"control1": 4., "control2": -2., "control3": 10.,
              "treated1": 5., "treated2": .5, "treated3": 11.5}
    scores = pd.DataFrame({"blue": values}).rename_axis("sample_id")
    full, selected, audit = matched_module_associations(scores, metadata, [manifest.contrasts[0].model_dump()])
    assert audit[0]["status"] == "completed"
    assert audit[0]["adapter_audit"]["n_independent_units"] == 3
    expected = ttest_rel([5., .5, 11.5], [4., -2., 10.])
    assert selected.iloc[0].pvalue == pytest.approx(expected.pvalue, rel=1e-8)
    assert selected.iloc[0].effect == pytest.approx(5/3)
    condition = full.loc[full.get("model_term", pd.Series(index=full.index, dtype=str)).eq('C(Q("declared_contrast_arm"), Sum)')]
    assert selected.iloc[0].FDR == condition.iloc[0].FDR
