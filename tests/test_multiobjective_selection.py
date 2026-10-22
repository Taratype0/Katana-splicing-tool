from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from katana_cosplice.config import (
    CoSpliceConfig,
    ExternalCohortConfig,
    MultiObjectiveSelectionConfig,
)
from katana_cosplice.consensus import _finalize_candidate_objectives
from katana_cosplice.multiobjective_selection import (
    build_multiobjective_candidate_scorer,
)


def _write_inputs(tmp_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, Path, Path]:
    rng = np.random.default_rng(44)
    pattern_a = np.linspace(-1.0, 1.0, 8)
    pattern_b = np.linspace(1.0, -1.0, 8)
    reference = pd.DataFrame(
        np.vstack(
            [
                pattern_a + rng.normal(0, 0.02, 8)
                for _ in range(6)
            ]
            + [
                pattern_b + rng.normal(0, 0.02, 8)
                for _ in range(6)
            ]
        ),
        index=[f"event_{index:02d}" for index in range(12)],
        columns=[f"R{index}" for index in range(8)],
    )
    target = pd.DataFrame(
        np.vstack(
            [
                pattern_a + rng.normal(0, 0.03, 8)
                for _ in range(6)
            ]
            + [
                pattern_b + rng.normal(0, 0.03, 8)
                for _ in range(6)
            ]
        ),
        index=reference.index,
        columns=[f"T{index}" for index in range(8)],
    )
    target_path = tmp_path / "tuning.tsv"
    target.to_csv(target_path, sep="\t")
    samples_path = tmp_path / "tuning_samples.tsv"
    pd.DataFrame(
        {
            "sample_id": target.columns,
            "condition": ["A"] * 4 + ["B"] * 4,
        }
    ).to_csv(samples_path, sep="\t", index=False)
    metadata = pd.DataFrame(
        {
            "event_key": reference.index,
            "gene_id": [f"ENSG{index:05d}" for index in range(12)],
            "gene_symbol": [f"Gene{index}" for index in range(12)],
            "event_type": ["SE"] * 12,
        }
    )
    gmt = tmp_path / "tuning.gmt"
    gmt.write_text(
        "PROGRAM_A\tdescription\t"
        + "\t".join(f"Gene{index}" for index in range(6))
        + "\nPROGRAM_B\tdescription\t"
        + "\t".join(f"Gene{index}" for index in range(6, 12))
        + "\n",
        encoding="utf-8",
    )
    return reference, metadata, target_path, samples_path


def test_multiobjective_scorer_rewards_preserved_functional_partition(
    tmp_path: Path,
) -> None:
    reference, metadata, target_path, samples_path = _write_inputs(tmp_path)
    settings = MultiObjectiveSelectionConfig(
        enabled=True,
        tuning_cohorts=["tuning"],
        holdout_cohorts=["holdout"],
        functional_gene_sets_gmt=tmp_path / "tuning.gmt",
        min_shared_events_per_module=4,
        min_module_genes=3,
    )
    scorer = build_multiobjective_candidate_scorer(
        reference,
        metadata,
        [
            ExternalCohortConfig("tuning", target_path, samples_path),
            # A missing holdout file must not be read during module selection.
            ExternalCohortConfig(
                "holdout",
                tmp_path / "missing_holdout.tsv",
                tmp_path / "missing_holdout_samples.tsv",
            ),
        ],
        settings,
    )
    coherent = np.asarray(["C01"] * 6 + ["C02"] * 6, dtype=object)
    mixed = np.asarray(["C01", "C02"] * 6, dtype=object)
    coherent_score = scorer(coherent)
    mixed_score = scorer(mixed)
    assert coherent_score["objective_external_preservation"] > 0
    assert mixed_score["objective_external_preservation"] > 0
    assert (
        coherent_score["objective_functional_coherence"]
        > mixed_score["objective_functional_coherence"]
    )
    audit = scorer.audit(coherent)
    assert set(audit["audit_type"]) == {
        "external_preservation",
        "functional_coherence",
    }
    assert set(
        audit.loc[
            audit["audit_type"].eq("external_preservation"), "cohort"
        ]
    ) == {"tuning"}


def test_multiobjective_config_rejects_tuning_holdout_leakage(
    tmp_path: Path,
) -> None:
    config = CoSpliceConfig(
        project_name="test",
        rmats_root=tmp_path / "rmats",
        sample_metadata=tmp_path / "samples.tsv",
        output_dir=tmp_path / "output",
        external_cohorts=[
            ExternalCohortConfig(
                "shared",
                tmp_path / "psi.tsv",
                tmp_path / "external_samples.tsv",
            )
        ],
        multiobjective_selection=MultiObjectiveSelectionConfig(
            enabled=True,
            tuning_cohorts=["shared"],
            holdout_cohorts=["shared"],
            functional_gene_sets_gmt=tmp_path / "terms.gmt",
        ),
    )
    with pytest.raises(ValueError, match="must be disjoint"):
        config.validate()


def test_multiobjective_selection_requires_all_three_objectives() -> None:
    candidates = pd.DataFrame(
        {
            "robust_bootstrap_stability": [0.95, 0.78, 0.70],
            "objective_external_preservation": [0.10, 0.72, 0.60],
            "objective_functional_coherence": [0.10, 0.70, 0.55],
            "objective_stability_weight": [0.40] * 3,
            "objective_preservation_weight": [0.35] * 3,
            "objective_functional_weight": [0.25] * 3,
            "objective_score_floor": [0.01] * 3,
            "n_preservation_module_cohort_tests": [2] * 3,
            "n_functional_modules_scored": [2] * 3,
        }
    )
    scored = _finalize_candidate_objectives(candidates, True)
    assert scored.loc[1, "multiobjective_selection_score"] > scored.loc[
        0, "multiobjective_selection_score"
    ]
    assert bool(scored.loc[0, "pareto_optimal"])
    assert bool(scored.loc[1, "pareto_optimal"])
    assert not bool(scored.loc[2, "pareto_optimal"])


def test_multiobjective_selection_disables_functional_weight_without_fdr_signal() -> None:
    candidates = pd.DataFrame(
        {
            "robust_bootstrap_stability": [0.80, 0.70],
            "objective_external_preservation": [0.60, 0.75],
            "objective_functional_coherence": [0.40, 0.60],
            "objective_stability_weight": [0.40, 0.40],
            "objective_preservation_weight": [0.35, 0.35],
            "objective_functional_weight": [0.25, 0.25],
            "functional_weight_mode": ["auto_significance"] * 2,
            "min_significant_functional_modules": [1, 1],
            "objective_score_floor": [0.01, 0.01],
            "n_preservation_module_cohort_tests": [2, 2],
            "n_functional_modules_scored": [2, 2],
            "n_significant_functional_modules": [0, 0],
            "n_modules": [3, 3],
            "grey_fraction": [0.10, 0.10],
            "module_balance_entropy": [0.80, 0.80],
            "mean_membership_confidence": [0.85, 0.85],
            "mean_confidence_margin": [0.50, 0.50],
        }
    )
    scored = _finalize_candidate_objectives(candidates, True)
    assert not scored["functional_objective_active"].any()
    assert (scored["effective_functional_weight"] == 0).all()
    assert scored["effective_stability_weight"].iloc[0] == pytest.approx(
        0.40 / 0.75
    )
    assert scored["effective_preservation_weight"].iloc[0] == pytest.approx(
        0.35 / 0.75
    )


def test_multiobjective_selection_falls_back_when_preservation_is_not_evaluable() -> None:
    candidates = pd.DataFrame(
        {
            "robust_bootstrap_stability": [0.80, 0.65],
            "objective_external_preservation": [0.0, 0.0],
            "objective_functional_coherence": [0.0, 0.0],
            "objective_stability_weight": [0.40, 0.40],
            "objective_preservation_weight": [0.35, 0.35],
            "objective_functional_weight": [0.25, 0.25],
            "functional_weight_mode": ["auto_significance"] * 2,
            "min_significant_functional_modules": [1, 1],
            "objective_score_floor": [0.01, 0.01],
            "n_preservation_module_cohort_tests": [0, 0],
            "n_functional_modules_scored": [0, 0],
            "n_significant_functional_modules": [0, 0],
            "n_modules": [1, 1],
            "grey_fraction": [0.40, 0.40],
            "module_balance_entropy": [0.0, 0.0],
            "mean_membership_confidence": [0.90, 0.90],
            "mean_confidence_margin": [0.90, 0.90],
        }
    )

    scored = _finalize_candidate_objectives(candidates, True)

    assert not scored["preservation_objective_active"].any()
    assert not scored["functional_objective_active"].any()
    assert (scored["effective_stability_weight"] == 1.0).all()
    assert (scored["effective_preservation_weight"] == 0.0).all()
