from __future__ import annotations

import numpy as np
import pandas as pd

from katana_cosplice.independent_validation import (
    fit_module_projection_models,
    project_module_scores,
    tumor_normal_module_association,
)


def test_locked_projection_and_paired_phenotype_association() -> None:
    rng = np.random.default_rng(19)
    events = [f"G{i}|ES|{i}.0" for i in range(8)]
    discovery_samples = [f"D{i}" for i in range(12)]
    discovery = pd.DataFrame(
        rng.normal(size=(8, 12)),
        index=events,
        columns=discovery_samples,
    )
    assignments = pd.DataFrame(
        {
            "event_key": events,
            "consensus_module": ["C01"] * 4 + ["C02"] * 4,
        }
    )
    models = fit_module_projection_models(discovery, assignments)
    assert set(models) == {"C01", "C02"}

    target_samples = [
        "P1",
        "P1_Norm",
        "P2",
        "P2_Norm",
        "P3",
        "P3_Norm",
        "P4",
        "P4_Norm",
    ]
    target = pd.DataFrame(
        rng.normal(size=(8, 8)),
        index=events,
        columns=target_samples,
    )
    target.loc[events[:4], ["P1", "P2", "P3", "P4"]] += 3.0
    scores, audit = project_module_scores(target, models)
    assert scores.columns.tolist() == ["C01", "C02"]
    assert audit["status"].eq("projected").all()

    metadata = pd.DataFrame(
        {
            "sample_id": target_samples,
            "phenotype": [
                "tumor",
                "normal",
                "tumor",
                "normal",
                "tumor",
                "normal",
                "tumor",
                "normal",
            ],
        }
    )
    result = tumor_normal_module_association(scores, metadata)
    assert set(result["module"]) == {"C01", "C02"}
    assert result["n_matched_pairs"].eq(4).all()
    assert {
        "welch_FDR",
        "paired_wilcoxon_FDR",
        "hedges_g",
        "association_evidence_class",
        "robust_all_and_paired_association",
    }.issubset(result)
