from __future__ import annotations

from katana_cosplice.method_extensions_benchmark import (
    run_synthetic_method_ablation,
)


def test_known_truth_method_ablation_has_guardrails() -> None:
    metrics = run_synthetic_method_ablation(seed=17)
    uncertainty = metrics.loc[
        metrics["algorithm"].eq("uncertainty_weighting")
        & metrics["metric"].eq("recoverable_top_edge_precision")
    ].set_index("configuration")["value"]
    assert uncertainty.between(0.0, 1.0).all()
    differential = metrics.loc[
        metrics["algorithm"].eq("differential_network")
    ].set_index("metric")["value"]
    assert differential["known_gained_module_detected"] == 1.0
    assert differential["stable_module_false_positive"] == 0.0
    multiview = metrics.loc[
        metrics["algorithm"].eq("multiview_AS_primary")
        & metrics["metric"].eq("true_edge_AUROC")
    ].set_index("configuration")["value"]
    assert multiview["informative_auxiliary_20pct"] > multiview["AS_only"]
    assert multiview["adversarial_auxiliary_20pct"] >= multiview["AS_only"] - 0.05
