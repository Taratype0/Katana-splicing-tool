from __future__ import annotations

import numpy as np
import pandas as pd

from katana_cosplice.bayesian_preservation import (
    bayesian_multi_cohort_preservation,
)
from katana_cosplice.config import UncertaintyConfig
from katana_cosplice.differential_network import compare_condition_networks
from katana_cosplice.multiview import fuse_event_edge_views
from katana_cosplice.uncertainty import estimate_beta_binomial_uncertainty


def test_uncertainty_is_a_required_production_layer_by_default() -> None:
    settings = UncertaintyConfig()

    assert settings.enabled is True
    assert settings.network_mode == "posterior"
    assert settings.filter_low_reliability is True
    assert settings.min_network_reliability == 0.15


def test_uncertainty_downweights_low_depth_inconsistent_events() -> None:
    samples = pd.DataFrame(
        {
            "sample_id": [f"s{index}" for index in range(8)],
            "condition": ["A"] * 4 + ["B"] * 4,
        }
    )
    ijc = pd.DataFrame(
        [
            [80, 82, 78, 81, 20, 18, 22, 19],
            [1, 4, 0, 3, 4, 0, 3, 1],
        ],
        index=["high_depth", "low_depth"],
        columns=samples["sample_id"],
    )
    sjc = pd.DataFrame(
        [
            [20, 18, 22, 19, 80, 82, 78, 81],
            [4, 1, 5, 2, 1, 5, 2, 4],
        ],
        index=ijc.index,
        columns=ijc.columns,
    )
    qc = pd.DataFrame(
        {"stability_score": [0.95, 0.35]},
        index=ijc.index,
    )
    result = estimate_beta_binomial_uncertainty(ijc, sjc, samples, qc)
    metrics = result.event_metrics
    assert metrics.loc["high_depth", "network_reliability"] > metrics.loc[
        "low_depth", "network_reliability"
    ]
    assert metrics.loc["high_depth", "median_posterior_sd"] < metrics.loc[
        "low_depth", "median_posterior_sd"
    ]
    assert result.posterior_psi.shape == ijc.shape


def test_differential_network_recovers_condition_specific_density() -> None:
    rng = np.random.default_rng(4)
    sample_ids = [f"A{index}" for index in range(8)] + [
        f"B{index}" for index in range(8)
    ]
    latent = np.linspace(-1.5, 1.5, 8)
    values_a = rng.normal(0.0, 1.0, size=(6, 8))
    values_b = np.vstack(
        [latent + rng.normal(0.0, 0.03, size=8) for _ in range(6)]
    )
    matrix = pd.DataFrame(
        np.hstack([values_a, values_b]),
        index=[f"e{index}" for index in range(6)],
        columns=sample_ids,
    )
    samples = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "condition": ["A"] * 8 + ["B"] * 8,
        }
    )
    assignments = pd.DataFrame(
        {
            "event_key": matrix.index,
            "consensus_module": ["C01"] * len(matrix),
        }
    )
    result = compare_condition_networks(
        matrix,
        samples,
        assignments,
        permutations=99,
        random_seed=12,
    )
    row = result.module_changes.iloc[0]
    assert row["density_delta_B_minus_A"] > 0.25
    assert row["differential_network_state"] in {
        "gained_in_B",
        "rewired",
        "density_changed_and_rewired",
    }
    assert bool(row["has_significant_density_change"])
    assert len(result.event_rewiring) == len(matrix)


def test_multiview_reweights_only_existing_as_edges() -> None:
    edges = pd.DataFrame(
        {
            "event_a": ["e1", "e1", "e2", "e3"],
            "event_b": ["e2", "e3", "e3", "e4"],
            "adjacency": [0.8, 0.7, 0.6, 0.75],
        }
    )
    features = pd.DataFrame(
        {
            "event_key": ["e1", "e2", "e3", "e4", "e5"],
            "rbp_score": [1.0, 1.1, 4.0, 4.1, 1.0],
            "dtu_score": [0.2, 0.3, 0.8, 0.9, 0.2],
        }
    )
    result = fuse_event_edge_views(
        edges,
        features,
        auxiliary_columns=["rbp_score", "dtu_score"],
        as_weight=0.8,
        edge_threshold=0.0,
        min_module_size=2,
    )
    observed_pairs = set(
        zip(result.fused_edges["event_a"], result.fused_edges["event_b"], strict=True)
    )
    assert observed_pairs == set(zip(edges["event_a"], edges["event_b"], strict=True))
    assert "e5" not in set(result.assignments["event_key"])
    assert result.fused_edges["AS_primary_weight"].ge(0.8).all()
    assert result.fused_edges["fused_adjacency"].between(0.0, 1.0).all()


def test_bayesian_preservation_combines_multiple_cohorts() -> None:
    preservation = pd.DataFrame(
        {
            "cohort": ["T1", "T2", "T3", "T1", "T2", "T3"],
            "consensus_module": ["C01"] * 3 + ["C02"] * 3,
            "Zdensity": [8.0, 7.0, 9.0, 0.2, -0.1, 0.4],
            "Zconnectivity": [7.0, 6.5, 8.0, 0.1, 0.3, -0.2],
            "density_empirical_pvalue": [0.002, 0.004, 0.001, 0.8, 0.7, 0.9],
            "connectivity_empirical_pvalue": [
                0.003,
                0.005,
                0.002,
                0.9,
                0.8,
                0.7,
            ],
            "program_score_preservation": [0.75, 0.70, 0.80, 0.02, -0.05, 0.01],
            "membership_profile_correlation": [
                0.70,
                0.65,
                0.72,
                0.01,
                0.03,
                -0.02,
            ],
        }
    )
    result = bayesian_multi_cohort_preservation(
        preservation,
        grid_size=101,
    )
    modules = result.module_posteriors.set_index("consensus_module")
    assert modules.loc[
        "C01", "conservative_probability_partial_preservation"
    ] > modules.loc["C02", "conservative_probability_partial_preservation"]
    assert set(result.metric_posteriors["metric"]) == {
        "density",
        "connectivity",
        "program",
        "membership",
    }
