from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control


@dataclass
class DifferentialNetworkResult:
    module_changes: pd.DataFrame
    event_rewiring: pd.DataFrame
    condition_edges: pd.DataFrame


def _correlation(values: np.ndarray, shrinkage: float) -> np.ndarray:
    centered = values - np.nanmean(values, axis=1, keepdims=True)
    centered = np.nan_to_num(centered, nan=0.0)
    norm = np.sqrt(np.sum(centered**2, axis=1, keepdims=True))
    standardized = centered / np.where(norm < 1e-12, 1.0, norm)
    correlation = np.clip(standardized @ standardized.T, -1.0, 1.0)
    n_samples = values.shape[1]
    weight = n_samples / (n_samples + shrinkage)
    correlation *= weight
    np.fill_diagonal(correlation, 1.0)
    return correlation


def _adjacency(values: np.ndarray, soft_power: int, shrinkage: float) -> np.ndarray:
    correlation = _correlation(values, shrinkage)
    adjacency = np.power(np.clip((1.0 + correlation) / 2.0, 0.0, 1.0), soft_power)
    np.fill_diagonal(adjacency, 0.0)
    return adjacency


def _module_statistics(adjacency: np.ndarray) -> tuple[float, np.ndarray]:
    if len(adjacency) < 2:
        return np.nan, np.full(len(adjacency), np.nan)
    upper = adjacency[np.triu_indices_from(adjacency, k=1)]
    return float(np.mean(upper)), np.mean(adjacency, axis=1)


def compare_condition_networks(
    matrix: pd.DataFrame,
    samples: pd.DataFrame,
    assignments: pd.DataFrame,
    *,
    condition_column: str = "condition",
    min_samples_per_condition: int = 3,
    soft_power: int = 6,
    shrinkage: float = 5.0,
    permutations: int = 200,
    max_output_edges_per_module_pair: int = 5_000,
    random_seed: int = 71,
) -> DifferentialNetworkResult:
    """Test condition-specific module density and topology changes.

    The global module labels define a fixed event universe. Separate
    shrinkage-correlation networks are built within each condition. Label
    permutations test density change and edge-pattern rewiring without
    pretending that three-replicate condition-specific community detection is
    stable enough to define new modules.
    """
    if not {"sample_id", condition_column}.issubset(samples.columns):
        raise ValueError(f"Sample metadata requires sample_id and {condition_column}.")
    module_column = "consensus_module" if "consensus_module" in assignments else "module"
    if "event_key" not in assignments or module_column not in assignments:
        raise ValueError("Assignments require event_key and module labels.")
    ordered = samples.set_index("sample_id").reindex(matrix.columns)
    if ordered[condition_column].isna().any():
        raise ValueError("Matrix samples are missing condition metadata.")
    labels = ordered[condition_column].astype(str)
    counts = labels.value_counts()
    conditions = sorted(counts[counts >= min_samples_per_condition].index)
    if len(conditions) < 2:
        raise ValueError("Need at least two conditions with sufficient replicates.")
    assignment_map = assignments.drop_duplicates("event_key").set_index("event_key")[module_column]
    common = matrix.index.intersection(assignment_map.index)
    matrix = matrix.loc[common]
    module_labels = assignment_map.reindex(common).astype(str)
    rng = np.random.default_rng(random_seed)
    module_rows: list[dict[str, object]] = []
    event_rows: list[dict[str, object]] = []
    edge_rows: list[dict[str, object]] = []

    for condition_a, condition_b in combinations(conditions, 2):
        sample_a = labels.index[labels.eq(condition_a)]
        sample_b = labels.index[labels.eq(condition_b)]
        combined = np.concatenate([sample_a.to_numpy(), sample_b.to_numpy()])
        n_a = len(sample_a)
        for module in sorted(set(module_labels).difference({"grey"})):
            keys = module_labels.index[module_labels.eq(module)]
            if len(keys) < 3:
                continue
            values = matrix.loc[keys, combined].to_numpy(dtype=float)
            adj_a = _adjacency(values[:, :n_a], soft_power, shrinkage)
            adj_b = _adjacency(values[:, n_a:], soft_power, shrinkage)
            density_a, connectivity_a = _module_statistics(adj_a)
            density_b, connectivity_b = _module_statistics(adj_b)
            density_delta = density_b - density_a
            edge_delta = adj_b - adj_a
            topology_delta = float(
                np.sqrt(np.mean(edge_delta[np.triu_indices_from(edge_delta, k=1)] ** 2))
            )
            null_density: list[float] = []
            null_topology: list[float] = []
            for _ in range(permutations):
                permutation = rng.permutation(values.shape[1])
                perm_a = _adjacency(values[:, permutation[:n_a]], soft_power, shrinkage)
                perm_b = _adjacency(values[:, permutation[n_a:]], soft_power, shrinkage)
                perm_density_a, _ = _module_statistics(perm_a)
                perm_density_b, _ = _module_statistics(perm_b)
                null_density.append(perm_density_b - perm_density_a)
                perm_delta = perm_b - perm_a
                null_topology.append(
                    float(
                        np.sqrt(
                            np.mean(
                                perm_delta[np.triu_indices_from(perm_delta, k=1)] ** 2
                            )
                        )
                    )
                )
            density_p = (1 + np.sum(np.abs(null_density) >= abs(density_delta))) / (
                permutations + 1
            )
            topology_p = (1 + np.sum(np.asarray(null_topology) >= topology_delta)) / (
                permutations + 1
            )
            pair = f"{condition_a}__vs__{condition_b}"
            module_rows.append(
                {
                    "comparison_pair": pair,
                    "condition_A": condition_a,
                    "condition_B": condition_b,
                    "module": module,
                    "n_events": len(keys),
                    "n_samples_A": len(sample_a),
                    "n_samples_B": len(sample_b),
                    "density_A": density_a,
                    "density_B": density_b,
                    "density_delta_B_minus_A": density_delta,
                    "topology_RMS_delta": topology_delta,
                    "density_permutation_pvalue": density_p,
                    "topology_permutation_pvalue": topology_p,
                    "network_method": "condition_specific_shrinkage_correlation",
                }
            )
            event_delta = np.mean(np.abs(edge_delta), axis=1)
            for index, event_key in enumerate(keys):
                event_rows.append(
                    {
                        "comparison_pair": pair,
                        "condition_A": condition_a,
                        "condition_B": condition_b,
                        "module": module,
                        "event_key": event_key,
                        "connectivity_A": connectivity_a[index],
                        "connectivity_B": connectivity_b[index],
                        "connectivity_delta_B_minus_A": (
                            connectivity_b[index] - connectivity_a[index]
                        ),
                        "mean_absolute_edge_rewiring": event_delta[index],
                    }
                )
            upper = np.triu_indices_from(adj_a, k=1)
            edge_candidates = pd.DataFrame(
                {
                    "row": upper[0],
                    "col": upper[1],
                    "absolute_delta": np.abs(edge_delta[upper]),
                }
            ).nlargest(max_output_edges_per_module_pair, "absolute_delta")
            for candidate in edge_candidates.itertuples(index=False):
                row = int(candidate.row)
                col = int(candidate.col)
                edge_rows.append(
                    {
                        "comparison_pair": pair,
                        "module": module,
                        "event_a": keys[row],
                        "event_b": keys[col],
                        "adjacency_A": adj_a[row, col],
                        "adjacency_B": adj_b[row, col],
                        "adjacency_delta_B_minus_A": edge_delta[row, col],
                    }
                )

    modules = pd.DataFrame(module_rows)
    if not modules.empty:
        modules["density_pairwise_FDR"] = modules.groupby(
            "comparison_pair", observed=True
        )[
            "density_permutation_pvalue"
        ].transform(lambda values: false_discovery_control(values.to_numpy()))
        modules["topology_pairwise_FDR"] = modules.groupby(
            "comparison_pair", observed=True
        )[
            "topology_permutation_pvalue"
        ].transform(lambda values: false_discovery_control(values.to_numpy()))
        modules["density_FDR"] = false_discovery_control(
            modules["density_permutation_pvalue"].to_numpy()
        )
        modules["topology_FDR"] = false_discovery_control(
            modules["topology_permutation_pvalue"].to_numpy()
        )
        significant_density = modules["density_FDR"].le(0.05)
        significant_topology = modules["topology_FDR"].le(0.05)
        modules["has_significant_density_change"] = significant_density
        modules["has_significant_topology_rewiring"] = significant_topology
        modules["density_change_direction"] = np.select(
            [
                significant_density & modules["density_delta_B_minus_A"].gt(0),
                significant_density & modules["density_delta_B_minus_A"].lt(0),
            ],
            ["gained_in_B", "lost_in_B"],
            default="not_significant",
        )
        modules["differential_network_state"] = np.select(
            [
                significant_density & significant_topology,
                significant_density & modules["density_delta_B_minus_A"].gt(0),
                significant_density & modules["density_delta_B_minus_A"].lt(0),
                significant_topology,
            ],
            [
                "density_changed_and_rewired",
                "gained_in_B",
                "lost_in_B",
                "rewired",
            ],
            default="stable_or_underpowered",
        )
    return DifferentialNetworkResult(
        modules,
        pd.DataFrame(event_rows),
        pd.DataFrame(edge_rows),
    )
