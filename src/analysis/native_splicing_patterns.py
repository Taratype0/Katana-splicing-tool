from __future__ import annotations

import math

import numpy as np
import pandas as pd


def sample_similarity_tables(
    psi_matrix: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return pairwise PSI correlations and the measured-event denominator."""
    if psi_matrix.empty or psi_matrix.shape[1] <= 2:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    event_column = psi_matrix.columns[0]
    numeric = psi_matrix.drop(columns=[event_column]).apply(pd.to_numeric, errors="coerce")
    samples = numeric.columns.astype(str).tolist()
    correlation = pd.DataFrame(np.nan, index=samples, columns=samples, dtype=float)
    shared_n = pd.DataFrame(0, index=samples, columns=samples, dtype=int)
    rows: list[dict[str, object]] = []

    for left_index, left in enumerate(samples):
        for right_index in range(left_index, len(samples)):
            right = samples[right_index]
            paired = numeric[[left, right]].dropna()
            n_shared = int(len(paired))
            value = float("nan")
            if left == right and n_shared:
                value = 1.0
            elif n_shared >= 3 and paired[left].nunique() > 1 and paired[right].nunique() > 1:
                value = float(paired[left].corr(paired[right], method="spearman"))
            correlation.loc[left, right] = value
            correlation.loc[right, left] = value
            shared_n.loc[left, right] = n_shared
            shared_n.loc[right, left] = n_shared
            if left != right:
                rows.append(
                    {
                        "sample_A": left,
                        "sample_B": right,
                        "spearman_PSI_similarity": value,
                        "shared_measured_events": n_shared,
                        "similarity_confidence": _overlap_confidence(n_shared),
                        "interpretation": (
                            "Interpret correlation together with shared_measured_events; "
                            "a high correlation based on few events is low-confidence."
                        ),
                    }
                )

    return (
        pd.DataFrame(rows),
        correlation.reset_index(names="sample_id"),
        shared_n.reset_index(names="sample_id"),
    )


def replicate_group_separation(psi_long: pd.DataFrame, similarity_long: pd.DataFrame) -> pd.DataFrame:
    """Summarize whether within-group replicates are more similar than between-group samples."""
    identity_column = "biological_sample_id" if "biological_sample_id" in psi_long.columns else "sample_id"
    required = {identity_column, "comparison_id", "sample_group"}
    if psi_long.empty or similarity_long.empty or not required.issubset(psi_long.columns):
        return pd.DataFrame()

    sample_meta = psi_long[[identity_column, "comparison_id", "sample_group"]].drop_duplicates()
    rows: list[dict[str, object]] = []
    for comparison_id, comparison_samples in sample_meta.groupby("comparison_id", sort=False):
        comparison_samples = comparison_samples.drop_duplicates(identity_column).set_index(identity_column)
        sample_ids = set(comparison_samples.index.astype(str))
        subset = similarity_long.loc[
            similarity_long["sample_A"].astype(str).isin(sample_ids)
            & similarity_long["sample_B"].astype(str).isin(sample_ids)
        ].copy()
        if subset.empty:
            continue
        subset["group_A"] = subset["sample_A"].map(comparison_samples["sample_group"])
        subset["group_B"] = subset["sample_B"].map(comparison_samples["sample_group"])
        subset["within_group"] = subset["group_A"].astype(str) == subset["group_B"].astype(str)
        values = pd.to_numeric(subset["spearman_PSI_similarity"], errors="coerce")
        within = values.loc[subset["within_group"]].dropna()
        between = values.loc[~subset["within_group"]].dropna()
        within_median = float(within.median()) if not within.empty else float("nan")
        between_median = float(between.median()) if not between.empty else float("nan")
        rows.append(
            {
                "comparison_id": str(comparison_id),
                "n_within_group_pairs": int(len(within)),
                "n_between_group_pairs": int(len(between)),
                "median_within_group_spearman": within_median,
                "median_between_group_spearman": between_median,
                "replicate_separation_delta": within_median - between_median,
                "replicate_qc_interpretation": (
                    "Positive values mean same-group replicates are more similar than between-group samples. "
                    "This is descriptive QC, not a differential-splicing significance test."
                ),
            }
        )
    return pd.DataFrame(rows)


def psi_pattern_clusters(
    psi_matrix: pd.DataFrame,
    event_metadata: pd.DataFrame,
    *,
    max_clusters: int = 8,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Cluster events by row-standardized sample PSI trajectories."""
    if psi_matrix.empty or psi_matrix.shape[1] <= 2:
        return pd.DataFrame(), pd.DataFrame()

    key_column = psi_matrix.columns[0]
    numeric = psi_matrix.drop(columns=[key_column]).apply(pd.to_numeric, errors="coerce")
    measured = numeric.notna().sum(axis=1)
    numeric = numeric.loc[measured >= 2].copy()
    keys = psi_matrix.loc[numeric.index, key_column].astype(str)
    if len(numeric) < 3:
        return pd.DataFrame(), pd.DataFrame()

    filled = numeric.T.fillna(numeric.mean(axis=1)).T
    row_mean = filled.mean(axis=1)
    row_std = filled.std(axis=1).replace(0, np.nan)
    standardized = filled.sub(row_mean, axis=0).div(row_std, axis=0)
    keep = standardized.notna().any(axis=1)
    standardized = standardized.loc[keep].fillna(0.0)
    keys = keys.loc[keep]
    if len(standardized) < 3:
        return pd.DataFrame(), pd.DataFrame()

    n_clusters = min(
        max(2, int(max_clusters)),
        max(2, min(8, int(round(math.sqrt(len(standardized)) / 4.0)))),
        len(standardized),
    )
    labels = _hierarchical_cluster_labels(standardized, n_clusters)
    output = pd.concat(
        [
            pd.DataFrame({key_column: keys.to_numpy(), "psi_pattern_cluster": labels}),
            standardized.reset_index(drop=True),
        ],
        axis=1,
    )
    if not event_metadata.empty and "event_matching_key" in event_metadata.columns:
        metadata_columns = [
            column
            for column in ["event_matching_key", "gene_symbol", "event_type", "event_id"]
            if column in event_metadata.columns
        ]
        metadata = event_metadata[metadata_columns].drop_duplicates("event_matching_key")
        output = output.merge(
            metadata,
            left_on=key_column,
            right_on="event_matching_key",
            how="left",
        )
        if key_column != "event_matching_key":
            output = output.drop(columns=["event_matching_key"])

    sample_columns = standardized.columns.astype(str).tolist()
    cluster_means = output.groupby("psi_pattern_cluster", dropna=False)[sample_columns].mean().reset_index()
    counts = output.groupby("psi_pattern_cluster", dropna=False).size().rename("n_events").reset_index()
    event_types = (
        output.groupby("psi_pattern_cluster", dropna=False)["event_type"]
        .agg(lambda values: "; ".join(f"{key}:{value}" for key, value in values.value_counts().items()))
        .rename("event_type_counts")
        .reset_index()
        if "event_type" in output.columns
        else pd.DataFrame()
    )
    summary = counts.merge(cluster_means, on="psi_pattern_cluster", how="left")
    if not event_types.empty:
        summary = summary.merge(event_types, on="psi_pattern_cluster", how="left")
    summary["cluster_method"] = (
        "Weighted-linkage hierarchical clustering on row-z-scored sample/replicate PSI trajectories; "
        f"k={n_clusters}. Descriptive pattern discovery, not a causal program model."
    )
    return output, summary


def row_zscore_matrix(psi_matrix: pd.DataFrame) -> pd.DataFrame:
    if psi_matrix.empty or psi_matrix.shape[1] <= 1:
        return pd.DataFrame()
    key_column = psi_matrix.columns[0]
    numeric = psi_matrix.drop(columns=[key_column]).apply(pd.to_numeric, errors="coerce")
    means = numeric.mean(axis=1)
    stds = numeric.std(axis=1).replace(0, np.nan)
    zscores = numeric.sub(means, axis=0).div(stds, axis=0)
    return pd.concat([psi_matrix[[key_column]].reset_index(drop=True), zscores.reset_index(drop=True)], axis=1)


def _hierarchical_cluster_labels(matrix: pd.DataFrame, n_clusters: int) -> list[int]:
    try:
        from scipy.cluster.hierarchy import fcluster, linkage

        linked = linkage(matrix.to_numpy(dtype=float), method="weighted", metric="cityblock")
        return fcluster(linked, t=n_clusters, criterion="maxclust").astype(int).tolist()
    except Exception:
        order = matrix.abs().mean(axis=1).rank(method="first").to_numpy(dtype=float)
        bins = np.linspace(0, len(order), n_clusters + 1)
        return np.clip(np.digitize(order - 1, bins[1:-1]) + 1, 1, n_clusters).astype(int).tolist()


def _overlap_confidence(n_shared: int) -> str:
    if n_shared >= 100:
        return "high"
    if n_shared >= 30:
        return "moderate"
    if n_shared >= 10:
        return "low"
    if n_shared > 0:
        return "very_low"
    return "not_available"
