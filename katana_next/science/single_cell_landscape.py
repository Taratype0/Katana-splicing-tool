from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


@dataclass(frozen=True)
class SingleCellLandscapeResult:
    event_qc: pd.DataFrame
    selected_psi: pd.DataFrame
    embedding: pd.DataFrame
    similarity: pd.DataFrame
    heatmap: pd.DataFrame
    diagnostics: pd.DataFrame


def build_single_cell_landscape(
    pseudobulk: pd.DataFrame,
    *,
    min_total_reads: float = 10.0,
    min_unit_fraction: float = 0.70,
    min_psi_variance: float = 0.0025,
    max_events: int = 500,
    heatmap_events: int = 60,
    random_seed: int = 17,
) -> SingleCellLandscapeResult:
    """Build a donor-aware PSI landscape without using cells as replicates."""

    required = {
        "event_id",
        "donor_id",
        "cell_type",
        "condition",
        "total",
        "psi_posterior",
        "measurement_se",
    }
    missing = sorted(required.difference(pseudobulk.columns))
    if missing:
        raise ValueError(f"Pseudobulk landscape input is missing columns: {missing}")
    if max_events < 2:
        raise ValueError("max_events must be at least 2.")

    frame = pseudobulk.copy()
    frame["unit_id"] = (
        frame["donor_id"].astype(str) + "::" + frame["cell_type"].astype(str)
        + "::" + frame["condition"].astype(str)
    )
    if frame.duplicated(["event_id", "unit_id"]).any():
        raise ValueError("Duplicate donor-condition-cell-type event observations; aggregate technical counts explicitly.")
    unit_metadata = (
        frame.groupby("unit_id", as_index=False)
        .agg(
            donor_id=("donor_id", "first"),
            cell_type=("cell_type", "first"),
            condition=("condition", "first"),
            batch=("batch", "first"),
            assay=("assay", "first"),
        )
        .set_index("unit_id")
    )
    frame["measured"] = frame["total"].ge(min_total_reads)
    psi = frame.pivot_table(
        index="event_id", columns="unit_id", values="psi_posterior", aggfunc="first"
    ).reindex(columns=unit_metadata.index)
    total = frame.pivot_table(
        index="event_id", columns="unit_id", values="total", aggfunc="sum"
    ).reindex(index=psi.index, columns=psi.columns)
    se = frame.pivot_table(
        index="event_id", columns="unit_id", values="measurement_se", aggfunc="first"
    ).reindex(index=psi.index, columns=psi.columns)
    measured = total.ge(min_total_reads)
    observed_psi = psi.where(measured)

    event_metadata = (
        frame.sort_values("event_id").drop_duplicates("event_id").set_index("event_id")
    )
    valid_fraction = measured.mean(axis=1)
    variance = observed_psi.var(axis=1, ddof=1)
    median_reads = total.where(measured).median(axis=1)
    median_se = se.where(measured).median(axis=1)
    read_reliability = median_reads / (median_reads + 20.0)
    reliability = (
        valid_fraction * read_reliability * (1.0 - 2.0 * median_se).clip(0.0, 1.0)
    )
    event_qc = pd.DataFrame(
        {
            "event_id": psi.index.astype(str),
            "gene_id": event_metadata.reindex(psi.index)["gene_id"]
            .fillna("")
            .to_numpy(),
            "gene_symbol": event_metadata.reindex(psi.index)["gene_symbol"]
            .fillna("")
            .to_numpy(),
            "event_type": event_metadata.reindex(psi.index)["event_type"]
            .fillna("junction_choice")
            .to_numpy(),
            "measured_units": measured.sum(axis=1).astype(int).to_numpy(),
            "total_units": len(psi.columns),
            "valid_fraction": valid_fraction.to_numpy(float),
            "psi_variance": variance.to_numpy(float),
            "median_total_reads": median_reads.to_numpy(float),
            "median_measurement_se": median_se.to_numpy(float),
            "network_reliability": reliability.to_numpy(float),
        }
    )
    event_qc["passes_landscape_gate"] = event_qc["valid_fraction"].ge(
        min_unit_fraction
    ) & event_qc["psi_variance"].ge(min_psi_variance)
    event_qc["selection_score"] = event_qc["psi_variance"].fillna(0.0) * event_qc[
        "network_reliability"
    ].fillna(0.0)
    event_qc = event_qc.sort_values(
        ["passes_landscape_gate", "selection_score", "event_id"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    passing = event_qc.loc[event_qc["passes_landscape_gate"], "event_id"].head(
        max_events
    )
    if len(passing) < 2:
        raise ValueError(
            "Fewer than two events passed the donor-unit coverage and variance gates."
        )
    selected = observed_psi.reindex(passing)
    imputed = selected.T.fillna(selected.median(axis=1)).T
    transformed = _logit(imputed)
    standardized = pd.DataFrame(
        StandardScaler().fit_transform(transformed.T).T,
        index=transformed.index,
        columns=transformed.columns,
    )

    embedding = _unit_embedding(
        standardized,
        unit_metadata,
        random_seed=random_seed,
    )
    similarity_matrix = selected.corr(method="spearman", min_periods=3)
    similarity_matrix = similarity_matrix.fillna(0.0).clip(-1.0, 1.0)
    np.fill_diagonal(similarity_matrix.values, 1.0)
    order = _similarity_order(similarity_matrix)
    order_lookup = {unit: index for index, unit in enumerate(order)}
    similarity = (
        similarity_matrix.rename_axis(index="unit_a", columns="unit_b")
        .stack(future_stack=True)
        .rename("spearman_similarity")
        .reset_index()
    )
    similarity["unit_a_order"] = similarity["unit_a"].map(order_lookup)
    similarity["unit_b_order"] = similarity["unit_b"].map(order_lookup)
    similarity["row_id"] = similarity["unit_a"] + "::" + similarity["unit_b"]

    selected_long = (
        selected.rename_axis(index="event_id", columns="unit_id")
        .stack(future_stack=True)
        .rename("psi_posterior")
        .reset_index()
        .merge(unit_metadata.reset_index(), on="unit_id", how="left")
        .merge(
            event_qc[["event_id", "gene_symbol", "network_reliability"]],
            on="event_id",
            how="left",
        )
    )
    selected_long["row_id"] = (
        selected_long["event_id"] + "::" + selected_long["unit_id"]
    )

    top_heatmap = event_qc.loc[event_qc["event_id"].isin(passing)].head(heatmap_events)
    heatmap_matrix = standardized.reindex(top_heatmap["event_id"])
    heatmap = (
        heatmap_matrix.rename_axis(index="event_id", columns="unit_id")
        .stack(future_stack=True)
        .rename("standardized_psi")
        .reset_index()
        .merge(unit_metadata.reset_index(), on="unit_id", how="left")
        .merge(top_heatmap[["event_id", "gene_symbol"]], on="event_id", how="left")
    )
    heatmap["event_label"] = heatmap["gene_symbol"].where(
        heatmap["gene_symbol"].ne(""), heatmap["event_id"]
    )
    heatmap["row_id"] = heatmap["event_id"] + "::" + heatmap["unit_id"]

    diagnostics = pd.DataFrame(
        [
            {
                "row_id": "landscape",
                "analysis_unit": "donor by cell type pseudobulk",
                "total_events": len(event_qc),
                "events_passing_gate": int(event_qc["passes_landscape_gate"].sum()),
                "events_selected": len(passing),
                "independent_donors": int(frame["donor_id"].nunique()),
                "descriptive_units": len(unit_metadata),
                "cell_types": int(frame["cell_type"].nunique()),
                "embedding_method": embedding["embedding_method"].iloc[0],
                "pca_variance_pc1": embedding["pca_variance_pc1"].iloc[0],
                "pca_variance_pc2": embedding["pca_variance_pc2"].iloc[0],
                "claim_boundary": "PCA/UMAP and similarity are descriptive; donors, not cells, define replication.",
            }
        ]
    )
    event_qc["row_id"] = event_qc["event_id"]
    return SingleCellLandscapeResult(
        event_qc=event_qc,
        selected_psi=selected_long,
        embedding=embedding,
        similarity=similarity,
        heatmap=heatmap,
        diagnostics=diagnostics,
    )


def _unit_embedding(
    standardized: pd.DataFrame,
    unit_metadata: pd.DataFrame,
    *,
    random_seed: int,
) -> pd.DataFrame:
    values = standardized.T.to_numpy(float)
    components = min(2, values.shape[0], values.shape[1])
    pca = PCA(n_components=components, random_state=random_seed)
    pca_values = pca.fit_transform(values)
    if components == 1:
        pca_values = np.column_stack([pca_values[:, 0], np.zeros(len(values))])
        explained = np.array([pca.explained_variance_ratio_[0], 0.0])
    else:
        explained = pca.explained_variance_ratio_
    umap_values = pca_values[:, :2]
    method = "PCA fallback"
    if len(values) >= 4:
        try:
            from umap import UMAP

            reducer = UMAP(
                n_components=2,
                n_neighbors=min(15, len(values) - 1),
                min_dist=0.20,
                metric="euclidean",
                random_state=random_seed,
                transform_seed=random_seed,
                n_jobs=1,
            )
            umap_values = reducer.fit_transform(values)
            method = "UMAP on donor-level posterior PSI"
        except (ImportError, ValueError, RuntimeError):
            pass
    output = unit_metadata.reset_index().copy()
    output["pca_1"] = pca_values[:, 0]
    output["pca_2"] = pca_values[:, 1]
    output["umap_1"] = umap_values[:, 0]
    output["umap_2"] = umap_values[:, 1]
    output["embedding_method"] = method
    output["pca_variance_pc1"] = float(explained[0])
    output["pca_variance_pc2"] = float(explained[1])
    output["row_id"] = output["unit_id"]
    return output


def _similarity_order(similarity: pd.DataFrame) -> list[str]:
    if len(similarity) < 3:
        return similarity.index.astype(str).tolist()
    distance = np.clip(1.0 - similarity.to_numpy(float), 0.0, 2.0)
    np.fill_diagonal(distance, 0.0)
    hierarchy = linkage(squareform(distance, checks=False), method="average")
    return similarity.index[leaves_list(hierarchy)].astype(str).tolist()


def _logit(frame: pd.DataFrame) -> pd.DataFrame:
    clipped = frame.clip(0.005, 0.995)
    return np.log(clipped / (1.0 - clipped))
