from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from katana_next.science.experimental_design import observation_units, bootstrap_blocks, donor_covariates
from scipy.sparse import csr_matrix, eye
from scipy.spatial.distance import squareform
from scipy.stats import f as f_distribution
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors

from katana_cosplice.config import NetworkConfig
from katana_cosplice.dynamic_tree import dynamic_tree_cut_labels
from katana_cosplice.network import (
    correlation_matrix,
    signed_soft_adjacency,
    topological_overlap_matrix,
)
from katana_next.science.statistics import bh_adjust
from katana_next.services.project_service import ProjectContext


@dataclass(frozen=True)
class SingleCellProgramResult:
    event_catalog: pd.DataFrame
    program_summary: pd.DataFrame
    cell_scores: pd.DataFrame
    donor_scores: pd.DataFrame
    network_edges: pd.DataFrame
    diagnostics: pd.DataFrame


def fit_single_cell_tom_labels(
    matrix: pd.DataFrame,
    *,
    reliability: np.ndarray,
    min_module_size: int,
    soft_power: int,
) -> tuple[np.ndarray, np.ndarray, str]:
    """Public single-cell adapter for the shared signed soft-TOM backbone."""

    return _fit_tom_labels(
        matrix,
        reliability=reliability,
        min_module_size=min_module_size,
        soft_power=soft_power,
        official=True,
    )


def smooth_cell_event_matrix(
    matrix: pd.DataFrame, coordinates: np.ndarray, neighbors: int
) -> pd.DataFrame:
    """Smooth measured PSI over a declared external cell manifold."""

    graph = _knn_graph(coordinates, neighbors)
    values = [
        _neighbor_average(matrix.loc[event_id].to_numpy(float), graph)
        for event_id in matrix.index
    ]
    return pd.DataFrame(values, index=matrix.index, columns=matrix.columns)


def load_declared_cell_embedding(
    context: ProjectContext, metadata: pd.DataFrame
) -> tuple[pd.DataFrame | None, str]:
    """Load an expression-derived cell manifold without deriving it from AS."""

    asset = next(
        (
            item
            for item in context.manifest.assets
            if item.kind == "single_cell_embedding" and item.path.is_file()
        ),
        None,
    )
    if asset is not None:
        frame = pd.read_csv(asset.path, sep=None, engine="python")
        source = f"declared asset: {asset.asset_id}"
    else:
        candidates = [
            column
            for column in metadata.columns
            if column.casefold().startswith(
                ("embedding_", "expression_umap", "expression_pca", "state_")
            )
        ]
        if len(candidates) < 2:
            return None, "unavailable: declare a single_cell_embedding asset"
        frame = metadata[["cell_id", *candidates]].copy()
        source = "numeric expression/state coordinates declared in single-cell metadata"
    if "cell_id" not in frame:
        raise ValueError("The single-cell embedding requires a cell_id column.")
    if frame["cell_id"].astype(str).duplicated().any():
        raise ValueError("The single-cell embedding contains duplicate cell IDs.")
    numeric = [
        column
        for column in frame.columns
        if column != "cell_id"
        and pd.to_numeric(frame[column], errors="coerce").notna().mean() >= 0.95
    ]
    if len(numeric) < 2:
        raise ValueError(
            "The single-cell embedding requires at least two numeric coordinates."
        )
    output = frame[["cell_id", *numeric]].copy()
    output["cell_id"] = output["cell_id"].astype(str)
    output[numeric] = output[numeric].apply(pd.to_numeric, errors="coerce")
    output = output.dropna(subset=numeric)
    output = output.rename(columns={column: f"state_{i + 1}" for i, column in enumerate(numeric)})
    return output, f"{source}; graph dimensions={len(numeric)}"


def build_single_cell_programs(
    metadata: pd.DataFrame,
    counts: pd.DataFrame,
    pseudobulk: pd.DataFrame,
    embedding: pd.DataFrame | None,
    *,
    embedding_source: str,
    min_cell_event_reads: float = 1.0,
    min_cell_fraction: float = 0.08,
    min_state_cells: int = 40,
    state_neighbors: int = 15,
    state_permutations: int = 100,
    state_fdr: float = 0.10,
    min_independent_donors: int = 15,
    min_unit_reads: float = 10.0,
    min_unit_fraction: float = 0.70,
    min_psi_variance: float = 0.0025,
    min_module_size: int = 5,
    max_events: int = 0,
    bootstrap_iterations: int = 20,
    soft_power: int = 6,
    random_seed: int = 17,
) -> SingleCellProgramResult:
    """Fit exploratory cell-state and replicate-aware donor-level AS discovery programs."""

    state = _state_programs(
        metadata,
        counts,
        embedding,
        embedding_source=embedding_source,
        min_cell_event_reads=min_cell_event_reads,
        min_cell_fraction=min_cell_fraction,
        min_state_cells=min_state_cells,
        state_neighbors=state_neighbors,
        state_permutations=state_permutations,
        state_fdr=state_fdr,
        min_module_size=min_module_size,
        max_events=max_events,
        soft_power=soft_power,
        random_seed=random_seed,
    )
    replicate = _replicate_programs(
        pseudobulk,
        min_independent_donors=min_independent_donors,
        min_unit_reads=min_unit_reads,
        min_unit_fraction=min_unit_fraction,
        min_psi_variance=min_psi_variance,
        min_module_size=min_module_size,
        max_events=max_events,
        bootstrap_iterations=bootstrap_iterations,
        soft_power=soft_power,
        random_seed=random_seed,
    )
    donor_scores = _concat_tables(state.donor_scores, replicate.donor_scores)
    program_summary = _program_condition_tests(
        _concat_tables(state.program_summary, replicate.program_summary),
        donor_scores,
    )
    return SingleCellProgramResult(
        event_catalog=_concat_tables(state.event_catalog, replicate.event_catalog),
        program_summary=program_summary,
        cell_scores=state.cell_scores,
        donor_scores=donor_scores,
        network_edges=_concat_tables(state.network_edges, replicate.network_edges),
        diagnostics=_concat_tables(state.diagnostics, replicate.diagnostics),
    )


def _state_programs(
    metadata: pd.DataFrame,
    counts: pd.DataFrame,
    embedding: pd.DataFrame | None,
    *,
    embedding_source: str,
    min_cell_event_reads: float,
    min_cell_fraction: float,
    min_state_cells: int,
    state_neighbors: int,
    state_permutations: int,
    state_fdr: float,
    min_module_size: int,
    max_events: int,
    soft_power: int,
    random_seed: int,
) -> SingleCellProgramResult:
    if embedding is None:
        return _empty_program_result(
            "state_association",
            "not_run",
            embedding_source,
            "Declare an expression-derived cell embedding to test cell-state association.",
        )
    cells = (
        metadata.merge(embedding, on="cell_id", how="inner", validate="one_to_one")
        .dropna(subset=["state_1", "state_2"])
        .reset_index(drop=True)
    )
    if len(cells) < min_state_cells:
        return _empty_program_result(
            "state_association",
            "not_run",
            embedding_source,
            f"Only {len(cells)} embedded cells; at least {min_state_cells} are required.",
        )
    cell_order = pd.Index(cells["cell_id"].astype(str), name="cell_id")
    frame = counts[counts["cell_id"].isin(cell_order)].copy()
    frame["measured"] = frame["total"].ge(min_cell_event_reads)
    psi = frame.pivot_table(
        index="event_id", columns="cell_id", values="psi_posterior", aggfunc="first"
    ).reindex(columns=cell_order)
    total = frame.pivot_table(
        index="event_id", columns="cell_id", values="total", aggfunc="sum"
    ).reindex(index=psi.index, columns=cell_order)
    observed = psi.where(total.ge(min_cell_event_reads))
    measured_fraction = observed.notna().mean(axis=1)
    variance = observed.var(axis=1, ddof=1)
    donors = cells.set_index("cell_id")["donor_id"]
    donor_count = observed.notna().T.groupby(donors, sort=False).any().sum(axis=0)
    event_metadata = frame.drop_duplicates("event_id").set_index("event_id")
    event_qc = pd.DataFrame(
        {
            "event_id": observed.index.astype(str),
            "gene_id": event_metadata.reindex(observed.index)["gene_id"]
            .fillna("")
            .to_numpy(),
            "gene_symbol": event_metadata.reindex(observed.index)["gene_symbol"]
            .fillna("")
            .to_numpy(),
            "event_type": event_metadata.reindex(observed.index)["event_type"]
            .fillna("junction_choice")
            .to_numpy(),
            "measured_cells": observed.notna().sum(axis=1).astype(int).to_numpy(),
            "measured_cell_fraction": measured_fraction.to_numpy(float),
            "measured_donors": donor_count.astype(int).to_numpy(),
            "psi_variance": variance.to_numpy(float),
            "median_reads": total.where(total.ge(min_cell_event_reads))
            .median(axis=1)
            .to_numpy(float),
        }
    )
    event_qc["eligible"] = (
        event_qc["measured_cells"].ge(min_state_cells)
        & event_qc["measured_cell_fraction"].ge(min_cell_fraction)
        & event_qc["measured_donors"].ge(4)
        & event_qc["psi_variance"].gt(0)
    )
    event_qc["selection_score"] = (
        event_qc["measured_cell_fraction"]
        * event_qc["psi_variance"].fillna(0.0)
        * event_qc["median_reads"]
        .fillna(0.0)
        .div(event_qc["median_reads"].fillna(0.0) + 5.0)
    )
    selected_ids = (
        event_qc.loc[event_qc["eligible"]]
        .sort_values(["selection_score", "event_id"], ascending=[False, True])[
            "event_id"
        ]
        .head(max_events if max_events > 0 else len(counts))
    )
    if len(selected_ids) < min_module_size:
        return _empty_program_result(
            "state_association",
            "not_run",
            embedding_source,
            f"Only {len(selected_ids)} events passed state-association coverage gates.",
        )

    coordinates = cells[[c for c in cells if c.startswith("state_") and c[6:].isdigit()]].to_numpy(float)
    graph = _knn_graph(coordinates, state_neighbors)
    permutation_strata = (
        cells["donor_id"].astype(str) + "::" + cells["cell_type"].astype(str)
    ).to_numpy()
    rng = np.random.default_rng(random_seed)
    smoothed_rows: list[np.ndarray] = []
    statistic_rows: list[dict[str, Any]] = []
    for event_index, event_id in enumerate(selected_ids):
        values = observed.loc[event_id].to_numpy(float)
        smooth = _neighbor_average(values, graph)
        smoothed_rows.append(smooth)
        score, null = moran_permutation_scores(values, graph, permutation_strata, rng, state_permutations)
        finite = null[np.isfinite(null)]
        pvalue = (
            float((1 + np.sum(finite >= score)) / (1 + len(finite)))
            if np.isfinite(score) and len(finite)
            else float("nan")
        )
        qc = event_qc.set_index("event_id").loc[event_id]
        statistic_rows.append(
            {
                "event_id": event_id,
                "gene_id": qc["gene_id"],
                "gene_symbol": qc["gene_symbol"],
                "event_type": qc["event_type"],
                "measured_cells": int(qc["measured_cells"]),
                "measured_cell_fraction": float(qc["measured_cell_fraction"]),
                "measured_donors": int(qc["measured_donors"]),
                "median_reads": float(qc["median_reads"]),
                "state_autocorrelation": score,
                "pvalue": pvalue,
                "permutations": state_permutations,
            }
        )
    statistics = pd.DataFrame(statistic_rows)
    statistics["fdr"] = bh_adjust(statistics["pvalue"])
    statistics["state_associated"] = statistics["fdr"].le(state_fdr) & statistics[
        "state_autocorrelation"
    ].gt(0)
    module_events = statistics.loc[statistics["state_associated"], "event_id"]
    smoothed = pd.DataFrame(
        np.vstack(smoothed_rows), index=selected_ids, columns=cell_order
    )
    smoothed = smoothed.reindex(module_events).dropna(thresh=min_state_cells, axis=0)
    if len(smoothed) >= min_module_size:
        imputed = smoothed.T.fillna(smoothed.median(axis=1)).T
        labels, tom, engine = _fit_tom_labels(
            imputed,
            reliability=np.ones(len(imputed), dtype=float),
            min_module_size=min_module_size,
            soft_power=soft_power,
        )
        label_series = pd.Series(labels, index=imputed.index)
    else:
        imputed = smoothed
        tom = np.empty((0, 0), dtype=float)
        engine = "not_run"
        label_series = pd.Series("grey", index=statistics["event_id"], dtype=str)
    statistics["module"] = statistics["event_id"].map(label_series).fillna("grey")
    statistics["program_track"] = "cell_state_exploratory"
    statistics["cell_type_scope"] = "all embedded cells"
    statistics["embedding_source"] = embedding_source

    cell_scores, donor_scores, summary, memberships = _score_state_modules(
        imputed,
        label_series,
        cells,
        statistics,
    )
    statistics = statistics.merge(memberships, on="event_id", how="left")
    statistics["row_id"] = statistics["program_track"] + "::" + statistics["event_id"]
    edges = _module_edges(
        imputed.index,
        label_series,
        tom,
        track="cell_state_exploratory",
        cell_type="all embedded cells",
    )
    diagnostics = pd.DataFrame(
        [
            {
                "row_id": "cell_state_exploratory",
                "program_track": "cell_state_exploratory",
                "status": "complete",
                "embedding_source": embedding_source,
                "cells": len(cells),
                "independent_donors": int(cells["donor_id"].nunique()),
                "events_tested": len(statistics),
                "events_fdr_significant": int(statistics["state_associated"].sum()),
                "modules": int(label_series[label_series.ne("grey")].nunique()),
                "network_engine": engine,
                "claim_boundary": "Donor-and-cell-type-stratified permutation supports cell-state association; modules remain exploratory and non-causal.",
            }
        ]
    )
    return SingleCellProgramResult(
        event_catalog=statistics,
        program_summary=summary,
        cell_scores=cell_scores,
        donor_scores=donor_scores,
        network_edges=edges,
        diagnostics=diagnostics,
    )


def _replicate_programs(
    pseudobulk: pd.DataFrame,
    *,
    min_independent_donors: int,
    min_unit_reads: float,
    min_unit_fraction: float,
    min_psi_variance: float,
    min_module_size: int,
    max_events: int,
    bootstrap_iterations: int,
    soft_power: int,
    random_seed: int,
) -> SingleCellProgramResult:
    event_tables: list[pd.DataFrame] = []
    summaries: list[pd.DataFrame] = []
    donor_scores: list[pd.DataFrame] = []
    edges: list[pd.DataFrame] = []
    diagnostics: list[dict[str, Any]] = []
    for cell_type, frame in pseudobulk.groupby("cell_type", sort=True):
        frame = observation_units(frame)
        donor_metadata = frame.drop_duplicates("observation_id").set_index("observation_id")
        donors = donor_metadata.index.astype(str)
        if donor_metadata["donor_id"].nunique() < min_independent_donors:
            diagnostics.append(
                {
                    "row_id": f"replicate::{cell_type}",
                    "program_track": "donor_replicate_discovery",
                    "cell_type": cell_type,
                    "status": "not_run",
                    "independent_donors": donor_metadata["donor_id"].nunique(),
                    "required_donors": min_independent_donors,
                    "message": "Insufficient independent donors for a replicate program network.",
                    "claim_boundary": "Cells are never substituted for missing biological replicates.",
                }
            )
            continue
        psi = frame.pivot_table(
            index="event_id",
            columns="observation_id",
            values="psi_posterior",
            aggfunc="first",
        ).reindex(columns=donors)
        total = frame.pivot_table(
            index="event_id", columns="observation_id", values="total", aggfunc="first"
        ).reindex(index=psi.index, columns=donors)
        se = frame.pivot_table(
            index="event_id",
            columns="observation_id",
            values="measurement_se",
            aggfunc="first",
        ).reindex(index=psi.index, columns=donors)
        observed = psi.where(total.ge(min_unit_reads))
        valid = observed.notna().mean(axis=1)
        variance = observed.var(axis=1, ddof=1)
        median_reads = total.where(total.ge(min_unit_reads)).median(axis=1)
        median_se = se.where(total.ge(min_unit_reads)).median(axis=1)
        reliability = (
            valid
            * median_reads.div(median_reads + 20.0)
            * (1.0 - 2.0 * median_se).clip(0.0, 1.0)
        )
        selected_ids = (
            pd.DataFrame(
                {
                    "valid": valid,
                    "variance": variance,
                    "reliability": reliability,
                }
            )
            .query("valid >= @min_unit_fraction and variance >= @min_psi_variance")
            .assign(score=lambda table: table["variance"] * table["reliability"])
            .sort_values(["score"], ascending=False)
            .head(max_events if max_events > 0 else len(psi))
            .index
        )
        if len(selected_ids) < min_module_size:
            diagnostics.append(
                {
                    "row_id": f"replicate::{cell_type}",
                    "program_track": "donor_replicate_discovery",
                    "cell_type": cell_type,
                    "status": "not_run",
                    "independent_donors": donor_metadata["donor_id"].nunique(),
                    "events_passing_gate": len(selected_ids),
                    "message": "Too few reliable variable events for program discovery.",
                    "claim_boundary": "No module is forced when event support is inadequate.",
                }
            )
            continue
        matrix = observed.reindex(selected_ids)
        matrix = matrix.T.fillna(matrix.median(axis=1)).T
        event_reliability = (
            reliability.reindex(selected_ids).fillna(0.0).to_numpy(float)
        )
        labels, tom, engine = _fit_tom_labels(
            matrix,
            reliability=event_reliability,
            min_module_size=min_module_size,
            soft_power=soft_power,
            official=True,
        )
        label_series = pd.Series(labels, index=matrix.index)
        conditions = donor_metadata.reindex(donors)["condition"].astype(str)
        confidence, bootstrap_ari = _bootstrap_membership_confidence(
            matrix,
            label_series,
            conditions,
            event_reliability,
            min_module_size=min_module_size,
            soft_power=soft_power,
            iterations=bootstrap_iterations,
            seed=random_seed,
            donor_metadata=donor_metadata,
            official=True,
        )
        event_metadata = frame.drop_duplicates("event_id").set_index("event_id")
        assignments = pd.DataFrame(
            {
                "event_id": matrix.index.astype(str),
                "gene_id": event_metadata.reindex(matrix.index)["gene_id"]
                .fillna("")
                .to_numpy(),
                "gene_symbol": event_metadata.reindex(matrix.index)["gene_symbol"]
                .fillna("")
                .to_numpy(),
                "event_type": event_metadata.reindex(matrix.index)["event_type"]
                .fillna("junction_choice")
                .to_numpy(),
                "module": label_series.to_numpy(),
                "network_reliability": event_reliability,
                "bootstrap_membership_confidence": confidence,
                "program_track": "donor_replicate_discovery",
                "cell_type_scope": str(cell_type),
                "state_autocorrelation": np.nan,
                "pvalue": np.nan,
                "fdr": np.nan,
                "state_associated": False,
            }
        )
        score_table, membership = _score_replicate_modules(
            matrix, label_series, donor_metadata.reindex(donors), str(cell_type)
        )
        assignments = assignments.merge(membership, on="event_id", how="left")
        assignments["row_id"] = (
            assignments["program_track"]
            + "::"
            + assignments["cell_type_scope"]
            + "::"
            + assignments["event_id"]
        )
        summary = _summarize_programs(assignments, score_table)
        event_tables.append(assignments)
        summaries.append(summary)
        donor_scores.append(score_table)
        edges.append(
            _module_edges(
                matrix.index,
                label_series,
                tom,
                track="donor_replicate_discovery",
                cell_type=str(cell_type),
            )
        )
        diagnostics.append(
            {
                "row_id": f"replicate::{cell_type}",
                "program_track": "donor_replicate_discovery",
                "cell_type": cell_type,
                "status": "complete",
                "independent_donors": donor_metadata["donor_id"].nunique(),
                "donor_condition_observations": len(donors),
                "events_selected": len(matrix),
                "modules": int(label_series[label_series.ne("grey")].nunique()),
                "grey_fraction": float(label_series.eq("grey").mean()),
                "mean_bootstrap_ari": float(np.nanmean(bootstrap_ari))
                if len(bootstrap_ari)
                else np.nan,
                "min_bootstrap_ari": float(np.nanmin(bootstrap_ari))
                if np.isfinite(bootstrap_ari).any()
                else np.nan,
                "network_engine": engine,
                "claim_boundary": "The independent donor is the replicate; modules are confidence-annotated rather than forced.",
            }
        )
    return SingleCellProgramResult(
        event_catalog=_concat_many(event_tables),
        program_summary=_concat_many(summaries),
        cell_scores=pd.DataFrame(),
        donor_scores=_concat_many(donor_scores),
        network_edges=_concat_many(edges),
        diagnostics=pd.DataFrame(diagnostics),
    )


def _fit_tom_labels(
    matrix: pd.DataFrame,
    *,
    reliability: np.ndarray,
    min_module_size: int,
    soft_power: int,
    official: bool = False,
) -> tuple[np.ndarray, np.ndarray, str]:
    if len(matrix) < 2:
        return np.repeat("grey", len(matrix)), np.eye(len(matrix)), "not_run"
    if len(matrix) ** 2 * 64 > 2 * 1024 ** 3:
        raise ValueError("This exploratory dense SC TOM requires over 2 GiB of working matrices. Use an explicit exploratory event cap or the blockwise bulk engine on donor pseudobulk. No events were silently dropped.")
    if official:
        from katana_cosplice.r_wgcna import exact_wgcna_labels
        labels, tom = exact_wgcna_labels(matrix, soft_power=soft_power, min_module_size=min_module_size,
                                         merge_cut_height=.25, deep_split=2, pam_respects_dendro=False,
                                         event_reliability=pd.Series(reliability, index=matrix.index),
                                         reliability_power=1.0, return_tom=True)
        return labels, tom, "official R WGCNA signed bicor/TOM + Dynamic Tree Cut; shared bulk production adapter"
    settings = NetworkConfig(
        method="bicor",
        signed=True,
        soft_power=soft_power,
        min_module_size=min_module_size,
        dynamic_deep_split=2,
        dynamic_cut_height=0.995,
        module_refinement=False,
    )
    correlation = correlation_matrix(matrix, settings.method)
    adjacency = signed_soft_adjacency(correlation, settings)
    if len(reliability) == len(matrix):
        adjacency *= np.sqrt(np.outer(reliability, reliability)).astype(np.float32)
    tom = topological_overlap_matrix(adjacency)
    distance = np.clip(1.0 - tom, 0.0, 1.0)
    np.fill_diagonal(distance, 0.0)
    hierarchy = linkage(squareform(distance, checks=False), method="average")
    try:
        raw = dynamic_tree_cut_labels(
            hierarchy,
            distance,
            min_module_size=min_module_size,
            deep_split=2,
            cut_height=0.995,
            pam_respects_dendro=False,
        )
        engine = "signed bicor soft adjacency + TOM + Dynamic Tree Cut"
    except (ValueError, RuntimeError, IndexError):
        raw = fcluster(hierarchy, t=0.80, criterion="distance")
        engine = "signed bicor soft adjacency + TOM + fixed-height fallback"
    return _stable_module_labels(raw, matrix.index, min_module_size), tom, engine


def _stable_module_labels(
    raw: np.ndarray, event_ids: pd.Index, min_module_size: int
) -> np.ndarray:
    series = pd.Series(np.asarray(raw), index=event_ids)
    groups = []
    for label, members in series.groupby(series, sort=False):
        if int(label) == 0 or len(members) < min_module_size:
            continue
        groups.append((str(members.index.min()), int(label), len(members)))
    groups.sort(key=lambda item: (-item[2], item[0]))
    mapping = {
        raw_label: f"M{index:02d}" for index, (_, raw_label, _) in enumerate(groups, 1)
    }
    return np.asarray([mapping.get(int(label), "grey") for label in raw], dtype=object)


def _knn_graph(coordinates: np.ndarray, neighbors: int) -> csr_matrix:
    count = len(coordinates)
    model = NearestNeighbors(
        n_neighbors=min(count, max(2, neighbors + 1)), metric="euclidean"
    ).fit(coordinates)
    distances, indices = model.kneighbors(coordinates)
    positive = distances[:, 1:][distances[:, 1:] > 0]
    scale = float(np.median(positive)) if len(positive) else 1.0
    rows: list[int] = []
    columns: list[int] = []
    weights: list[float] = []
    for source in range(count):
        for distance, target in zip(
            distances[source, 1:], indices[source, 1:], strict=True
        ):
            rows.append(source)
            columns.append(int(target))
            weights.append(float(np.exp(-(float(distance) ** 2) / (2.0 * scale**2))))
    graph = csr_matrix((weights, (rows, columns)), shape=(count, count))
    graph = graph.maximum(graph.T)
    graph.setdiag(0.0)
    graph.eliminate_zeros()
    return graph


def _neighbor_average(values: np.ndarray, graph: csr_matrix) -> np.ndarray:
    observed = np.isfinite(values)
    augmented = graph + eye(graph.shape[0], format="csr")
    numerator = augmented @ np.where(observed, values, 0.0)
    denominator = augmented @ observed.astype(float)
    return np.divide(
        numerator,
        denominator,
        out=np.full(len(values), np.nan, dtype=float),
        where=denominator > 0,
    )


def _morans_i(values: np.ndarray, graph: csr_matrix) -> float:
    observed = np.flatnonzero(np.isfinite(values))
    if len(observed) < 4:
        return float("nan")
    subgraph = graph[observed][:, observed]
    weight_sum = float(subgraph.sum())
    centered = values[observed] - np.nanmean(values[observed])
    denominator = float(centered @ centered)
    if weight_sum <= 0 or denominator <= 1e-12:
        return float("nan")
    numerator = float(centered @ (subgraph @ centered))
    return float(len(observed) / weight_sum * numerator / denominator)


def moran_permutation_scores(values, graph, strata, rng, iterations, batch_size=128):
    """Same stratified null as the scalar kernel, with a fixed missingness mask.

    Cache the observed subgraph and multiply a batch of permutations at once.
    The RNG sequence matches _shuffle_within_groups; no cell or donor is added.
    """
    score = _morans_i(values, graph)
    if not np.isfinite(score):
        return score, np.full(iterations, np.nan)
    observed = np.flatnonzero(np.isfinite(values))
    subgraph = graph[observed][:, observed]
    centered = np.asarray(values)[observed] - np.mean(np.asarray(values)[observed])
    observed_strata = np.asarray(strata)[observed]
    groups = [np.flatnonzero(observed_strata == group) for group in pd.unique(strata)]
    groups = [group for group in groups if len(group) > 1]
    scale = len(observed) / float(subgraph.sum()) / float(centered @ centered)
    null = np.empty(iterations)
    for start in range(0, iterations, batch_size):
        count = min(batch_size, iterations-start)
        permutations = np.tile(centered[:, None], (1, count))
        for column in range(count):
            for group in groups:
                permutations[group, column] = centered[rng.permutation(group)]
        null[start:start+count] = np.sum(permutations * (subgraph @ permutations), axis=0) * scale
    return score, null


def _shuffle_within_groups(
    values: np.ndarray, groups: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    shuffled = values.copy()
    for group in pd.unique(groups):
        indices = np.flatnonzero((groups == group) & np.isfinite(values))
        if len(indices) > 1:
            shuffled[indices] = values[rng.permutation(indices)]
    return shuffled


def _score_state_modules(
    matrix: pd.DataFrame,
    labels: pd.Series,
    cells: pd.DataFrame,
    statistics: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cell_rows: list[pd.DataFrame] = []
    membership_rows: list[dict[str, Any]] = []
    for module in sorted(set(labels).difference({"grey"})):
        event_ids = labels.index[labels.eq(module)]
        values = matrix.loc[event_ids]
        score = _first_component_score(values)
        for event_id in event_ids:
            membership_rows.append(
                {
                    "event_id": event_id,
                    "module_kme": _safe_correlation(values.loc[event_id], score),
                }
            )
        table = cells[
            [
                "cell_id",
                "sample_id",
                "donor_id",
                "cell_type",
                "condition",
                "batch",
                "state_1",
                "state_2",
            ]
        ].copy()
        table["module"] = module
        table["program_score"] = score
        table["program_track"] = "cell_state_exploratory"
        table["cell_type_scope"] = "all embedded cells"
        table["row_id"] = table["cell_id"] + "::" + module
        cell_rows.append(table)
    cell_scores = _concat_many(cell_rows)
    donor_scores = (
        cell_scores.groupby(
            [
                "program_track",
                "cell_type_scope",
                "module",
                "donor_id",
                "condition",
                "batch",
            ],
            as_index=False,
        ).agg(program_score=("program_score", "mean"), cells=("cell_id", "nunique"))
        if not cell_scores.empty
        else pd.DataFrame()
    )
    if not donor_scores.empty:
        donor_scores["row_id"] = (
            donor_scores["program_track"]
            + "::"
            + donor_scores["module"]
            + "::"
            + donor_scores["donor_id"]
            + "::" + donor_scores["condition"].astype(str)
            + "::" + donor_scores["batch"].astype(str)
        )
    membership = pd.DataFrame(membership_rows, columns=["event_id", "module_kme"])
    catalog = statistics.copy()
    if not membership.empty:
        catalog = catalog.merge(membership, on="event_id", how="left")
    summary = _summarize_programs(catalog, donor_scores)
    return cell_scores, donor_scores, summary, membership


def _score_replicate_modules(
    matrix: pd.DataFrame,
    labels: pd.Series,
    metadata: pd.DataFrame,
    cell_type: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    score_rows: list[pd.DataFrame] = []
    membership_rows: list[dict[str, Any]] = []
    for module in sorted(set(labels).difference({"grey"})):
        event_ids = labels.index[labels.eq(module)]
        values = matrix.loc[event_ids]
        score = _first_component_score(values)
        for event_id in event_ids:
            membership_rows.append(
                {
                    "event_id": event_id,
                    "module_kme": _safe_correlation(values.loc[event_id], score),
                }
            )
        table = metadata.reset_index()[["donor_id", "condition", "batch"]].copy()
        table["module"] = module
        table["program_score"] = score
        table["program_track"] = "donor_replicate_discovery"
        table["cell_type_scope"] = cell_type
        table["cells"] = np.nan
        table["row_id"] = (
            table["program_track"]
            + "::"
            + cell_type
            + "::"
            + module
            + "::"
            + table["donor_id"].astype(str)
            + "::" + table["condition"].astype(str)
        )
        score_rows.append(table)
    return _concat_many(score_rows), pd.DataFrame(
        membership_rows, columns=["event_id", "module_kme"]
    )


def _summarize_programs(
    assignments: pd.DataFrame, scores: pd.DataFrame
) -> pd.DataFrame:
    if assignments.empty or "module" not in assignments:
        return pd.DataFrame()
    retained = assignments[assignments["module"].ne("grey")].copy()
    rows: list[dict[str, Any]] = []
    for keys, group in retained.groupby(
        ["program_track", "cell_type_scope", "module"], sort=True
    ):
        track, cell_type, module = keys
        hub = group.sort_values(
            "module_kme", key=lambda values: values.abs(), ascending=False
        ).iloc[0]
        score_group = (
            scores[
                scores["program_track"].eq(track)
                & scores["cell_type_scope"].eq(cell_type)
                & scores["module"].eq(module)
            ]
            if not scores.empty
            else pd.DataFrame()
        )
        rows.append(
            {
                "row_id": f"{track}::{cell_type}::{module}",
                "program_track": track,
                "cell_type_scope": cell_type,
                "module": module,
                "events": len(group),
                "genes": int(group["gene_symbol"].replace("", np.nan).nunique()),
                "hub_event": hub["event_id"],
                "hub_gene": hub.get("gene_symbol", ""),
                "hub_kme": hub.get("module_kme", np.nan),
                "median_membership_confidence": pd.to_numeric(
                    group.get(
                        "bootstrap_membership_confidence", pd.Series(dtype=float)
                    ),
                    errors="coerce",
                ).median(),
                "median_state_fdr": pd.to_numeric(
                    group.get("fdr", pd.Series(dtype=float)), errors="coerce"
                ).median(),
                "score_units": int(score_group["donor_id"].nunique())
                if not score_group.empty
                else 0,
                "scientific_role": "exploratory donor-level program"
                if track == "donor_replicate_discovery"
                else "exploratory cell-state-associated program",
            }
        )
    return pd.DataFrame(rows)


def _program_condition_tests(
    summary: pd.DataFrame, donor_scores: pd.DataFrame
) -> pd.DataFrame:
    """Attach donor-level condition omnibus tests to discovered programs."""

    if summary.empty:
        return summary
    tests: list[dict[str, Any]] = []
    keys = ["program_track", "cell_type_scope", "module"]
    for key_values, group in donor_scores.groupby(keys, sort=True):
        track, cell_type, module = key_values
        values = group.dropna(subset=["program_score", "condition"]).copy()
        condition_counts = values.groupby("condition")["donor_id"].nunique()
        row: dict[str, Any] = {
            "program_track": track,
            "cell_type_scope": cell_type,
            "module": module,
            "condition_test_status": "insufficient_conditions",
            "condition_omnibus_f": np.nan,
            "condition_omnibus_pvalue": np.nan,
            "condition_df": 0,
            "residual_df": 0,
            "highest_score_condition": "",
            "lowest_score_condition": "",
            "max_condition_mean_difference": np.nan,
            "condition_test_batch_adjusted": False,
        }
        if len(condition_counts) < 2 or int(condition_counts.min()) < 2:
            tests.append(row)
            continue
        condition_means = values.groupby("condition")["program_score"].mean()
        row["highest_score_condition"] = str(condition_means.idxmax())
        row["lowest_score_condition"] = str(condition_means.idxmin())
        row["max_condition_mean_difference"] = float(
            condition_means.max() - condition_means.min()
        )
        condition_design = pd.get_dummies(
            values["condition"].astype(str), drop_first=True, dtype=float
        ).to_numpy(float)
        covariates, paired, batch_added = donor_covariates(values)
        null_design = np.column_stack([np.ones(len(values)), covariates]) if covariates is not None else np.ones((len(values), 1))
        full_design = np.column_stack([null_design, condition_design])
        row["condition_test_batch_adjusted"] = batch_added
        row["condition_test_paired_donor_adjusted"] = paired
        if (
            np.linalg.matrix_rank(null_design) < null_design.shape[1]
            or np.linalg.matrix_rank(full_design) < full_design.shape[1]
        ):
            row["condition_test_status"] = "condition_batch_confounded"
            tests.append(row)
            continue
        outcome = values["program_score"].to_numpy(float)
        null_residual = (
            outcome - null_design @ np.linalg.lstsq(null_design, outcome, rcond=None)[0]
        )
        full_residual = (
            outcome - full_design @ np.linalg.lstsq(full_design, outcome, rcond=None)[0]
        )
        null_rss = float(null_residual @ null_residual)
        full_rss = float(full_residual @ full_residual)
        numerator_df = full_design.shape[1] - null_design.shape[1]
        residual_df = len(values) - full_design.shape[1]
        if numerator_df <= 0 or residual_df <= 0 or full_rss <= 1e-12:
            row["condition_test_status"] = "insufficient_residual_df"
            tests.append(row)
            continue
        statistic = max(
            0.0,
            ((null_rss - full_rss) / numerator_df) / (full_rss / residual_df),
        )
        row.update(
            {
                "condition_test_status": "fit",
                "condition_omnibus_f": statistic,
                "condition_omnibus_pvalue": float(
                    f_distribution.sf(statistic, numerator_df, residual_df)
                ),
                "condition_df": numerator_df,
                "residual_df": residual_df,
            }
        )
        tests.append(row)
    test_frame = pd.DataFrame(tests)
    if test_frame.empty:
        return summary
    test_frame["condition_omnibus_fdr"] = bh_adjust(test_frame["condition_omnibus_pvalue"])
    test_frame["condition_test_family"] = "all module-condition omnibus tests in this run; exploratory post-discovery associations"
    return summary.merge(test_frame, on=keys, how="left", validate="one_to_one")


def _bootstrap_membership_confidence(
    matrix: pd.DataFrame,
    reference: pd.Series,
    conditions: pd.Series,
    reliability: np.ndarray,
    *,
    min_module_size: int,
    soft_power: int,
    iterations: int,
    seed: int,
    donor_metadata: pd.DataFrame | None = None,
    official: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    counts = np.zeros((len(matrix), len(matrix)), dtype=float)
    aris: list[float] = []
    rng = np.random.default_rng(seed)
    for iteration in range(iterations):
        indices = bootstrap_blocks(donor_metadata, rng) if donor_metadata is not None else np.concatenate(
            [
                rng.choice(group, size=len(group), replace=True)
                for condition in pd.unique(conditions)
                for group in [np.flatnonzero(conditions.to_numpy() == condition)]
            ]
        )
        bootstrap = matrix.iloc[:, indices].copy()
        bootstrap.columns = [
            f"boot_{iteration}_{index}" for index in range(len(indices))
        ]
        labels, _tom, _engine = _fit_tom_labels(
            bootstrap,
            reliability=reliability,
            min_module_size=min_module_size,
            soft_power=soft_power,
            official=official,
        )
        aris.append(adjusted_rand_score(reference.to_numpy(), labels))
        same = labels[:, None] == labels[None, :]
        retained = (labels[:, None] != "grey") & (labels[None, :] != "grey")
        counts += same & retained
    coassignment = counts / max(1, iterations)
    confidence = np.zeros(len(reference), dtype=float)
    for index, module in enumerate(reference):
        members = np.flatnonzero(reference.to_numpy() == module)
        members = members[members != index]
        confidence[index] = (
            float(coassignment[index, members].mean())
            if module != "grey" and len(members)
            else 0.0
        )
    return confidence, np.asarray(aris, dtype=float)


def _module_edges(
    event_ids: pd.Index,
    labels: pd.Series,
    tom: np.ndarray,
    *,
    track: str,
    cell_type: str,
    max_edges: int = 2_000,
) -> pd.DataFrame:
    if len(tom) != len(event_ids):
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    for left in range(len(event_ids)):
        for right in range(left + 1, len(event_ids)):
            module = labels.iloc[left]
            if module == "grey" or module != labels.iloc[right]:
                continue
            rows.append(
                {
                    "source_event": str(event_ids[left]),
                    "target_event": str(event_ids[right]),
                    "module": module,
                    "tom_similarity": float(tom[left, right]),
                    "program_track": track,
                    "cell_type_scope": cell_type,
                }
            )
    if not rows:
        return pd.DataFrame()
    output = (
        pd.DataFrame(rows)
        .sort_values("tom_similarity", ascending=False)
        .head(max_edges)
    )
    output["row_id"] = output["source_event"] + "::" + output["target_event"]
    return output


def _first_component_score(matrix: pd.DataFrame) -> np.ndarray:
    values = matrix.to_numpy(float)
    centered = values - values.mean(axis=1, keepdims=True)
    scale = centered.std(axis=1, ddof=1, keepdims=True)
    standardized = centered / np.where(scale < 1e-9, 1.0, scale)
    _u, _s, right = np.linalg.svd(standardized, full_matrices=False)
    score = right[0]
    if np.corrcoef(score, standardized.mean(axis=0))[0, 1] < 0:
        score = -score
    return score


def _safe_correlation(values: pd.Series, score: np.ndarray) -> float:
    array = values.to_numpy(float)
    if len(array) < 3 or np.std(array) < 1e-12 or np.std(score) < 1e-12:
        return float("nan")
    return float(np.corrcoef(array, score)[0, 1])


def _empty_program_result(
    track: str, status: str, source: str, message: str
) -> SingleCellProgramResult:
    diagnostics = pd.DataFrame(
        [
            {
                "row_id": track,
                "program_track": track,
                "status": status,
                "embedding_source": source,
                "message": message,
            }
        ]
    )
    return SingleCellProgramResult(
        event_catalog=pd.DataFrame(),
        program_summary=pd.DataFrame(),
        cell_scores=pd.DataFrame(),
        donor_scores=pd.DataFrame(),
        network_edges=pd.DataFrame(),
        diagnostics=diagnostics,
    )


def _concat_tables(left: pd.DataFrame, right: pd.DataFrame) -> pd.DataFrame:
    return _concat_many([left, right])


def _concat_many(frames: list[pd.DataFrame]) -> pd.DataFrame:
    available = [frame for frame in frames if frame is not None and not frame.empty]
    return (
        pd.concat(available, ignore_index=True, sort=False)
        if available
        else pd.DataFrame()
    )
