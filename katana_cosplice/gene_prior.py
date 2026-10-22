from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from igraph import Graph
import leidenalg

from katana_cosplice.config import GenePriorConfig
from katana_cosplice.mantel_cosplicex import (
    canberra_distance_vectors,
    mantel_similarity,
    unsigned_adjacency,
)


@dataclass
class GeneMantelPrior:
    gene_adjacency: pd.DataFrame
    gene_modules: pd.Series
    event_genes: pd.Series
    gene_audit: pd.DataFrame
    event_audit: pd.DataFrame
    prior_mode: str = "continuous_adjacency"


def canonical_gene_id(value: object) -> str:
    return str(value).strip().split(".", 1)[0]


def event_gene_id(event_key: object) -> str:
    return canonical_gene_id(str(event_key).split("|", 1)[0])


def _read_table(path: Path, *, index_col: int | None = None) -> pd.DataFrame:
    suffixes = "".join(path.suffixes).lower()
    sep = "," if suffixes.endswith(".csv") or suffixes.endswith(".csv.gz") else "\t"
    return pd.read_csv(path, sep=sep, index_col=index_col, low_memory=False)


def _feature_gene_ids(
    feature_ids: pd.Index,
    settings: GenePriorConfig,
) -> pd.Series:
    features = pd.Series(feature_ids.astype(str), index=feature_ids.astype(str), dtype="object")
    if settings.feature_gene_map is not None:
        mapping = _read_table(settings.feature_gene_map)
        required = {settings.feature_id_column, settings.gene_id_column}
        missing = required.difference(mapping.columns)
        if missing:
            raise ValueError(
                f"gene_prior.feature_gene_map is missing columns: {sorted(missing)}"
            )
        lookup = (
            mapping.loc[:, [settings.feature_id_column, settings.gene_id_column]]
            .dropna()
            .drop_duplicates(settings.feature_id_column)
            .set_index(settings.feature_id_column)[settings.gene_id_column]
        )
        genes = features.map(lookup)
    elif settings.encoded_gene_rule == "exact_feature_id":
        genes = features
    else:
        genes = features.str.rsplit("_", n=1).str[0]
    return genes.map(canonical_gene_id)


def build_gene_mantel_prior(
    event_keys: pd.Index,
    sample_ids: pd.Index,
    settings: GenePriorConfig,
) -> GeneMantelPrior:
    if not settings.enabled or settings.exon_count_matrix is None:
        raise ValueError("gene_prior is not enabled or exon_count_matrix is missing")
    counts = _read_table(settings.exon_count_matrix, index_col=0)
    counts.index = counts.index.astype(str)
    counts.columns = counts.columns.astype(str)
    missing_samples = sorted(set(sample_ids.astype(str)).difference(counts.columns))
    if missing_samples:
        raise ValueError(
            "gene_prior exon count matrix is missing network samples: "
            f"{missing_samples[:10]}"
        )
    counts = counts.loc[:, sample_ids.astype(str)].apply(pd.to_numeric, errors="coerce")
    counts = counts.fillna(0.0)
    if (counts.to_numpy(dtype=float) < 0).any():
        raise ValueError("gene_prior exon count matrix cannot contain negative values")

    event_genes = pd.Series(
        [event_gene_id(value) for value in event_keys],
        index=event_keys.astype(str),
        name="gene_id",
        dtype="object",
    )
    target_genes = set(event_genes)
    feature_genes = _feature_gene_ids(counts.index, settings)
    matched_features = feature_genes.isin(target_genes)
    counts = counts.loc[matched_features]
    feature_genes = feature_genes.loc[matched_features]

    profiles: dict[str, pd.DataFrame] = {}
    audit_rows: list[dict[str, object]] = []
    for gene_id, feature_ids in feature_genes.groupby(feature_genes, sort=False).groups.items():
        profile = counts.loc[list(feature_ids)]
        variable_features = int(profile.std(axis=1, ddof=0).gt(0).sum())
        if len(profile) < settings.min_exons or variable_features < 2:
            continue
        median_total = float(profile.sum(axis=0).median())
        if median_total <= 0:
            continue
        profiles[str(gene_id)] = profile
        audit_rows.append(
            {
                "gene_id": str(gene_id),
                "n_features": int(len(profile)),
                "n_variable_features": variable_features,
                "median_feature_count_sum": median_total,
                "n_network_events": int(event_genes.eq(str(gene_id)).sum()),
            }
        )
    gene_audit = pd.DataFrame(audit_rows)
    if gene_audit.empty:
        raise ValueError(
            "gene_prior found no network genes with at least "
            f"{settings.min_exons} mapped variable exon features"
        )
    gene_audit = gene_audit.sort_values(
        ["n_network_events", "median_feature_count_sum", "gene_id"],
        ascending=[False, False, True],
        kind="mergesort",
    )
    if settings.max_genes is not None:
        gene_audit = gene_audit.head(settings.max_genes)
    retained_genes = gene_audit["gene_id"].astype(str).tolist()
    profiles = {gene: profiles[gene] for gene in retained_genes}
    if len(profiles) < 2:
        raise ValueError("gene_prior requires at least two eligible network genes")

    vectors = canberra_distance_vectors(profiles)
    correlation = mantel_similarity(vectors)
    adjacency = unsigned_adjacency(correlation, settings.soft_power)
    modules = _gene_program_scaffold(adjacency, settings)
    gene_audit["gene_prior_module"] = gene_audit["gene_id"].map(modules).fillna("grey")
    available = event_genes.isin(adjacency.index)
    event_modules = event_genes.map(modules).fillna("grey")
    event_audit = pd.DataFrame(
        {
            "event_key": event_genes.index,
            "gene_id": event_genes.to_numpy(),
            "gene_prior_available": available.to_numpy(),
            "gene_prior_module": event_modules.to_numpy(),
        }
    )
    return GeneMantelPrior(
        gene_adjacency=adjacency,
        gene_modules=modules,
        event_genes=event_genes,
        gene_audit=gene_audit.reset_index(drop=True),
        event_audit=event_audit,
        prior_mode=settings.prior_mode,
    )


def _gene_program_scaffold(
    adjacency: pd.DataFrame,
    settings: GenePriorConfig,
) -> pd.Series:
    genes = adjacency.index.astype(str)
    if len(genes) < settings.scaffold_min_module_size:
        return pd.Series("grey", index=genes, name="gene_prior_module")
    values = adjacency.to_numpy(dtype=float).copy()
    np.fill_diagonal(values, -np.inf)
    neighbors = min(settings.scaffold_knn_neighbors, len(genes) - 1)
    directed = np.zeros((len(genes), len(genes)), dtype=bool)
    for row in range(len(genes)):
        selected = np.argpartition(values[row], -neighbors)[-neighbors:]
        directed[row, selected] = True
    mutual = directed & directed.T
    rows, cols = np.where(np.triu(mutual, k=1))
    graph = Graph(n=len(genes), directed=False)
    if len(rows):
        graph.add_edges(list(zip(rows.tolist(), cols.tolist(), strict=True)))
        graph.es["weight"] = adjacency.to_numpy(dtype=float)[rows, cols].tolist()
        partition = leidenalg.find_partition(
            graph,
            leidenalg.RBConfigurationVertexPartition,
            weights="weight",
            resolution_parameter=settings.scaffold_resolution,
            seed=settings.random_seed,
        )
        raw = np.asarray(partition.membership, dtype=int)
    else:
        raw = np.arange(len(genes), dtype=int)
    counts = pd.Series(raw).value_counts()
    retained = [
        int(label)
        for label, count in counts.items()
        if count >= settings.scaffold_min_module_size
    ]
    retained.sort(key=lambda label: (-int(counts[label]), label))
    names = {label: f"GP{index:02d}" for index, label in enumerate(retained, start=1)}
    return pd.Series(
        [names.get(int(label), "grey") for label in raw],
        index=genes,
        name="gene_prior_module",
        dtype="object",
    )


def edge_gene_prior(
    prior: GeneMantelPrior,
    event_keys: pd.Index,
    rows: np.ndarray,
    cols: np.ndarray,
) -> np.ndarray:
    genes = prior.event_genes.reindex(event_keys.astype(str))
    lookup = {gene: index for index, gene in enumerate(prior.gene_adjacency.index.astype(str))}
    gene_indices = np.asarray([lookup.get(str(gene), -1) for gene in genes], dtype=int)
    left = gene_indices[rows]
    right = gene_indices[cols]
    values = np.zeros(len(rows), dtype=float)
    valid = (left >= 0) & (right >= 0)
    if valid.any():
        adjacency = prior.gene_adjacency.to_numpy(dtype=float)
        values[valid] = adjacency[left[valid], right[valid]]
        if prior.prior_mode == "scaffold_reinforcement":
            modules = prior.gene_modules.reindex(prior.gene_adjacency.index).fillna("grey")
            module_values = modules.to_numpy(dtype=object)
            same_program = (
                valid
                & (module_values[left.clip(min=0)] != "grey")
                & (module_values[left.clip(min=0)] == module_values[right.clip(min=0)])
            )
            values[same_program] = 1.0
    return np.clip(values, 0.0, 1.0)


def save_gene_mantel_prior(
    prior: GeneMantelPrior,
    matrix_path: Path,
    gene_audit_path: Path,
    event_audit_path: Path,
) -> None:
    matrix_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        matrix_path,
        gene_ids=prior.gene_adjacency.index.astype(str).to_numpy(dtype=str),
        adjacency=prior.gene_adjacency.to_numpy(dtype=np.float32),
        gene_modules=prior.gene_modules.reindex(
            prior.gene_adjacency.index
        ).fillna("grey").astype(str).to_numpy(dtype=str),
        prior_mode=np.asarray([prior.prior_mode], dtype=str),
    )
    prior.gene_audit.to_csv(gene_audit_path, sep="\t", index=False)
    prior.event_audit.to_csv(event_audit_path, sep="\t", index=False)


def load_gene_mantel_prior(
    matrix_path: Path,
    event_keys: pd.Index,
) -> GeneMantelPrior:
    payload = np.load(matrix_path, allow_pickle=False)
    gene_ids = pd.Index(payload["gene_ids"].astype(str), name="gene_id")
    adjacency = pd.DataFrame(payload["adjacency"], index=gene_ids, columns=gene_ids)
    gene_modules = pd.Series(
        payload["gene_modules"].astype(str),
        index=gene_ids,
        name="gene_prior_module",
        dtype="object",
    )
    prior_mode = str(payload["prior_mode"].astype(str)[0])
    event_genes = pd.Series(
        [event_gene_id(value) for value in event_keys],
        index=event_keys.astype(str),
        name="gene_id",
        dtype="object",
    )
    event_audit = pd.DataFrame(
        {
            "event_key": event_genes.index,
            "gene_id": event_genes.to_numpy(),
            "gene_prior_available": event_genes.isin(gene_ids).to_numpy(),
            "gene_prior_module": event_genes.map(gene_modules).fillna("grey").to_numpy(),
        }
    )
    return GeneMantelPrior(
        adjacency,
        gene_modules,
        event_genes,
        pd.DataFrame(),
        event_audit,
        prior_mode,
    )
