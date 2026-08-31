from __future__ import annotations

import subprocess
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from sklearn.metrics import adjusted_rand_score

from katana_cosplice.benchmark import _r_subprocess_env


Progress = Callable[[str], None] | None


@dataclass(frozen=True)
class MantelCoSplicExConfig:
    """Published CoSplicEx network settings from Iancu et al. (2015)."""

    soft_power: float = 6.0
    min_module_size: int = 20
    cut_height: float = 0.99995
    max_core_scatter: float = 0.995
    deep_split: int = 4
    bootstrap_iterations: int = 20
    bootstrap_fraction: float = 0.9
    random_seed: int = 73


@dataclass
class MantelCoSplicExResult:
    assignments: pd.DataFrame
    gene_audit: pd.DataFrame
    bootstrap_stability: pd.DataFrame
    method_summary: pd.DataFrame
    similarity: pd.DataFrame


def exon_gene_id(exon_id: object) -> str:
    """Extract the gene identifier used by the published GSE62669 exon table."""
    return str(exon_id).rsplit("_", 1)[0].strip()


def load_normalized_exon_counts(path: Path) -> pd.DataFrame:
    counts = pd.read_csv(path, sep="\t", index_col=0)
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    if counts.empty or counts.shape[1] < 4:
        raise ValueError("CoSplicEx requires a non-empty exon-by-sample matrix with >=4 samples.")
    if (counts.to_numpy(float) < 0).any():
        raise ValueError("Normalized exon counts must be non-negative.")
    return counts


def select_gene_profiles(
    exon_counts: pd.DataFrame,
    *,
    max_genes: int | None = None,
    min_exons: int = 2,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """Select multi-exon genes without using phenotype or validation information."""
    gene_ids = exon_counts.index.to_series().map(exon_gene_id)
    profiles: dict[str, pd.DataFrame] = {}
    rows: list[dict[str, object]] = []
    for gene_id, positions in gene_ids.groupby(gene_ids, sort=False).groups.items():
        profile = exon_counts.loc[list(positions)]
        if len(profile) < min_exons:
            continue
        sample_totals = profile.sum(axis=0)
        median_total = float(sample_totals.median())
        variable_exons = int(profile.std(axis=1, ddof=0).gt(0).sum())
        if median_total <= 0 or variable_exons < 2:
            continue
        rows.append(
            {
                "gene_id": str(gene_id),
                "n_exons": int(len(profile)),
                "n_variable_exons": variable_exons,
                "median_gene_exon_count": median_total,
                "mean_gene_exon_count": float(sample_totals.mean()),
            }
        )
        profiles[str(gene_id)] = profile
    audit = pd.DataFrame(rows)
    if audit.empty:
        raise ValueError("No eligible multi-exon genes were found.")
    audit = audit.sort_values(
        ["median_gene_exon_count", "n_variable_exons", "gene_id"],
        ascending=[False, False, True],
        kind="mergesort",
    )
    if max_genes is not None and max_genes > 0:
        audit = audit.head(max_genes)
    selected = audit["gene_id"].tolist()
    return {gene: profiles[gene] for gene in selected}, audit.reset_index(drop=True)


def canberra_distance_vectors(
    profiles: dict[str, pd.DataFrame],
    *,
    sample_indices: np.ndarray | None = None,
) -> pd.DataFrame:
    """Return one lower-triangle sample-distance vector per gene."""
    if not profiles:
        raise ValueError("At least one gene profile is required.")
    vectors: dict[str, np.ndarray] = {}
    pair_count: int | None = None
    for gene_id, profile in profiles.items():
        values = profile.to_numpy(dtype=float, copy=True).T
        if sample_indices is not None:
            values = values[np.asarray(sample_indices, dtype=int)]
        vector = pdist(values, metric="canberra")
        vector[~np.isfinite(vector)] = 0.0
        if pair_count is None:
            pair_count = len(vector)
        elif len(vector) != pair_count:
            raise ValueError("All gene profiles must contain the same sample set.")
        vectors[gene_id] = vector
    frame = pd.DataFrame(vectors)
    varying = frame.std(axis=0, ddof=0).gt(0)
    return frame.loc[:, varying]


def mantel_similarity(distance_vectors: pd.DataFrame) -> pd.DataFrame:
    """Compute the CoSplicEx Mantel correlation and unsigned power-6 similarity."""
    if distance_vectors.shape[1] < 2:
        raise ValueError("At least two non-constant gene distance profiles are required.")
    values = distance_vectors.to_numpy(dtype=float, copy=True)
    values -= values.mean(axis=0, keepdims=True)
    scale = np.sqrt(np.square(values).sum(axis=0, keepdims=True))
    scale[scale == 0] = 1.0
    standardized = values / scale
    correlation = np.clip(standardized.T @ standardized, -1.0, 1.0)
    np.fill_diagonal(correlation, 1.0)
    genes = distance_vectors.columns.astype(str)
    return pd.DataFrame(correlation, index=genes, columns=genes)


def unsigned_adjacency(
    mantel_correlation: pd.DataFrame,
    soft_power: float = 6.0,
) -> pd.DataFrame:
    values = np.power(np.abs(mantel_correlation.to_numpy(dtype=float)), soft_power)
    values = np.clip((values + values.T) / 2.0, 0.0, 1.0)
    np.fill_diagonal(values, 1.0)
    return pd.DataFrame(
        values,
        index=mantel_correlation.index,
        columns=mantel_correlation.columns,
    )


def _canonical_labels(raw: pd.Series, min_module_size: int) -> pd.Series:
    text = raw.astype(str).replace({"0": "grey", "nan": "grey"})
    counts = text.loc[text.ne("grey")].value_counts()
    retained = counts.loc[counts.ge(min_module_size)].index.tolist()
    retained.sort(key=lambda value: (-int(counts[value]), str(value)))
    names = {value: f"MCX{index:02d}" for index, value in enumerate(retained, start=1)}
    return text.map(names).fillna("grey")


def cluster_published_cosplicex(
    adjacency: pd.DataFrame,
    config: MantelCoSplicExConfig,
    *,
    rscript: str | Path,
) -> pd.Series:
    """Run the paper's average-linkage plus dynamic hybrid tree cut in R."""
    with tempfile.TemporaryDirectory(prefix="katana_mantel_cosplicex_") as temp_dir:
        temp = Path(temp_dir)
        similarity_path = temp / "adjacency.csv"
        output_path = temp / "modules.csv"
        script_path = temp / "cluster.R"
        adjacency.to_csv(similarity_path, index=True, index_label="gene_id")
        script_path.write_text(
            """
args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(dynamicTreeCut))
adj <- as.matrix(read.csv(args[1], row.names=1, check.names=FALSE))
adj[is.na(adj)] <- 0
adj <- (adj + t(adj)) / 2
diag(adj) <- 1
diss <- 1 - adj
diag(diss) <- 0
tree <- hclust(as.dist(diss), method='average')
modules <- cutreeHybrid(
  dendro=tree,
  distM=diss,
  cutHeight=as.numeric(args[3]),
  maxCoreScatter=as.numeric(args[4]),
  minClusterSize=as.integer(args[5]),
  deepSplit=as.integer(args[6]),
  pamStage=TRUE,
  maxPamDist=1,
  pamRespectsDendro=FALSE,
  useMedoids=FALSE,
  respectSmallClusters=FALSE,
  verbose=0
)
write.csv(
  data.frame(gene_id=rownames(adj), raw_module=modules$labels),
  args[2],
  row.names=FALSE,
  quote=TRUE
)
""".strip(),
            encoding="utf-8",
        )
        completed = subprocess.run(
            [
                str(rscript),
                str(script_path),
                str(similarity_path),
                str(output_path),
                str(config.cut_height),
                str(config.max_core_scatter),
                str(config.min_module_size),
                str(config.deep_split),
            ],
            capture_output=True,
            text=True,
            timeout=1800,
            check=False,
            env=_r_subprocess_env(rscript),
        )
        if completed.returncode != 0 or not output_path.exists():
            raise RuntimeError(
                "Mantel-CoSplicEx dynamic tree cut failed: "
                + (completed.stderr.strip() or completed.stdout.strip())
            )
        result = pd.read_csv(output_path).set_index("gene_id").reindex(adjacency.index)
        if result["raw_module"].isna().any():
            raise RuntimeError("Mantel-CoSplicEx output did not preserve the selected gene universe.")
        return _canonical_labels(result["raw_module"], config.min_module_size)


def fit_mantel_cosplicex(
    profiles: dict[str, pd.DataFrame],
    config: MantelCoSplicExConfig,
    *,
    rscript: str | Path,
    progress: Progress = None,
) -> tuple[pd.Series, pd.DataFrame]:
    if progress:
        progress("Computing per-gene Canberra sample-distance profiles")
    vectors = canberra_distance_vectors(profiles)
    if progress:
        progress("Computing gene-gene Mantel correlations")
    correlation = mantel_similarity(vectors)
    adjacency = unsigned_adjacency(correlation, config.soft_power)
    if progress:
        progress("Running published dynamic hybrid tree cut")
    labels = cluster_published_cosplicex(adjacency, config, rscript=rscript)
    return labels, correlation


def benchmark_mantel_cosplicex(
    exon_counts: pd.DataFrame,
    config: MantelCoSplicExConfig,
    *,
    rscript: str | Path,
    max_genes: int | None,
    progress: Progress = None,
) -> MantelCoSplicExResult:
    profiles, audit = select_gene_profiles(exon_counts, max_genes=max_genes)
    started = time.perf_counter()
    labels, similarity = fit_mantel_cosplicex(
        profiles,
        config,
        rscript=rscript,
        progress=progress,
    )
    initial_seconds = time.perf_counter() - started
    rng = np.random.default_rng(config.random_seed)
    n_samples = exon_counts.shape[1]
    sample_size = max(4, int(round(config.bootstrap_fraction * n_samples)))
    rows: list[dict[str, object]] = []
    for iteration in range(1, config.bootstrap_iterations + 1):
        indices = np.sort(rng.choice(n_samples, size=sample_size, replace=False))
        subset = {
            gene: profile.iloc[:, indices]
            for gene, profile in profiles.items()
        }
        current, _ = fit_mantel_cosplicex(
            subset,
            config,
            rscript=rscript,
        )
        aligned = current.reindex(labels.index).fillna("grey")
        rows.append(
            {
                "iteration": iteration,
                "sample_fraction": config.bootstrap_fraction,
                "adjusted_rand_index": adjusted_rand_score(labels, aligned),
            }
        )
        if progress:
            progress(
                f"Mantel-CoSplicEx bootstrap {iteration}/{config.bootstrap_iterations}"
            )
    stability = pd.DataFrame(rows)
    assignments = pd.DataFrame(
        {
            "gene_id": labels.index.astype(str),
            "mantel_cosplicex_module": labels.to_numpy(object),
        }
    )
    counts = labels.loc[labels.ne("grey")].value_counts()
    non_grey = int(counts.sum())
    n_modules = int(len(counts))
    summary = pd.DataFrame(
        [
            {
                "method": "Mantel-CoSplicEx (published gene-level comparator)",
                "biological_unit": "gene-level multi-exon profile",
                "n_features": int(len(labels)),
                "n_samples": int(n_samples),
                "n_modules": n_modules,
                "grey_fraction": float(labels.eq("grey").mean()),
                "largest_module_fraction": (
                    float(counts.max() / max(non_grey, 1)) if n_modules else 0.0
                ),
                "valid_partition": bool(
                    n_modules >= 2
                    and float(counts.max() / max(non_grey, 1)) <= 0.75
                ),
                "mean_bootstrap_ARI": (
                    float(stability["adjusted_rand_index"].mean())
                    if not stability.empty
                    else np.nan
                ),
                "median_bootstrap_ARI": (
                    float(stability["adjusted_rand_index"].median())
                    if not stability.empty
                    else np.nan
                ),
                "min_bootstrap_ARI": (
                    float(stability["adjusted_rand_index"].min())
                    if not stability.empty
                    else np.nan
                ),
                "runtime_seconds_initial_fit": initial_seconds,
                "soft_power": config.soft_power,
                "distance": "Canberra",
                "network_similarity": "absolute Mantel/Pearson correlation",
                "clustering": "average linkage + dynamicTreeCut::cutreeHybrid",
            }
        ]
    )
    return MantelCoSplicExResult(
        assignments=assignments,
        gene_audit=audit,
        bootstrap_stability=stability,
        method_summary=summary,
        similarity=similarity,
    )
