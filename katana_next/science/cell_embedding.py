"""Label-blind cell embeddings and held-out-donor checks, not differential tests."""
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.metrics import balanced_accuracy_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler


@dataclass(frozen=True)
class CellEmbeddingResult:
    cells: pd.DataFrame
    features: pd.DataFrame
    coordinates: pd.DataFrame
    donor_validation: pd.DataFrame
    diagnostics: pd.DataFrame
    distributions: pd.DataFrame
    latent_coordinates: pd.DataFrame


def cell_feature_values(counts, *, mode="splicing", min_reads=10):
    if mode == "splicing":
        required = {"cell_id", "event_id", "inclusion", "exclusion"}
        if not required <= set(counts):
            raise ValueError(f"Splicing embedding needs {sorted(required)}")
        frame = counts.copy()
        values = frame[["inclusion", "exclusion"]].apply(pd.to_numeric, errors="raise").to_numpy()
        if not np.isfinite(values).all() or (values < 0).any() or (values != np.floor(values)).any():
            raise ValueError("Cell-event inputs must be nonnegative integer evidence counts")
        frame["depth"] = values.sum(axis=1)
        frame["value"] = np.where(frame.depth.ge(min_reads), (values[:, 0] + .5) / (frame.depth + 1), np.nan)
        frame["feature_id"] = frame.event_id
    elif mode == "expression":
        if not {"cell_id", "gene_id", "count"} <= set(counts):
            raise ValueError("Expression embedding needs long-form cell_id/gene_id/count; TPM is not read counts")
        frame = counts.copy()
        values = pd.to_numeric(frame["count"], errors="raise").to_numpy()
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError("Invalid expression count values")
        frame["depth"] = values
        totals = frame.groupby("cell_id").depth.transform("sum")
        frame["value"] = np.log1p(frame.depth / totals.replace(0, np.nan) * 10000)
        frame["feature_id"] = frame.gene_id
    else:
        raise ValueError("mode must be splicing or expression")
    for col in ("cell_id", "feature_id"):
        if frame[col].isna().any() or frame[col].astype(str).str.strip().eq("").any():
            raise ValueError("Every measurement needs cell and feature identity")
        frame[col] = frame[col].astype(str)
    if frame.duplicated(["cell_id", "feature_id"]).any():
        raise ValueError("Duplicate cell-feature observations; aggregate technical counts before import")
    return frame


def _feature_panel(values, training_cells, *, min_fraction, max_features, expression=False):
    observed = values.loc[values.cell_id.isin(training_cells)]
    stats = observed.groupby("feature_id").value.agg(["count", "var"])
    if expression:
        # Sparse expression omissions are zeros; missing splicing measurements are not.
        n = len(training_cells)
        sums = observed.groupby("feature_id").value.sum()
        squares = observed.assign(square=observed.value**2).groupby("feature_id").square.sum()
        stats["var"] = ((squares - sums**2 / n) / max(n - 1, 1)).clip(lower=0)
        stats["count"] = observed.loc[observed.value.gt(0)].groupby("feature_id").size().reindex(stats.index, fill_value=0)
    stats["measured_fraction"] = stats["count"] / len(training_cells)
    stats["retained"] = stats.measured_fraction.ge(min_fraction) & stats["var"].gt(0)
    selected = stats.loc[stats.retained].sort_values(["var", "feature_id"], ascending=[False, True]).head(max_features).index
    stats["selected_for_embedding"] = stats.index.isin(selected)
    return stats.reset_index(), list(selected)


def _matrix(values, cells, features, *, expression=False):
    if len(cells) * len(features) * 8 > 256 * 1024**2:
        raise ValueError("Selected embedding matrix exceeds 256 MiB; explicitly reduce cells/features")
    matrix = values.loc[values.feature_id.isin(features)].pivot(index="cell_id", columns="feature_id", values="value")
    matrix = matrix.reindex(index=cells, columns=features)
    return matrix.fillna(0) if expression else matrix


def build_cell_embedding(counts, metadata, *, mode="splicing", min_reads=10,
                         min_feature_fraction=.2, min_observed_features=5, max_features=500,
                         n_neighbors=15, min_dist=.3, pca_components=30, random_seed=17,
                         gene_filter="", evaluate_cell_types=True):
    required = {"cell_id", "donor_id", "cell_type", "condition"}
    if not required <= set(metadata) or metadata.cell_id.duplicated().any():
        raise ValueError("Unique cell metadata with donor, cell type and condition is required")
    if metadata[list(required)].isna().any().any() or metadata[list(required)].astype(str).apply(lambda s: s.str.strip().eq("")).any().any():
        raise ValueError("Embedding metadata must contain nonempty identifiers and labels")
    if not 0 < min_feature_fraction <= 1 or max_features < 2 or n_neighbors < 2 or min_observed_features < 1 or min_reads < 1:
        raise ValueError("Invalid embedding coverage or dimension parameters")
    if not np.isfinite(min_dist) or not 0 <= min_dist <= 1 or pca_components < 2:
        raise ValueError("UMAP min_dist must be between 0 and 1; PCA components >= 2")
    metadata = metadata.copy()
    for column in required:
        metadata[column] = metadata[column].astype(str)
    values = cell_feature_values(counts, mode=mode, min_reads=min_reads)
    unknown = set(values.cell_id) - set(metadata.cell_id)
    if unknown:
        raise ValueError(f"Measurement cells absent from metadata: {sorted(unknown)[:5]}")
    full_depth = values.groupby("cell_id").depth.sum()
    if gene_filter:
        keep = values.get("gene_id", pd.Series("", index=values.index)).eq(gene_filter)
        keep |= values.get("gene_symbol", pd.Series("", index=values.index)).eq(gene_filter)
        values = values.loc[keep].copy()
        if values.feature_id.nunique() < 2:
            raise ValueError("A gene-specific embedding needs at least two measured variable features; do not force a 2-D map of one event")
    measurable = values.value.notna() if mode == "splicing" else values.value.gt(0)
    counts_measured = values.loc[measurable].groupby("cell_id").feature_id.nunique()
    cells = metadata.copy()
    cells["observed_features"] = cells.cell_id.map(counts_measured).fillna(0).astype(int)
    cells["total_evidence_counts"] = cells.cell_id.map(full_depth).fillna(0)
    cells["included"] = cells.observed_features.ge(min_observed_features)
    cells["exclusion_reason"] = np.where(cells.included, "", "too_few_measured_features")
    included = cells.loc[cells.included].copy()
    if len(included) < 5:
        raise ValueError("Fewer than five cells pass the declared measurement-coverage gate")
    features, panel = _feature_panel(values, set(included.cell_id), min_fraction=min_feature_fraction, max_features=max_features, expression=mode == "expression")
    if len(panel) < 2:
        raise ValueError("Fewer than two variable features meet the label-blind coverage gate")
    raw = _matrix(values, list(included.cell_id), panel, expression=mode == "expression")
    included["missing_feature_fraction"] = raw.isna().mean(axis=1).to_numpy()
    imputer = SimpleImputer(strategy="median")
    scaler = StandardScaler()
    transformed = scaler.fit_transform(imputer.fit_transform(raw))
    components = min(int(pca_components), len(included) - 1, len(panel))
    pca = PCA(n_components=components, random_state=random_seed)
    latent = pca.fit_transform(transformed)
    from umap import UMAP
    embedding = UMAP(n_components=2, n_neighbors=min(n_neighbors, len(included) - 1),
                     min_dist=min_dist, random_state=random_seed, n_jobs=1).fit_transform(latent)
    coordinates = included.copy()
    coordinates["umap_1"], coordinates["umap_2"] = embedding[:, 0], embedding[:, 1]
    coordinates["pc1"], coordinates["pc2"] = latent[:, 0], latent[:, 1]
    latent_coordinates = pd.DataFrame(latent, columns=[f"pc_{i + 1}" for i in range(components)])
    latent_coordinates.insert(0, "cell_id", included.cell_id.to_numpy())
    latent_coordinates["row_id"] = latent_coordinates.cell_id
    coordinates["embedding_source"] = mode
    coordinates["imputation"] = "median of measured cells for visualization" if mode == "splicing" else "sparse count zeros; log1p(CP10K)"
    rows = []
    if evaluate_cell_types:
        for donor, heldout in included.groupby("donor_id"):
            train = included.loc[included.donor_id.ne(donor)]
            record = {"donor_id": donor, "test_cells": len(heldout), "training_cells": len(train),
                      "mode": mode, "balanced_accuracy": np.nan, "majority_baseline": np.nan,
                      "depth_missingness_baseline": np.nan, "features": 0,
                      "status": "not_evaluable"}
            if train.donor_id.nunique() < 2 or train.cell_type.nunique() < 2:
                record["reason"] = "need at least two training donors and cell types"
            else:
                _, train_panel = _feature_panel(values, set(train.cell_id), min_fraction=min_feature_fraction, max_features=max_features, expression=mode == "expression")
                record["features"] = len(train_panel)
                if len(train_panel) < 2:
                    record["reason"] = "too_few_training_features"
                else:
                    a = _matrix(values, list(train.cell_id), train_panel, expression=mode == "expression")
                    b = _matrix(values, list(heldout.cell_id), train_panel, expression=mode == "expression")
                    nuisance_a = np.column_stack((np.log1p(train.total_evidence_counts), a.isna().mean(axis=1)))
                    nuisance_b = np.column_stack((np.log1p(heldout.total_evidence_counts), b.isna().mean(axis=1)))
                    nuisance_scale = StandardScaler()
                    nuisance_model = KNeighborsClassifier(n_neighbors=min(n_neighbors, len(a)))
                    nuisance_predictions = nuisance_model.fit(nuisance_scale.fit_transform(nuisance_a), train.cell_type).predict(nuisance_scale.transform(nuisance_b))
                    imp, scale = SimpleImputer(strategy="median"), StandardScaler()
                    a = scale.fit_transform(imp.fit_transform(a))
                    b = scale.transform(imp.transform(b))
                    model = KNeighborsClassifier(n_neighbors=min(n_neighbors, len(a)))
                    predictions = model.fit(a, train.cell_type).predict(b)
                    majority = train.cell_type.value_counts().index[0]
                    record.update(status="evaluated", reason="", balanced_accuracy=balanced_accuracy_score(heldout.cell_type, predictions),
                                  depth_missingness_baseline=balanced_accuracy_score(heldout.cell_type, nuisance_predictions),
                                  majority_baseline=balanced_accuracy_score(heldout.cell_type, np.repeat(majority, len(heldout))))
            rows.append(record)
    diagnostics = [{"metric": "cells_input", "value": len(cells)}, {"metric": "cells_retained", "value": len(included)},
                   {"metric": "variable_features_retained", "value": len(panel)},
                   {"metric": "PCA_variance_explained", "value": float(pca.explained_variance_ratio_.sum())}]
    for axis in ("umap_1", "umap_2"):
        for nuisance in ("total_evidence_counts", "missing_feature_fraction"):
            x, y = coordinates[axis], coordinates[nuisance]
            rho = float(spearmanr(x, y).statistic) if x.nunique() > 1 and y.nunique() > 1 else np.nan
            diagnostics.append({"metric": f"{axis}_vs_{nuisance}_spearman", "value": rho})
    distributions = values.loc[values.cell_id.isin(included.cell_id) & values.feature_id.isin(panel[:8]), ["cell_id", "feature_id", "value", "depth"]]
    distributions = distributions.merge(included[["cell_id", "cell_type", "donor_id", "condition"]], on="cell_id", validate="many_to_one")
    for frame in (cells, features, coordinates, distributions):
        frame["row_id"] = frame.index.astype(str)
    validation = pd.DataFrame(rows, columns=["donor_id", "test_cells", "training_cells", "mode", "balanced_accuracy", "majority_baseline", "depth_missingness_baseline", "features", "status", "reason"])
    validation["row_id"] = validation.index.astype(str)
    return CellEmbeddingResult(cells, features, coordinates, validation, pd.DataFrame(diagnostics), distributions, latent_coordinates)
