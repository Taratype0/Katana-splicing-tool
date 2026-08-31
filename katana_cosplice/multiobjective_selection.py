from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from katana_cosplice.config import (
    ExternalCohortConfig,
    MultiObjectiveSelectionConfig,
)
from katana_cosplice.utils import bh_adjust


def _read_matrix(path: Path) -> pd.DataFrame:
    matrix = pd.read_csv(path, sep=None, engine="python", index_col=0)
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    if not matrix.index.is_unique or not matrix.columns.is_unique:
        raise ValueError(f"Tuning PSI matrix labels must be unique: {path}")
    if matrix.shape[1] < 3:
        raise ValueError(
            f"Tuning PSI matrix requires at least three samples: {path}"
        )
    return matrix


def _validate_sample_metadata(path: Path, matrix_columns: pd.Index) -> None:
    samples = pd.read_csv(path, sep=None, engine="python")
    sample_column = (
        "sample_id"
        if "sample_id" in samples
        else "sample"
        if "sample" in samples
        else ""
    )
    if not sample_column:
        raise ValueError(
            f"Tuning sample metadata requires sample_id or sample: {path}"
        )
    available = set(samples[sample_column].astype(str))
    missing = [
        sample
        for sample in matrix_columns.astype(str)
        if sample not in available
    ]
    if missing:
        raise ValueError(
            "Tuning sample metadata does not cover matrix columns "
            f"for {path}: {missing[:5]}"
        )


def _load_gmt(path: Path) -> dict[str, set[str]]:
    gene_sets: dict[str, set[str]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            members = {
                value.strip().casefold()
                for value in fields[2:]
                if value.strip()
            }
            if members:
                gene_sets[fields[0]] = members
    if not gene_sets:
        raise ValueError(f"No usable gene sets were loaded from tuning GMT: {path}")
    return gene_sets


def _standardize_rows(frame: pd.DataFrame) -> np.ndarray:
    values = frame.to_numpy(dtype=float)
    means = np.nanmean(values, axis=1, keepdims=True)
    means = np.where(np.isfinite(means), means, 0.0)
    values = np.where(np.isfinite(values), values, means)
    centered = values - values.mean(axis=1, keepdims=True)
    scale = centered.std(axis=1, ddof=1, keepdims=True)
    return centered / np.where(scale < 1e-8, 1.0, scale)


def _absolute_correlation(standardized: np.ndarray) -> np.ndarray:
    denominator = max(1, standardized.shape[1] - 1)
    correlation = (standardized @ standardized.T) / denominator
    correlation = np.clip(np.abs(correlation), 0.0, 1.0)
    np.fill_diagonal(correlation, 1.0)
    return correlation


def _safe_correlation(left: np.ndarray, right: np.ndarray) -> float:
    if (
        len(left) < 3
        or np.std(left) < 1e-8
        or np.std(right) < 1e-8
    ):
        return 0.0
    value = float(np.corrcoef(left, right)[0, 1])
    return value if np.isfinite(value) else 0.0


def _preservation_components(
    reference: pd.DataFrame,
    target: pd.DataFrame,
) -> dict[str, float]:
    reference_values = _standardize_rows(reference)
    target_values = _standardize_rows(target)
    reference_corr = _absolute_correlation(reference_values)
    target_corr = _absolute_correlation(target_values)
    upper = np.triu_indices(len(reference), k=1)
    reference_density = (
        float(reference_corr[upper].mean()) if len(upper[0]) else 0.0
    )
    target_density = (
        float(target_corr[upper].mean()) if len(upper[0]) else 0.0
    )
    epsilon = 1e-6
    density_similarity = float(
        np.exp(
            -abs(
                np.log(
                    (target_density + epsilon)
                    / (reference_density + epsilon)
                )
            )
        )
    )
    reference_connectivity = reference_corr.sum(axis=1) - 1.0
    target_connectivity = target_corr.sum(axis=1) - 1.0
    connectivity = float(
        np.clip(
            _safe_correlation(reference_connectivity, target_connectivity),
            0.0,
            1.0,
        )
    )

    _ref_u, reference_s, reference_vt = np.linalg.svd(
        reference_values, full_matrices=False
    )
    _target_u, target_s, target_vt = np.linalg.svd(
        target_values, full_matrices=False
    )
    reference_score = reference_vt[0]
    target_score = target_vt[0]
    reference_membership = np.asarray(
        [_safe_correlation(row, reference_score) for row in reference_values]
    )
    target_membership = np.asarray(
        [_safe_correlation(row, target_score) for row in target_values]
    )
    membership = abs(
        _safe_correlation(reference_membership, target_membership)
    )
    reference_total = float(np.square(reference_s).sum())
    target_total = float(np.square(target_s).sum())
    reference_variance = (
        float(reference_s[0] ** 2 / reference_total)
        if reference_total
        else 0.0
    )
    target_variance = (
        float(target_s[0] ** 2 / target_total) if target_total else 0.0
    )
    variance_similarity = (
        min(
            reference_variance / target_variance,
            target_variance / reference_variance,
        )
        if reference_variance > 0 and target_variance > 0
        else 0.0
    )
    program_preservation = float(
        np.clip(membership * variance_similarity, 0.0, 1.0)
    )
    module_score = float(
        np.prod(
            np.clip(
                [density_similarity, connectivity, program_preservation],
                1e-6,
                1.0,
            )
        )
        ** (1.0 / 3.0)
    )
    return {
        "density_similarity": density_similarity,
        "connectivity_preservation": connectivity,
        "membership_profile_correlation": membership,
        "program_variance_similarity": variance_similarity,
        "program_score_preservation": program_preservation,
        "module_preservation_score": module_score,
    }


@dataclass
class MultiObjectiveCandidateScorer:
    """Score candidate partitions using tuning data only.

    Holdout cohorts are intentionally not accepted by this object. They remain
    available for confirmatory evaluation after the partition is locked.
    """

    reference_matrix: pd.DataFrame
    event_metadata: pd.DataFrame
    tuning_matrices: dict[str, pd.DataFrame]
    gene_sets: dict[str, set[str]]
    settings: MultiObjectiveSelectionConfig

    def __post_init__(self) -> None:
        self.reference_matrix = self.reference_matrix.copy()
        self.reference_matrix.index = self.reference_matrix.index.astype(str)
        metadata = self.event_metadata.copy()
        metadata["event_key"] = metadata["event_key"].astype(str)
        metadata = metadata.drop_duplicates("event_key").set_index("event_key")
        self._event_gene = pd.Series(index=self.reference_matrix.index, dtype=object)
        for event_key in self.reference_matrix.index:
            if event_key in metadata.index:
                row = metadata.loc[event_key]
                symbol = str(row.get("gene_symbol", "") or "").strip()
                gene_id = str(row.get("gene_id", "") or "").strip().split(".", 1)[0]
                self._event_gene.loc[event_key] = (
                    symbol.casefold() if symbol else gene_id.casefold()
                )
            else:
                self._event_gene.loc[event_key] = (
                    str(event_key).split("|", 1)[0].split(".", 1)[0].casefold()
                )
        self._annotated_genes = set().union(*self.gene_sets.values())
        self._cache: dict[
            tuple[str, ...], tuple[dict[str, float | int], pd.DataFrame]
        ] = {}

    def __call__(self, labels: np.ndarray) -> dict[str, float | int]:
        metrics, _audit = self._evaluate(labels)
        return metrics.copy()

    def audit(self, labels: np.ndarray) -> pd.DataFrame:
        _metrics, audit = self._evaluate(labels)
        return audit.copy()

    def _evaluate(
        self, labels: np.ndarray
    ) -> tuple[dict[str, float | int], pd.DataFrame]:
        key = tuple(map(str, labels))
        if key in self._cache:
            return self._cache[key]
        label_series = pd.Series(
            np.asarray(labels, dtype=object),
            index=self.reference_matrix.index,
            dtype=object,
        )
        label_series = label_series.loc[label_series.ne("grey")]
        preservation_scores: list[float] = []
        audit_rows: list[dict[str, object]] = []
        eligible_preservation_modules = 0
        for cohort, target in self.tuning_matrices.items():
            shared = self.reference_matrix.index.intersection(target.index)
            for module, module_events in label_series.groupby(
                label_series, observed=True
            ):
                event_keys = module_events.index.intersection(shared)
                if (
                    len(event_keys)
                    < self.settings.min_shared_events_per_module
                ):
                    continue
                eligible_preservation_modules += 1
                components = _preservation_components(
                    self.reference_matrix.loc[event_keys],
                    target.loc[event_keys],
                )
                shared_fraction = float(
                    len(event_keys) / max(1, len(module_events))
                )
                score = float(
                    np.sqrt(
                        components["module_preservation_score"]
                        * shared_fraction
                    )
                )
                preservation_scores.append(score)
                audit_rows.append(
                    {
                        "audit_type": "external_preservation",
                        "cohort": cohort,
                        "module": str(module),
                        "shared_events": len(event_keys),
                        "reference_module_events": len(module_events),
                        "shared_event_fraction": shared_fraction,
                        **components,
                        "objective_contribution": score,
                    }
                )
        preservation_score = (
            float(np.median(preservation_scores))
            if preservation_scores
            else 0.0
        )

        module_genes: dict[str, set[str]] = {}
        for module, module_events in label_series.groupby(
            label_series, observed=True
        ):
            genes = {
                str(value)
                for value in self._event_gene.reindex(module_events.index)
                if pd.notna(value) and str(value)
            }
            if len(genes) >= self.settings.min_module_genes:
                module_genes[str(module)] = genes
        universe = set().union(*module_genes.values()) if module_genes else set()
        enrichment_rows: list[dict[str, object]] = []
        for module, genes in module_genes.items():
            for term, members in self.gene_sets.items():
                selected = members & universe
                overlap = genes & selected
                pvalue = (
                    float(
                        hypergeom.sf(
                            len(overlap) - 1,
                            len(universe),
                            len(selected),
                            len(genes),
                        )
                    )
                    if overlap and selected and universe
                    else 1.0
                )
                enrichment_rows.append(
                    {
                        "module": module,
                        "term": term,
                        "module_genes": len(genes),
                        "overlap_count": len(overlap),
                        "pvalue": pvalue,
                    }
                )
        functional_scores: list[float] = []
        significant_functional_modules = 0
        enrichment = pd.DataFrame(enrichment_rows)
        if not enrichment.empty:
            # Correct terms within each module so candidates with more modules
            # are not penalized solely for having a larger total test count.
            enrichment["FDR"] = enrichment.groupby(
                "module", observed=True
            )["pvalue"].transform(lambda values: bh_adjust(values))
            for module, rows in enrichment.groupby("module", observed=True):
                best = rows.sort_values(["FDR", "pvalue"]).iloc[0]
                genes = module_genes[str(module)]
                annotation_coverage = float(
                    len(genes & self._annotated_genes) / max(1, len(genes))
                )
                fdr = float(np.clip(best["FDR"], 1e-300, 1.0))
                if fdr <= self.settings.functional_fdr_cutoff:
                    significant_functional_modules += 1
                pvalue = float(np.clip(best["pvalue"], 1e-300, 1.0))
                fdr_scale = max(
                    1e-12,
                    -np.log10(self.settings.functional_fdr_cutoff),
                )
                fdr_strength = float(
                    np.clip(-np.log10(fdr) / fdr_scale, 0.0, 1.0)
                )
                # Formal claims still use the reported FDR. Candidate selection
                # needs a continuous objective, however: with thousands of GMT
                # terms, every candidate can otherwise receive exactly zero
                # before any partition is frozen. A small raw-p contribution
                # preserves ranking information without turning it into a
                # significance claim; untouched holdout replication remains the
                # confirmatory endpoint.
                raw_strength = float(
                    np.clip(-np.log10(pvalue) / -np.log10(0.05), 0.0, 1.0)
                )
                enrichment_strength = 0.75 * fdr_strength + 0.25 * raw_strength
                score = annotation_coverage * enrichment_strength
                functional_scores.append(score)
                audit_rows.append(
                    {
                        "audit_type": "functional_coherence",
                        "cohort": "",
                        "module": str(module),
                        "module_genes": len(genes),
                        "annotated_gene_fraction": annotation_coverage,
                        "best_term": best["term"],
                        "best_term_overlap": int(best["overlap_count"]),
                        "best_term_pvalue": pvalue,
                        "best_term_FDR": fdr,
                        "fdr_enrichment_strength": fdr_strength,
                        "raw_pvalue_enrichment_strength": raw_strength,
                        "enrichment_strength": enrichment_strength,
                        "objective_contribution": score,
                    }
                )
        functional_score = (
            float(np.median(functional_scores))
            if functional_scores
            else 0.0
        )
        metrics: dict[str, float | int] = {
            "objective_external_preservation": preservation_score,
            "objective_functional_coherence": functional_score,
            "objective_stability_weight": self.settings.stability_weight,
            "objective_preservation_weight": self.settings.preservation_weight,
            "objective_functional_weight": self.settings.functional_weight,
            "functional_weight_mode": self.settings.functional_weight_mode,
            "min_significant_functional_modules": (
                self.settings.min_significant_functional_modules
            ),
            "objective_score_floor": self.settings.score_floor,
            "n_tuning_cohorts_scored": len(self.tuning_matrices),
            "n_preservation_module_cohort_tests": eligible_preservation_modules,
            "n_functional_modules_scored": len(functional_scores),
            "n_significant_functional_modules": significant_functional_modules,
        }
        audit = pd.DataFrame(audit_rows)
        self._cache[key] = (metrics, audit)
        return metrics, audit


def build_multiobjective_candidate_scorer(
    reference_matrix: pd.DataFrame,
    event_metadata: pd.DataFrame,
    external_cohorts: list[ExternalCohortConfig],
    settings: MultiObjectiveSelectionConfig,
) -> MultiObjectiveCandidateScorer:
    if not settings.enabled:
        raise ValueError("Multi-objective scorer requested while configuration is disabled")
    cohort_lookup = {cohort.name: cohort for cohort in external_cohorts}
    tuning_matrices: dict[str, pd.DataFrame] = {}
    for name in settings.tuning_cohorts:
        cohort = cohort_lookup[name]
        if not cohort.psi_matrix.exists():
            raise FileNotFoundError(
                f"Missing tuning PSI matrix for {name}: {cohort.psi_matrix}"
            )
        if not cohort.sample_metadata.exists():
            raise FileNotFoundError(
                f"Missing tuning sample metadata for {name}: "
                f"{cohort.sample_metadata}"
            )
        tuning_matrices[name] = _read_matrix(cohort.psi_matrix)
        _validate_sample_metadata(
            cohort.sample_metadata, tuning_matrices[name].columns
        )
    gmt = settings.functional_gene_sets_gmt
    if gmt is None or not gmt.exists():
        raise FileNotFoundError(f"Missing tuning functional GMT: {gmt}")
    return MultiObjectiveCandidateScorer(
        reference_matrix=reference_matrix,
        event_metadata=event_metadata,
        tuning_matrices=tuning_matrices,
        gene_sets=_load_gmt(gmt),
        settings=settings,
    )
