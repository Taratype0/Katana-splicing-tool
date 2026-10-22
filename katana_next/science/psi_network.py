"""Public PSI-only network contract, distinct from read-supported rMATS QC."""
from __future__ import annotations
from dataclasses import dataclass, asdict
import numpy as np
import pandas as pd
from katana_cosplice.config import NetworkConfig, ConsensusConfig
from katana_cosplice.network import fit_network
from katana_cosplice.consensus import build_bootstrap_coassignment_consensus


@dataclass(frozen=True)
class PsiFilterConfig:
    min_valid_fraction: float = .8
    min_variance: float = .0025
    min_range: float = .5
    min_alternative_fraction: float = .1

    def validate(self):
        for key, value in asdict(self).items():
            if not np.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f"{key} must be finite and within [0,1]")


def psi_filter_audit(matrix, filters=None):
    """Evaluate every gate independently, retaining overlapping failure reasons."""
    filters = filters or PsiFilterConfig()
    filters.validate()
    values = matrix.to_numpy(dtype=float)
    if not matrix.index.is_unique or not matrix.columns.is_unique:
        raise ValueError("Unique event and sample IDs are required")
    if np.isinf(values).any() or np.any((values[np.isfinite(values)] < 0) | (values[np.isfinite(values)] > 1)):
        raise ValueError("Measured PSI must lie in [0,1]; missing values remain missing")
    metrics = pd.DataFrame(index=matrix.index)
    metrics["valid_fraction"] = matrix.notna().mean(axis=1)
    metrics["psi_variance"] = matrix.var(axis=1)
    metrics["psi_range"] = matrix.max(axis=1)-matrix.min(axis=1)
    # Keep the existing denominator: all discovery samples, including missing PSI.
    metrics["alternative_fraction"] = matrix.ge(.2).mul(matrix.le(.8)).mean(axis=1)
    gates = {
        "valid_fraction": ("valid_fraction", filters.min_valid_fraction),
        "variance": ("psi_variance", filters.min_variance),
        "range": ("psi_range", filters.min_range),
        "alternative_fraction": ("alternative_fraction", filters.min_alternative_fraction),
    }
    for name, (column, threshold) in gates.items():
        metrics["threshold_"+name] = threshold
        metrics["failed_"+name] = ~metrics[column].ge(threshold)
    failures = ["failed_"+name for name in gates]
    metrics["failure_count"] = metrics[failures].sum(axis=1)
    metrics["filter_reasons"] = [";".join(name for name in gates if row["failed_"+name]) or "passed"
                                 for _, row in metrics.iterrows()]
    metrics["passed"] = metrics.failure_count.eq(0)
    metrics["stability_score"] = metrics.valid_fraction
    metrics["read_support_status"] = "not_assessed_PSI_only"
    metrics["repeat_support_status"] = "not_assessed_PSI_only"
    return metrics


def prepare_psi_network(matrix, samples, *, filters=None, filter_profile='legacy', max_features=2000):
    if not matrix.index.is_unique or not matrix.columns.is_unique or not samples.sample_id.is_unique:
        raise ValueError("Unique event, sample and metadata IDs are required")
    design = samples.set_index("sample_id").loc[matrix.columns]
    if not {"condition", "biological_sample_id"}.issubset(design):
        raise ValueError("Declare condition and biological_sample_id; cells are not donor replicates")
    if design.biological_sample_id.isna().any():
        raise ValueError("Missing biological sample identity")
    if design.biological_sample_id.nunique() < 15:
        raise ValueError("At least 15 independent biological samples are required")
    if filter_profile not in {'legacy','coverage_variable'}:
        raise ValueError('Unknown PSI filter_profile')
    if filter_profile=='coverage_variable':
        from katana_next.science.psi_robust import coverage_variable_audit
        metrics=coverage_variable_audit(matrix,filters or PsiFilterConfig(),max_features=max_features)
    else:
        metrics = psi_filter_audit(matrix, filters)
        metrics['filter_profile']='legacy'
    selected = matrix.loc[metrics.passed].copy()
    if len(selected) < 40 and filter_profile=='legacy':
        raise ValueError(f"Only {len(selected)} events qualify for PSI-only network fitting (minimum 40)")
    # Learn imputation means from discovery only, never a holdout or phenotype.
    means = selected.mean(axis=1)
    complete = selected.T.fillna(means).T
    return selected, complete, metrics, means


def fit_psi_backbone(complete, qc, network, samples):
    if len(complete)>=40:
        return fit_network(complete,qc.loc[complete.index],network,samples=samples)
    from katana_cosplice.network import NetworkResult
    return NetworkResult(
        pd.DataFrame({'event_key':complete.index,'module':'grey'}),
        pd.DataFrame(columns=['source','target','weight']),pd.DataFrame(index=complete.columns),
        pd.DataFrame([dict(status='no_estimable_network',reason='fewer_than_40_eligible_events',
                           eligible_events=len(complete),n_final_modules=0)]),
        pd.DataFrame(),complete,None)


def fit_psi_network(matrix, samples, *, rscript=None, progress=None, filters=None, network=None,
                    filter_profile='legacy',max_features=2000):
    from katana_next.science.frozen_psi import _samples
    _samples(matrix, samples,allow_empty_events=filter_profile=='coverage_variable')
    measured, complete, qc, means = prepare_psi_network(matrix, samples, filters=filters,
        filter_profile=filter_profile,max_features=max_features)
    network = network or NetworkConfig(wgcna_rscript=rscript)
    baseline = fit_psi_backbone(complete,qc,network,samples)
    if baseline.program_scores.shape[1]==0:
        from katana_cosplice.consensus import ConsensusResult
        confidence=ConsensusResult(baseline.assignments.copy(),pd.DataFrame(),baseline.program_scores.copy(),
            pd.DataFrame([dict(status='not_run_no_estimable_modules',iterations=0)]),
            pd.DataFrame(),np.empty((0,0)),np.array([],dtype=str),baseline.diagnostics.copy(),pd.DataFrame(),pd.DataFrame())
        return measured,qc,means,baseline,confidence
    if progress:
        progress("Discovery WGCNA fit complete; starting 50 module-confidence refits")
    confidence = build_bootstrap_coassignment_consensus(complete, baseline.assignments, samples, network,
        ConsensusConfig(bootstrap_iterations=50, clustering_method="backbone_confidence_only"), progress=progress)
    return measured, qc, means, baseline, confidence
