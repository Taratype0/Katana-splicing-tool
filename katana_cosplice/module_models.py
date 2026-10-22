from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from scipy.optimize import minimize
from scipy.special import betaln, expit, gammaln, logit, logsumexp
from scipy.stats import chi2, spearmanr

from katana_cosplice.config import ModuleStatisticsConfig
from katana_cosplice.utils import bh_adjust


@dataclass
class ModuleStatisticsResult:
    event_models: pd.DataFrame
    module_models: pd.DataFrame
    variance_partition: pd.DataFrame
    simple_delta_benchmark: pd.DataFrame
    model_diagnostics: pd.DataFrame


def _design(metadata: pd.DataFrame, condition_column: str, batch_column: str | None) -> tuple[np.ndarray, np.ndarray, list[str], list[str]]:
    condition = pd.get_dummies(metadata[condition_column].astype(str), prefix="condition", drop_first=True, dtype=float)
    batch = (
        pd.get_dummies(metadata[batch_column].astype(str), prefix="batch", drop_first=True, dtype=float)
        if batch_column and batch_column in metadata
        else pd.DataFrame(index=metadata.index)
    )
    intercept = np.ones((len(metadata), 1), dtype=float)
    null = np.column_stack([intercept, batch.to_numpy(dtype=float)])
    full = np.column_stack([null, condition.to_numpy(dtype=float)])
    null_names = ["Intercept", *batch.columns.tolist()]
    full_names = [*null_names, *condition.columns.tolist()]
    return null, full, null_names, full_names


def _beta_binomial_nll(parameters: np.ndarray, design: np.ndarray, successes: np.ndarray, totals: np.ndarray) -> float:
    beta = parameters[:-1]
    rho = 1e-7 + (0.95 - 1e-7) * expit(parameters[-1])
    mu = np.clip(expit(design @ beta), 1e-7, 1.0 - 1e-7)
    concentration = (1.0 - rho) / rho
    alpha = mu * concentration
    beta_shape = (1.0 - mu) * concentration
    failures = totals - successes
    log_choose = gammaln(totals + 1.0) - gammaln(successes + 1.0) - gammaln(failures + 1.0)
    log_likelihood = log_choose + betaln(successes + alpha, failures + beta_shape) - betaln(alpha, beta_shape)
    if not np.isfinite(log_likelihood).all():
        return 1e100
    return float(-log_likelihood.sum())


def _fit_model(design: np.ndarray, successes: np.ndarray, totals: np.ndarray, max_iterations: int) -> tuple[np.ndarray, float, bool, str]:
    pooled = float((successes.sum() + 0.5) / (totals.sum() + 1.0))
    initial = np.zeros(design.shape[1] + 1, dtype=float)
    initial[0] = float(logit(np.clip(pooled, 1e-5, 1 - 1e-5)))
    initial[-1] = float(logit(0.05 / 0.95))
    bounds = [(-20.0, 20.0)] * design.shape[1] + [(-14.0, 8.0)]
    fit = minimize(
        _beta_binomial_nll,
        initial,
        args=(design, successes, totals),
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": max_iterations, "ftol": 1e-9},
    )
    return fit.x, float(fit.fun), bool(fit.success), str(fit.message)


def _beta_binomial_mixed_nll(
    parameters: np.ndarray,
    design: np.ndarray,
    successes: np.ndarray,
    totals: np.ndarray,
    groups: np.ndarray,
    nodes: np.ndarray,
    weights: np.ndarray,
) -> float:
    beta = parameters[:-2]
    rho = 1e-7 + (0.95 - 1e-7) * expit(parameters[-2])
    sigma = np.exp(parameters[-1])
    concentration = (1.0 - rho) / rho
    failures = totals - successes
    log_choose = gammaln(totals + 1.0) - gammaln(successes + 1.0) - gammaln(failures + 1.0)
    total_log_likelihood = 0.0
    for group in pd.unique(groups):
        mask = groups == group
        quadrature_log_likelihood = []
        for node, weight in zip(nodes, weights, strict=False):
            random_intercept = np.sqrt(2.0) * sigma * node
            mu = np.clip(expit(design[mask] @ beta + random_intercept), 1e-7, 1.0 - 1e-7)
            alpha = mu * concentration
            beta_shape = (1.0 - mu) * concentration
            ll = log_choose[mask] + betaln(successes[mask] + alpha, failures[mask] + beta_shape) - betaln(alpha, beta_shape)
            quadrature_log_likelihood.append(np.log(weight) + float(ll.sum()))
        total_log_likelihood += float(logsumexp(quadrature_log_likelihood) - 0.5 * np.log(np.pi))
    return -total_log_likelihood if np.isfinite(total_log_likelihood) else 1e100


def _fit_mixed_model(
    design: np.ndarray,
    successes: np.ndarray,
    totals: np.ndarray,
    groups: np.ndarray,
    max_iterations: int,
    quadrature_points: int,
) -> tuple[np.ndarray, float, bool, str]:
    pooled = float((successes.sum() + 0.5) / (totals.sum() + 1.0))
    initial = np.zeros(design.shape[1] + 2, dtype=float)
    initial[0] = float(logit(np.clip(pooled, 1e-5, 1 - 1e-5)))
    initial[-2] = float(logit(0.05 / 0.95))
    initial[-1] = np.log(0.25)
    nodes, weights = np.polynomial.hermite.hermgauss(quadrature_points)
    bounds = [(-20.0, 20.0)] * design.shape[1] + [(-14.0, 8.0), (-8.0, 3.0)]
    fit = minimize(
        _beta_binomial_mixed_nll,
        initial,
        args=(design, successes, totals, groups, nodes, weights),
        method="L-BFGS-B", bounds=bounds,
        options={"maxiter": max_iterations, "ftol": 1e-8},
    )
    return fit.x, float(fit.fun), bool(fit.success), str(fit.message)


def _condition_means(values: pd.Series, metadata: pd.DataFrame, condition_column: str) -> dict[str, float]:
    aligned = pd.DataFrame({"value": pd.to_numeric(values, errors="coerce"), "condition": metadata[condition_column].astype(str)}, index=metadata.index)
    return aligned.groupby("condition", observed=True)["value"].mean().to_dict()


def fit_replicate_aware_models(
    psi: pd.DataFrame,
    ijc: pd.DataFrame,
    sjc: pd.DataFrame,
    samples: pd.DataFrame,
    assignments: pd.DataFrame,
    program_scores: pd.DataFrame,
    settings: ModuleStatisticsConfig,
    *,
    batch_column: str | None = None,
) -> ModuleStatisticsResult:
    if "sample_id" not in samples or settings.condition_column not in samples:
        raise ValueError(f"Sample metadata requires sample_id and {settings.condition_column}.")
    metadata = samples.set_index("sample_id")
    common_samples = psi.columns.intersection(metadata.index).intersection(ijc.columns).intersection(sjc.columns)
    metadata = metadata.loc[common_samples]
    if metadata[settings.condition_column].nunique() < 2:
        raise ValueError("Replicate-aware statistics require at least two conditions.")
    counts = metadata[settings.condition_column].value_counts()
    if (counts < settings.min_samples_per_condition).any():
        raise ValueError(
            "Each condition needs at least "
            f"{settings.min_samples_per_condition} samples for replicate-aware modeling: {counts.to_dict()}"
        )
    if batch_column and batch_column not in metadata:
        raise ValueError(f"Configured optional batch column is missing: {batch_column}")
    random_effect_column = settings.random_effect_column
    if random_effect_column and random_effect_column not in metadata:
        raise ValueError(f"Configured random-effect column is missing: {random_effect_column}")
    if random_effect_column and metadata[random_effect_column].value_counts().max() < 2:
        raise ValueError("Random-effect modeling requires at least one donor/block with repeated observations.")
    null_design, full_design, _, full_names = _design(metadata, settings.condition_column, batch_column)
    if np.linalg.matrix_rank(full_design) < full_design.shape[1]:
        raise ValueError(
            "Condition and optional batch design is rank deficient. Batch is confounded with condition; "
            "remove the batch column or provide samples spanning condition and batch levels."
        )
    assignment_events = (
        pd.Index(assignments["event_key"].dropna().astype(str).unique())
        if "event_key" in assignments
        else pd.Index([], dtype=str)
    )
    if assignment_events.empty:
        raise ValueError("Replicate-aware module statistics require non-empty event module assignments.")
    events = psi.index.intersection(ijc.index).intersection(sjc.index).intersection(assignment_events)
    total_reads = (ijc.loc[events, common_samples] + sjc.loc[events, common_samples]).median(axis=1)
    if settings.max_events is not None and len(events) > settings.max_events:
        events = total_reads.nlargest(settings.max_events).index

    rows: list[dict[str, object]] = []
    diagnostics: list[dict[str, object]] = []
    for event_key in events:
        successes = pd.to_numeric(ijc.loc[event_key, common_samples], errors="coerce").to_numpy(dtype=float)
        totals = successes + pd.to_numeric(sjc.loc[event_key, common_samples], errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(successes) & np.isfinite(totals) & (totals >= settings.min_total_reads)
        valid_conditions = metadata.iloc[np.flatnonzero(valid)][settings.condition_column].value_counts()
        if len(valid_conditions) < 2 or (valid_conditions < settings.min_samples_per_condition).any():
            diagnostics.append({"event_key": event_key, "status": "insufficient_read_supported_replicates", "message": str(valid_conditions.to_dict())})
            continue
        y = successes[valid]
        n = totals[valid]
        null = null_design[valid]
        full = full_design[valid]
        try:
            if random_effect_column:
                groups = metadata.iloc[np.flatnonzero(valid)][random_effect_column].astype(str).to_numpy()
                null_parameters, null_nll, null_ok, null_message = _fit_mixed_model(
                    null, y, n, groups, settings.optimizer_max_iterations, settings.quadrature_points
                )
                full_parameters, full_nll, full_ok, full_message = _fit_mixed_model(
                    full, y, n, groups, settings.optimizer_max_iterations, settings.quadrature_points
                )
                random_effect_sd = float(np.exp(full_parameters[-1]))
                coefficient_parameters = full_parameters[:-2]
                model_name = "beta-binomial random-intercept LRT: condition + optional batch versus optional batch"
            else:
                null_parameters, null_nll, null_ok, null_message = _fit_model(null, y, n, settings.optimizer_max_iterations)
                full_parameters, full_nll, full_ok, full_message = _fit_model(full, y, n, settings.optimizer_max_iterations)
                random_effect_sd = np.nan
                coefficient_parameters = full_parameters[:-1]
                model_name = "beta-binomial LRT: condition + optional batch versus optional batch"
            likelihood_ratio = max(0.0, 2.0 * (null_nll - full_nll))
            df = full.shape[1] - null.shape[1]
            pvalue = float(chi2.sf(likelihood_ratio, df)) if df > 0 else np.nan
            rho_parameter = full_parameters[-2] if random_effect_column else full_parameters[-1]
            rho = 1e-7 + (0.95 - 1e-7) * expit(rho_parameter)
            condition_means = _condition_means(psi.loc[event_key, common_samples], metadata, settings.condition_column)
            means = np.asarray(list(condition_means.values()), dtype=float)
            delta = float(np.nanmax(means) - np.nanmin(means)) if np.isfinite(means).any() else np.nan
            row: dict[str, object] = {
                "event_key": event_key,
                "n_samples": int(valid.sum()),
                "n_conditions": int(len(valid_conditions)),
                "median_total_junction_reads": float(np.median(n)),
                "beta_binomial_rho": float(rho),
                "likelihood_ratio": likelihood_ratio,
                "model_df": df,
                "condition_effect_pvalue": pvalue,
                "simple_max_delta_PSI": delta,
                "condition_mean_PSI": ";".join(f"{key}={value:.6g}" for key, value in condition_means.items()),
                "converged": bool(null_ok and full_ok),
                "random_effect_column": random_effect_column or "",
                "random_intercept_sd": random_effect_sd,
                "model": model_name,
            }
            for name, coefficient in zip(full_names, coefficient_parameters, strict=False):
                row[f"coefficient_{name}"] = float(coefficient)
            rows.append(row)
            diagnostics.append({
                "event_key": event_key, "status": "finished" if null_ok and full_ok else "optimizer_warning",
                "message": f"null={null_message}; full={full_message}", "null_nll": null_nll, "full_nll": full_nll,
            })
        except Exception as exc:
            diagnostics.append({"event_key": event_key, "status": "failed", "message": str(exc)})
    event_models = pd.DataFrame(rows)
    if not event_models.empty:
        module_column = "consensus_module" if "consensus_module" in assignments else "module"
        if {"event_key", module_column}.issubset(assignments):
            event_models = event_models.merge(
                assignments[["event_key", module_column]].drop_duplicates("event_key").rename(columns={module_column: "module"}),
                on="event_key",
                how="left",
                validate="one_to_one",
            )
        event_models["condition_effect_FDR"] = bh_adjust(event_models["condition_effect_pvalue"])
        event_models["replicate_supported_condition_effect"] = event_models["condition_effect_FDR"].le(settings.fdr_cutoff)
        event_models = event_models.sort_values(["condition_effect_FDR", "condition_effect_pvalue", "simple_max_delta_PSI"], ascending=[True, True, False], na_position="last")

    scores = program_scores.copy()
    scores.index = scores.index.astype(str)
    score_samples = scores.index.intersection(metadata.index)
    module_rows: list[dict[str, object]] = []
    variance_rows: list[dict[str, object]] = []
    condition_dummies = pd.get_dummies(metadata.loc[score_samples, settings.condition_column].astype(str), drop_first=True, dtype=float)
    batch_dummies = (
        pd.get_dummies(metadata.loc[score_samples, batch_column].astype(str), drop_first=True, dtype=float)
        if batch_column else pd.DataFrame(index=score_samples)
    )
    null_x = sm.add_constant(batch_dummies, has_constant="add")
    full_x = sm.add_constant(pd.concat([batch_dummies, condition_dummies], axis=1), has_constant="add")
    for module in scores.columns:
        y = pd.to_numeric(scores.loc[score_samples, module], errors="coerce")
        valid = y.notna()
        if valid.sum() <= full_x.shape[1] + 1:
            diagnostics.append({
                "event_key": f"module:{module}", "status": "insufficient_module_samples",
                "message": f"{int(valid.sum())} usable samples for {full_x.shape[1]} model coefficients",
            })
            continue
        try:
            if random_effect_column:
                groups = metadata.loc[score_samples].loc[valid, random_effect_column].astype(str)
                null_fit = sm.MixedLM(y.loc[valid], null_x.loc[valid], groups=groups).fit(reml=False, method="lbfgs", disp=False)
                full_fit = sm.MixedLM(y.loc[valid], full_x.loc[valid], groups=groups).fit(reml=False, method="lbfgs", disp=False)
                likelihood_ratio = max(0.0, 2.0 * (full_fit.llf - null_fit.llf))
                pvalue = float(chi2.sf(likelihood_ratio, max(1, len(condition_dummies.columns))))
                model_label = "linear mixed model program score with donor/block random intercept"
                total_sum_squares = float(np.square(y.loc[valid] - y.loc[valid].mean()).sum())
                null_r2 = max(0.0, 1.0 - float(np.square(y.loc[valid] - null_fit.fittedvalues).sum()) / total_sum_squares) if total_sum_squares > 0 else 0.0
                full_r2 = max(0.0, 1.0 - float(np.square(y.loc[valid] - full_fit.fittedvalues).sum()) / total_sum_squares) if total_sum_squares > 0 else 0.0
            else:
                null_fit = sm.OLS(y.loc[valid], null_x.loc[valid]).fit(cov_type="HC3")
                full_fit = sm.OLS(y.loc[valid], full_x.loc[valid]).fit(cov_type="HC3")
                condition_columns = condition_dummies.columns.tolist()
                restriction = np.zeros((len(condition_columns), len(full_fit.params)))
                for row_index, name in enumerate(condition_columns):
                    restriction[row_index, full_fit.params.index.get_loc(name)] = 1.0
                test = full_fit.f_test(restriction) if condition_columns else None
                pvalue = float(np.asarray(test.pvalue).item()) if test is not None else np.nan
                model_label = "OLS program score with HC3 robust covariance"
                null_r2 = float(null_fit.rsquared)
                full_r2 = float(full_fit.rsquared)
        except Exception as exc:
            diagnostics.append({
                "event_key": f"module:{module}", "status": "module_model_failed",
                "message": str(exc),
            })
            continue
        condition_r2 = max(0.0, full_r2 - null_r2)
        module_rows.append({
            "module": module, "n_samples": int(valid.sum()), "condition_effect_pvalue": pvalue,
            "condition_partial_R2": condition_r2, "full_model_R2": full_r2,
            "batch_included": bool(batch_column), "random_effect_column": random_effect_column or "",
            "model": model_label,
        })
        variance_rows.append({
            "module": module, "condition_variance_fraction": condition_r2,
            "batch_or_baseline_variance_fraction": max(0.0, null_r2),
            "unexplained_residual_fraction": max(0.0, 1.0 - full_r2),
        })
    module_models = pd.DataFrame(module_rows)
    if not module_models.empty:
        module_models["condition_effect_FDR"] = bh_adjust(module_models["condition_effect_pvalue"])
        module_models["condition_driven"] = module_models["condition_effect_FDR"].le(settings.fdr_cutoff)
        module_models = module_models.sort_values(["condition_effect_FDR", "condition_partial_R2"], ascending=[True, False], na_position="last")
    variance_partition = pd.DataFrame(variance_rows)
    if not variance_partition.empty and not module_models.empty:
        variance_partition = variance_partition.merge(
            module_models[["module", "condition_effect_FDR", "condition_driven"]],
            on="module",
            how="left",
            validate="one_to_one",
        )
        variance_partition["module_class"] = np.where(
            variance_partition["condition_driven"].fillna(False),
            "condition_driven",
            "unresolved_or_non_condition_driven",
        )

    benchmark_rows: list[dict[str, object]] = []
    if not event_models.empty:
        model_score = -np.log10(event_models["condition_effect_FDR"].clip(lower=1e-300))
        delta_score = event_models["simple_max_delta_PSI"]
        if model_score.nunique(dropna=True) >= 2 and delta_score.nunique(dropna=True) >= 2:
            correlation, correlation_p = spearmanr(model_score, delta_score, nan_policy="omit")
        else:
            correlation, correlation_p = np.nan, np.nan
        for top_n in sorted(set([min(50, len(event_models)), min(100, len(event_models)), min(500, len(event_models))])):
            if top_n <= 0:
                continue
            model_top = set(event_models.nsmallest(top_n, "condition_effect_FDR")["event_key"])
            delta_top = set(event_models.nlargest(top_n, "simple_max_delta_PSI")["event_key"])
            union = model_top | delta_top
            benchmark_rows.append({
                "top_n": top_n, "overlap_events": len(model_top & delta_top),
                "jaccard": len(model_top & delta_top) / len(union) if union else np.nan,
                "spearman_model_vs_delta": float(correlation), "spearman_pvalue": float(correlation_p),
                "model_only_events": len(model_top - delta_top), "delta_only_events": len(delta_top - model_top),
                "interpretation": "beta-binomial uses replicate read-count uncertainty; delta PSI is an effect-size-only baseline",
            })
    return ModuleStatisticsResult(
        event_models=event_models.reset_index(drop=True),
        module_models=module_models.reset_index(drop=True),
        variance_partition=variance_partition.reset_index(drop=True),
        simple_delta_benchmark=pd.DataFrame(benchmark_rows),
        model_diagnostics=pd.DataFrame(diagnostics),
    )


def render_module_statistics_plots(result: ModuleStatisticsResult, output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    if not result.event_models.empty:
        frame = result.event_models.copy()
        frame["-log10 FDR"] = -np.log10(frame["condition_effect_FDR"].clip(lower=1e-300))
        fig, ax = plt.subplots(figsize=(9, 6))
        sns.scatterplot(data=frame, x="simple_max_delta_PSI", y="-log10 FDR", hue="beta_binomial_rho", palette="viridis", s=32, ax=ax)
        ax.set_title("Replicate-aware beta-binomial evidence versus simple delta PSI")
        ax.set_xlabel("Simple maximum condition mean delta PSI")
        ax.set_ylabel("-log10 beta-binomial condition-effect FDR")
        fig.tight_layout()
        path = output_dir / "beta_binomial_vs_simple_delta_dpsi.png"
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
    if not result.variance_partition.empty:
        frame = result.variance_partition.set_index("module")[["condition_variance_fraction", "batch_or_baseline_variance_fraction", "unexplained_residual_fraction"]]
        fig, ax = plt.subplots(figsize=(max(9, len(frame) * 0.45), 6))
        frame.plot(kind="bar", stacked=True, ax=ax, color=["#0f766e", "#94a3b8", "#f59e0b"])
        ax.set_title("CoSplice module program-score variance partition")
        ax.set_ylabel("Variance fraction")
        ax.set_xlabel("Module")
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
        fig.tight_layout()
        path = output_dir / "module_variance_partition.png"
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
    return paths
