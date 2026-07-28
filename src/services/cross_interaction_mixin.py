from __future__ import annotations

import math
import re

import numpy as np
import pandas as pd

from src.analysis.common import EVENT_TYPES, event_uid_from_row
from src.analysis.stable_as import apply_stable_as_filter
from src.domain.models.project import ComparisonPairDefinition


class CrossInteractionMixin:
    """Replicate-aware cross-comparison AS interaction calculations."""

    @staticmethod
    def _empty_cross_as_interaction_tables() -> dict[str, pd.DataFrame]:
        base_columns = [
            "gene",
            "event_type",
            "event_id",
            "event_matching_key",
            "comparison_A",
            "comparison_B",
            "experiment_A_display_name",
            "experiment_B_display_name",
            "dPSI_A",
            "dPSI_B",
            "FDR_A",
            "FDR_B",
            "significant_A",
            "significant_B",
            "direction_A",
            "direction_B",
            "delta_dPSI",
            "abs_delta_dPSI",
            "same_direction",
            "opposite_direction",
            "same_direction_strength_delta",
            "interaction_log_odds",
            "interaction_se",
            "interaction_z",
            "interaction_pvalue",
            "interaction_FDR",
            "interaction_model",
            "interaction_model_status",
            "interaction_overdispersion_scale",
            "interaction_model_n_observations",
            "interaction_model_df_resid",
            "interaction_model_pearson_chi2",
            "interaction_model_deviance",
            "interaction_model_total_reads",
            "quasi_binomial_dispersion",
            "quasi_binomial_overdispersed",
            "beta_binomial_rho_estimate",
            "beta_binomial_interaction_log_odds",
            "beta_binomial_likelihood_ratio",
            "beta_binomial_pvalue",
            "beta_binomial_FDR",
            "beta_binomial_rho",
            "beta_binomial_full_loglik",
            "beta_binomial_null_loglik",
            "beta_binomial_model_status",
            "beta_binomial_warning",
            "model_valid",
            "model_diagnostic_class",
            "model_warning",
            "eb_prior_variance",
            "shrinkage_weight",
            "shrunken_interaction_log_odds",
            "shrunken_interaction_z",
            "shrunken_interaction_pvalue",
            "shrunken_interaction_FDR",
            "interaction_priority_score",
            "simple_delta_score",
            "simple_delta_priority_score",
            "delta_dPSI_rank",
            "interaction_model_rank",
            "delta_vs_model_rank_shift",
            "simple_delta_top20",
            "model_top20",
            "benchmark_class",
            "benchmark_score_delta",
            "interaction_class",
            "read_support_class",
            "replicate_consistency_score",
            "mean_total_reads_A",
            "mean_total_reads_B",
            "replicate_count_A",
            "replicate_count_B",
            "exclusive_side",
            "counterpart_event_found",
            "significant_side_mean_total_reads",
            "non_significant_side_mean_total_reads",
            "significant_side_replicate_count",
            "non_significant_side_replicate_count",
            "non_significant_side_dPSI",
            "non_significant_side_FDR",
            "non_significant_side_read_support",
            "event_class",
            "unified_result_type",
            "unified_display_class",
            "support_class",
            "abs_dPSI_significance_cutoff_used",
            "large_delta_dPSI_cutoff_used",
            "min_total_reads_used",
        ]
        empty = pd.DataFrame(columns=base_columns)
        return {
            "unified_event_table": empty.copy(),
            "shared_interaction": empty.copy(),
            "opposite_direction": empty.copy(),
            "large_delta": empty.copy(),
            "exclusive": empty.copy(),
            "benchmark": empty.copy(),
            "diagnostics": pd.DataFrame(
                columns=[
                    "pair_id",
                    "pair_name",
                    "comparison_A",
                    "comparison_B",
                    "n_events_A",
                    "n_events_B",
                    "n_shared_significant",
                    "n_opposite_direction",
                    "n_large_delta",
                    "n_exclusive_supported",
                    "n_exclusive_low_power",
                    "n_exclusive_counterpart_not_found",
                    "n_exclusive_not_shared_significant",
                    "n_model_converged",
                    "n_model_fallback",
                    "n_model_failed",
                    "median_quasi_binomial_dispersion",
                    "max_quasi_binomial_dispersion",
                    "median_beta_binomial_rho_estimate",
                    "median_beta_binomial_rho",
                    "n_beta_binomial_converged",
                    "n_beta_binomial_failed",
                    "n_beta_binomial_FDR_0_05",
                    "n_overdispersed_events",
                    "delta_top20_model_top20_overlap",
                    "delta_top20_model_top20_overlap_fraction",
                    "delta_vs_model_spearman_rank_correlation",
                    "n_delta_only_top20",
                    "n_model_only_top20",
                    "abs_dPSI_significance_cutoff_used",
                    "large_delta_dPSI_cutoff_used",
                    "min_total_reads_used",
                ]
            ),
        }

    @staticmethod
    def _empty_cross_as_interaction_batch_tables() -> dict[str, pd.DataFrame]:
        empty_tables = CrossInteractionMixin._empty_cross_as_interaction_tables()
        return {
            "batch_summary": pd.DataFrame(
                columns=[
                    "batch_pair_order",
                    "pair_id",
                    "pair_name",
                    "comparison_A",
                    "comparison_B",
                    "n_shared_significant",
                    "n_high_confidence_interaction",
                    "n_moderate_interaction",
                    "n_opposite_direction",
                    "n_same_direction_large_delta",
                    "n_exclusive",
                    "n_exclusive_A",
                    "n_exclusive_B",
                    "n_exclusive_supported",
                    "n_exclusive_A_supported",
                    "n_exclusive_B_supported",
                    "n_exclusive_low_power",
                    "n_exclusive_counterpart_not_found",
                    "top_gene",
                    "top_event_type",
                    "top_abs_delta_dPSI",
                    "top_interaction_priority_score",
                    "delta_top20_model_top20_overlap",
                    "delta_vs_model_spearman_rank_correlation",
                    "abs_dPSI_significance_cutoff_used",
                    "large_delta_dPSI_cutoff_used",
                    "min_total_reads_used",
                ]
            ),
            "batch_top_interactions": empty_tables["shared_interaction"].copy(),
            "batch_exclusive": empty_tables["exclusive"].copy(),
            "batch_benchmark": empty_tables["benchmark"].copy(),
            "batch_diagnostics": empty_tables["diagnostics"].copy(),
        }

    @staticmethod
    def _parse_float_list(value: object) -> list[float]:
        if value is None:
            return []
        try:
            if pd.isna(value):
                return []
        except (TypeError, ValueError):
            pass
        values: list[float] = []
        for token in str(value).replace(";", ",").split(","):
            text = token.strip()
            if not text or text.upper() in {"NA", "NAN", "NONE"}:
                continue
            try:
                values.append(float(text))
            except ValueError:
                continue
        return values

    @staticmethod
    def _paired_read_totals(ijc_values: list[float], sjc_values: list[float]) -> list[float]:
        n_items = max(len(ijc_values), len(sjc_values))
        totals: list[float] = []
        for idx in range(n_items):
            ijc = ijc_values[idx] if idx < len(ijc_values) and np.isfinite(ijc_values[idx]) else 0.0
            sjc = sjc_values[idx] if idx < len(sjc_values) and np.isfinite(sjc_values[idx]) else 0.0
            totals.append(float(ijc) + float(sjc))
        return totals

    def _apply_stable_as_columns(self, frame: pd.DataFrame, *, dpsi_column: str = "dPSI", reverse_splicing: bool = False) -> pd.DataFrame:
        working = frame.copy()
        fdr = pd.to_numeric(working.get("FDR"), errors="coerce")
        standardized_dpsi = pd.to_numeric(working.get(dpsi_column), errors="coerce")
        return apply_stable_as_filter(
            working,
            fdr=fdr,
            standardized_dpsi=standardized_dpsi,
            reverse_splicing=reverse_splicing,
            thresholds=self.current_thresholds,
        )

    @staticmethod
    def _safe_logit(value: float) -> float:
        if not np.isfinite(value):
            return float("nan")
        clipped = min(max(float(value), 1e-4), 1.0 - 1e-4)
        return float(np.log(clipped / (1.0 - clipped)))

    @staticmethod
    def _normal_two_sided_pvalue(z_value: float) -> float:
        if not np.isfinite(z_value):
            return float("nan")
        return float(math.erfc(abs(float(z_value)) / math.sqrt(2.0)))

    @staticmethod
    def _bh_adjust(pvalues: pd.Series) -> pd.Series:
        numeric = pd.to_numeric(pvalues, errors="coerce")
        adjusted = pd.Series(np.nan, index=numeric.index, dtype="float64")
        valid = numeric.dropna()
        if valid.empty:
            return adjusted
        ordered = valid.sort_values()
        m = float(len(ordered))
        previous = 1.0
        for rank in range(len(ordered), 0, -1):
            item_index = ordered.index[rank - 1]
            pvalue = float(ordered.iloc[rank - 1])
            estimate = min(previous, pvalue * m / float(rank))
            adjusted.loc[item_index] = min(max(estimate, 0.0), 1.0)
            previous = estimate
        return adjusted

    def _interaction_group_stats(self, psi_values: list[float], ijc_values: list[float], sjc_values: list[float]) -> dict[str, float]:
        count_psi: list[float] = []
        totals: list[float] = []
        max_len = max(len(ijc_values), len(sjc_values))
        for idx in range(max_len):
            ijc = ijc_values[idx] if idx < len(ijc_values) and np.isfinite(ijc_values[idx]) else 0.0
            sjc = sjc_values[idx] if idx < len(sjc_values) and np.isfinite(sjc_values[idx]) else 0.0
            total = ijc + sjc
            totals.append(total)
            count_psi.append((ijc + 0.5) / (total + 1.0) if total > 0 else float("nan"))
        usable_psi = [float(value) for value in psi_values if np.isfinite(value)]
        if not usable_psi:
            usable_psi = [float(value) for value in count_psi if np.isfinite(value)]
        logits = [self._safe_logit(value) for value in usable_psi if np.isfinite(value)]
        n = len(logits)
        mean_logit = float(np.nanmean(logits)) if logits else float("nan")
        var_logit = float(np.nanvar(logits, ddof=1)) if len(logits) > 1 else float("nan")
        mean_psi = float(np.nanmean(usable_psi)) if usable_psi else float("nan")
        finite_totals = [float(value) for value in totals if np.isfinite(value)]
        mean_total = float(np.nanmean(finite_totals)) if finite_totals else float("nan")
        total_sum = float(np.nansum(finite_totals)) if finite_totals else 0.0
        return {
            "mean_psi": mean_psi,
            "mean_logit": mean_logit,
            "var_logit": var_logit,
            "n": float(n),
            "mean_total_reads": mean_total,
            "total_reads_sum": total_sum,
        }

    def _replicate_direction_consistency(self, psi_group1: list[float], psi_group2: list[float], dpsi_value: float) -> float:
        return pairwise_direction_consistency(psi_group1, psi_group2, dpsi_value)

    def _fit_binomial_glm_interaction(self, row: pd.Series) -> dict[str, object]:
        """Fit IJC/SJC ~ comparison + group + comparison:group.

        interaction_log_odds is reported as comparison A group effect minus
        comparison B group effect so its sign matches delta_dPSI = dPSI_A -
        dPSI_B. The GLM coefficient itself is B-minus-A, so we flip it.
        """

        def as_list(key: str) -> list[float]:
            value = row.get(key, [])
            if isinstance(value, list):
                return [float(item) for item in value if np.isfinite(item)]
            return self._parse_float_list(value)

        records: list[dict[str, float]] = []
        specs = [
            ("A", 0.0, 1.0, as_list("ijc_group1_values_A"), as_list("sjc_group1_values_A")),
            ("A", 0.0, 0.0, as_list("ijc_group2_values_A"), as_list("sjc_group2_values_A")),
            ("B", 1.0, 1.0, as_list("ijc_group1_values_B"), as_list("sjc_group1_values_B")),
            ("B", 1.0, 0.0, as_list("ijc_group2_values_B"), as_list("sjc_group2_values_B")),
        ]
        for _label, comparison_b, group1, ijc_values, sjc_values in specs:
            for ijc, sjc in zip(ijc_values, sjc_values, strict=False):
                if not np.isfinite(ijc) or not np.isfinite(sjc):
                    continue
                total = float(ijc) + float(sjc)
                if total <= 0:
                    continue
                records.append(
                    {
                        "ijc": max(float(ijc), 0.0),
                        "sjc": max(float(sjc), 0.0),
                        "comparison_b": comparison_b,
                        "group1": group1,
                        "interaction": comparison_b * group1,
                    }
                )
        if len(records) < 4:
            return {
                "interaction_log_odds": float("nan"),
                "interaction_se": float("nan"),
                "interaction_z": float("nan"),
                "interaction_pvalue": float("nan"),
                "interaction_model": "binomial_glm",
                "interaction_model_status": "insufficient_replicate_count",
                "interaction_overdispersion_scale": float("nan"),
                "interaction_model_n_observations": len(records),
                "interaction_model_df_resid": float("nan"),
                "interaction_model_pearson_chi2": float("nan"),
                "interaction_model_deviance": float("nan"),
                "interaction_model_total_reads": float("nan"),
                "quasi_binomial_dispersion": float("nan"),
                "quasi_binomial_overdispersed": False,
                "beta_binomial_rho_estimate": float("nan"),
                "model_valid": False,
                "model_diagnostic_class": "insufficient_replicates",
                "model_warning": "Fewer than four replicate group observations with IJC/SJC support.",
            }
        try:
            import statsmodels.api as sm

            frame = pd.DataFrame(records)
            endog = frame[["ijc", "sjc"]].to_numpy(dtype=float)
            exog = frame[["comparison_b", "group1", "interaction"]].to_numpy(dtype=float)
            exog = sm.add_constant(exog, has_constant="add")
            model = sm.GLM(endog, exog, family=sm.families.Binomial())
            result = model.fit(maxiter=100, disp=0)
            raw_coef = float(result.params[3])
            raw_se = float(result.bse[3])
            scale = 1.0
            pearson_chi2 = float(getattr(result, "pearson_chi2", float("nan")))
            df_resid = float(getattr(result, "df_resid", float("nan")))
            deviance = float(getattr(result, "deviance", float("nan")))
            try:
                if np.isfinite(df_resid) and df_resid > 0 and np.isfinite(pearson_chi2):
                    scale = max(1.0, pearson_chi2 / df_resid)
            except Exception:
                scale = 1.0
            adjusted_se = raw_se * math.sqrt(scale) if np.isfinite(raw_se) and raw_se > 0 else float("nan")
            interaction_log_odds = -raw_coef
            interaction_z = (-raw_coef / adjusted_se) if np.isfinite(adjusted_se) and adjusted_se > 0 else float("nan")
            interaction_pvalue = self._normal_two_sided_pvalue(interaction_z)
            total_reads = float((frame["ijc"] + frame["sjc"]).sum())
            mean_total = float((frame["ijc"] + frame["sjc"]).mean()) if not frame.empty else float("nan")
            rho = (scale - 1.0) / max(mean_total - 1.0, 1.0) if np.isfinite(mean_total) and mean_total > 1 else 0.0
            rho = float(np.clip(rho, 0.0, 0.99)) if np.isfinite(rho) else float("nan")
            overdispersed = bool(np.isfinite(scale) and scale >= 1.5)
            model_valid = bool(np.isfinite(interaction_pvalue) and np.isfinite(adjusted_se) and adjusted_se > 0)
            return {
                "interaction_log_odds": interaction_log_odds,
                "interaction_se": adjusted_se,
                "interaction_z": interaction_z,
                "interaction_pvalue": interaction_pvalue,
                "interaction_model": "binomial_glm_ijc_sjc_with_pearson_overdispersion",
                "interaction_model_status": "converged" if bool(getattr(result, "converged", True)) else "not_converged",
                "interaction_overdispersion_scale": scale,
                "interaction_model_n_observations": len(records),
                "interaction_model_df_resid": df_resid,
                "interaction_model_pearson_chi2": pearson_chi2,
                "interaction_model_deviance": deviance,
                "interaction_model_total_reads": total_reads,
                "quasi_binomial_dispersion": scale,
                "quasi_binomial_overdispersed": overdispersed,
                "beta_binomial_rho_estimate": rho,
                "model_valid": model_valid,
                "model_diagnostic_class": "quasi_binomial_overdispersed" if overdispersed else "quasi_binomial_ok",
                "model_warning": "" if model_valid else "Model converged but interaction coefficient or SE was not finite.",
            }
        except Exception as exc:
            return {
                "interaction_log_odds": float("nan"),
                "interaction_se": float("nan"),
                "interaction_z": float("nan"),
                "interaction_pvalue": float("nan"),
                "interaction_model": "binomial_glm_ijc_sjc_with_pearson_overdispersion",
                "interaction_model_status": f"failed: {exc}",
                "interaction_overdispersion_scale": float("nan"),
                "interaction_model_n_observations": len(records),
                "interaction_model_df_resid": float("nan"),
                "interaction_model_pearson_chi2": float("nan"),
                "interaction_model_deviance": float("nan"),
                "interaction_model_total_reads": float("nan"),
                "quasi_binomial_dispersion": float("nan"),
                "quasi_binomial_overdispersed": False,
                "beta_binomial_rho_estimate": float("nan"),
                "model_valid": False,
                "model_diagnostic_class": "model_failed",
                "model_warning": str(exc),
            }

    def _fit_beta_binomial_interaction(self, row: pd.Series) -> dict[str, object]:
        """Fit a beta-binomial interaction model on replicate IJC/SJC counts.

        The full model is inclusion ~ comparison + group + comparison:group.
        The null model removes the interaction term. The likelihood-ratio
        p-value tests whether the AS effect differs between comparison A and B.
        """

        def as_list(key: str) -> list[float]:
            value = row.get(key, [])
            if isinstance(value, list):
                return [float(item) for item in value if np.isfinite(item)]
            return self._parse_float_list(value)

        records: list[tuple[float, float, float, float, float]] = []
        specs = [
            (0.0, 1.0, as_list("ijc_group1_values_A"), as_list("sjc_group1_values_A")),
            (0.0, 0.0, as_list("ijc_group2_values_A"), as_list("sjc_group2_values_A")),
            (1.0, 1.0, as_list("ijc_group1_values_B"), as_list("sjc_group1_values_B")),
            (1.0, 0.0, as_list("ijc_group2_values_B"), as_list("sjc_group2_values_B")),
        ]
        for comparison_b, group1, ijc_values, sjc_values in specs:
            for ijc, sjc in zip(ijc_values, sjc_values, strict=False):
                if not np.isfinite(ijc) or not np.isfinite(sjc):
                    continue
                success = max(float(ijc), 0.0)
                failure = max(float(sjc), 0.0)
                total = success + failure
                if total <= 0:
                    continue
                records.append((success, total, comparison_b, group1, comparison_b * group1))
        if len(records) < 4:
            return {
                "beta_binomial_interaction_log_odds": float("nan"),
                "beta_binomial_likelihood_ratio": float("nan"),
                "beta_binomial_pvalue": float("nan"),
                "beta_binomial_rho": float("nan"),
                "beta_binomial_full_loglik": float("nan"),
                "beta_binomial_null_loglik": float("nan"),
                "beta_binomial_model_status": "insufficient_replicate_count",
                "beta_binomial_warning": "Fewer than four replicate group observations with IJC/SJC support.",
            }
        try:
            from scipy.optimize import minimize
            from scipy.special import betaln, expit, gammaln
            from scipy.stats import chi2

            frame = pd.DataFrame(records, columns=["success", "total", "comparison_b", "group1", "interaction"])
            y = frame["success"].to_numpy(dtype=float)
            n = frame["total"].to_numpy(dtype=float)
            x_full = np.column_stack(
                [
                    np.ones(len(frame), dtype=float),
                    frame["comparison_b"].to_numpy(dtype=float),
                    frame["group1"].to_numpy(dtype=float),
                    frame["interaction"].to_numpy(dtype=float),
                ]
            )
            x_null = x_full[:, :3]

            def beta_binomial_loglik(params: np.ndarray, design: np.ndarray) -> float:
                beta = params[:-1]
                rho = float(expit(params[-1]) * 0.95)
                rho = min(max(rho, 1e-6), 0.95)
                mu = np.clip(expit(design @ beta), 1e-6, 1.0 - 1e-6)
                concentration = (1.0 - rho) / rho
                alpha = np.clip(mu * concentration, 1e-6, 1e8)
                beta_param = np.clip((1.0 - mu) * concentration, 1e-6, 1e8)
                log_comb = gammaln(n + 1.0) - gammaln(y + 1.0) - gammaln(n - y + 1.0)
                value = log_comb + betaln(y + alpha, n - y + beta_param) - betaln(alpha, beta_param)
                return float(np.sum(value))

            def negative_loglik(params: np.ndarray, design: np.ndarray) -> float:
                value = beta_binomial_loglik(params, design)
                if not np.isfinite(value):
                    return 1e100
                return -value

            observed_rate = np.clip((y.sum() + 0.5) / (n.sum() + 1.0), 1e-4, 1.0 - 1e-4)
            intercept = float(np.log(observed_rate / (1.0 - observed_rate)))
            rho_start = float(pd.to_numeric(row.get("beta_binomial_rho_estimate", np.nan), errors="coerce"))
            if not np.isfinite(rho_start) or rho_start <= 0:
                rho_start = 0.02
            rho_start = min(max(rho_start, 1e-5), 0.80)
            rho_param_start = float(np.log((rho_start / 0.95) / max(1.0 - (rho_start / 0.95), 1e-6)))
            full_start = np.array([intercept, 0.0, 0.0, 0.0, rho_param_start], dtype=float)
            null_start = np.array([intercept, 0.0, 0.0, rho_param_start], dtype=float)
            bounds_full = [(-20.0, 20.0), (-20.0, 20.0), (-20.0, 20.0), (-20.0, 20.0), (-12.0, 3.0)]
            bounds_null = [(-20.0, 20.0), (-20.0, 20.0), (-20.0, 20.0), (-12.0, 3.0)]
            full_fit = minimize(
                negative_loglik,
                full_start,
                args=(x_full,),
                method="L-BFGS-B",
                bounds=bounds_full,
                options={"maxiter": 120, "ftol": 1e-7},
            )
            null_fit = minimize(
                negative_loglik,
                null_start,
                args=(x_null,),
                method="L-BFGS-B",
                bounds=bounds_null,
                options={"maxiter": 120, "ftol": 1e-7},
            )
            full_ll = beta_binomial_loglik(full_fit.x, x_full)
            null_ll = beta_binomial_loglik(null_fit.x, x_null)
            lr_stat = max(0.0, 2.0 * (full_ll - null_ll))
            pvalue = float(chi2.sf(lr_stat, df=1))
            rho = float(expit(full_fit.x[-1]) * 0.95)
            warning_parts = []
            if not full_fit.success:
                warning_parts.append(f"full model: {full_fit.message}")
            if not null_fit.success:
                warning_parts.append(f"null model: {null_fit.message}")
            return {
                "beta_binomial_interaction_log_odds": -float(full_fit.x[3]),
                "beta_binomial_likelihood_ratio": lr_stat,
                "beta_binomial_pvalue": pvalue,
                "beta_binomial_rho": rho,
                "beta_binomial_full_loglik": full_ll,
                "beta_binomial_null_loglik": null_ll,
                "beta_binomial_model_status": "converged" if full_fit.success and null_fit.success else "fit_warning",
                "beta_binomial_warning": "; ".join(warning_parts),
            }
        except Exception as exc:
            return {
                "beta_binomial_interaction_log_odds": float("nan"),
                "beta_binomial_likelihood_ratio": float("nan"),
                "beta_binomial_pvalue": float("nan"),
                "beta_binomial_rho": float("nan"),
                "beta_binomial_full_loglik": float("nan"),
                "beta_binomial_null_loglik": float("nan"),
                "beta_binomial_model_status": "failed",
                "beta_binomial_warning": str(exc),
            }

    def _cross_as_interaction_source_frame(
        self,
        comparison_id: str | None,
        abs_dpsi_cutoff: float,
        min_total_reads: int,
        *,
        significant_only: bool = True,
        candidate_keys: set[str] | None = None,
    ) -> pd.DataFrame:
        if not comparison_id or self.current_project is None:
            return pd.DataFrame()
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None:
            return pd.DataFrame()
        mode = self.current_project.selected_rmats_mode or "JC"
        annotation_lookup = self._annotation_gene_symbol_lookup()
        left_label, right_label = self._comparison_splicing_direction_labels(comparison)
        rows: list[dict[str, object]] = []
        for event_type in EVENT_TYPES:
            source_path = self._resolve_rmats_event_table_path(comparison, event_type, mode)
            if source_path is None or not source_path.exists():
                continue
            try:
                frame = pd.read_csv(source_path, sep="\t", low_memory=False)
            except Exception:
                continue
            if frame.empty:
                continue
            working = frame.copy()
            fdr_values = pd.to_numeric(working.get("FDR"), errors="coerce")
            dpsi_values = pd.to_numeric(working.get("IncLevelDifference"), errors="coerce")
            if comparison.reverse_splicing:
                dpsi_values = -dpsi_values
            working["FDR"] = fdr_values
            working["dPSI"] = dpsi_values
            working = self._apply_stable_as_columns(working, dpsi_column="dPSI", reverse_splicing=bool(comparison.reverse_splicing))
            fdr_values = pd.to_numeric(working.get("FDR"), errors="coerce")
            dpsi_values = pd.to_numeric(working.get("dPSI"), errors="coerce")
            normalized_gene_ids = self._normalized_gene_id_series(working)
            if "geneSymbol" in working.columns:
                gene_symbols = working["geneSymbol"].astype("object")
            elif "gene_symbol" in working.columns:
                gene_symbols = working["gene_symbol"].astype("object")
            else:
                gene_symbols = pd.Series(pd.NA, index=working.index, dtype="object")
            missing_symbols = gene_symbols.isna() | gene_symbols.astype(str).str.strip().eq("")
            if missing_symbols.any() and annotation_lookup:
                gene_symbols.loc[missing_symbols] = normalized_gene_ids.loc[missing_symbols].map(annotation_lookup)
            gene_symbols = gene_symbols.combine_first(normalized_gene_ids)

            event_ids = working.apply(lambda row: event_uid_from_row(row, event_type), axis=1).astype(str)
            coordinates = working.apply(lambda row: self._rmats_event_coordinates(row, event_type), axis=1).astype(str)
            matching_ids = event_ids.where(event_ids.str.strip().ne(""), coordinates)
            key_genes = normalized_gene_ids.fillna("").astype(str)
            key_genes = key_genes.where(key_genes.str.strip().ne(""), gene_symbols.fillna("").astype(str))
            event_keys = (
                key_genes
                + "|"
                + event_type
                + "|"
                + matching_ids.fillna("").astype(str)
            )
            significant_mask = working["stable_as_event"].fillna(False).map(self._safe_bool) & dpsi_values.abs().ge(abs_dpsi_cutoff)
            if candidate_keys is not None:
                keep_mask = event_keys.isin(candidate_keys)
            elif significant_only:
                keep_mask = significant_mask
            else:
                keep_mask = pd.Series(True, index=working.index)
            working = working.loc[keep_mask].copy()
            if working.empty:
                continue
            fdr_values = fdr_values.loc[working.index]
            dpsi_values = dpsi_values.loc[working.index]
            gene_symbols = gene_symbols.loc[working.index]
            normalized_gene_ids = normalized_gene_ids.loc[working.index]
            event_ids = event_ids.loc[working.index]
            coordinates = coordinates.loc[working.index]
            event_keys = event_keys.loc[working.index]
            significant_mask = significant_mask.loc[working.index]
            stable_as_event = working["stable_as_event"].loc[working.index]
            rmats_sig_before_stability = working["rmats_sig_before_stability"].loc[working.index]

            for idx, row in working.iterrows():
                gene_id = str(normalized_gene_ids.loc[idx] or "").strip()
                gene_symbol = str(gene_symbols.loc[idx] or gene_id or "").strip()
                event_id = event_ids.loc[idx]
                coordinate_text = coordinates.loc[idx]
                event_matching_key = event_keys.loc[idx]
                fdr = fdr_values.loc[idx]
                dpsi = dpsi_values.loc[idx]

                psi_1 = self._parse_float_list(row.get("IncLevel1"))
                psi_2 = self._parse_float_list(row.get("IncLevel2"))
                ijc_1 = self._parse_float_list(row.get("IJC_SAMPLE_1"))
                ijc_2 = self._parse_float_list(row.get("IJC_SAMPLE_2"))
                sjc_1 = self._parse_float_list(row.get("SJC_SAMPLE_1"))
                sjc_2 = self._parse_float_list(row.get("SJC_SAMPLE_2"))
                if comparison.reverse_splicing:
                    psi_1, psi_2 = psi_2, psi_1
                    ijc_1, ijc_2 = ijc_2, ijc_1
                    sjc_1, sjc_2 = sjc_2, sjc_1

                stats_1 = self._interaction_group_stats(psi_1, ijc_1, sjc_1)
                stats_2 = self._interaction_group_stats(psi_2, ijc_2, sjc_2)
                effect_logit = stats_1["mean_logit"] - stats_2["mean_logit"]
                var_terms = []
                if stats_1["n"] > 1 and np.isfinite(stats_1["var_logit"]):
                    var_terms.append(stats_1["var_logit"] / max(stats_1["n"], 1.0))
                if stats_2["n"] > 1 and np.isfinite(stats_2["var_logit"]):
                    var_terms.append(stats_2["var_logit"] / max(stats_2["n"], 1.0))
                effect_se = float(np.sqrt(np.nansum(var_terms))) if var_terms else float("nan")
                if not np.isfinite(effect_se) or effect_se <= 0:
                    effect_se = float(
                        np.sqrt(
                            1.0 / max(stats_1["total_reads_sum"], 1.0)
                            + 1.0 / max(stats_2["total_reads_sum"], 1.0)
                        )
                    )
                mean_total_reads = float(np.nanmean([stats_1["mean_total_reads"], stats_2["mean_total_reads"]]))
                replicate_count = int(max(stats_1["n"], stats_2["n"]))
                consistency = self._replicate_direction_consistency(psi_1, psi_2, float(dpsi) if pd.notna(dpsi) else float("nan"))
                significant = bool(significant_mask.loc[idx])
                support_ok = bool(np.isfinite(mean_total_reads) and mean_total_reads >= min_total_reads and replicate_count >= 2)
                rows.append(
                    {
                        "comparison_id": comparison.comparison_id,
                        "comparison_display_name": comparison.display_resolved_name,
                        "gene": gene_symbol or gene_id,
                        "gene_symbol": gene_symbol,
                        "gene_id": gene_id,
                        "event_type": event_type,
                        "event_id": event_id,
                        "event_matching_key": event_matching_key,
                        "coordinates": coordinate_text,
                        "dPSI": float(dpsi) if pd.notna(dpsi) else float("nan"),
                        "FDR": float(fdr) if pd.notna(fdr) else float("nan"),
                        "significant": significant,
                        "rmats_sig_before_stability": bool(rmats_sig_before_stability.loc[idx]),
                        "stable_as_event": bool(stable_as_event.loc[idx]),
                        "direction": (
                            f"inclusion higher in {left_label} than {right_label}"
                            if pd.notna(dpsi) and float(dpsi) > 0
                            else f"inclusion lower in {left_label} than {right_label}"
                            if pd.notna(dpsi) and float(dpsi) < 0
                            else "no directional shift"
                        ),
                        "group1_label": left_label,
                        "group2_label": right_label,
                        "psi_group1_mean": stats_1["mean_psi"],
                        "psi_group2_mean": stats_2["mean_psi"],
                        "psi_group1_values": psi_1,
                        "psi_group2_values": psi_2,
                        "ijc_group1_values": ijc_1,
                        "ijc_group2_values": ijc_2,
                        "sjc_group1_values": sjc_1,
                        "sjc_group2_values": sjc_2,
                        "mean_total_reads": mean_total_reads,
                        "replicate_count": replicate_count,
                        "effect_logit": effect_logit,
                        "effect_se": effect_se,
                        "replicate_consistency_score": consistency,
                        "read_support_ok": support_ok,
                        "abs_dPSI_significance_cutoff_used": abs_dpsi_cutoff,
                        "min_total_reads_used": min_total_reads,
                    }
                )
        if not rows:
            return pd.DataFrame()
        output = pd.DataFrame(rows)
        output = output.sort_values(["significant", "FDR", "event_type", "gene"], ascending=[False, True, True, True], na_position="last")
        return output.drop_duplicates(subset=["event_matching_key"], keep="first").reset_index(drop=True)

    def _build_cross_comparison_as_interactions(
        self,
        pair: ComparisonPairDefinition,
        *,
        abs_dpsi_cutoff: float,
        large_delta_dpsi_cutoff: float,
        min_total_reads: int,
    ) -> dict[str, pd.DataFrame]:
        comparison_a = self._cross_as_interaction_source_frame(pair.comparison_a, abs_dpsi_cutoff, min_total_reads)
        comparison_b = self._cross_as_interaction_source_frame(pair.comparison_b, abs_dpsi_cutoff, min_total_reads)
        if comparison_a.empty and comparison_b.empty:
            return self._empty_cross_as_interaction_tables()
        n_significant_events_a = len(comparison_a)
        n_significant_events_b = len(comparison_b)

        keys_a = set(comparison_a.get("event_matching_key", pd.Series(dtype="object")).dropna().astype(str))
        keys_b = set(comparison_b.get("event_matching_key", pd.Series(dtype="object")).dropna().astype(str))
        a_only_keys = keys_a - keys_b
        b_only_keys = keys_b - keys_a
        if a_only_keys and pair.comparison_b:
            counterpart_b = self._cross_as_interaction_source_frame(
                pair.comparison_b,
                abs_dpsi_cutoff,
                min_total_reads,
                significant_only=False,
                candidate_keys=a_only_keys,
            )
            if not counterpart_b.empty:
                comparison_b = (
                    pd.concat([comparison_b, counterpart_b], axis=0, ignore_index=True)
                    .sort_values(["significant", "FDR"], ascending=[False, True], na_position="last")
                    .drop_duplicates(subset=["event_matching_key"], keep="first")
                    .reset_index(drop=True)
                )
        if b_only_keys and pair.comparison_a:
            counterpart_a = self._cross_as_interaction_source_frame(
                pair.comparison_a,
                abs_dpsi_cutoff,
                min_total_reads,
                significant_only=False,
                candidate_keys=b_only_keys,
            )
            if not counterpart_a.empty:
                comparison_a = (
                    pd.concat([comparison_a, counterpart_a], axis=0, ignore_index=True)
                    .sort_values(["significant", "FDR"], ascending=[False, True], na_position="last")
                    .drop_duplicates(subset=["event_matching_key"], keep="first")
                    .reset_index(drop=True)
                )

        merged = comparison_a.merge(comparison_b, on="event_matching_key", how="outer", suffixes=("_A", "_B"))
        if merged.empty:
            return self._empty_cross_as_interaction_tables()

        def coalesce(left: str, right: str, default: object = "") -> pd.Series:
            return (
                merged.get(left, pd.Series(pd.NA, index=merged.index))
                .combine_first(merged.get(right, pd.Series(pd.NA, index=merged.index)))
                .fillna(default)
            )

        merged["gene"] = coalesce("gene_A", "gene_B", "").astype(str)
        merged["event_type"] = coalesce("event_type_A", "event_type_B", "").astype(str)
        merged["event_id"] = coalesce("event_id_A", "event_id_B", "").astype(str)
        merged["comparison_A"] = pair.comparison_a or ""
        merged["comparison_B"] = pair.comparison_b or ""
        merged["experiment_A_display_name"] = (
            comparison_a.get("comparison_display_name", pd.Series(dtype="object")).iloc[0] if not comparison_a.empty else str(pair.comparison_a or "")
        )
        merged["experiment_B_display_name"] = (
            comparison_b.get("comparison_display_name", pd.Series(dtype="object")).iloc[0] if not comparison_b.empty else str(pair.comparison_b or "")
        )
        for column in ("dPSI_A", "dPSI_B", "FDR_A", "FDR_B", "effect_logit_A", "effect_logit_B", "effect_se_A", "effect_se_B", "mean_total_reads_A", "mean_total_reads_B", "replicate_consistency_score_A", "replicate_consistency_score_B"):
            merged[column] = pd.to_numeric(merged.get(column), errors="coerce")
        merged["significant_A"] = merged.get("significant_A", pd.Series(False, index=merged.index)).map(
            lambda value: bool(value) if pd.notna(value) else False
        )
        merged["significant_B"] = merged.get("significant_B", pd.Series(False, index=merged.index)).map(
            lambda value: bool(value) if pd.notna(value) else False
        )
        merged["direction_A"] = merged.get("direction_A", pd.Series("", index=merged.index)).fillna("").astype(str)
        merged["direction_B"] = merged.get("direction_B", pd.Series("", index=merged.index)).fillna("").astype(str)
        merged["replicate_count_A"] = pd.to_numeric(merged.get("replicate_count_A"), errors="coerce").fillna(0).astype(int)
        merged["replicate_count_B"] = pd.to_numeric(merged.get("replicate_count_B"), errors="coerce").fillna(0).astype(int)
        merged["delta_dPSI"] = merged["dPSI_A"] - merged["dPSI_B"]
        merged["abs_delta_dPSI"] = merged["delta_dPSI"].abs()
        sign_a = np.sign(merged["dPSI_A"].fillna(0))
        sign_b = np.sign(merged["dPSI_B"].fillna(0))
        merged["same_direction"] = (sign_a == sign_b) & sign_a.ne(0) & sign_b.ne(0)
        merged["opposite_direction"] = sign_a.ne(sign_b) & sign_a.ne(0) & sign_b.ne(0)
        merged["same_direction_strength_delta"] = (merged["dPSI_A"].abs() - merged["dPSI_B"].abs()).abs()
        merged["fallback_interaction_log_odds"] = merged["effect_logit_A"] - merged["effect_logit_B"]
        merged["fallback_interaction_se"] = np.sqrt(
            pd.to_numeric(merged["effect_se_A"], errors="coerce").pow(2)
            + pd.to_numeric(merged["effect_se_B"], errors="coerce").pow(2)
        )
        merged["fallback_interaction_z"] = merged["fallback_interaction_log_odds"] / merged["fallback_interaction_se"].replace(0, np.nan)
        merged["fallback_interaction_pvalue"] = merged["fallback_interaction_z"].map(self._normal_two_sided_pvalue)
        merged["interaction_log_odds"] = merged["fallback_interaction_log_odds"]
        merged["interaction_se"] = merged["fallback_interaction_se"]
        merged["interaction_z"] = merged["fallback_interaction_z"]
        merged["interaction_pvalue"] = merged["fallback_interaction_pvalue"]
        merged["interaction_model"] = "fallback_logit_difference"
        merged["interaction_model_status"] = "fallback_before_binomial_glm"
        merged["interaction_overdispersion_scale"] = np.nan
        merged["replicate_consistency_score"] = pd.concat(
            [
                merged["replicate_consistency_score_A"],
                merged["replicate_consistency_score_B"],
            ],
            axis=1,
        ).min(axis=1, skipna=True)
        merged["read_support_class"] = np.where(
            merged["mean_total_reads_A"].ge(min_total_reads) & merged["mean_total_reads_B"].ge(min_total_reads),
            "adequate_read_support",
            "low_power_or_missing_read_support",
        )
        merged["abs_dPSI_significance_cutoff_used"] = abs_dpsi_cutoff
        merged["large_delta_dPSI_cutoff_used"] = large_delta_dpsi_cutoff
        merged["min_total_reads_used"] = min_total_reads

        shared = merged.loc[merged["significant_A"] & merged["significant_B"]].copy()
        if not shared.empty:
            glm_results = shared.apply(lambda row: pd.Series(self._fit_binomial_glm_interaction(row)), axis=1)
            for column in glm_results.columns:
                shared[column] = glm_results[column]
            fallback_mask = shared["interaction_pvalue"].isna()
            if fallback_mask.any():
                shared.loc[fallback_mask, "interaction_log_odds"] = shared.loc[fallback_mask, "fallback_interaction_log_odds"]
                shared.loc[fallback_mask, "interaction_se"] = shared.loc[fallback_mask, "fallback_interaction_se"]
                shared.loc[fallback_mask, "interaction_z"] = shared.loc[fallback_mask, "fallback_interaction_z"]
                shared.loc[fallback_mask, "interaction_pvalue"] = shared.loc[fallback_mask, "fallback_interaction_pvalue"]
                shared.loc[fallback_mask, "interaction_model"] = "fallback_logit_difference"
                shared.loc[fallback_mask, "interaction_model_status"] = (
                    "binomial_glm_unavailable; fallback to replicate logit difference"
                )
                shared.loc[fallback_mask, "interaction_overdispersion_scale"] = np.nan
                shared.loc[fallback_mask, "quasi_binomial_dispersion"] = np.nan
                shared.loc[fallback_mask, "quasi_binomial_overdispersed"] = False
                shared.loc[fallback_mask, "beta_binomial_rho_estimate"] = np.nan
                shared.loc[fallback_mask, "model_valid"] = False
                shared.loc[fallback_mask, "model_diagnostic_class"] = "fallback_logit_difference"
                shared.loc[fallback_mask, "model_warning"] = "Binomial GLM unavailable; using replicate logit-difference fallback."
            beta_results = shared.apply(lambda row: pd.Series(self._fit_beta_binomial_interaction(row)), axis=1)
            for column in beta_results.columns:
                shared[column] = beta_results[column]
        shared["interaction_FDR"] = self._bh_adjust(shared["interaction_pvalue"]) if not shared.empty else pd.Series(dtype="float64")
        shared["beta_binomial_FDR"] = self._bh_adjust(shared["beta_binomial_pvalue"]) if not shared.empty and "beta_binomial_pvalue" in shared.columns else pd.Series(dtype="float64")
        if not shared.empty:
            effect = pd.to_numeric(shared["interaction_log_odds"], errors="coerce")
            se = pd.to_numeric(shared["interaction_se"], errors="coerce")
            valid = effect.notna() & se.notna() & se.gt(0)
            if valid.sum() > 1:
                observed_variance = float(effect.loc[valid].var(ddof=1))
                # A few near-separated GLMs can have enormous SEs; the median gives a
                # robust empirical-Bayes prior scale without letting those dominate.
                sampling_variance = float(se.loc[valid].pow(2).median())
                prior_variance = max(observed_variance - sampling_variance, 0.0)
            else:
                prior_variance = 0.0
            if prior_variance > 0:
                shrinkage_weight = prior_variance / (prior_variance + se.pow(2))
            else:
                shrinkage_weight = pd.Series(0.0, index=shared.index, dtype="float64")
            shared["eb_prior_variance"] = prior_variance
            shared["shrinkage_weight"] = shrinkage_weight.clip(lower=0.0, upper=1.0)
            shared["shrunken_interaction_log_odds"] = effect * shared["shrinkage_weight"]
            shared["shrunken_interaction_z"] = shared["shrunken_interaction_log_odds"] / se.replace(0, np.nan)
            shared["shrunken_interaction_pvalue"] = shared["shrunken_interaction_z"].map(self._normal_two_sided_pvalue)
            shared["shrunken_interaction_FDR"] = self._bh_adjust(shared["shrunken_interaction_pvalue"])
            effect_component = (shared["abs_delta_dPSI"] / max(large_delta_dpsi_cutoff, 1e-6)).clip(lower=0.0, upper=3.0) / 3.0
            confidence_component = (
                -np.log10(pd.to_numeric(shared["shrunken_interaction_pvalue"], errors="coerce").clip(lower=1e-300))
            ).clip(lower=0.0, upper=5.0) / 5.0
            consistency_component = pd.to_numeric(shared["replicate_consistency_score"], errors="coerce").fillna(0.0).clip(lower=0.0, upper=1.0)
            read_support_component = (
                pd.concat(
                    [
                        pd.to_numeric(shared["mean_total_reads_A"], errors="coerce"),
                        pd.to_numeric(shared["mean_total_reads_B"], errors="coerce"),
                    ],
                    axis=1,
                )
                .min(axis=1, skipna=True)
                .fillna(0.0)
                .div(max(float(min_total_reads), 1.0) * 5.0)
                .clip(lower=0.0, upper=1.0)
            )
            shared["interaction_priority_score"] = (
                100.0
                * (
                    0.35 * effect_component
                    + 0.25 * confidence_component
                    + 0.20 * consistency_component
                    + 0.20 * read_support_component
                )
            ).round(3)
            simple_delta_score = pd.to_numeric(shared["abs_delta_dPSI"], errors="coerce")
            shared["simple_delta_score"] = simple_delta_score
            shared["simple_delta_priority_score"] = (
                100.0
                * (simple_delta_score / max(float(large_delta_dpsi_cutoff), 1e-6)).clip(lower=0.0, upper=3.0)
                / 3.0
            ).round(3)
            shared["delta_dPSI_rank"] = simple_delta_score.rank(method="min", ascending=False, na_option="bottom").astype("Int64")
            shared["interaction_model_rank"] = pd.to_numeric(
                shared["interaction_priority_score"], errors="coerce"
            ).rank(method="min", ascending=False, na_option="bottom").astype("Int64")
            shared["delta_vs_model_rank_shift"] = (
                pd.to_numeric(shared["delta_dPSI_rank"], errors="coerce")
                - pd.to_numeric(shared["interaction_model_rank"], errors="coerce")
            )
            top_n = min(20, len(shared))
            shared["simple_delta_top20"] = pd.to_numeric(shared["delta_dPSI_rank"], errors="coerce").le(top_n)
            shared["model_top20"] = pd.to_numeric(shared["interaction_model_rank"], errors="coerce").le(top_n)
            shared["benchmark_score_delta"] = shared["interaction_priority_score"] - shared["simple_delta_priority_score"]
            shared["benchmark_class"] = np.select(
                [
                    shared["simple_delta_top20"] & shared["model_top20"],
                    shared["simple_delta_top20"] & ~shared["model_top20"],
                    ~shared["simple_delta_top20"] & shared["model_top20"],
                ],
                [
                    "top_by_both_delta_and_model",
                    "delta_only_top20",
                    "model_only_top20",
                ],
                default="not_top20",
            )
        else:
            for column in (
                "eb_prior_variance",
                "shrinkage_weight",
                "shrunken_interaction_log_odds",
                "shrunken_interaction_z",
                "shrunken_interaction_pvalue",
                "shrunken_interaction_FDR",
                "interaction_priority_score",
                "simple_delta_score",
                "simple_delta_priority_score",
                "delta_dPSI_rank",
                "interaction_model_rank",
                "delta_vs_model_rank_shift",
                "benchmark_score_delta",
            ):
                shared[column] = pd.Series(dtype="float64")
            for column in ("simple_delta_top20", "model_top20"):
                shared[column] = pd.Series(dtype="bool")
            shared["benchmark_class"] = pd.Series(dtype="object")
        shared["interaction_class"] = np.where(
            shared["shrunken_interaction_FDR"].le(0.05) & shared["abs_delta_dPSI"].ge(large_delta_dpsi_cutoff) & shared["read_support_class"].eq("adequate_read_support"),
            "high_confidence_interaction_candidate",
            np.where(
                shared["shrunken_interaction_pvalue"].le(0.05) & shared["abs_delta_dPSI"].ge(large_delta_dpsi_cutoff),
                "moderate_interaction_candidate",
                "pattern_followup_only",
            ),
        )
        shared["event_class"] = "shared_significant_interaction_tested"
        opposite = shared.loc[shared["opposite_direction"]].copy()
        opposite["event_class"] = "opposite_direction_shared_significant"
        large_delta = shared.loc[
            shared["same_direction"] & shared["same_direction_strength_delta"].ge(large_delta_dpsi_cutoff)
        ].copy()
        large_delta["event_class"] = "same_direction_large_dpsi_delta"
        benchmark = shared.sort_values(
            ["benchmark_class", "interaction_model_rank", "delta_dPSI_rank"],
            ascending=[True, True, True],
            na_position="last",
        ).copy()

        exclusive = merged.loc[
            (merged["significant_A"] & ~merged["significant_B"])
            | (~merged["significant_A"] & merged["significant_B"])
        ].copy()
        if not exclusive.empty:
            exclusive["interaction_FDR"] = np.nan
            exclusive["interaction_class"] = "not_interaction_tested_exclusive_event"
            exclusive["exclusive_side"] = np.where(exclusive["significant_A"], "comparison_A_only", "comparison_B_only")
            exclusive["counterpart_event_found"] = np.where(
                exclusive["significant_A"],
                exclusive["dPSI_B"].notna(),
                exclusive["dPSI_A"].notna(),
            )
            exclusive["significant_side_mean_total_reads"] = np.where(
                exclusive["significant_A"],
                exclusive["mean_total_reads_A"],
                exclusive["mean_total_reads_B"],
            )
            exclusive["non_significant_side_mean_total_reads"] = np.where(
                exclusive["significant_A"],
                exclusive["mean_total_reads_B"],
                exclusive["mean_total_reads_A"],
            )
            exclusive["significant_side_replicate_count"] = np.where(
                exclusive["significant_A"],
                exclusive["replicate_count_A"],
                exclusive["replicate_count_B"],
            )
            exclusive["non_significant_side_replicate_count"] = np.where(
                exclusive["significant_A"],
                exclusive["replicate_count_B"],
                exclusive["replicate_count_A"],
            )
            exclusive["non_significant_side_dPSI"] = np.where(
                exclusive["significant_A"],
                exclusive["dPSI_B"],
                exclusive["dPSI_A"],
            )
            exclusive["non_significant_side_FDR"] = np.where(
                exclusive["significant_A"],
                exclusive["FDR_B"],
                exclusive["FDR_A"],
            )
            exclusive["non_significant_side_read_support"] = np.select(
                [
                    exclusive["significant_A"] & ~exclusive["counterpart_event_found"],
                    ~exclusive["significant_A"] & ~exclusive["counterpart_event_found"],
                    exclusive["significant_A"] & exclusive["non_significant_side_mean_total_reads"].astype(float).ge(min_total_reads),
                    ~exclusive["significant_A"] & exclusive["non_significant_side_mean_total_reads"].astype(float).ge(min_total_reads),
                    exclusive["significant_A"],
                    ~exclusive["significant_A"],
                ],
                [
                    "counterpart_event_not_found_in_B",
                    "counterpart_event_not_found_in_A",
                    "supported_non_significant_in_B",
                    "supported_non_significant_in_A",
                    "low_power_or_absent_in_B",
                    "low_power_or_absent_in_A",
                ],
                default="not_classified",
            )
            exclusive["event_class"] = exclusive["exclusive_side"] + "|" + exclusive["non_significant_side_read_support"]

        output_columns = self._empty_cross_as_interaction_tables()["shared_interaction"].columns.tolist()
        model_status = shared.get("interaction_model_status", pd.Series(dtype="object")).fillna("").astype(str) if not shared.empty else pd.Series(dtype="object")
        diagnostic_class = shared.get("model_diagnostic_class", pd.Series(dtype="object")).fillna("").astype(str) if not shared.empty else pd.Series(dtype="object")
        top20_overlap = int((shared.get("simple_delta_top20", pd.Series(dtype=bool)) & shared.get("model_top20", pd.Series(dtype=bool))).sum()) if not shared.empty else 0
        top20_union = int((shared.get("simple_delta_top20", pd.Series(dtype=bool)) | shared.get("model_top20", pd.Series(dtype=bool))).sum()) if not shared.empty else 0
        rank_corr = float("nan")
        if not shared.empty and len(shared) > 1:
            rank_frame = shared[["delta_dPSI_rank", "interaction_model_rank"]].apply(pd.to_numeric, errors="coerce").dropna()
            if len(rank_frame) > 1:
                rank_corr = float(rank_frame["delta_dPSI_rank"].corr(rank_frame["interaction_model_rank"], method="spearman"))
        diagnostics = pd.DataFrame(
            [
                {
                    "pair_id": pair.pair_id,
                    "pair_name": pair.resolved_name,
                    "comparison_A": pair.comparison_a,
                    "comparison_B": pair.comparison_b,
                    "n_events_A": n_significant_events_a,
                    "n_events_B": n_significant_events_b,
                    "n_events_A_with_counterparts": len(comparison_a),
                    "n_events_B_with_counterparts": len(comparison_b),
                    "n_shared_significant": len(shared),
                    "n_opposite_direction": len(opposite),
                    "n_large_delta": len(large_delta),
                    "n_exclusive_supported": int(exclusive.get("non_significant_side_read_support", pd.Series(dtype="object")).astype(str).str.startswith("supported_").sum()) if not exclusive.empty else 0,
                    "n_exclusive_low_power": int(exclusive.get("non_significant_side_read_support", pd.Series(dtype="object")).astype(str).str.startswith("low_power").sum()) if not exclusive.empty else 0,
                    "n_exclusive_counterpart_not_found": int(exclusive.get("non_significant_side_read_support", pd.Series(dtype="object")).astype(str).str.startswith("counterpart_event_not_found").sum()) if not exclusive.empty else 0,
                    "n_exclusive_not_shared_significant": len(exclusive),
                    "n_model_converged": int(model_status.eq("converged").sum()),
                    "n_model_fallback": int(diagnostic_class.eq("fallback_logit_difference").sum()),
                    "n_model_failed": int(diagnostic_class.eq("model_failed").sum()),
                    "median_quasi_binomial_dispersion": float(pd.to_numeric(shared.get("quasi_binomial_dispersion", pd.Series(dtype="float64")), errors="coerce").median()) if not shared.empty else float("nan"),
                    "max_quasi_binomial_dispersion": float(pd.to_numeric(shared.get("quasi_binomial_dispersion", pd.Series(dtype="float64")), errors="coerce").max()) if not shared.empty else float("nan"),
                    "median_beta_binomial_rho_estimate": float(pd.to_numeric(shared.get("beta_binomial_rho_estimate", pd.Series(dtype="float64")), errors="coerce").median()) if not shared.empty else float("nan"),
                    "median_beta_binomial_rho": float(pd.to_numeric(shared.get("beta_binomial_rho", pd.Series(dtype="float64")), errors="coerce").median()) if not shared.empty else float("nan"),
                    "n_beta_binomial_converged": int(shared.get("beta_binomial_model_status", pd.Series(dtype="object")).astype(str).eq("converged").sum()) if not shared.empty else 0,
                    "n_beta_binomial_failed": int(shared.get("beta_binomial_model_status", pd.Series(dtype="object")).astype(str).isin(["failed", "fit_warning"]).sum()) if not shared.empty else 0,
                    "n_beta_binomial_FDR_0_05": int(pd.to_numeric(shared.get("beta_binomial_FDR", pd.Series(dtype="float64")), errors="coerce").le(0.05).sum()) if not shared.empty else 0,
                    "n_overdispersed_events": int(shared.get("quasi_binomial_overdispersed", pd.Series(dtype=bool)).map(self._safe_bool).sum()) if not shared.empty else 0,
                    "delta_top20_model_top20_overlap": top20_overlap,
                    "delta_top20_model_top20_overlap_fraction": float(top20_overlap / top20_union) if top20_union else float("nan"),
                    "delta_vs_model_spearman_rank_correlation": rank_corr,
                    "n_delta_only_top20": int(shared.get("benchmark_class", pd.Series(dtype="object")).astype(str).eq("delta_only_top20").sum()) if not shared.empty else 0,
                    "n_model_only_top20": int(shared.get("benchmark_class", pd.Series(dtype="object")).astype(str).eq("model_only_top20").sum()) if not shared.empty else 0,
                    "abs_dPSI_significance_cutoff_used": abs_dpsi_cutoff,
                    "large_delta_dPSI_cutoff_used": large_delta_dpsi_cutoff,
                    "min_total_reads_used": min_total_reads,
                }
            ]
        )
        unified_parts: list[pd.DataFrame] = []
        if not shared.empty:
            shared_unified = shared.copy()
            shared_unified["unified_result_type"] = "shared_significant_modeled"
            shared_unified["unified_display_class"] = "Shared significant modeled"
            shared_unified.loc[shared_unified.get("opposite_direction", pd.Series(False, index=shared_unified.index)).map(self._safe_bool), "unified_display_class"] = "Opposite direction"
            shared_unified.loc[
                shared_unified.get("same_direction", pd.Series(False, index=shared_unified.index)).map(self._safe_bool)
                & pd.to_numeric(shared_unified.get("same_direction_strength_delta", pd.Series(np.nan, index=shared_unified.index)), errors="coerce").ge(large_delta_dpsi_cutoff),
                "unified_display_class",
            ] = "Same direction + large dPSI delta"
            shared_unified["support_class"] = shared_unified.get("read_support_class", "shared_significant")
            unified_parts.append(shared_unified)
        if not exclusive.empty:
            exclusive_unified = exclusive.copy()
            side = exclusive_unified.get("exclusive_side", pd.Series("", index=exclusive_unified.index)).fillna("").astype(str)
            support = exclusive_unified.get("non_significant_side_read_support", pd.Series("", index=exclusive_unified.index)).fillna("").astype(str)
            exclusive_unified["unified_result_type"] = np.where(
                side.eq("comparison_A_only"),
                "A_specific_significant",
                np.where(side.eq("comparison_B_only"), "B_specific_significant", "exclusive_significant"),
            )
            exclusive_unified["unified_display_class"] = np.where(
                side.eq("comparison_A_only"),
                "A-specific significant",
                np.where(side.eq("comparison_B_only"), "B-specific significant", "Specific significant"),
            )
            exclusive_unified["support_class"] = support
            unified_parts.append(exclusive_unified)
        if unified_parts:
            unified = pd.concat(unified_parts, ignore_index=True, sort=False)
            unified["_sort_priority"] = unified["unified_result_type"].map(
                {
                    "shared_significant_modeled": 0,
                    "A_specific_significant": 1,
                    "B_specific_significant": 2,
                    "exclusive_significant": 3,
                }
            ).fillna(9)
            unified = unified.sort_values(
                ["_sort_priority", "interaction_priority_score", "abs_delta_dPSI", "gene"],
                ascending=[True, False, False, True],
                na_position="last",
            ).drop(columns=["_sort_priority"], errors="ignore")
        else:
            unified = pd.DataFrame(columns=output_columns)
        return {
            "unified_event_table": unified.reindex(columns=output_columns).reset_index(drop=True),
            "shared_interaction": shared.reindex(columns=output_columns).reset_index(drop=True),
            "opposite_direction": opposite.reindex(columns=output_columns).reset_index(drop=True),
            "large_delta": large_delta.reindex(columns=output_columns).reset_index(drop=True),
            "exclusive": exclusive.reindex(columns=output_columns).reset_index(drop=True),
            "benchmark": benchmark.reindex(columns=output_columns).reset_index(drop=True),
            "diagnostics": diagnostics.reset_index(drop=True),
        }
