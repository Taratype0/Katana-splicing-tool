from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

from katana_cosplice.utils import bh_adjust


ASSOCIATION_COLUMNS = [
    "module",
    "factor",
    "contrast",
    "level",
    "effect",
    "statistic",
    "pvalue",
    "test",
    "test_warning",
]
FACTORIAL_COLUMNS = [
    "module",
    "model_term",
    "sum_sq",
    "df",
    "F",
    "pvalue",
    "partial_eta_squared",
    "n_samples",
    "residual_df",
    "model_r_squared",
]
COEFFICIENT_COLUMNS = [
    "module",
    "coefficient",
    "estimate",
    "std_error",
    "t",
    "pvalue",
    "ci_low",
    "ci_high",
    "n_samples",
]


def _aligned_metadata(program_scores: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    if not program_scores.index.is_unique:
        raise ValueError("Program score sample index must be unique.")
    if "sample_id" not in samples.columns or not samples["sample_id"].is_unique:
        raise ValueError("Association metadata requires a unique sample_id column.")
    missing = sorted(set(program_scores.index).difference(samples["sample_id"]))
    if missing:
        raise ValueError(f"Program score samples are absent from metadata: {missing}")
    return samples.set_index("sample_id").reindex(program_scores.index)


def summarize_program_scores(
    program_scores: pd.DataFrame,
    samples: pd.DataFrame,
    factors: list[str],
) -> pd.DataFrame:
    metadata = _aligned_metadata(program_scores, samples)
    rows: list[dict[str, object]] = []
    for factor in factors:
        if factor not in metadata.columns:
            continue
        for level, sample_index in metadata.groupby(factor, sort=False, observed=True).groups.items():
            subset = program_scores.loc[list(sample_index)]
            for module in program_scores.columns:
                values = pd.to_numeric(subset[module], errors="coerce").dropna()
                rows.append(
                    {
                        "factor": factor,
                        "level": str(level),
                        "module": module,
                        "n_samples": len(values),
                        "mean_program_score": float(values.mean()) if len(values) else np.nan,
                        "sd_program_score": float(values.std(ddof=1)) if len(values) > 1 else np.nan,
                        "median_program_score": float(values.median()) if len(values) else np.nan,
                    }
                )
    return pd.DataFrame(rows)


def associate_program_scores(
    program_scores: pd.DataFrame,
    samples: pd.DataFrame,
    factors: list[str] | None = None,
) -> pd.DataFrame:
    metadata = _aligned_metadata(program_scores, samples)
    if factors is None:
        factors = ["condition"]
    rows: list[dict[str, object]] = []
    for factor in factors:
        if factor not in metadata.columns:
            continue
        labels = metadata[factor].astype(str)
        levels = [level for level in labels.dropna().unique() if level not in {"", "nan", "None"}]
        if len(levels) < 2:
            continue
        for module in program_scores.columns:
            values = pd.to_numeric(program_scores[module], errors="coerce")
            groups = [values.loc[labels.eq(level)].dropna().to_numpy() for level in levels]
            groups = [group for group in groups if len(group) >= 2]
            if len(groups) >= 2:
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    statistic, pvalue = stats.f_oneway(*groups)
                test_warning = "; ".join(dict.fromkeys(str(item.message) for item in caught))
                grand_mean = float(values.mean())
                ss_between = sum(len(group) * (float(np.mean(group)) - grand_mean) ** 2 for group in groups)
                ss_total = float(np.nansum((values.to_numpy() - grand_mean) ** 2))
                eta_squared = ss_between / ss_total if ss_total > 0 else 0.0
            else:
                statistic, pvalue, eta_squared = np.nan, np.nan, np.nan
                test_warning = "insufficient groups with at least two samples"
            rows.append(
                {
                    "module": module,
                    "factor": factor,
                    "contrast": "omnibus",
                    "level": "all",
                    "effect": eta_squared,
                    "statistic": statistic,
                    "pvalue": pvalue,
                    "test": "one-way ANOVA; effect=eta_squared",
                    "test_warning": test_warning,
                }
            )
            for level in levels:
                in_level = values.loc[labels.eq(level)].dropna()
                outside = values.loc[~labels.eq(level)].dropna()
                if len(in_level) < 2 or len(outside) < 2:
                    statistic, pvalue = np.nan, np.nan
                    test_warning = "insufficient samples for Welch t-test"
                else:
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter("always")
                        statistic, pvalue = stats.ttest_ind(in_level, outside, equal_var=False, nan_policy="omit")
                    test_warning = "; ".join(dict.fromkeys(str(item.message) for item in caught))
                rows.append(
                    {
                        "module": module,
                        "factor": factor,
                        "contrast": f"{level}_vs_rest",
                        "level": level,
                        "effect": float(in_level.mean() - outside.mean()) if len(in_level) and len(outside) else np.nan,
                        "statistic": statistic,
                        "pvalue": pvalue,
                        "test": "Welch t-test; effect=mean program-score difference",
                        "test_warning": test_warning,
                    }
                )
    result = pd.DataFrame(rows, columns=ASSOCIATION_COLUMNS)
    if not result.empty:
        result["FDR"] = bh_adjust(result["pvalue"])
        result = result.sort_values(["FDR", "pvalue", "module"], na_position="last").reset_index(drop=True)
    return result


def render_association_plot(table: pd.DataFrame, output_path: Path) -> Path | None:
    subset = table.loc[table["contrast"].ne("omnibus")].copy() if not table.empty else pd.DataFrame()
    if subset.empty:
        return None
    matrix = subset.pivot_table(index="module", columns="contrast", values="effect", aggfunc="first")
    if matrix.empty:
        return None
    width = max(11.0, matrix.shape[1] * 0.42)
    height = max(5.5, matrix.shape[0] * 0.5)
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(matrix, cmap="vlag", center=0.0, ax=ax, cbar_kws={"label": "Program-score difference"})
    ax.set_title("Co-splicing module association with experimental groups", fontweight="bold", pad=14)
    ax.set_xlabel("Experimental level vs all other samples")
    ax.set_ylabel("Co-splicing module")
    ax.tick_params(axis="x", rotation=90)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def factorial_program_association(
    program_scores: pd.DataFrame,
    samples: pd.DataFrame,
    factors: list[str],
    batch_column: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not factors:
        return pd.DataFrame(columns=FACTORIAL_COLUMNS), pd.DataFrame(columns=COEFFICIENT_COLUMNS)
    metadata = _aligned_metadata(program_scores, samples).copy()
    required = [*factors, *([batch_column] if batch_column else [])]
    missing = [column for column in required if column not in metadata.columns]
    if missing:
        raise ValueError(f"Factorial association metadata columns are missing: {missing}")
    encoded = [f'C(Q("{column}"), Sum)' for column in factors]
    rhs = " * ".join(encoded)
    if batch_column:
        rhs += f' + C(Q("{batch_column}"), Sum)'
    anova_rows: list[dict[str, object]] = []
    coefficient_rows: list[dict[str, object]] = []
    for module in program_scores.columns:
        data = metadata.loc[:, required].copy()
        data["score"] = pd.to_numeric(program_scores[module], errors="coerce")
        data = data.dropna()
        if len(data) <= 4:
            continue
        model = smf.ols(f"score ~ {rhs}", data=data).fit()
        table = anova_lm(model, typ=3)
        residual_ss = float(table.loc["Residual", "sum_sq"]) if "Residual" in table.index else np.nan
        for term, row in table.drop(index=["Intercept", "Residual"], errors="ignore").iterrows():
            sum_sq = float(row.get("sum_sq", np.nan))
            anova_rows.append(
                {
                    "module": module,
                    "model_term": str(term),
                    "sum_sq": sum_sq,
                    "df": float(row.get("df", np.nan)),
                    "F": float(row.get("F", np.nan)),
                    "pvalue": float(row.get("PR(>F)", np.nan)),
                    "partial_eta_squared": sum_sq / (sum_sq + residual_ss) if np.isfinite(residual_ss) and sum_sq + residual_ss > 0 else np.nan,
                    "n_samples": int(model.nobs),
                    "residual_df": float(model.df_resid),
                    "model_r_squared": float(model.rsquared),
                }
            )
        confidence = model.conf_int()
        for term in model.params.index:
            if term == "Intercept":
                continue
            coefficient_rows.append(
                {
                    "module": module,
                    "coefficient": str(term),
                    "estimate": float(model.params[term]),
                    "std_error": float(model.bse[term]),
                    "t": float(model.tvalues[term]),
                    "pvalue": float(model.pvalues[term]),
                    "ci_low": float(confidence.loc[term, 0]),
                    "ci_high": float(confidence.loc[term, 1]),
                    "n_samples": int(model.nobs),
                }
            )
    anova = pd.DataFrame(anova_rows, columns=FACTORIAL_COLUMNS)
    coefficients = pd.DataFrame(coefficient_rows, columns=COEFFICIENT_COLUMNS)
    if not anova.empty:
        anova["FDR"] = bh_adjust(anova["pvalue"])
        anova = anova.sort_values(["FDR", "module", "model_term"], na_position="last").reset_index(drop=True)
    if not coefficients.empty:
        coefficients["FDR"] = bh_adjust(coefficients["pvalue"])
        coefficients = coefficients.sort_values(["FDR", "module", "coefficient"], na_position="last").reset_index(drop=True)
    return anova, coefficients


def render_factorial_association_plot(table: pd.DataFrame, output_path: Path) -> Path | None:
    if table.empty:
        return None
    display = table.copy()
    display["display_term"] = display["model_term"].map(_display_model_term)
    matrix = display.pivot_table(index="module", columns="display_term", values="partial_eta_squared", aggfunc="first")
    if matrix.empty:
        return None
    preferred_order = [
        "genotype",
        "DP status",
        "treatment",
        "genotype x DP status",
        "genotype x treatment",
        "DP status x treatment",
        "genotype x DP status x treatment",
        "batch",
    ]
    matrix = matrix.reindex(columns=[term for term in preferred_order if term in matrix.columns])
    width = max(10.5, matrix.shape[1] * 1.15)
    height = max(6.5, matrix.shape[0] * 0.55 + 2.0)
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(matrix, cmap="YlGnBu", vmin=0.0, vmax=1.0, ax=ax, cbar_kws={"label": "Partial eta-squared"})
    network_name = output_path.parent.name.title()
    ax.set_title(f"{network_name} factorial effects on co-splicing program scores", fontweight="bold", pad=14)
    ax.set_xlabel("Factorial model term")
    ax.set_ylabel("Co-splicing module")
    ax.tick_params(axis="x", rotation=35)
    ax.tick_params(axis="y", rotation=0)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def _display_model_term(term: str) -> str:
    labels = {
        'C(Q("genotype"), Sum)': "genotype",
        'C(Q("DP_status"), Sum)': "DP status",
        'C(Q("treatment"), Sum)': "treatment",
        'C(Q("batch"), Sum)': "batch",
    }
    parts = [labels.get(part, part) for part in str(term).split(":")]
    return " x ".join(parts)
