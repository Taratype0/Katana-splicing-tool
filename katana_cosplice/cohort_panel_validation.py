from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import binomtest, wilcoxon


@dataclass(frozen=True)
class CohortPanelValidation:
    cohort_registry: pd.DataFrame
    internal_comparison: pd.DataFrame
    preservation_comparison: pd.DataFrame
    scope_summary: pd.DataFrame
    functional_comparison: pd.DataFrame
    claim_audit: pd.DataFrame


def _read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def _cohort_registry(config_path: Path, panel_manifest_path: Path) -> pd.DataFrame:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    manifest = json.loads(panel_manifest_path.read_text(encoding="utf-8"))
    config_rows = pd.DataFrame(config["cohorts"])
    manifest_rows = pd.DataFrame(manifest["cohorts"])
    return config_rows.merge(manifest_rows, on="cohort", how="left", validate="one_to_one")


def _internal_comparison(internal: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "min_bootstrap_ARI",
        "median_bootstrap_ARI",
        "weighted_mean_within_module_coassignment",
        "n_modules",
        "grey_fraction",
        "valid_partition",
    ]
    available = [column for column in metrics if column in internal.columns]
    wide = internal.pivot(index="cohort", columns="method", values=available)
    wide.columns = [f"{metric}__{method}" for metric, method in wide.columns]
    wide = wide.reset_index()
    for metric in (
        "min_bootstrap_ARI",
        "median_bootstrap_ARI",
        "weighted_mean_within_module_coassignment",
    ):
        left = f"{metric}__CoSplice"
        right = f"{metric}__Exact WGCNA"
        if left in wide and right in wide:
            wide[f"{metric}__CoSplice_minus_WGCNA"] = (
                pd.to_numeric(wide[left], errors="coerce")
                - pd.to_numeric(wide[right], errors="coerce")
            )
    valid_left = pd.Series(
        wide.get("valid_partition__CoSplice", False),
        index=wide.index,
        dtype="boolean",
    ).fillna(False)
    valid_right = pd.Series(
        wide.get("valid_partition__Exact WGCNA", False),
        index=wide.index,
        dtype="boolean",
    ).fillna(False)
    wide["both_partitions_valid"] = (valid_left & valid_right).astype(bool)
    return wide


def _preservation_comparison(
    external: pd.DataFrame,
    registry: pd.DataFrame,
) -> pd.DataFrame:
    external = external.copy()
    def numeric_column(name: str) -> pd.Series:
        if name not in external.columns:
            return pd.Series(np.nan, index=external.index, dtype=float)
        return pd.to_numeric(external[name], errors="coerce")

    tested = numeric_column("tested_modules")
    retained = (
        numeric_column("preserved_modules").fillna(0)
        + numeric_column("partially_preserved_modules").fillna(0)
    )
    external["preserved_or_partial_module_fraction"] = np.where(
        tested.gt(0),
        retained / tested,
        np.nan,
    )
    family = registry.set_index("cohort")["family"].astype(str)
    metrics = [
        "best_module_Zsummary",
        "median_Zsummary",
        "median_density_ratio",
        "median_connectivity_correlation",
        "median_membership_profile_correlation",
        "median_program_score_preservation",
        "mapped_events",
        "tested_modules",
        "preserved_modules",
        "partially_preserved_modules",
        "preserved_or_partial_module_fraction",
    ]
    available = [column for column in metrics if column in external.columns]
    wide = external.pivot(
        index=["discovery_cohort", "external_cohort"],
        columns="method",
        values=available,
    )
    wide.columns = [f"{metric}__{method}" for metric, method in wide.columns]
    wide = wide.reset_index()
    wide["discovery_family"] = wide["discovery_cohort"].map(family)
    wide["external_family"] = wide["external_cohort"].map(family)
    wide["validation_scope"] = np.where(
        wide["discovery_family"].eq(wide["external_family"])
        & wide["discovery_family"].ne("other"),
        "organ_matched_primary",
        "cross_organ_stress",
    )
    for metric in metrics:
        left = f"{metric}__CoSplice"
        right = f"{metric}__Exact WGCNA"
        if left in wide and right in wide:
            wide[f"{metric}__CoSplice_minus_WGCNA"] = (
                pd.to_numeric(wide[left], errors="coerce")
                - pd.to_numeric(wide[right], errors="coerce")
            )
    return wide


def _scope_summary(paired: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for scope, group in paired.groupby("validation_scope", sort=False):
        z_delta = pd.to_numeric(
            group.get("best_module_Zsummary__CoSplice_minus_WGCNA"),
            errors="coerce",
        ).dropna()
        program_delta = pd.to_numeric(
            group.get("median_program_score_preservation__CoSplice_minus_WGCNA"),
            errors="coerce",
        ).dropna()
        median_z_delta = pd.to_numeric(
            group.get("median_Zsummary__CoSplice_minus_WGCNA"),
            errors="coerce",
        ).dropna()
        retained_fraction_delta = pd.to_numeric(
            group.get(
                "preserved_or_partial_module_fraction__CoSplice_minus_WGCNA"
            ),
            errors="coerce",
        ).dropna()
        nonzero_z = z_delta.loc[z_delta.ne(0)]
        z_wins = int((nonzero_z > 0).sum())
        sign_pvalue = (
            float(
                binomtest(
                    z_wins,
                    n=len(nonzero_z),
                    p=0.5,
                    alternative="greater",
                ).pvalue
            )
            if len(nonzero_z)
            else np.nan
        )
        try:
            wilcoxon_pvalue = (
                float(
                    wilcoxon(
                        z_delta,
                        alternative="greater",
                        zero_method="wilcox",
                    ).pvalue
                )
                if len(z_delta)
                else np.nan
            )
        except ValueError:
            wilcoxon_pvalue = np.nan
        if len(z_delta):
            rng = np.random.default_rng(20260730)
            bootstrap_medians = np.median(
                rng.choice(
                    z_delta.to_numpy(dtype=float),
                    size=(10_000, len(z_delta)),
                    replace=True,
                ),
                axis=1,
            )
            median_ci_low, median_ci_high = np.quantile(
                bootstrap_medians,
                [0.025, 0.975],
            )
        else:
            median_ci_low, median_ci_high = np.nan, np.nan
        rows.append(
            {
                "validation_scope": scope,
                "directed_projection_count": len(group),
                "paired_projection_count": len(z_delta),
                "fraction_CoSplice_higher_best_Zsummary": float((z_delta > 0).mean()),
                "median_best_Zsummary_difference": float(z_delta.median()),
                "median_best_Zsummary_difference_CI95_low": float(median_ci_low),
                "median_best_Zsummary_difference_CI95_high": float(median_ci_high),
                "best_Zsummary_one_sided_sign_test_pvalue": sign_pvalue,
                "best_Zsummary_one_sided_wilcoxon_pvalue": wilcoxon_pvalue,
                "fraction_CoSplice_higher_median_Zsummary": float(
                    (median_z_delta > 0).mean()
                ),
                "median_median_Zsummary_difference": float(median_z_delta.median()),
                "fraction_CoSplice_higher_program_score_preservation": float(
                    (program_delta > 0).mean()
                ),
                "median_program_score_preservation_difference": float(
                    program_delta.median()
                ),
                "fraction_CoSplice_higher_retained_module_fraction": float(
                    (retained_fraction_delta > 0).mean()
                ),
                "median_retained_module_fraction_difference": float(
                    retained_fraction_delta.median()
                ),
            }
        )
    return pd.DataFrame(rows)


def _functional_comparison(
    benchmark_dir: Path,
    registry: pd.DataFrame,
) -> pd.DataFrame:
    summary = _read_table(benchmark_dir / "biological_interpretability_summary.tsv")
    if summary.empty:
        return summary
    family = registry.set_index("cohort")["family"].astype(str)
    summary = summary.sort_values("method").reset_index(drop=True)
    pairwise_inputs = (
        (
            "cross_cohort_functional_theme_replication.tsv",
            "best_term_jaccard",
            "organ_matched_best_term_jaccard",
        ),
        (
            "cross_cohort_module_gene_replication.tsv",
            "best_gene_jaccard",
            "organ_matched_best_gene_jaccard",
        ),
    )
    for filename, metric, prefix in pairwise_inputs:
        frame = _read_table(benchmark_dir / filename)
        if frame.empty or metric not in frame:
            continue
        frame["discovery_family"] = frame["discovery_cohort"].map(family)
        frame["validation_family"] = frame["validation_cohort"].map(family)
        frame = frame.loc[
            frame["discovery_family"].eq(frame["validation_family"])
            & frame["discovery_family"].ne("other")
        ].copy()
        if frame.empty:
            continue
        pairwise = (
            frame.groupby(
                ["method", "discovery_cohort", "validation_cohort"],
                as_index=False,
            )[metric]
            .max()
            .groupby("method")[metric]
            .agg(["count", "median", "mean"])
            .reset_index()
            .rename(
                columns={
                    "count": f"{prefix}_projection_count",
                    "median": f"{prefix}_median",
                    "mean": f"{prefix}_mean",
                }
            )
        )
        summary = summary.merge(pairwise, on="method", how="left")
    return summary


def _claim_audit(
    internal: pd.DataFrame,
    scope: pd.DataFrame,
    functional: pd.DataFrame,
) -> pd.DataFrame:
    both_valid = internal.loc[internal["both_partitions_valid"]]
    stable_delta = pd.to_numeric(
        both_valid.get("min_bootstrap_ARI__CoSplice_minus_WGCNA"),
        errors="coerce",
    )
    primary = scope.loc[scope["validation_scope"].eq("organ_matched_primary")]
    best_module_preservation_pass = bool(
        not primary.empty
        and float(primary.iloc[0]["fraction_CoSplice_higher_best_Zsummary"]) >= 0.70
        and float(primary.iloc[0]["median_best_Zsummary_difference"]) > 0
        and float(
            primary.iloc[0]["best_Zsummary_one_sided_wilcoxon_pvalue"]
        )
        < 0.05
    )
    overall_preservation_pass = bool(
        not primary.empty
        and float(primary.iloc[0]["fraction_CoSplice_higher_median_Zsummary"])
        >= 0.70
        and float(primary.iloc[0]["median_median_Zsummary_difference"]) > 0
        and float(
            primary.iloc[0][
                "median_program_score_preservation_difference"
            ]
        )
        > 0
    )
    stability_pass = bool(
        len(stable_delta) >= 3
        and float((stable_delta > 0).mean()) >= 0.70
        and float(stable_delta.median()) > 0
    )
    cosplice_valid = int(
        pd.Series(
            internal.get("valid_partition__CoSplice", False),
            index=internal.index,
            dtype="boolean",
        )
        .fillna(False)
        .sum()
    )
    wgcna_valid = int(
        pd.Series(
            internal.get("valid_partition__Exact WGCNA", False),
            index=internal.index,
            dtype="boolean",
        )
        .fillna(False)
        .sum()
    )
    usable_partition_pass = bool(
        len(internal) >= 3
        and cosplice_valid > wgcna_valid
        and cosplice_valid / len(internal) >= 0.70
    )

    function_pass = False
    function_evidence = "Functional comparison unavailable."
    if not functional.empty and {"CoSplice", "Exact WGCNA"}.issubset(
        set(functional["method"].astype(str))
    ):
        indexed = functional.set_index("method")
        metric = "organ_matched_best_term_jaccard_median"
        count_metric = "organ_matched_best_term_jaccard_projection_count"
        if metric in indexed and count_metric in indexed:
            left = pd.to_numeric(indexed.at["CoSplice", metric], errors="coerce")
            right = pd.to_numeric(indexed.at["Exact WGCNA", metric], errors="coerce")
            left_n = pd.to_numeric(
                indexed.at["CoSplice", count_metric],
                errors="coerce",
            )
            right_n = pd.to_numeric(
                indexed.at["Exact WGCNA", count_metric],
                errors="coerce",
            )
            function_pass = bool(
                pd.notna(left)
                and pd.notna(right)
                and left_n >= 3
                and right_n >= 3
                and left > right
            )
            function_evidence = (
                "Organ-matched median best term Jaccard: "
                f"CoSplice={left} across {left_n} projections; "
                f"WGCNA={right} across {right_n} projections."
            )

    return pd.DataFrame(
        [
            {
                "claim": "CoSplice modules are more internally stable than exact WGCNA.",
                "status": "supported" if stability_pass else "not_supported",
                "evidence": (
                    f"{len(stable_delta)} cohorts had valid partitions for both methods; "
                    f"CoSplice won {int((stable_delta > 0).sum())}."
                ),
            },
            {
                "claim": "CoSplice produces usable multi-module partitions more consistently.",
                "status": "supported" if usable_partition_pass else "not_supported",
                "evidence": (
                    f"Valid partitions: CoSplice={cosplice_valid}/{len(internal)}; "
                    f"exact WGCNA={wgcna_valid}/{len(internal)}."
                ),
            },
            {
                "claim": (
                    "At least one CoSplice module has stronger organ-matched "
                    "external preservation in each paired projection."
                ),
                "status": (
                    "supported"
                    if best_module_preservation_pass
                    else "not_supported"
                ),
                "evidence": (
                    primary.to_dict(orient="records")[0]
                    if not primary.empty
                    else "No organ-matched primary projections."
                ),
            },
            {
                "claim": (
                    "The complete CoSplice module set has stronger organ-matched "
                    "external preservation than exact WGCNA."
                ),
                "status": (
                    "supported" if overall_preservation_pass else "not_supported"
                ),
                "evidence": (
                    primary.to_dict(orient="records")[0]
                    if not primary.empty
                    else "No organ-matched primary projections."
                ),
            },
            {
                "claim": "CoSplice modules have stronger functional-theme replication.",
                "status": "supported" if function_pass else "not_supported",
                "evidence": function_evidence,
            },
            {
                "claim": "Ten independent studies were validated.",
                "status": "not_tested",
                "evidence": (
                    "The panel contains ten independent TCGA disease cohorts processed "
                    "within one TCGA/SpliceSeq ecosystem, not ten independent studies."
                ),
            },
        ]
    )


def evaluate_cohort_panel(
    benchmark_dir: Path,
    config_path: Path,
    panel_manifest_path: Path,
) -> CohortPanelValidation:
    benchmark_dir = benchmark_dir.resolve()
    registry = _cohort_registry(config_path.resolve(), panel_manifest_path.resolve())
    internal = _read_table(benchmark_dir / "cohort_internal_benchmark.tsv")
    external = _read_table(
        benchmark_dir / "directed_projection_preservation_summary.tsv"
    )
    if internal.empty or external.empty:
        raise ValueError("The formal cohort benchmark is incomplete.")
    internal_comparison = _internal_comparison(internal)
    preservation = _preservation_comparison(external, registry)
    scope = _scope_summary(preservation)
    functional = _functional_comparison(benchmark_dir, registry)
    claims = _claim_audit(internal_comparison, scope, functional)
    return CohortPanelValidation(
        registry,
        internal_comparison,
        preservation,
        scope,
        functional,
        claims,
    )


def render_cohort_panel(
    result: CohortPanelValidation,
    output: Path,
) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(19, 6))
    internal = result.internal_comparison.copy()
    ari_columns = [
        "min_bootstrap_ARI__CoSplice",
        "min_bootstrap_ARI__Exact WGCNA",
    ]
    ari = internal.melt(
        id_vars="cohort",
        value_vars=[column for column in ari_columns if column in internal],
        var_name="method",
        value_name="minimum_bootstrap_ARI",
    )
    ari["method"] = ari["method"].str.split("__").str[-1]
    sns.barplot(
        data=ari,
        x="cohort",
        y="minimum_bootstrap_ARI",
        hue="method",
        hue_order=["CoSplice", "Exact WGCNA"],
        ax=axes[0],
    )
    method_offsets = {"CoSplice": -0.2, "Exact WGCNA": 0.2}
    for method in ("CoSplice", "Exact WGCNA"):
        valid_column = f"valid_partition__{method}"
        ari_column = f"min_bootstrap_ARI__{method}"
        if valid_column not in internal or ari_column not in internal:
            continue
        for cohort_index, row in internal.reset_index(drop=True).iterrows():
            valid_value = row[valid_column]
            if pd.notna(valid_value) and bool(valid_value):
                continue
            y_value = pd.to_numeric(row[ari_column], errors="coerce")
            if pd.isna(y_value):
                continue
            axes[0].scatter(
                cohort_index + method_offsets[method],
                float(y_value),
                marker="x",
                s=90,
                linewidth=2.2,
                color="#B91C1C",
                zorder=5,
            )
    axes[0].set_ylim(0, 1)
    axes[0].set_title("Internal stability")
    axes[0].tick_params(axis="x", rotation=90)
    axes[0].text(
        0.01,
        0.98,
        "Red X: invalid single/zero-module partition",
        transform=axes[0].transAxes,
        va="top",
        fontsize=9,
        color="#991B1B",
    )

    primary = result.preservation_comparison.loc[
        result.preservation_comparison["validation_scope"].eq(
            "organ_matched_primary"
        )
    ].copy()
    primary["pair"] = (
        primary["discovery_cohort"] + " -> " + primary["external_cohort"]
    )
    z = primary.melt(
        id_vars="pair",
        value_vars=[
            column
            for column in (
                "best_module_Zsummary__CoSplice",
                "best_module_Zsummary__Exact WGCNA",
            )
            if column in primary
        ],
        var_name="method",
        value_name="best_Zsummary",
    )
    z["method"] = z["method"].str.split("__").str[-1]
    sns.barplot(data=z, x="pair", y="best_Zsummary", hue="method", ax=axes[1])
    axes[1].axhline(2, color="#D97706", linestyle="--", linewidth=1)
    axes[1].set_title("Organ-matched external preservation")
    axes[1].tick_params(axis="x", rotation=90)

    scope = result.scope_summary.melt(
        id_vars="validation_scope",
        value_vars=[
            "fraction_CoSplice_higher_best_Zsummary",
            "fraction_CoSplice_higher_median_Zsummary",
            "fraction_CoSplice_higher_program_score_preservation",
            "fraction_CoSplice_higher_retained_module_fraction",
        ],
        var_name="metric",
        value_name="fraction",
    )
    scope["metric"] = scope["metric"].map(
        {
            "fraction_CoSplice_higher_best_Zsummary": "best module Zsummary",
            "fraction_CoSplice_higher_median_Zsummary": "median module Zsummary",
            "fraction_CoSplice_higher_program_score_preservation": (
                "program-score preservation"
            ),
            "fraction_CoSplice_higher_retained_module_fraction": (
                "retained-module fraction"
            ),
        }
    )
    sns.barplot(
        data=scope,
        x="validation_scope",
        y="fraction",
        hue="metric",
        ax=axes[2],
    )
    axes[2].axhline(0.5, color="#64748B", linestyle="--", linewidth=1)
    axes[2].set_ylim(0, 1)
    axes[2].set_title("Fraction of projections favoring CoSplice")
    axes[2].tick_params(axis="x", rotation=20)
    fig.suptitle(
        "Locked 10-cohort CoSplice validation",
        fontsize=17,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)
    return output


def write_cohort_panel_validation(
    result: CohortPanelValidation,
    output_dir: Path,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    tables = {
        "cohort_registry.tsv": result.cohort_registry,
        "internal_method_comparison.tsv": result.internal_comparison,
        "paired_preservation_comparison.tsv": result.preservation_comparison,
        "validation_scope_summary.tsv": result.scope_summary,
        "functional_comparison.tsv": result.functional_comparison,
        "claim_audit.tsv": result.claim_audit,
    }
    paths: list[Path] = []
    for name, frame in tables.items():
        path = output_dir / name
        frame.to_csv(path, sep="\t", index=False)
        paths.append(path)
    paths.append(
        render_cohort_panel(result, output_dir / "ten_cohort_validation_summary.png")
    )
    report_path = output_dir / "TEN_COHORT_VALIDATION_REPORT.md"
    valid_columns = [
        "valid_partition__CoSplice",
        "valid_partition__Exact WGCNA",
    ]
    valid_counts = {
        column.removeprefix("valid_partition__"): int(
            pd.Series(
                result.internal_comparison.get(
                    column,
                    pd.Series(index=result.internal_comparison.index, dtype="boolean"),
                ),
                index=result.internal_comparison.index,
                dtype="boolean",
            )
            .fillna(False)
            .astype(bool)
            .sum()
        )
        for column in valid_columns
    }
    primary = result.scope_summary.loc[
        result.scope_summary["validation_scope"].eq("organ_matched_primary")
    ]
    stress = result.scope_summary.loc[
        result.scope_summary["validation_scope"].eq("cross_organ_stress")
    ]

    def scope_line(frame: pd.DataFrame) -> str:
        if frame.empty:
            return "not available"
        row = frame.iloc[0]
        return (
            f"{int(row['directed_projection_count'])} directed projections; "
            f"CoSplice higher best Zsummary in "
            f"{float(row['fraction_CoSplice_higher_best_Zsummary']):.1%}; "
            f"median best-Z delta "
            f"{float(row['median_best_Zsummary_difference']):+.3f}; "
            f"median module-set Z delta "
            f"{float(row['median_median_Zsummary_difference']):+.3f}; "
            f"median program-score delta "
            f"{float(row['median_program_score_preservation_difference']):+.3f}"
        )

    claim_lines = [
        f"- **{row.claim}** `{row.status}`: {row.evidence}"
        for row in result.claim_audit.itertuples(index=False)
    ]
    cohort_lines = [
        (
            f"- {row.cohort}: family `{row.family}`, "
            f"{int(row.n_samples)} samples, {int(row.n_events)} eligible events"
        )
        for row in result.cohort_registry.itertuples(index=False)
    ]
    report = "\n".join(
        [
            "# Locked 10-cohort CoSplice validation",
            "",
            "## Scope",
            "",
            (
                "This benchmark compares CoSplice with exact R WGCNA on the same "
                "discovery-event universe. It retains negative and inconclusive results."
            ),
            (
                "The panel contains independent TCGA disease cohorts, not ten "
                "independent sequencing studies. Organ-matched projections are the "
                "primary preservation tests; cross-organ projections are stress tests."
            ),
            "",
            "## Cohorts",
            "",
            *cohort_lines,
            "",
            "## Locked design",
            "",
            "- 500 discovery events selected by variance within the discovery cohort only.",
            "- The external cohort retains every eligible event.",
            "- 20 matched bootstraps for internal stability.",
            "- 500 permutations for external module preservation.",
            "- Exact event-key matching; no fuzzy coordinate matching.",
            "- Functional annotation is evaluated after module discovery.",
            "",
            "## Result summary",
            "",
            (
                f"- Valid multi-module partitions: CoSplice "
                f"{valid_counts.get('CoSplice', 0)}/"
                f"{len(result.internal_comparison)}; exact WGCNA "
                f"{valid_counts.get('Exact WGCNA', 0)}/"
                f"{len(result.internal_comparison)}."
            ),
            f"- Organ-matched primary tests: {scope_line(primary)}.",
            f"- Cross-organ stress tests: {scope_line(stress)}.",
            "",
            "## Claim audit",
            "",
            *claim_lines,
            "",
            "## Interpretation boundary",
            "",
            (
                "A high bootstrap ARI from a single-module or zero-module partition is "
                "not treated as evidence of a useful stable program structure. "
                "Functional replication is support for biological interpretability, "
                "not proof of causal regulation."
            ),
            "",
        ]
    )
    report_path.write_text(report, encoding="utf-8")
    paths.append(report_path)
    return paths
