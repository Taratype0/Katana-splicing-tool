from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats


@dataclass(frozen=True)
class MatchedBackgroundCoverageResult:
    event_catalog: pd.DataFrame
    inhibitor_summary: pd.DataFrame
    coverage_class_summary: pd.DataFrame


@dataclass(frozen=True)
class SignificantEventVennResult:
    membership_catalog: pd.DataFrame
    region_summary: pd.DataFrame


def _prepare_event_table(
    frame: pd.DataFrame,
    *,
    fdr_cutoff: float,
    abs_dpsi_cutoff: float,
) -> pd.DataFrame:
    required = {"event_key", "FDR", "dPSI"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Event table is missing required columns: {missing}")

    output = frame.copy()
    output["event_key"] = output["event_key"].astype(str)
    output["FDR"] = pd.to_numeric(output["FDR"], errors="coerce")
    output["dPSI"] = pd.to_numeric(output["dPSI"], errors="coerce")
    output = (
        output.sort_values(["event_key", "FDR"], na_position="last")
        .drop_duplicates("event_key", keep="first")
        .reset_index(drop=True)
    )
    output["significant"] = output["FDR"].le(fdr_cutoff) & output["dPSI"].abs().ge(abs_dpsi_cutoff)
    return output


def build_matched_background_mut_coverage(
    mut_events: pd.DataFrame,
    inhibitor_events: dict[str, pd.DataFrame],
    *,
    context: str,
    fdr_cutoff: float = 0.05,
    abs_dpsi_cutoff: float = 0.10,
) -> MatchedBackgroundCoverageResult:
    """Measure how much of a MUT AS signature is covered by matched-background inhibitors.

    The MUT significant-event set is the fixed denominator. Each inhibitor table must
    already use the same canonical direction as MUT, for example treatment minus
    matched WT. Event identity is exact and supplied by ``event_key``.
    """

    if not inhibitor_events:
        raise ValueError("At least one inhibitor event table is required.")

    mut = _prepare_event_table(
        mut_events,
        fdr_cutoff=fdr_cutoff,
        abs_dpsi_cutoff=abs_dpsi_cutoff,
    )
    mut_significant = mut.loc[mut["significant"]].copy()
    if mut_significant.empty:
        empty_summary = pd.DataFrame(
            columns=[
                "context",
                "inhibitor",
                "n_MUT_significant",
                "n_inhibitor_significant",
                "n_MUT_events_tested_by_inhibitor",
                "n_MUT_events_covered",
                "MUT_coverage_pct",
                "n_same_direction",
                "same_direction_pct",
                "n_opposite_direction",
                "median_abs_delta_dPSI",
                "same_direction_binomial_pvalue",
                "dPSI_spearman_rho",
                "dPSI_spearman_pvalue",
            ]
        )
        return MatchedBackgroundCoverageResult(
            event_catalog=mut_significant,
            inhibitor_summary=empty_summary,
            coverage_class_summary=pd.DataFrame(columns=["context", "coverage_class", "n_events", "pct_MUT_signature"]),
        )

    metadata_columns = [
        column
        for column in [
            "event_key",
            "gene_id",
            "gene_symbol",
            "event_type",
            "chr",
            "strand",
            "event_coordinates",
        ]
        if column in mut_significant.columns
    ]
    catalog = mut_significant.loc[:, metadata_columns + ["FDR", "dPSI"]].rename(
        columns={"FDR": "MUT_FDR", "dPSI": "MUT_dPSI"}
    )
    catalog.insert(0, "context", str(context))

    summaries: list[dict[str, object]] = []
    inhibitor_names = list(inhibitor_events)
    for inhibitor_name, inhibitor_frame in inhibitor_events.items():
        inhibitor = _prepare_event_table(
            inhibitor_frame,
            fdr_cutoff=fdr_cutoff,
            abs_dpsi_cutoff=abs_dpsi_cutoff,
        )
        prefix = str(inhibitor_name)
        values = inhibitor.loc[:, ["event_key", "FDR", "dPSI", "significant"]].rename(
            columns={
                "FDR": f"{prefix}_FDR",
                "dPSI": f"{prefix}_dPSI",
                "significant": f"{prefix}_significant",
            }
        )
        catalog = catalog.merge(values, on="event_key", how="left", validate="one_to_one")
        catalog[f"{prefix}_tested"] = catalog[f"{prefix}_dPSI"].notna()
        catalog[f"{prefix}_significant"] = catalog[f"{prefix}_significant"].eq(True)
        catalog[f"{prefix}_covers_MUT"] = catalog[f"{prefix}_significant"]
        catalog[f"{prefix}_same_direction"] = (
            catalog["MUT_dPSI"].mul(catalog[f"{prefix}_dPSI"]).gt(0)
            & catalog[f"{prefix}_covers_MUT"]
        )
        catalog[f"{prefix}_opposite_direction"] = (
            catalog["MUT_dPSI"].mul(catalog[f"{prefix}_dPSI"]).lt(0)
            & catalog[f"{prefix}_covers_MUT"]
        )
        catalog[f"{prefix}_delta_dPSI"] = catalog[f"{prefix}_dPSI"] - catalog["MUT_dPSI"]
        catalog[f"{prefix}_abs_delta_dPSI"] = catalog[f"{prefix}_delta_dPSI"].abs()

        covered = catalog.loc[catalog[f"{prefix}_covers_MUT"]]
        same_direction = int(covered[f"{prefix}_same_direction"].sum())
        opposite_direction = int(covered[f"{prefix}_opposite_direction"].sum())
        n_covered = len(covered)
        direction_pvalue = (
            float(stats.binomtest(same_direction, n_covered, p=0.5, alternative="greater").pvalue)
            if n_covered
            else np.nan
        )
        if n_covered >= 3:
            correlation = stats.spearmanr(
                covered["MUT_dPSI"],
                covered[f"{prefix}_dPSI"],
                nan_policy="omit",
            )
            spearman_rho = float(correlation.statistic)
            spearman_pvalue = float(correlation.pvalue)
        else:
            spearman_rho = np.nan
            spearman_pvalue = np.nan
        summaries.append(
            {
                "context": str(context),
                "inhibitor": prefix,
                "n_MUT_significant": len(catalog),
                "n_inhibitor_significant": int(inhibitor["significant"].sum()),
                "n_MUT_events_tested_by_inhibitor": int(catalog[f"{prefix}_tested"].sum()),
                "n_MUT_events_covered": n_covered,
                "MUT_coverage_pct": 100.0 * n_covered / len(catalog),
                "n_same_direction": same_direction,
                "same_direction_pct": 100.0 * same_direction / n_covered if n_covered else np.nan,
                "n_opposite_direction": opposite_direction,
                "median_abs_delta_dPSI": (
                    float(covered[f"{prefix}_abs_delta_dPSI"].median()) if n_covered else np.nan
                ),
                "same_direction_binomial_pvalue": direction_pvalue,
                "dPSI_spearman_rho": spearman_rho,
                "dPSI_spearman_pvalue": spearman_pvalue,
            }
        )

    coverage_flags = [f"{name}_covers_MUT" for name in inhibitor_names]
    if len(inhibitor_names) == 2:
        first, second = inhibitor_names
        catalog["coverage_class"] = np.select(
            [
                catalog[f"{first}_covers_MUT"] & catalog[f"{second}_covers_MUT"],
                catalog[f"{first}_covers_MUT"],
                catalog[f"{second}_covers_MUT"],
            ],
            ["covered_by_both", f"{first}_only", f"{second}_only"],
            default="MUT_only",
        )
        both = catalog.loc[catalog["coverage_class"].eq("covered_by_both")]
        all_same = (
            np.sign(both["MUT_dPSI"]).eq(np.sign(both[f"{first}_dPSI"]))
            & np.sign(both["MUT_dPSI"]).eq(np.sign(both[f"{second}_dPSI"]))
        )
        summaries.append(
            {
                "context": str(context),
                "inhibitor": f"{first}+{second}",
                "n_MUT_significant": len(catalog),
                "n_inhibitor_significant": np.nan,
                "n_MUT_events_tested_by_inhibitor": int(
                    (catalog[f"{first}_tested"] & catalog[f"{second}_tested"]).sum()
                ),
                "n_MUT_events_covered": len(both),
                "MUT_coverage_pct": 100.0 * len(both) / len(catalog),
                "n_same_direction": int(all_same.sum()),
                "same_direction_pct": 100.0 * float(all_same.mean()) if len(both) else np.nan,
                "n_opposite_direction": int((~all_same).sum()),
                "median_abs_delta_dPSI": np.nan,
                "same_direction_binomial_pvalue": np.nan,
                "dPSI_spearman_rho": np.nan,
                "dPSI_spearman_pvalue": np.nan,
            }
        )
    else:
        catalog["coverage_class"] = catalog[coverage_flags].apply(
            lambda row: "+".join(
                name for name, covered in zip(inhibitor_names, row, strict=True) if bool(covered)
            )
            or "MUT_only",
            axis=1,
        )

    class_summary = (
        catalog.groupby("coverage_class", dropna=False)
        .size()
        .rename("n_events")
        .reset_index()
    )
    class_summary.insert(0, "context", str(context))
    class_summary["pct_MUT_signature"] = 100.0 * class_summary["n_events"] / len(catalog)

    catalog["fdr_cutoff_used"] = float(fdr_cutoff)
    catalog["abs_dPSI_cutoff_used"] = float(abs_dpsi_cutoff)
    return MatchedBackgroundCoverageResult(
        event_catalog=catalog,
        inhibitor_summary=pd.DataFrame(summaries),
        coverage_class_summary=class_summary,
    )


def build_significant_event_venn_summary(
    event_tables: dict[str, pd.DataFrame],
    *,
    context: str,
    fdr_cutoff: float = 0.05,
    abs_dpsi_cutoff: float = 0.10,
) -> SignificantEventVennResult:
    """Build exact membership regions for three significant rMATS event sets."""

    expected = ["MUT", "CLK3", "Pla"]
    if list(event_tables) != expected:
        raise ValueError(f"event_tables must be ordered as {expected}.")

    significant_keys: dict[str, set[str]] = {}
    metadata_frames: list[pd.DataFrame] = []
    for label, frame in event_tables.items():
        prepared = _prepare_event_table(
            frame,
            fdr_cutoff=fdr_cutoff,
            abs_dpsi_cutoff=abs_dpsi_cutoff,
        )
        significant = prepared.loc[prepared["significant"]].copy()
        significant_keys[label] = set(significant["event_key"])
        metadata_columns = [
            column
            for column in [
                "event_key",
                "gene_id",
                "gene_symbol",
                "event_type",
                "chr",
                "strand",
                "event_coordinates",
            ]
            if column in significant.columns
        ]
        metadata_frames.append(significant.loc[:, metadata_columns].copy())

    union_keys = set().union(*significant_keys.values())
    metadata_catalog = (
        pd.concat(metadata_frames, ignore_index=True)
        .drop_duplicates("event_key", keep="first")
        .set_index("event_key")
        .reindex(sorted(union_keys))
        .reset_index()
    )
    membership = metadata_catalog.copy()
    membership.insert(0, "context", str(context))
    for label in expected:
        membership[f"{label}_significant"] = membership["event_key"].isin(
            significant_keys[label]
        )

    membership["membership_region"] = np.select(
        [
            membership["MUT_significant"]
            & ~membership["CLK3_significant"]
            & ~membership["Pla_significant"],
            ~membership["MUT_significant"]
            & membership["CLK3_significant"]
            & ~membership["Pla_significant"],
            ~membership["MUT_significant"]
            & ~membership["CLK3_significant"]
            & membership["Pla_significant"],
            membership["MUT_significant"]
            & membership["CLK3_significant"]
            & ~membership["Pla_significant"],
            membership["MUT_significant"]
            & ~membership["CLK3_significant"]
            & membership["Pla_significant"],
            ~membership["MUT_significant"]
            & membership["CLK3_significant"]
            & membership["Pla_significant"],
            membership["MUT_significant"]
            & membership["CLK3_significant"]
            & membership["Pla_significant"],
        ],
        [
            "MUT_only",
            "CLK3_only",
            "PlaB_only",
            "MUT_CLK3_only",
            "MUT_PlaB_only",
            "CLK3_PlaB_only",
            "MUT_CLK3_PlaB",
        ],
        default="unclassified",
    )

    region_order = [
        "MUT_only",
        "CLK3_only",
        "PlaB_only",
        "MUT_CLK3_only",
        "MUT_PlaB_only",
        "CLK3_PlaB_only",
        "MUT_CLK3_PlaB",
    ]
    counts = membership["membership_region"].value_counts()
    totals = {label: len(significant_keys[label]) for label in expected}
    summary = pd.DataFrame(
        {
            "context": str(context),
            "region": region_order,
            "n_events": [int(counts.get(region, 0)) for region in region_order],
            "n_MUT_significant_total": totals["MUT"],
            "n_CLK3_significant_total": totals["CLK3"],
            "n_PlaB_significant_total": totals["Pla"],
            "fdr_cutoff_used": float(fdr_cutoff),
            "abs_dPSI_cutoff_used": float(abs_dpsi_cutoff),
        }
    )
    return SignificantEventVennResult(
        membership_catalog=membership,
        region_summary=summary,
    )
