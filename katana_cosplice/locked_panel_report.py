from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def _read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def collect_locked_runs(root: Path) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for summary_path in sorted(root.glob("*/locked_three_arm_summary.tsv")):
        run_dir = summary_path.parent
        plan_path = run_dir / "analysis_plan.json"
        freeze_path = run_dir / "freeze_manifest.json"
        search_path = (
            run_dir
            / "01_prelock"
            / "three_objective_cosplice"
            / "parameter_search.tsv"
        )
        if not plan_path.exists() or not freeze_path.exists():
            continue
        plan = _read_json(plan_path)
        freeze = _read_json(freeze_path)
        frame = pd.read_csv(summary_path, sep="\t")
        frame.insert(0, "run_id", run_dir.name)
        frame.insert(1, "event_type", str(plan.get("event_type", "ES")))
        frame.insert(2, "discovery", str(plan.get("discovery", "")))
        frame.insert(3, "tuning", ";".join(map(str, plan.get("tuning_cohorts", []))))
        frame.insert(4, "holdout", str(plan.get("holdout", "")))
        frame.insert(5, "n_frozen_events", int(freeze.get("n_frozen_events", 0)))
        frame["functional_objective_active"] = False
        if search_path.exists():
            search = pd.read_csv(search_path, sep="\t")
            selected = search.loc[search["selected"].astype(str).str.casefold().eq("true")]
            if not selected.empty and "functional_objective_active" in selected:
                frame["functional_objective_active"] = bool(
                    selected.iloc[0]["functional_objective_active"]
                )
        frame["valid_multi_module_partition"] = (
            pd.to_numeric(frame["n_modules"], errors="coerce").ge(2)
            & pd.to_numeric(
                frame["assigned_event_fraction"], errors="coerce"
            ).ge(0.50)
            & pd.to_numeric(
                frame["largest_module_fraction"], errors="coerce"
            ).le(0.60)
        )
        rows.append(frame)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def _render_panel(frame: pd.DataFrame, output: Path) -> Path:
    plot = frame.copy()
    method_labels = {
        "stability_only": "CoSplice: stability only",
        "three_objective_cosplice": "CoSplice: adaptive multi-objective",
        "event_level_wgcna": "Event-level WGCNA",
    }
    method_order = list(method_labels.values())
    plot["method_display"] = plot["method"].map(method_labels).fillna(plot["method"])
    plot["label"] = (
        plot["event_type"].astype(str)
        + " / n="
        + plot["n_frozen_events"].astype(str)
        + " / "
        + plot["holdout"].astype(str)
    )
    order = list(dict.fromkeys(plot["label"]))
    fig, axes = plt.subplots(2, 2, figsize=(18, max(10, 0.55 * len(order) + 5)))
    palette = ["#2878B5", "#5A9F68", "#E07A1F"]
    panels = [
        ("mean_bootstrap_ARI", "Internal bootstrap stability", "Mean adjusted Rand index"),
        ("median_holdout_Zsummary", "External module preservation", "Median Zsummary-like"),
        ("n_modules", "Recovered module granularity", "Non-grey modules"),
        (
            "median_holdout_program_score_preservation",
            "Program-score preservation",
            "Median preservation score",
        ),
    ]
    for axis, (column, title, xlabel) in zip(axes.flat, panels, strict=True):
        sns.barplot(
            data=plot,
            y="label",
            x=column,
            hue="method_display",
            hue_order=method_order,
            palette=palette,
            order=order,
            ax=axis,
        )
        axis.set_title(title, loc="left", fontweight="bold")
        axis.set_xlabel(xlabel)
        axis.set_ylabel("")
        legend = axis.get_legend()
        if legend is not None:
            legend.remove()
        if column == "median_holdout_Zsummary":
            axis.axvline(2, color="#D97706", linestyle="--", linewidth=1)
            axis.axvline(10, color="#0F766E", linestyle="--", linewidth=1)
        sns.despine(ax=axis)
    fig.suptitle(
        "Locked CoSplice validation across event universes and holdouts",
        fontsize=16,
        fontweight="bold",
    )
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.005),
    )
    fig.tight_layout(rect=[0, 0.055, 1, 0.97])
    path = output / "locked_ablation_panel.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def _markdown_table(frame: pd.DataFrame) -> str:
    display = frame.copy()
    for column in display.select_dtypes(include=["number"]).columns:
        display[column] = display[column].map(
            lambda value: "" if pd.isna(value) else f"{float(value):.3f}"
        )
    rows = [
        "| " + " | ".join(map(str, display.columns)) + " |",
        "| " + " | ".join(["---"] * len(display.columns)) + " |",
    ]
    rows.extend(
        "| " + " | ".join(map(str, row)) + " |"
        for row in display.itertuples(index=False, name=None)
    )
    return "\n".join(rows)


def build_locked_panel_report(root: Path) -> list[Path]:
    root = root.resolve()
    frame = collect_locked_runs(root)
    if frame.empty:
        raise ValueError(f"No completed locked ablation runs found under {root}")
    summary_path = root / "LOCKED_ABLATION_PANEL_SUMMARY.tsv"
    frame.to_csv(summary_path, sep="\t", index=False)
    plot_path = _render_panel(frame, root)
    compact_columns = [
        "run_id",
        "event_type",
        "n_frozen_events",
        "holdout",
        "method",
        "n_modules",
        "mean_bootstrap_ARI",
        "median_holdout_Zsummary",
        "median_holdout_program_score_preservation",
        "valid_multi_module_partition",
        "functional_objective_active",
    ]
    compact = frame[compact_columns]
    cosplice = frame.loc[frame["method"].eq("stability_only")]
    wgcna = frame.loc[frame["method"].eq("event_level_wgcna")]
    functional_active_runs = int(
        frame.groupby("run_id")["functional_objective_active"].max().sum()
    )
    cosplice_mean_ari = float(cosplice["mean_bootstrap_ARI"].mean())
    wgcna_mean_ari = float(wgcna["mean_bootstrap_ARI"].mean())
    cosplice_mean_program = float(
        cosplice["median_holdout_program_score_preservation"].mean()
    )
    wgcna_mean_program = float(
        wgcna["median_holdout_program_score_preservation"].mean()
    )
    report_path = root / "LOCKED_ABLATION_PANEL_REPORT.md"
    report_path.write_text(
        "\n".join(
            [
                "# Locked CoSplice validation panel",
                "",
                "## Scope",
                "",
                f"- Completed locked runs: {frame['run_id'].nunique()}",
                f"- Event types: {', '.join(sorted(frame['event_type'].unique()))}",
                f"- Holdouts: {', '.join(sorted(frame['holdout'].unique()))}",
                "- Exact WGCNA is run on the identical discovery-selected event universe.",
                "- KICH is the same-organ kidney holdout; LUAD/LUSC are transportability challenges.",
                "",
                "## Primary conclusions",
                "",
                f"- CoSplice valid multi-module partitions: {int(cosplice['valid_multi_module_partition'].sum())}/{len(cosplice)} runs.",
                f"- Exact WGCNA valid multi-module partitions: {int(wgcna['valid_multi_module_partition'].sum())}/{len(wgcna)} runs.",
                f"- Runs in which formal functional evidence activated the optimization weight: {functional_active_runs}/{frame['run_id'].nunique()}.",
                f"- Mean bootstrap ARI: CoSplice {cosplice_mean_ari:.3f}; exact WGCNA {wgcna_mean_ari:.3f}.",
                f"- Mean holdout program-score preservation: CoSplice {cosplice_mean_program:.3f}; exact WGCNA {wgcna_mean_program:.3f}.",
                "- Internal stability, module granularity and external preservation are reported separately; no fused score is used to claim universal superiority.",
                "- The adaptive multi-objective arm equals stability-only whenever no module passes the predeclared functional FDR gate; this is an intentional negative-control result, not a missing analysis.",
                "- Zero GO/KEGG replication is retained as a negative result and blocks a functional-superiority claim.",
                "",
                "## Run-level results",
                "",
                _markdown_table(compact),
                "",
                "## Interpretation boundary",
                "",
                "The TCGA cohorts are independent patient cohorts but share one processing ecosystem. Cross-cancer preservation tests transportability, not same-disease replication. Event-type sensitivity shows whether the architecture generalizes, but each event class may require separately frozen network hyperparameters. A method is not called superior solely because a coarse or single-module partition has a high ARI.",
                "",
                "## Reproducibility",
                "",
                "Every child run contains a pre-holdout freeze manifest, input and assignment SHA256 fingerprints, holdout access audit, method-specific assignments, local enrichment tables and a standalone report.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    parameters_path = root / "locked_panel_parameters.json"
    parameters_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "input_root": str(root),
                "completed_run_ids": sorted(frame["run_id"].unique()),
                "generated_files": [
                    str(summary_path),
                    str(plot_path),
                    str(report_path),
                ],
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return [summary_path, plot_path, report_path, parameters_path]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Aggregate completed frozen CoSplice ablations into one release report."
    )
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args(argv)
    for path in build_locked_panel_report(args.root):
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
