from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import adjusted_rand_score

matplotlib.use("Agg")
import matplotlib.pyplot as plt


@dataclass
class HybridBenchmarkResult:
    event_assignments: pd.DataFrame
    branch_summary: pd.DataFrame
    tight_diffuse_overlap: pd.DataFrame
    tight_incremental_programs: pd.DataFrame
    architecture_summary: pd.DataFrame
    files: list[Path]


def _program_sets(labels: pd.Series, prefix: str) -> dict[str, set[str]]:
    values = labels.fillna("grey").astype(str)
    return {
        f"{prefix}{module}": set(values.index[values.eq(module)].astype(str))
        for module in sorted(set(values).difference({"grey"}))
    }


def _overlap_table(
    tight_labels: pd.Series,
    diffuse_labels: pd.Series,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    tight_sets = _program_sets(tight_labels, "T-")
    diffuse_sets = _program_sets(diffuse_labels, "D-")
    for tight_program, tight_events in tight_sets.items():
        for diffuse_program, diffuse_events in diffuse_sets.items():
            shared = len(tight_events & diffuse_events)
            union = len(tight_events | diffuse_events)
            rows.append(
                {
                    "tight_program": tight_program,
                    "diffuse_program": diffuse_program,
                    "shared_events": shared,
                    "jaccard": shared / union if union else 0.0,
                    "tight_fraction_shared": shared / len(tight_events),
                    "diffuse_fraction_shared": shared / len(diffuse_events),
                }
            )
    return pd.DataFrame(rows)


def _method_row(summary: pd.DataFrame, method: str) -> pd.Series:
    rows = summary.loc[summary["method"].astype(str).eq(method)]
    if rows.empty:
        return pd.Series(dtype=object)
    return rows.iloc[0]


def build_hybrid_benchmark(
    assignments: pd.DataFrame,
    method_summary: pd.DataFrame,
) -> HybridBenchmarkResult:
    """Audit the production tight/diffuse architecture without inventing one fused partition."""
    required = {"event_key", "module_katana_sparse", "module_wgcna_exact"}
    missing = required.difference(assignments.columns)
    if missing:
        raise ValueError(
            "Hybrid benchmark requires same-universe tight and exact-WGCNA assignments: "
            + ", ".join(sorted(missing))
        )

    indexed = assignments.loc[:, sorted(required)].copy()
    indexed["event_key"] = indexed["event_key"].astype(str)
    indexed = indexed.drop_duplicates("event_key").set_index("event_key")
    tight = indexed["module_katana_sparse"].fillna("grey").astype(str)
    diffuse = indexed["module_wgcna_exact"].fillna("grey").astype(str)
    event_assignments = pd.DataFrame(
        {
            "event_key": indexed.index,
            "tight_program": np.where(tight.ne("grey"), "T-" + tight, "grey"),
            "diffuse_program": np.where(diffuse.ne("grey"), "D-" + diffuse, "grey"),
        }
    )
    tight_retained = tight.ne("grey").to_numpy()
    diffuse_retained = diffuse.ne("grey").to_numpy()
    event_assignments["program_scope"] = np.select(
        [
            tight_retained & diffuse_retained,
            tight_retained,
            diffuse_retained,
        ],
        ["both", "tight_only", "diffuse_only"],
        default="unassigned",
    )

    overlap = _overlap_table(tight, diffuse).reindex(
        columns=[
            "tight_program",
            "diffuse_program",
            "shared_events",
            "jaccard",
            "tight_fraction_shared",
            "diffuse_fraction_shared",
        ]
    )
    tight_sets = _program_sets(tight, "T-")
    incremental_rows: list[dict[str, object]] = []
    for tight_program, tight_events in tight_sets.items():
        module_overlap = overlap.loc[overlap["tight_program"].eq(tight_program)]
        best = (
            module_overlap.sort_values(
                ["jaccard", "shared_events"],
                ascending=False,
            ).iloc[0]
            if not module_overlap.empty
            else pd.Series(dtype=object)
        )
        tight_only = set(
            event_assignments.loc[
                event_assignments["tight_program"].eq(tight_program)
                & event_assignments["program_scope"].eq("tight_only"),
                "event_key",
            ]
        )
        incremental_rows.append(
            {
                "tight_program": tight_program,
                "n_events": len(tight_events),
                "n_tight_only_events": len(tight_only),
                "tight_only_fraction": len(tight_only) / len(tight_events),
                "best_matching_diffuse_program": best.get("diffuse_program", ""),
                "best_diffuse_jaccard": best.get("jaccard", np.nan),
                "best_diffuse_shared_events": best.get("shared_events", 0),
                "evidence_scope": (
                    "Same-input membership increment. Biological value requires bootstrap "
                    "confidence and independent-cohort preservation."
                ),
            }
        )
    incremental = pd.DataFrame(incremental_rows).reindex(
        columns=[
            "tight_program",
            "n_events",
            "n_tight_only_events",
            "tight_only_fraction",
            "best_matching_diffuse_program",
            "best_diffuse_jaccard",
            "best_diffuse_shared_events",
            "evidence_scope",
        ]
    )

    n_events = len(event_assignments)
    scope_counts = event_assignments["program_scope"].value_counts()
    tight_row = _method_row(method_summary, "katana_sparse")
    diffuse_row = _method_row(method_summary, "wgcna_exact")
    branch_rows: list[dict[str, object]] = []
    for branch, row, labels, role in (
        (
            "tight",
            tight_row,
            tight,
            "Sparse signed CoSplice branch for compact local programs.",
        ),
        (
            "diffuse",
            diffuse_row,
            diffuse,
            "Exact WGCNA soft-adjacency/TOM branch for broad programs.",
        ),
    ):
        branch_rows.append(
            {
                "branch": branch,
                "method": row.get("method", ""),
                "role": role,
                "n_modules": int(len(set(labels).difference({"grey"}))),
                "retained_events": int(labels.ne("grey").sum()),
                "retained_fraction": float(labels.ne("grey").mean()),
                "grey_fraction": float(labels.eq("grey").mean()),
                "mean_bootstrap_ARI": row.get("mean_bootstrap_ARI", np.nan),
                "min_bootstrap_ARI": row.get("min_bootstrap_ARI", np.nan),
                "silhouette_cosine": row.get("silhouette_cosine", np.nan),
                "within_between_separation": row.get(
                    "within_between_separation",
                    np.nan,
                ),
                "runtime_seconds": row.get("runtime_seconds", np.nan),
            }
        )
    branch_summary = pd.DataFrame(branch_rows)

    architecture_summary = pd.DataFrame(
        [
            {
                "metric": "n_events",
                "value": n_events,
                "interpretation": "Shared event universe used by both branches.",
            },
            {
                "metric": "both_fraction",
                "value": scope_counts.get("both", 0) / n_events,
                "interpretation": "Events retained by both tight and diffuse branches.",
            },
            {
                "metric": "tight_only_fraction",
                "value": scope_counts.get("tight_only", 0) / n_events,
                "interpretation": (
                    "Events added by the sparse branch beyond exact WGCNA membership."
                ),
            },
            {
                "metric": "diffuse_only_fraction",
                "value": scope_counts.get("diffuse_only", 0) / n_events,
                "interpretation": "Broad-program events retained only by exact WGCNA/TOM.",
            },
            {
                "metric": "unassigned_fraction",
                "value": scope_counts.get("unassigned", 0) / n_events,
                "interpretation": "Events assigned to neither branch.",
            },
            {
                "metric": "tight_vs_diffuse_adjusted_rand_index",
                "value": adjusted_rand_score(tight, diffuse),
                "interpretation": (
                    "Agreement between the two complementary partitions; lower values may reflect "
                    "different granularity and are not automatically a defect."
                ),
            },
            {
                "metric": "diffuse_exact_wgcna_equivalence",
                "value": 1.0,
                "interpretation": (
                    "The benchmark diffuse labels are the exact-WGCNA labels by construction. "
                    "This is an implementation contract, not evidence of superiority."
                ),
            },
        ]
    )
    return HybridBenchmarkResult(
        event_assignments=event_assignments,
        branch_summary=branch_summary,
        tight_diffuse_overlap=overlap,
        tight_incremental_programs=incremental,
        architecture_summary=architecture_summary,
        files=[],
    )


def render_hybrid_benchmark_plots(
    result: HybridBenchmarkResult,
    output_dir: Path,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    generated: list[Path] = []

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5))
    order = ["both", "tight_only", "diffuse_only", "unassigned"]
    colors = {
        "both": "#0F766E",
        "tight_only": "#2563EB",
        "diffuse_only": "#D97706",
        "unassigned": "#CBD5E1",
    }
    counts = (
        result.event_assignments["program_scope"]
        .value_counts()
        .reindex(order, fill_value=0)
    )
    axes[0].bar(counts.index, counts.values, color=[colors[item] for item in counts.index])
    axes[0].set_title("Hybrid event coverage", fontweight="bold")
    axes[0].set_ylabel("Events")
    axes[0].tick_params(axis="x", rotation=25)

    branches = result.branch_summary.set_index("branch")
    metric_columns = [
        column
        for column in ["mean_bootstrap_ARI", "silhouette_cosine", "within_between_separation"]
        if column in branches and branches[column].notna().any()
    ]
    if metric_columns:
        branch_metrics = branches[metric_columns].apply(pd.to_numeric, errors="coerce")
        sns.heatmap(
            branch_metrics,
            annot=True,
            fmt=".3f",
            cmap="YlGnBu",
            ax=axes[1],
            cbar_kws={"label": "Raw metric"},
        )
        axes[1].set_title("Branch-level quality", fontweight="bold")
        axes[1].set_xlabel("")
        axes[1].set_ylabel("")
    else:
        axes[1].axis("off")

    if not result.tight_diffuse_overlap.empty:
        overlap = result.tight_diffuse_overlap.pivot(
            index="tight_program",
            columns="diffuse_program",
            values="jaccard",
        )
        sns.heatmap(
            overlap,
            cmap="mako",
            vmin=0,
            vmax=1,
            ax=axes[2],
            cbar_kws={"label": "Event-set Jaccard"},
        )
        axes[2].set_title("Tight-to-diffuse overlap", fontweight="bold")
        axes[2].set_xlabel("Diffuse WGCNA program")
        axes[2].set_ylabel("Tight CoSplice program")
    else:
        axes[2].axis("off")
    fig.suptitle(
        "Katana-CoSplice hybrid architecture audit",
        fontsize=16,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    path = output_dir / "hybrid_architecture_audit.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    generated.append(path)
    return generated


def write_hybrid_benchmark(
    result: HybridBenchmarkResult,
    output_dir: Path,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    tables = {
        "hybrid_event_program_assignments.tsv": result.event_assignments,
        "hybrid_branch_summary.tsv": result.branch_summary,
        "hybrid_tight_diffuse_overlap.tsv": result.tight_diffuse_overlap,
        "hybrid_tight_incremental_programs.tsv": result.tight_incremental_programs,
        "hybrid_architecture_summary.tsv": result.architecture_summary,
    }
    files: list[Path] = []
    for name, frame in tables.items():
        path = output_dir / name
        frame.to_csv(path, sep="\t", index=False)
        files.append(path)
    files.extend(render_hybrid_benchmark_plots(result, output_dir / "plots"))
    result.files = files
    return files
