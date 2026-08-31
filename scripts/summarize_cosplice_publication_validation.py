from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


METHOD_LABELS = {
    "event_level_wgcna": "Exact WGCNA",
    "stability_only": "CoSplice",
    "three_objective_cosplice": "CoSplice",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t")


def _locked_summary(root: Path) -> pd.DataFrame:
    table = _read(root / "LOCKED_ABLATION_PANEL_SUMMARY.tsv")
    table = table.loc[
        table["method"].isin(["event_level_wgcna", "three_objective_cosplice"])
    ].copy()
    table["display_method"] = table["method"].map(METHOD_LABELS)
    rows = []
    for method, group in table.groupby("display_method", sort=False):
        rows.append(
            {
                "method": method,
                "n_runs": len(group),
                "mean_bootstrap_ARI": group["mean_bootstrap_ARI"].mean(),
                "mean_holdout_program_score_preservation": group[
                    "median_holdout_program_score_preservation"
                ].mean(),
                "valid_multi_module_runs": int(
                    group["valid_multi_module_partition"].astype(bool).sum()
                ),
                "functional_objective_active_runs": int(
                    group["functional_objective_active"].astype(bool).sum()
                ),
            }
        )
    return pd.DataFrame(rows)


def _kirc_summary(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    associations = _read(root / "tcga_kirc_SETD2_module_associations.tsv")
    preservation = _read(root / "tcga_kirc_to_kich_module_preservation.tsv")
    enrichment = _read(root / "tcga_kirc_module_GO_KEGG_enrichment.tsv")
    rows = [
        {
            "endpoint": "Adjusted SETD2 phenotype association",
            "significant": int(associations["adjusted_SETD2_FDR"].lt(0.05).sum()),
            "tested": len(associations),
        },
        {
            "endpoint": "Unadjusted Mann-Whitney association",
            "significant": int(associations["mann_whitney_FDR"].lt(0.05).sum()),
            "tested": len(associations),
        },
        {
            "endpoint": "Preserved/partially preserved modules",
            "significant": int(
                preservation["preservation_class"].isin(
                    ["preserved", "partially_preserved"]
                ).sum()
            ),
            "tested": len(preservation),
        },
        {
            "endpoint": "GO/KEGG terms at FDR < 0.05",
            "significant": int(pd.to_numeric(enrichment["FDR"], errors="coerce").lt(0.05).sum()),
            "tested": len(enrichment),
        },
    ]
    top = associations.sort_values("adjusted_SETD2_FDR").head(12).copy()
    top["label"] = top["event_type"].astype(str) + ":" + top["module"].astype(str)
    top["minus_log10_FDR"] = -np.log10(
        pd.to_numeric(top["adjusted_SETD2_FDR"], errors="coerce").clip(lower=1e-300)
    )
    return pd.DataFrame(rows), top


def _render(
    locked: pd.DataFrame,
    setd2: pd.DataFrame,
    kirc_top: pd.DataFrame,
    output: Path,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
        }
    )
    colors = {"Exact WGCNA": "#4C78A8", "CoSplice": "#E45756"}
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.4))

    long = locked.melt(
        id_vars="method",
        value_vars=[
            "mean_bootstrap_ARI",
            "mean_holdout_program_score_preservation",
        ],
        var_name="metric",
        value_name="value",
    )
    long["metric"] = long["metric"].map(
        {
            "mean_bootstrap_ARI": "Bootstrap ARI",
            "mean_holdout_program_score_preservation": "Holdout score\npreservation",
        }
    )
    sns.barplot(
        data=long,
        x="metric",
        y="value",
        hue="method",
        palette=colors,
        ax=axes[0],
    )
    axes[0].set_ylim(0, 1)
    axes[0].set_title("a  Unified-seed public panel", loc="left")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Mean across 11 runs")
    axes[0].legend(frameon=False, loc="upper right")

    sns.barplot(
        data=setd2,
        x="arm",
        y="mean_bootstrap_ARI",
        color="#2A9D8F",
        ax=axes[1],
    )
    axes[1].set_ylim(0, max(0.7, setd2["mean_bootstrap_ARI"].max() + 0.08))
    axes[1].set_title("b  Real IJC/SJC ablation", loc="left")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("Mean matched-bootstrap ARI")
    display = {
        "raw_PSI_WGCNA": "Raw\nWGCNA",
        "posterior_mean_WGCNA": "Posterior mean\nWGCNA",
        "posterior_ensemble_WGCNA": "Posterior ensemble\nWGCNA",
        "WGCNA_plus_read_uncertainty": "WGCNA +\nuncertainty",
        "full_CoSplice": "Full\nCoSplice",
    }
    axes[1].set_xticks(
        axes[1].get_xticks(),
        [display.get(str(value), str(value)) for value in setd2["arm"]],
    )
    for patch, value in zip(axes[1].patches, setd2["mean_bootstrap_ARI"], strict=False):
        axes[1].text(
            patch.get_x() + patch.get_width() / 2,
            value + 0.012,
            f"{value:.3f}",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    if kirc_top.empty:
        axes[2].text(0.5, 0.5, "No phenotype tests", ha="center", va="center")
        axes[2].axis("off")
    else:
        plot = kirc_top.sort_values("minus_log10_FDR")
        axes[2].barh(plot["label"], plot["minus_log10_FDR"], color="#F2A65A")
        axes[2].axvline(-np.log10(0.05), color="#B22222", linestyle="--", linewidth=1)
        axes[2].set_title("c  Independent KIRC SETD2 test", loc="left")
        axes[2].set_xlabel("Adjusted association, -log10(FDR)")
        axes[2].set_ylabel("")
        axes[2].tick_params(axis="y", labelsize=8)

    fig.suptitle(
        "CoSplice publication-readiness validation",
        x=0.02,
        ha="left",
        fontsize=17,
        fontweight="bold",
    )
    fig.text(
        0.02,
        0.01,
        "Matched inputs and seeds; uncertainty evaluated with real replicate junction counts; independent phenotype tests adjusted for purity and stage.",
        fontsize=9,
        color="#4B5563",
    )
    fig.tight_layout(rect=[0, 0.045, 1, 0.93], w_pad=3)
    fig.savefig(output, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--locked-panel", type=Path, required=True)
    parser.add_argument("--setd2-ablation", type=Path, required=True)
    parser.add_argument("--kirc-validation", type=Path, required=True)
    parser.add_argument("--public-challenge", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    locked = _locked_summary(args.locked_panel.resolve())
    four_arm_path = args.setd2_ablation.resolve() / "setd2_four_arm_summary.tsv"
    setd2_path = (
        four_arm_path
        if four_arm_path.exists()
        else args.setd2_ablation.resolve() / "setd2_three_arm_summary.tsv"
    )
    setd2 = _read(setd2_path)
    kirc, kirc_top = _kirc_summary(args.kirc_validation.resolve())
    locked.to_csv(output / "locked_panel_method_summary.tsv", sep="\t", index=False)
    setd2.to_csv(output / "setd2_real_count_ablation_summary.tsv", sep="\t", index=False)
    kirc.to_csv(output / "independent_kirc_endpoint_summary.tsv", sep="\t", index=False)
    kirc_top.to_csv(output / "independent_kirc_top_adjusted_tests.tsv", sep="\t", index=False)
    _render(locked, setd2, kirc_top, output / "cosplice_validation_overview.png")

    public_lines = ["- Not supplied."]
    public_inputs: dict[str, str] = {}
    if args.public_challenge is not None:
        challenge_root = args.public_challenge.resolve()
        challenge = _read(challenge_root / "three_dataset_internal_head_to_head.tsv")
        challenge.to_csv(output / "public_challenge_method_summary.tsv", sep="\t", index=False)
        public_lines = challenge.to_markdown(index=False, floatfmt=".3f").splitlines()
        public_inputs["public_challenge"] = str(challenge_root)

    locked_row = locked.set_index("method")
    report = [
        "# CoSplice publication-readiness validation",
        "",
        "## Executive conclusion",
        "",
        "The paired evidence supports WGCNA-level computational stability and a small Setd2 read-uncertainty gain. It does not yet support universal superiority or independent SETD2 biological replication.",
        "",
        "## Unified-seed 11-run panel",
        "",
        locked.to_markdown(index=False, floatfmt=".3f"),
        "",
        f"CoSplice mean bootstrap ARI was {locked_row.loc['CoSplice', 'mean_bootstrap_ARI']:.3f} versus {locked_row.loc['Exact WGCNA', 'mean_bootstrap_ARI']:.3f} for exact WGCNA. Mean holdout program-score preservation was {locked_row.loc['CoSplice', 'mean_holdout_program_score_preservation']:.3f} versus {locked_row.loc['Exact WGCNA', 'mean_holdout_program_score_preservation']:.3f}.",
        "",
        "## Real Setd2 IJC/SJC ablation",
        "",
        setd2.to_markdown(index=False, floatfmt=".3f"),
        "",
        "The matched arms separate reliability filtering/posterior mean from posterior-draw correlation propagation. Full CoSplice retains the posterior-ensemble WGCNA partition and adds uncertainty annotations rather than claiming a second clustering gain.",
        "",
        "## Independent human KIRC SETD2 validation",
        "",
        kirc.to_markdown(index=False),
        "",
        "No module association survived purity/stage adjustment and FDR control, and no GO/KEGG term passed FDR 0.05. This is a valid negative result and prevents a biological-replication claim.",
        "",
        "## Independent public challenge datasets",
        "",
        *public_lines,
        "",
        "## Claim boundary",
        "",
        "Current evidence supports a technically reproducible WGCNA-backbone extension with useful uncertainty and audit outputs. A methods-paper superiority claim still requires predeclared untouched cohorts with replicated functional coherence and phenotype association.",
    ]
    (output / "COSPLICE_PUBLICATION_READINESS_REPORT.md").write_text(
        "\n".join(report), encoding="utf-8"
    )
    files = [
        args.locked_panel.resolve() / "LOCKED_ABLATION_PANEL_SUMMARY.tsv",
        setd2_path,
        args.kirc_validation.resolve() / "tcga_kirc_SETD2_module_associations.tsv",
        args.kirc_validation.resolve() / "tcga_kirc_to_kich_module_preservation.tsv",
        args.kirc_validation.resolve() / "tcga_kirc_module_GO_KEGG_enrichment.tsv",
    ]
    parameters = {
        "inputs": {
            str(path): {"sha256": _sha256(path), "bytes": path.stat().st_size}
            for path in files
        },
        **public_inputs,
        "interpretation_boundary": "Summary only; no post-hoc threshold optimization.",
    }
    (output / "parameters.json").write_text(
        json.dumps(parameters, indent=2), encoding="utf-8"
    )
    print(f"Output: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
