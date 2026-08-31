from __future__ import annotations

import argparse
import hashlib
import html
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from katana_cosplice.biological_interpretability import (
    evaluate_biological_interpretability,
)
from katana_cosplice.config import (
    ConsensusConfig,
    ExternalCohortConfig,
    NetworkConfig,
    PreservationConfig,
)
from katana_cosplice.consensus import build_bootstrap_coassignment_consensus
from katana_cosplice.independent_validation import (
    fit_module_projection_models,
    project_module_scores,
    tumor_normal_module_association,
)
from katana_cosplice.network import discover_module_labels
from katana_cosplice.preservation import (
    evaluate_external_module_preservation,
    render_preservation_plots,
)


def _read_matrix(path: Path) -> pd.DataFrame:
    matrix = pd.read_csv(path, sep="\t", index_col=0)
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    return matrix.apply(pd.to_numeric, errors="coerce")


def _fingerprint(paths: list[Path], settings: dict[str, object]) -> str:
    digest = hashlib.sha256(json.dumps(settings, sort_keys=True).encode())
    for path in paths:
        stat = path.stat()
        digest.update(str(path.resolve()).encode())
        digest.update(str(stat.st_size).encode())
        digest.update(str(stat.st_mtime_ns).encode())
    return digest.hexdigest()


def _assignments(
    event_keys: pd.Index,
    labels: np.ndarray,
    *,
    module_column: str = "consensus_module",
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "event_key": event_keys.astype(str),
            module_column: labels.astype(str),
        }
    )


def _fit_locked_validation_modules(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    locked_events: set[str],
    network: NetworkConfig,
    consensus: ConsensusConfig,
    output: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected = pd.Index(
        [event for event in locked_events if event in matrix.index],
        dtype=object,
    )
    locked = matrix.reindex(selected)
    quality = (
        locked.notna().mean(axis=1).ge(0.75)
        & locked.std(axis=1, skipna=True).gt(1e-8)
    )
    locked = locked.loc[quality]
    locked = locked.apply(lambda row: row.fillna(row.mean()), axis=1)
    sample_metadata = metadata[["sample_id", "phenotype"]].rename(
        columns={"phenotype": "condition"}
    )
    sample_metadata = sample_metadata.loc[
        sample_metadata["sample_id"].astype(str).isin(locked.columns)
    ]
    sample_metadata["sample_id"] = sample_metadata["sample_id"].astype(str)
    sample_metadata = sample_metadata.set_index("sample_id").reindex(
        locked.columns
    ).rename_axis("sample_id").reset_index()
    partition_settings = NetworkConfig(
        **{
            **asdict(network),
            "consensus_bootstrap_iterations": 0,
            "community_consensus_runs": 1,
            "bootstrap_iterations": 0,
        }
    )
    observed = _assignments(
        locked.index,
        discover_module_labels(locked, partition_settings),
        module_column="module",
    )
    result = build_bootstrap_coassignment_consensus(
        locked,
        observed,
        sample_metadata,
        network,
        consensus,
        progress=print,
    )
    result.assignments.to_csv(
        output / "locked_validation_module_assignments.tsv",
        sep="\t",
        index=False,
    )
    result.parameter_search.to_csv(
        output / "locked_validation_parameter_search.tsv",
        sep="\t",
        index=False,
    )
    audit = pd.DataFrame(
        [
            {
                "locked_discovery_event_union": len(locked_events),
                "exact_events_present_in_validation": len(selected),
                "events_passing_validation_missingness_and_variance_QC": len(
                    locked
                ),
                "validation_event_ranking_or_selection": "none",
                "validation_QC": ">=75% finite PSI and non-zero variance",
                "status": "finished",
            }
        ]
    )
    audit.to_csv(
        output / "locked_validation_event_audit.tsv",
        sep="\t",
        index=False,
    )
    return result.assignments, audit


def _render_phenotype(
    association: pd.DataFrame,
    output: Path,
) -> None:
    if association.empty:
        return
    plot = association.copy()
    plot["label"] = (
        plot["discovery_cohort"].astype(str)
        + " / "
        + plot["module"].astype(str)
    )
    plot = plot.sort_values("hedges_g")
    colors = np.where(
        plot["robust_all_and_paired_association"],
        "#0F766E",
        np.where(
            plot["all_sample_association_FDR_0.05"],
            "#D97706",
            "#94A3B8",
        ),
    )
    fig, ax = plt.subplots(
        figsize=(10, max(5.5, 0.45 * len(plot) + 1.5))
    )
    ax.barh(plot["label"], plot["hedges_g"], color=colors)
    ax.axvline(0, color="#111827", linewidth=0.8)
    lower = min(-1.0, float(plot["hedges_g"].min()) - 1.0)
    upper = max(1.0, float(plot["hedges_g"].max()) + 1.0)
    ax.set_xlim(lower, upper)
    for y, (_, row) in enumerate(plot.iterrows()):
        ax.text(
            row["hedges_g"],
            y,
            f"  g={row['hedges_g']:.2f}; q={row['welch_FDR']:.2g}",
            va="center",
            ha="left" if row["hedges_g"] >= 0 else "right",
            fontsize=8,
        )
    ax.set_xlabel(
        "Projected module-score effect (Hedges g, tumor vs adjacent normal)"
    )
    ax.set_ylabel("")
    ax.set_title(
        "Independent phenotype association of locked discovery modules",
        fontweight="bold",
        pad=14,
    )
    ax.text(
        0,
        1.01,
        "Teal: significant in all samples and matched pairs; amber: all-sample "
        "tests only. PC1 sign is module-specific.",
        transform=ax.transAxes,
        fontsize=8.5,
        color="#475569",
    )
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _render_replication(
    gene_replication: pd.DataFrame,
    theme_replication: pd.DataFrame,
    validation_name: str,
    output: Path,
) -> None:
    genes = gene_replication.loc[
        gene_replication["validation_cohort"].eq(validation_name)
    ].copy()
    themes = theme_replication.loc[
        theme_replication["validation_cohort"].eq(validation_name)
    ].copy()
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.8))
    if genes.empty:
        axes[0].text(0.5, 0.5, "No exact gene-level module match", ha="center")
        axes[0].axis("off")
    else:
        genes["label"] = (
            genes["discovery_cohort"].astype(str)
            + " / "
            + genes["discovery_module"].astype(str)
        )
        sns.barplot(
            data=genes.sort_values("best_gene_jaccard", ascending=False),
            x="best_gene_jaccard",
            y="label",
            color="#2563EB",
            ax=axes[0],
        )
        axes[0].set_xlim(0, 1)
        axes[0].set_title("Best independent module gene overlap")
        axes[0].set_xlabel("Unique-gene Jaccard")
        axes[0].set_ylabel("")
    if themes.empty:
        axes[1].text(
            0.5,
            0.5,
            "No independently significant GO/KEGG themes to match",
            ha="center",
            wrap=True,
        )
        axes[1].axis("off")
    else:
        themes["label"] = (
            themes["discovery_cohort"].astype(str)
            + " / "
            + themes["discovery_module"].astype(str)
        )
        sns.barplot(
            data=themes.sort_values("best_term_jaccard", ascending=False),
            x="best_term_jaccard",
            y="label",
            color="#D97706",
            ax=axes[1],
        )
        axes[1].set_xlim(0, 1)
        axes[1].set_title("Best independently enriched GO/KEGG overlap")
        axes[1].set_xlabel("Significant-term Jaccard")
        axes[1].set_ylabel("")
    fig.suptitle(
        f"Functional replication in untouched {validation_name}",
        fontsize=15,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _write_report(
    output: Path,
    validation_name: str,
    preservation: pd.DataFrame,
    phenotype: pd.DataFrame,
    gene_replication: pd.DataFrame,
    theme_replication: pd.DataFrame,
    generated: list[Path],
) -> Path:
    preserved = preservation["preservation_class"].value_counts()
    significant = phenotype.loc[
        phenotype["robust_all_and_paired_association"]
    ]
    all_sample_significant = phenotype.loc[
        phenotype["all_sample_association_FDR_0.05"]
    ]
    genes = gene_replication.loc[
        gene_replication["validation_cohort"].eq(validation_name)
    ]
    themes = theme_replication.loc[
        theme_replication["validation_cohort"].eq(validation_name)
    ]
    best_gene = (
        float(genes["best_gene_jaccard"].max()) if not genes.empty else np.nan
    )
    best_theme = (
        float(themes["best_term_jaccard"].max())
        if not themes.empty
        else np.nan
    )
    rows = "".join(
        f"<tr><td>{html.escape(str(row.discovery_cohort))}</td>"
        f"<td>{html.escape(str(row.module))}</td>"
        f"<td>{row.hedges_g:.3f}</td><td>{row.welch_FDR:.3g}</td>"
        f"<td>{int(row.n_matched_pairs)}</td></tr>"
        for row in significant.itertuples()
    )
    report = output / "LOCKED_EXTERNAL_VALIDATION_REPORT.html"
    report.write_text(
        f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Locked CoSplice external validation</title>
<style>
body{{font-family:Arial,sans-serif;max-width:1180px;margin:36px auto;color:#172033}}
h1,h2{{color:#0f4c5c}} .card{{border:1px solid #d9e2e8;border-radius:12px;padding:18px;margin:14px 0}}
.metric{{display:inline-block;margin:8px 24px 8px 0;font-size:1.1rem}}
img{{max-width:100%;border:1px solid #e5e7eb;border-radius:8px}}
table{{border-collapse:collapse;width:100%}}th,td{{padding:8px;border-bottom:1px solid #e5e7eb;text-align:left}}
.limit{{background:#fff7ed;border-left:4px solid #d97706;padding:12px}}
</style></head><body>
<h1>Locked external CoSplice validation: {html.escape(validation_name)}</h1>
<p>This validation uses exact event keys and parameters fixed before the validation
cohort was analyzed. The validation cohort was not used for event ranking or tuning.</p>
<div class="card">
<span class="metric"><b>Preserved:</b> {int(preserved.get("preserved", 0))}</span>
<span class="metric"><b>Partially preserved:</b> {int(preserved.get("partially_preserved", 0))}</span>
<span class="metric"><b>Not preserved:</b> {int(preserved.get("not_preserved", 0))}</span>
<span class="metric"><b>Robust phenotype-associated modules:</b> {len(significant)}</span>
<span class="metric"><b>All-sample associated modules:</b> {len(all_sample_significant)}</span>
<span class="metric"><b>Best independent gene Jaccard:</b> {best_gene:.3f}</span>
<span class="metric"><b>Best independent GO/KEGG Jaccard:</b> {best_theme:.3f}</span>
</div>
<h2>Module preservation</h2>
<img src="plots/external_module_preservation_zsummary.png">
<h2>Independent phenotype association</h2>
<img src="plots/locked_module_tumor_normal_association.png">
<table><thead><tr><th>Discovery cohort</th><th>Module</th><th>Hedges g</th>
<th>Welch FDR</th><th>Matched pairs</th></tr></thead><tbody>{rows or '<tr><td colspan="5">No module passed both all-sample and paired validation.</td></tr>'}</tbody></table>
<h2>Independent biological replication</h2>
<img src="plots/independent_functional_replication.png">
<p class="limit"><b>Interpretation boundary:</b> preservation, independently
recovered GO/KEGG themes and tumor-normal association support reproducibility and
phenotype relevance. They do not establish causal regulators or disease mechanism.
Cross-cancer non-preservation is informative and is not treated as pipeline failure.</p>
<h2>Reproducibility</h2><p>All tables, exact-event mappings, parameter settings,
input fingerprints and generated files are stored beside this report.</p>
</body></html>""",
        encoding="utf-8",
    )
    generated.append(report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Validate locked CoSplice discovery modules in one untouched cohort "
            "without validation-set event ranking or parameter tuning."
        )
    )
    parser.add_argument(
        "--discovery",
        action="append",
        nargs=3,
        metavar=("NAME", "PSI_MATRIX", "ASSIGNMENTS"),
        required=True,
    )
    parser.add_argument("--training-parameters", type=Path, required=True)
    parser.add_argument("--validation-name", required=True)
    parser.add_argument("--validation-matrix", type=Path, required=True)
    parser.add_argument("--validation-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--organism", default="hsapiens")
    parser.add_argument("--skip-functional-enrichment", action="store_true")
    args = parser.parse_args()

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    plots = output / "plots"
    plots.mkdir(exist_ok=True)
    parameter_source = json.loads(
        args.training_parameters.read_text(encoding="utf-8")
    )
    network = NetworkConfig(**parameter_source["network"])
    consensus = ConsensusConfig(**parameter_source["consensus"])
    preservation_settings = PreservationConfig(
        **parameter_source["preservation"]
    )
    validation_matrix = _read_matrix(args.validation_matrix)
    validation_metadata = pd.read_csv(
        args.validation_metadata,
        sep="\t",
    )
    discovery_inputs = [
        (str(name), Path(matrix), Path(assignments))
        for name, matrix, assignments in args.discovery
    ]
    input_paths = [
        args.training_parameters,
        args.validation_matrix,
        args.validation_metadata,
        *[path for _, matrix, assignments in discovery_inputs for path in (matrix, assignments)],
    ]
    fingerprint_settings = {
        "validation_name": args.validation_name,
        "network": asdict(network),
        "consensus": asdict(consensus),
        "preservation": asdict(preservation_settings),
        "organism": args.organism,
    }
    fingerprint = _fingerprint(input_paths, fingerprint_settings)

    external = ExternalCohortConfig(
        args.validation_name,
        args.validation_matrix,
        args.validation_metadata,
    )
    preservation_frames: list[pd.DataFrame] = []
    mapping_frames: list[pd.DataFrame] = []
    association_frames: list[pd.DataFrame] = []
    projection_audits: list[pd.DataFrame] = []
    assignment_sets: dict[str, dict[str, pd.DataFrame]] = {}
    locked_events: set[str] = set()
    for name, matrix_path, assignments_path in discovery_inputs:
        print(f"Locked external validation: {name} -> {args.validation_name}")
        matrix = _read_matrix(matrix_path)
        assignments = pd.read_csv(assignments_path, sep="\t")
        assignments["event_key"] = assignments["event_key"].astype(str)
        locked_events.update(assignments["event_key"])
        assignment_sets[name] = {"CoSplice": assignments}
        result = evaluate_external_module_preservation(
            matrix,
            assignments,
            [external],
            preservation_settings,
            progress=print,
        )
        current = result.module_preservation.copy()
        current.insert(0, "discovery_cohort", name)
        preservation_frames.append(current)
        mapping = result.event_mapping.copy()
        mapping.insert(0, "discovery_cohort", name)
        mapping_frames.append(mapping)
        models = fit_module_projection_models(matrix, assignments)
        scores, audit = project_module_scores(validation_matrix, models)
        audit.insert(0, "validation_cohort", args.validation_name)
        audit.insert(0, "discovery_cohort", name)
        projection_audits.append(audit)
        association = tumor_normal_module_association(
            scores,
            validation_metadata,
        )
        association.insert(0, "validation_cohort", args.validation_name)
        association.insert(0, "discovery_cohort", name)
        association_frames.append(association)

    validation_assignments, event_audit = _fit_locked_validation_modules(
        validation_matrix,
        validation_metadata,
        locked_events,
        network,
        consensus,
        output,
    )
    assignment_sets[args.validation_name] = {
        "CoSplice": validation_assignments
    }
    biology = evaluate_biological_interpretability(
        assignment_sets,
        organism=args.organism,
        run_functional_enrichment=not args.skip_functional_enrichment,
    )

    preservation = pd.concat(preservation_frames, ignore_index=True)
    mappings = pd.concat(mapping_frames, ignore_index=True)
    association = pd.concat(association_frames, ignore_index=True)
    projection_audit = pd.concat(projection_audits, ignore_index=True)
    preservation.to_csv(
        output / "locked_external_module_preservation.tsv",
        sep="\t",
        index=False,
    )
    mappings.to_csv(
        output / "locked_exact_event_mapping.tsv",
        sep="\t",
        index=False,
    )
    association.to_csv(
        output / "locked_module_phenotype_association.tsv",
        sep="\t",
        index=False,
    )
    projection_audit.to_csv(
        output / "locked_module_projection_audit.tsv",
        sep="\t",
        index=False,
    )
    biology.module_gene_audit.to_csv(
        output / "module_biological_interpretability.tsv",
        sep="\t",
        index=False,
    )
    biology.functional_enrichment.to_csv(
        output / "module_go_kegg_enrichment.tsv",
        sep="\t",
        index=False,
    )
    biology.enrichment_sources.to_csv(
        output / "functional_enrichment_sources.tsv",
        sep="\t",
        index=False,
    )
    biology.module_gene_replication.to_csv(
        output / "independent_module_gene_replication.tsv",
        sep="\t",
        index=False,
    )
    biology.functional_theme_replication.to_csv(
        output / "independent_functional_theme_replication.tsv",
        sep="\t",
        index=False,
    )
    biology.method_summary.to_csv(
        output / "biological_interpretability_summary.tsv",
        sep="\t",
        index=False,
    )
    generated: list[Path] = []
    preservation_plot_table = preservation.copy()
    preservation_plot_table["consensus_module"] = (
        preservation_plot_table["discovery_cohort"].astype(str)
        + " / "
        + preservation_plot_table["consensus_module"].astype(str)
    )
    preservation_result = type(
        "CombinedPreservation",
        (),
        {"module_preservation": preservation_plot_table},
    )()
    generated.extend(render_preservation_plots(preservation_result, plots))
    phenotype_plot = plots / "locked_module_tumor_normal_association.png"
    _render_phenotype(association, phenotype_plot)
    if phenotype_plot.exists():
        generated.append(phenotype_plot)
    replication_plot = plots / "independent_functional_replication.png"
    _render_replication(
        biology.module_gene_replication,
        biology.functional_theme_replication,
        args.validation_name,
        replication_plot,
    )
    if replication_plot.exists():
        generated.append(replication_plot)
    _write_report(
        output,
        args.validation_name,
        preservation,
        association,
        biology.module_gene_replication,
        biology.functional_theme_replication,
        generated,
    )
    parameters = {
        **fingerprint_settings,
        "fingerprint": fingerprint,
        "design": (
            "Untouched-cohort confirmation. Discovery modules and all network, "
            "consensus and preservation parameters were fixed before validation."
        ),
        "event_matching": "exact gene|event_type|TCGA_SpliceSeq_AS_ID",
        "validation_event_selection": (
            "No validation variance ranking. Independent validation modules use "
            "the union of discovery-selected events followed only by >=75% finite "
            "PSI and non-zero variance QC."
        ),
        "phenotype_test": (
            "Projected locked module PC1 scores; Welch and Mann-Whitney across all "
            "tumor/normal samples plus paired Wilcoxon sensitivity analysis for "
            "patient-matched tumor-normal pairs; BH-FDR per test family."
        ),
        "input_files": [str(path.resolve()) for path in input_paths],
        "locked_discovery_event_union": len(locked_events),
        "validation_event_audit": event_audit.iloc[0].to_dict(),
        "generated_files": [
            str(path)
            for path in sorted(
                [*output.glob("*.tsv"), output / "parameters.json", *generated]
            )
        ],
        "interpretation_boundary": (
            "Preservation, independent functional-theme recovery and phenotype "
            "association support reproducibility and relevance, not causality."
        ),
    }
    (output / "parameters.json").write_text(
        json.dumps(parameters, indent=2, default=str),
        encoding="utf-8",
    )
    print(json.dumps(
        {
            "output": str(output),
            "fingerprint": fingerprint,
            "preservation_rows": len(preservation),
            "phenotype_rows": len(association),
            "validation_modules": int(
                validation_assignments.loc[
                    validation_assignments["consensus_module"].ne("grey"),
                    "consensus_module",
                ].nunique()
            ),
        },
        indent=2,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
