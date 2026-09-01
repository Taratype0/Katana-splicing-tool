from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import hypergeom

from katana_cosplice.config import ModuleStatisticsConfig
from katana_cosplice.module_models import fit_replicate_aware_models
from katana_cosplice.utils import bh_adjust


RUN_MAP = {
    "ES": "primary_ES500_KICH",
    "AA": "sensitivity_AA500_KICH",
    "AD": "sensitivity_AD500_KICH",
    "AP": "sensitivity_AP500_KICH",
    "AT": "sensitivity_AT500_KICH",
    "ME": "sensitivity_ME68_KICH",
    "RI": "sensitivity_RI500_KICH",
}
METHODS = {
    "CoSplice": ("three_objective_cosplice", "consensus_module"),
    "event-level WGCNA": ("event_level_wgcna", "wgcna_module"),
}


def _event_gene(value: object) -> str:
    return str(value).split("|", 1)[0].strip().upper()


def _split_genes(value: object) -> set[str]:
    return {
        item.strip().upper()
        for item in str(value).replace(",", ";").split(";")
        if item.strip() and item.strip().lower() != "nan"
    }


def _condition_means(value: object) -> dict[str, float]:
    result: dict[str, float] = {}
    for item in str(value).split(";"):
        if "=" not in item:
            continue
        key, raw = item.split("=", 1)
        try:
            result[key] = float(raw)
        except ValueError:
            continue
    return result


def _boolean(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False)
    return values.astype(str).str.strip().str.casefold().isin({"true", "1", "yes"})


def _fit_event_chunk(
    payload: tuple[
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame,
        ModuleStatisticsConfig,
        str | None,
    ],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    psi, ijc, sjc, samples, model_settings, batch_column = payload
    assignments = pd.DataFrame(
        {"event_key": psi.index.astype(str), "consensus_module": "external_test_set"}
    )
    result = fit_replicate_aware_models(
        psi,
        ijc,
        sjc,
        samples,
        assignments,
        pd.DataFrame(index=samples["sample_id"].astype(str)),
        model_settings,
        batch_column=batch_column,
    )
    return result.event_models, result.model_diagnostics


def _load_frozen_modules(locked_root: Path) -> pd.DataFrame:
    expected = [
        locked_root / run_id / "01_prelock" / directory / "discovery_assignments.tsv"
        for run_id in RUN_MAP.values()
        for directory, _column in METHODS.values()
    ]
    missing = [path for path in expected if not path.exists()]
    if missing:
        preview = "\n".join(f"- {path}" for path in missing[:10])
        suffix = f"\n- ... and {len(missing) - 10} more" if len(missing) > 10 else ""
        raise FileNotFoundError(
            "The frozen module root is incomplete. Missing assignment files:\n"
            f"{preview}{suffix}"
        )
    parts: list[pd.DataFrame] = []
    for event_type, run_id in RUN_MAP.items():
        for display, (directory, column) in METHODS.items():
            path = locked_root / run_id / "01_prelock" / directory / "discovery_assignments.tsv"
            frame = pd.read_csv(path, sep="\t", dtype=str)
            frame = frame.loc[frame[column].fillna("grey").ne("grey")].copy()
            frame["gene_symbol"] = frame["event_key"].map(_event_gene)
            frame["event_type"] = event_type
            frame["run_id"] = run_id
            frame["method"] = display
            frame["module"] = frame[column]
            parts.append(frame[["run_id", "event_type", "method", "module", "event_key", "gene_symbol"]])
    return pd.concat(parts, ignore_index=True)


def _event_models(cosplice_output: Path, cache_dir: Path) -> pd.DataFrame:
    cache = cache_dir / "GSE297044_measurement_qc_beta_binomial_event_models.tsv"
    cache_parameters = cache_dir / "GSE297044_measurement_qc_beta_binomial_parameters.json"
    matrix_dir = cosplice_output / "02_matrix"
    inputs = [
        matrix_dir / "sample_level_psi_matrix.tsv",
        matrix_dir / "sample_level_ijc_matrix.tsv",
        matrix_dir / "sample_level_sjc_matrix.tsv",
        matrix_dir / "sample_metadata_resolved.tsv",
    ]
    settings = {
        "condition_column": "condition",
        "min_total_reads": 10,
        "min_samples_per_condition": 3,
        "fdr_cutoff": 0.05,
        "max_events": None,
        "batch_column": "batch",
        "model_workers": 6,
        "eligibility": (
            "all competing-junction events with >=10 total reads in all three "
            "replicates of both conditions; no effect-size or dynamic-range prefilter"
        ),
    }
    digest = hashlib.sha256(json.dumps(settings, sort_keys=True).encode())
    for path in inputs:
        stat = path.stat()
        digest.update(str(path.resolve()).encode())
        digest.update(str(stat.st_size).encode())
        digest.update(str(stat.st_mtime_ns).encode())
    fingerprint = digest.hexdigest()
    cached_parameters = (
        json.loads(cache_parameters.read_text(encoding="utf-8"))
        if cache_parameters.exists()
        else {}
    )
    if cache.exists() and cached_parameters.get("fingerprint") == fingerprint:
        return pd.read_csv(cache, sep="\t", low_memory=False)
    psi = pd.read_csv(matrix_dir / "sample_level_psi_matrix.tsv", sep="\t", index_col=0)
    ijc = pd.read_csv(matrix_dir / "sample_level_ijc_matrix.tsv", sep="\t", index_col=0).reindex(psi.index)
    sjc = pd.read_csv(matrix_dir / "sample_level_sjc_matrix.tsv", sep="\t", index_col=0).reindex(psi.index)
    samples = pd.read_csv(matrix_dir / "sample_metadata_resolved.tsv", sep="\t")
    model_settings = ModuleStatisticsConfig(
        **{
            key: value
            for key, value in settings.items()
            if key not in {"batch_column", "eligibility", "model_workers"}
        }
    )
    worker_count = max(1, min(int(settings["model_workers"]), os.cpu_count() or 1))
    event_chunks = [
        pd.Index(chunk)
        for chunk in np.array_split(psi.index.to_numpy(), worker_count)
        if len(chunk)
    ]
    payloads = [
        (
            psi.loc[index],
            ijc.loc[index],
            sjc.loc[index],
            samples,
            model_settings,
            str(settings["batch_column"]),
        )
        for index in event_chunks
    ]
    if worker_count == 1:
        fitted = [_fit_event_chunk(payloads[0])]
    else:
        with ProcessPoolExecutor(max_workers=worker_count) as executor:
            fitted = list(executor.map(_fit_event_chunk, payloads))
    event_models = pd.concat([item[0] for item in fitted], ignore_index=True)
    diagnostics = pd.concat([item[1] for item in fitted], ignore_index=True)
    event_models["condition_effect_FDR"] = bh_adjust(
        pd.to_numeric(event_models["condition_effect_pvalue"], errors="coerce")
    )
    event_models["replicate_supported_condition_effect"] = event_models[
        "condition_effect_FDR"
    ].le(model_settings.fdr_cutoff)
    event_models = event_models.sort_values(
        ["condition_effect_FDR", "condition_effect_pvalue", "simple_max_delta_PSI"],
        ascending=[True, True, False],
        na_position="last",
    )
    event_models.to_csv(cache, sep="\t", index=False)
    diagnostics.to_csv(
        cache_dir / "GSE297044_measurement_qc_beta_binomial_model_diagnostics.tsv",
        sep="\t",
        index=False,
    )
    cache_parameters.write_text(
        json.dumps(
            {
                "fingerprint": fingerprint,
                "settings": settings,
                "inputs": [str(path.resolve()) for path in inputs],
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return event_models


def _load_gse(
    cosplice_output: Path,
    cache_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    models = _event_models(cosplice_output, cache_dir)
    events = pd.read_csv(cosplice_output / "02_matrix" / "event_metadata.tsv", sep="\t", low_memory=False)
    qc = pd.read_csv(cosplice_output / "03_stable_as_qc" / "stable_as_event_qc.tsv", sep="\t", low_memory=False)
    models = models.merge(
        events[["event_key", "gene_symbol", "gene_id", "event_coordinates", "cluster_id"]],
        on="event_key",
        how="left",
        validate="one_to_one",
    ).merge(qc, on="event_key", how="left", suffixes=("", "_qc"))
    deltas: list[float] = []
    for value in models["condition_mean_PSI"]:
        means = _condition_means(value)
        deltas.append(means.get("KO_DMSO", np.nan) - means.get("WT_DMSO", np.nan))
    models["KO_minus_WT_delta_PSI"] = deltas
    models["mapped_genes"] = models["gene_symbol"].map(_split_genes)
    return models, events, qc


def evaluate(
    frozen: pd.DataFrame,
    models: pd.DataFrame,
    *,
    fdr_cutoff: float,
    dpsi_cutoff: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    modeled = models.loc[models["mapped_genes"].map(bool)].copy()
    universe = set().union(*modeled["mapped_genes"].tolist()) if len(modeled) else set()
    responsive_mask = (
        pd.to_numeric(modeled["condition_effect_FDR"], errors="coerce").le(fdr_cutoff)
        & modeled["KO_minus_WT_delta_PSI"].abs().ge(dpsi_cutoff)
        & _boolean(modeled["converged"])
    )
    responsive = modeled.loc[responsive_mask].copy()
    responsive_genes = set().union(*responsive["mapped_genes"].tolist()) if len(responsive) else set()

    rows: list[dict[str, object]] = []
    membership_rows: list[dict[str, object]] = []
    for keys, group in frozen.groupby(["run_id", "event_type", "method", "module"], observed=True):
        run_id, event_type, method, module = keys
        genes = set(group["gene_symbol"]) & universe
        overlap = genes & responsive_genes
        pvalue = (
            float(hypergeom.sf(len(overlap) - 1, len(universe), len(responsive_genes), len(genes)))
            if universe and responsive_genes and genes and overlap
            else 1.0
        )
        module_events = modeled.loc[modeled["mapped_genes"].map(lambda values: bool(values & genes))]
        responsive_events = module_events.loc[responsive_mask.reindex(module_events.index, fill_value=False)]
        rows.append(
            {
                "run_id": run_id,
                "event_type": event_type,
                "method": method,
                "module": module,
                "frozen_module_events": group["event_key"].nunique(),
                "frozen_module_unique_genes": group["gene_symbol"].nunique(),
                "module_genes_testable_in_GSE297044": len(genes),
                "GSE297044_responsive_genes_in_module": len(overlap),
                "responsive_gene_fraction": len(overlap) / len(genes) if genes else np.nan,
                "responsive_gene_enrichment_pvalue": pvalue,
                "GSE297044_tested_junctions_in_module_genes": len(module_events),
                "GSE297044_responsive_junctions_in_module_genes": len(responsive_events),
                "median_abs_delta_PSI_in_module_genes": float(module_events["KO_minus_WT_delta_PSI"].abs().median()) if len(module_events) else np.nan,
                "median_total_junction_reads_in_module_genes": float(pd.to_numeric(module_events["median_total_junction_reads"], errors="coerce").median()) if len(module_events) else np.nan,
                "responsive_genes": ";".join(sorted(overlap)),
            }
        )
        for event in responsive_events.itertuples(index=False):
            membership_rows.append(
                {
                    "run_id": run_id,
                    "event_type": event_type,
                    "method": method,
                    "module": module,
                    "GSE297044_event_key": event.event_key,
                    "gene_symbol": event.gene_symbol,
                    "event_coordinates": event.event_coordinates,
                    "KO_minus_WT_delta_PSI": event.KO_minus_WT_delta_PSI,
                    "beta_binomial_FDR": event.condition_effect_FDR,
                    "median_total_junction_reads": event.median_total_junction_reads,
                }
            )
    summary = pd.DataFrame(rows)
    if not summary.empty:
        summary["responsive_gene_enrichment_FDR"] = bh_adjust(summary["responsive_gene_enrichment_pvalue"])
        summary = summary.sort_values(
            ["responsive_gene_enrichment_FDR", "responsive_gene_fraction"],
            ascending=[True, False],
        )
    audit = pd.DataFrame(
        [
            {"metric": "GSE297044_all_beta_binomial_models", "value": len(models)},
            {"metric": "GSE297044_converged_beta_binomial_models", "value": int(_boolean(models["converged"]).sum())},
            {"metric": "GSE297044_beta_binomial_tested_junctions", "value": len(modeled)},
            {"metric": "GSE297044_responsive_junctions", "value": len(responsive)},
            {"metric": "GSE297044_testable_genes", "value": len(universe)},
            {"metric": "GSE297044_responsive_genes", "value": len(responsive_genes)},
            {"metric": "responsive_definition", "value": f"beta-binomial FDR <= {fdr_cutoff}; abs(delta PSI) >= {dpsi_cutoff}; converged"},
            {"metric": "cross_platform_exact_event_matches", "value": 0},
            {"metric": "cross_platform_mapping_scope", "value": "gene/program level; TCGA SpliceSeq event IDs lack genomic coordinates"},
        ]
    )
    return summary, pd.DataFrame(membership_rows), audit


def _render(
    summary: pd.DataFrame,
    models: pd.DataFrame,
    output: Path,
    *,
    fdr_cutoff: float,
    dpsi_cutoff: float,
) -> None:
    sns.set_theme(style="ticks", context="notebook", font="Arial")
    fig, axes = plt.subplots(2, 2, figsize=(15.5, 10.5), constrained_layout=True)
    mapped = models.loc[models["mapped_genes"].map(bool)].copy()
    mapped_fdr = pd.to_numeric(mapped["condition_effect_FDR"], errors="coerce")
    mapped_delta = mapped["KO_minus_WT_delta_PSI"].abs()
    mapped_converged = _boolean(mapped["converged"])
    responsive = mapped.loc[
        mapped_fdr.le(fdr_cutoff)
        & mapped_delta.ge(dpsi_cutoff)
        & mapped_converged
    ].copy()

    funnel = pd.DataFrame(
        {
            "stage": [
                "Modelled + gene annotated",
                f"FDR <= {fdr_cutoff:g}",
                f"FDR + |delta PSI| >= {dpsi_cutoff:g}",
            ],
            "events": [len(mapped), int(mapped_fdr.le(fdr_cutoff).sum()), len(responsive)],
        }
    )
    sns.barplot(
        data=funnel,
        y="stage",
        x="events",
        hue="stage",
        palette=["#4C566A", "#2A6FBB", "#D1495B"],
        legend=False,
        ax=axes[0, 0],
    )
    axes[0, 0].set_xscale("log")
    for patch, value in zip(axes[0, 0].patches, funnel["events"], strict=False):
        axes[0, 0].text(
            value * 1.08,
            patch.get_y() + patch.get_height() / 2,
            f"{value:,}",
            va="center",
            fontsize=10,
            fontweight="bold",
        )
    axes[0, 0].set_title("a  Replicate/read-supported event testing", loc="left", fontweight="bold")
    axes[0, 0].set_xlabel("Number of junction events (log scale)")
    axes[0, 0].set_ylabel("")

    if len(responsive):
        responsive["direction"] = np.where(
            responsive["KO_minus_WT_delta_PSI"].ge(0),
            "Higher usage in SETD2 KO",
            "Lower usage in SETD2 KO",
        )
        sns.histplot(
            data=responsive,
            x="KO_minus_WT_delta_PSI",
            hue="direction",
            bins=24,
            palette={
                "Lower usage in SETD2 KO": "#2A6FBB",
                "Higher usage in SETD2 KO": "#D1495B",
            },
            ax=axes[0, 1],
        )
    axes[0, 1].axvline(0, color="#1F2937", linewidth=0.9)
    axes[0, 1].set_title("b  Independent SETD2-responsive junctions", loc="left", fontweight="bold")
    axes[0, 1].set_xlabel("SETD2 KO minus WT junction-usage delta PSI")
    axes[0, 1].set_ylabel("Junction events")

    plot = summary.copy()
    plot["label"] = plot["event_type"].astype(str) + " " + plot["module"].astype(str)
    plot["-log10 FDR"] = -np.log10(plot["responsive_gene_enrichment_FDR"].clip(lower=1e-300))
    top = plot.sort_values(["responsive_gene_enrichment_FDR", "responsive_gene_fraction"]).head(14)
    sns.scatterplot(
        data=top,
        x="-log10 FDR",
        y="label",
        size="GSE297044_responsive_genes_in_module",
        hue="method",
        style="method",
        palette={"CoSplice": "#0072B2", "event-level WGCNA": "#D55E00"},
        sizes=(45, 260),
        legend=False,
        ax=axes[1, 0],
    )
    for _index, row in top.iterrows():
        axes[1, 0].annotate(
            f"n={int(row['GSE297044_responsive_genes_in_module'])}",
            (row["-log10 FDR"], row["label"]),
            xytext=(7, 0),
            textcoords="offset points",
            va="center",
            fontsize=7.5,
            color="#374151",
        )
    axes[1, 0].axvline(-np.log10(0.05), color="#111827", linestyle="--", linewidth=1)
    axes[1, 0].set_title("c  Frozen KIRC program enrichment", loc="left", fontweight="bold")
    axes[1, 0].set_xlabel("-log10 enrichment FDR (dashed line = 0.05)")
    axes[1, 0].set_ylabel("")
    axes[1, 0].text(
        0.98,
        0.02,
        "blue circle: CoSplice\norange cross: event-level WGCNA\nlabel: responsive genes",
        transform=axes[1, 0].transAxes,
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="#4B5563",
    )

    paired = plot.copy()
    paired["module_number"] = paired["module"].astype(str).str.extract(r"(\d+)$")[0]
    paired = paired.pivot_table(
        index=["event_type", "module_number"],
        columns="method",
        values="-log10 FDR",
        aggfunc="first",
    ).dropna()
    if len(paired):
        axes[1, 1].scatter(
            paired["event-level WGCNA"],
            paired["CoSplice"],
            color="#217A70",
            edgecolor="white",
            linewidth=0.6,
            s=65,
        )
        limit = max(1.45, float(paired.to_numpy().max()) * 1.15)
        axes[1, 1].plot([0, limit], [0, limit], color="#9CA3AF", linewidth=1)
        axes[1, 1].axhline(-np.log10(0.05), color="#111827", linestyle="--", linewidth=0.8)
        axes[1, 1].axvline(-np.log10(0.05), color="#111827", linestyle="--", linewidth=0.8)
        axes[1, 1].set_xlim(0, limit)
        axes[1, 1].set_ylim(0, limit)
        axes[1, 1].text(
            0.04,
            0.96,
            f"Matched frozen programs: {len(paired)}",
            transform=axes[1, 1].transAxes,
            va="top",
            fontsize=9,
        )
    axes[1, 1].set_title("d  CoSplice versus event-level WGCNA", loc="left", fontweight="bold")
    axes[1, 1].set_xlabel("event-level WGCNA: -log10 enrichment FDR")
    axes[1, 1].set_ylabel("CoSplice: -log10 enrichment FDR")
    sns.despine(fig=fig)
    fig.suptitle(
        "Independent SETD2 ccRCC validation: event response is robust, frozen programs are not preserved",
        x=0.01,
        y=1.035,
        ha="left",
        fontsize=18,
        fontweight="bold",
    )
    fig.text(
        0.01,
        1.005,
        f"GSE297044, n=3 SETD2 KO and n=3 WT; beta-binomial FDR <= {fdr_cutoff:g}; |delta PSI| >= {dpsi_cutoff:g}",
        ha="left",
        va="bottom",
        fontsize=10,
        color="#4B5563",
    )
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate frozen KIRC programs in GSE297044 without retuning."
    )
    parser.add_argument("--locked-root", type=Path, required=True)
    parser.add_argument("--cosplice-output", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--abs-dpsi", type=float, default=0.10)
    args = parser.parse_args(argv)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    frozen = _load_frozen_modules(args.locked_root.resolve())
    models, _events, _qc = _load_gse(args.cosplice_output.resolve(), output)
    summary, membership, audit = evaluate(
        frozen,
        models,
        fdr_cutoff=args.fdr,
        dpsi_cutoff=args.abs_dpsi,
    )
    frozen.to_csv(output / "frozen_KIRC_module_membership.tsv", sep="\t", index=False)
    summary.to_csv(output / "GSE297044_frozen_program_replication.tsv", sep="\t", index=False)
    membership.to_csv(output / "GSE297044_responsive_junctions_in_frozen_programs.tsv", sep="\t", index=False)
    audit.to_csv(output / "validation_scope_audit.tsv", sep="\t", index=False)
    _render(
        summary,
        models,
        output / "GSE297044_SETD2_program_validation.png",
        fdr_cutoff=args.fdr,
        dpsi_cutoff=args.abs_dpsi,
    )

    significant = summary.loc[summary["responsive_gene_enrichment_FDR"].le(0.05)]
    best = summary.iloc[0]
    method_summary = (
        summary.groupby("method", observed=True)
        .agg(
            frozen_programs=("module", "size"),
            minimum_enrichment_FDR=("responsive_gene_enrichment_FDR", "min"),
            significant_programs=("responsive_gene_enrichment_FDR", lambda values: int(values.le(0.05).sum())),
        )
        .reset_index()
    )
    audit_values = dict(zip(audit["metric"], audit["value"], strict=False))
    report = [
        "# GSE297044 independent SETD2 ccRCC program validation",
        "",
        "## Primary result",
        "",
        f"- All replicate/read-eligible beta-binomial models: {int(audit_values['GSE297044_all_beta_binomial_models']):,}",
        f"- Converged beta-binomial models: {int(audit_values['GSE297044_converged_beta_binomial_models']):,}",
        f"- Gene-annotated tested junctions: {int(audit_values['GSE297044_beta_binomial_tested_junctions']):,}",
        f"- SETD2-responsive junctions after global FDR and effect-size filtering: {int(audit_values['GSE297044_responsive_junctions']):,}",
        f"- Responsive genes: {int(audit_values['GSE297044_responsive_genes']):,}",
        f"- Frozen KIRC modules with gene-program enrichment FDR <= 0.05: {len(significant)}",
        "",
        "GSE297044 contains a robust replicate/read-supported SETD2 junction-usage response, but none of the frozen TCGA KIRC gene programs is externally enriched at FDR <= 0.05. The best nominal result is "
        f"{best['event_type']} {best['module']} ({best['method']}): p={best['responsive_gene_enrichment_pvalue']:.4g}, FDR={best['responsive_gene_enrichment_FDR']:.4g}. This validation therefore does not demonstrate frozen-program preservation or a CoSplice advantage over event-level WGCNA.",
        "",
        "## Method comparison",
        "",
        method_summary.to_markdown(index=False, floatfmt=".4g"),
        "",
        "## Scientific boundary",
        "",
        "This analysis is frozen before examining GSE297044 results. It tests whether genes from locked KIRC splicing programs contain more replicate/read-supported SETD2-responsive junctions than expected. TCGA SpliceSeq event IDs do not provide harmonizable genomic coordinates here, so this is gene/program-level replication, not exact event preservation. Direction concordance is not claimed without an exact event map.",
        "",
        "## Significant replicated programs",
        "",
        significant.head(30).to_markdown(index=False, floatfmt=".4g") if len(significant) else "None passed FDR 0.05.",
    ]
    (output / "GSE297044_SETD2_PROGRAM_VALIDATION_REPORT.md").write_text(
        "\n".join(report), encoding="utf-8"
    )
    parameters = {
        "locked_root": str(args.locked_root.resolve()),
        "cosplice_output": str(args.cosplice_output.resolve()),
        "fdr_cutoff": args.fdr,
        "abs_delta_psi_cutoff": args.abs_dpsi,
        "mapping_scope": "gene/program-level only unless exact coordinates become available",
        "parameter_tuning_on_GSE297044": False,
        "event_model_cache": str((output / "GSE297044_measurement_qc_beta_binomial_event_models.tsv").resolve()),
        "event_model_parameter_file": str((output / "GSE297044_measurement_qc_beta_binomial_parameters.json").resolve()),
        "preregistration_amendment": str(
            (
                args.cosplice_output.resolve().parent
                / "inputs"
                / "preregistration_amendment_01.json"
            ).resolve()
        ),
        "generated_files": sorted(
            path.name
            for path in output.iterdir()
            if path.is_file() and path.name != "parameters.json"
        ),
        "generated_plot_files": sorted(
            path.name for path in output.glob("*.png") if path.is_file()
        ),
    }
    (output / "parameters.json").write_text(json.dumps(parameters, indent=2), encoding="utf-8")
    print(report[4])
    print(report[5])
    print(report[6])
    print(f"Output: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
