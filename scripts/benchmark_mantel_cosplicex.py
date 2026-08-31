from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from katana_cosplice.functional_enrichment import gprofiler_enrichment
from katana_cosplice.mantel_cosplicex import (
    MantelCoSplicExConfig,
    benchmark_mantel_cosplicex,
    load_normalized_exon_counts,
)
from katana_cosplice.module_annotation import BUILTIN_RBP_SYMBOLS
from katana_cosplice.utils import bh_adjust, write_json


ATTRIBUTE_RE = re.compile(r'(\S+)\s+"([^"]+)"')


def _gtf_gene_symbols(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attributes = dict(ATTRIBUTE_RE.findall(fields[8]))
            gene_id = re.sub(r"\.\d+$", "", attributes.get("gene_id", ""))
            gene_name = attributes.get("gene_name", "")
            if gene_id and gene_name:
                mapping[gene_id] = gene_name
    return mapping


def _gene_module_biology(
    assignments: pd.DataFrame,
    *,
    gene_symbols: dict[str, str],
) -> pd.DataFrame:
    current = assignments.copy()
    current["gene_symbol"] = current["gene_id"].map(gene_symbols).fillna("")
    current = current.loc[current["mantel_cosplicex_module"].ne("grey")]
    universe = set(current["gene_symbol"].str.casefold()).difference({""})
    rbp = {symbol.casefold() for symbol in BUILTIN_RBP_SYMBOLS}
    rbp_universe = universe & rbp
    rows: list[dict[str, object]] = []
    for module, group in current.groupby("mantel_cosplicex_module", sort=True):
        symbols = sorted(set(group["gene_symbol"]).difference({""}))
        overlap = sorted(symbol for symbol in symbols if symbol.casefold() in rbp_universe)
        pvalue = (
            float(
                hypergeom.sf(
                    len(overlap) - 1,
                    len(universe),
                    len(rbp_universe),
                    len(symbols),
                )
            )
            if universe and rbp_universe and symbols and overlap
            else np.nan
        )
        rows.append(
            {
                "method": "Mantel-CoSplicEx",
                "module": str(module),
                "n_genes": int(len(group)),
                "n_mapped_gene_symbols": len(symbols),
                "rbp_splicing_factor_gene_count": len(overlap),
                "rbp_splicing_factor_genes": ";".join(overlap),
                "rbp_splicing_factor_enrichment_pvalue": pvalue,
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["rbp_splicing_factor_enrichment_FDR"] = bh_adjust(
            result["rbp_splicing_factor_enrichment_pvalue"]
        )
    return result


def _plot_head_to_head(summary: pd.DataFrame, output: Path) -> None:
    display = summary.copy()
    display["method_short"] = display["method"].replace(
        {
            "Mantel-CoSplicEx (published gene-level comparator)": "Mantel-CoSplicEx",
            "CoSplice": "CoSplice",
        }
    )
    metrics = [
        ("n_modules", "Non-grey modules", None),
        ("grey_fraction", "Grey fraction", (0, 1)),
        ("median_bootstrap_ARI", "Median bootstrap ARI", (0, 1)),
        ("functional_annotated_module_fraction", "Functionally annotated modules", (0, 1)),
    ]
    fig, axes = plt.subplots(1, len(metrics), figsize=(17, 4.6))
    palette = {"CoSplice": "#007C91", "Mantel-CoSplicEx": "#D97904"}
    for axis, (column, title, ylim) in zip(axes, metrics, strict=True):
        values = pd.to_numeric(display.get(column), errors="coerce")
        bars = axis.bar(
            display["method_short"],
            values.fillna(0),
            color=[palette.get(value, "#777777") for value in display["method_short"]],
            width=0.62,
        )
        axis.set_title(title, loc="left", fontweight="bold")
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(axis="x", rotation=20)
        if ylim:
            axis.set_ylim(*ylim)
        for bar, value in zip(bars, values, strict=True):
            label = "NA" if not np.isfinite(value) else f"{value:.3g}"
            axis.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
                label,
                ha="center",
                va="bottom",
                fontsize=9,
            )
    fig.suptitle(
        "Same-source public comparison: event-level CoSplice vs published gene-level Mantel-CoSplicEx",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.01,
        "GSE62669, identical 60 samples and exon-count source; node universes differ by design, so assignment ARI is not reported.",
        ha="center",
        fontsize=9,
        color="#4B5563",
    )
    fig.tight_layout(rect=[0, 0.06, 1, 0.92])
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)


def _write_report(
    output: Path,
    summary: pd.DataFrame,
    biology: pd.DataFrame,
    enrichment: pd.DataFrame,
) -> None:
    indexed = summary.set_index("method")
    mantel = indexed.loc["Mantel-CoSplicEx (published gene-level comparator)"]
    cosplice = indexed.loc["CoSplice"]
    mantel_terms = int(
        pd.to_numeric(enrichment.get("FDR"), errors="coerce").le(0.05).sum()
    ) if not enrichment.empty else 0
    rbp_modules = int(
        pd.to_numeric(
            biology.get("rbp_splicing_factor_enrichment_FDR"), errors="coerce"
        ).le(0.05).sum()
    ) if not biology.empty else 0
    text = f"""# Strict Mantel-CoSplicEx comparator

## Bottom line

This is a same-source, same-sample public comparison on GSE62669. CoSplice uses event-level
within-gene exon-usage features; the published Mantel-CoSplicEx method uses one multi-exon
distance profile per gene. The comparison is therefore a capability and validation comparison,
not an event-to-gene assignment contest.

- CoSplice: {int(cosplice['n_modules'])} modules, grey fraction
  {float(cosplice['grey_fraction']):.3f}, median bootstrap ARI
  {float(cosplice['median_bootstrap_ARI']):.3f}.
- Mantel-CoSplicEx: {int(mantel['n_modules'])} modules, grey fraction
  {float(mantel['grey_fraction']):.3f}, median bootstrap ARI
  {float(mantel['median_bootstrap_ARI']):.3f}.
- Mantel-CoSplicEx functional output: {mantel_terms} GO/KEGG terms at FDR <= 0.05 and
  {rbp_modules} RBP/splicing-factor-enriched modules at FDR <= 0.05.

## Biological interpretability

Mantel-CoSplicEx modules are directly gene-centred and are naturally suitable for GO/KEGG and
hub-gene interpretation. CoSplice retains event identity and therefore answers the more specific
question of which AS events form a reproducible program. A biologically useful CoSplice result
requires both a stable event partition and module-level functional or external-cohort support.

## Scientific boundary

No claim of universal superiority is made. Mantel-CoSplicEx is the closest published formal peer,
but it does not identify the specific event or isoform driving a gene-level module. CoSplice adds
that event resolution, uncertainty/read-support weighting and external preservation; these are
useful only when validated rather than treated as automatic biological truth.
"""
    (output / "STRICT_MANTEL_COSPLICEX_REPORT.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the published gene-level Mantel-CoSplicEx comparator on GSE62669."
    )
    parser.add_argument("--exon-counts", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rscript", type=Path, required=True)
    parser.add_argument("--genes", type=int, default=500)
    parser.add_argument("--bootstraps", type=int, default=20)
    parser.add_argument("--cosplice-summary", type=Path, required=True)
    parser.add_argument("--cosplice-biology", type=Path, required=True)
    parser.add_argument("--skip-go", action="store_true")
    args = parser.parse_args()

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    counts = load_normalized_exon_counts(args.exon_counts.resolve())
    config = MantelCoSplicExConfig(bootstrap_iterations=args.bootstraps)
    result = benchmark_mantel_cosplicex(
        counts,
        config,
        rscript=args.rscript.resolve(),
        max_genes=args.genes,
        progress=print,
    )
    result.assignments.to_csv(
        output / "mantel_cosplicex_gene_modules.tsv", sep="\t", index=False
    )
    result.gene_audit.to_csv(
        output / "mantel_cosplicex_gene_selection_audit.tsv", sep="\t", index=False
    )
    result.bootstrap_stability.to_csv(
        output / "mantel_cosplicex_bootstrap_stability.tsv", sep="\t", index=False
    )
    result.similarity.to_csv(
        output / "mantel_cosplicex_mantel_correlation.tsv", sep="\t", index=True
    )

    symbol_map = _gtf_gene_symbols(args.gtf.resolve())
    biology = _gene_module_biology(result.assignments, gene_symbols=symbol_map)
    biology.to_csv(
        output / "mantel_cosplicex_biological_audit.tsv", sep="\t", index=False
    )
    module_genes = {
        str(module): group["gene_id"].astype(str).tolist()
        for module, group in result.assignments.loc[
            result.assignments["mantel_cosplicex_module"].ne("grey")
        ].groupby("mantel_cosplicex_module", sort=True)
    }
    enrichment = pd.DataFrame()
    provider: dict[str, object] = {"status": "disabled"}
    if not args.skip_go and module_genes:
        try:
            enrichment, provider = gprofiler_enrichment(
                module_genes,
                organism="mmusculus",
                background=result.assignments["gene_id"].astype(str).tolist(),
            )
            provider["status"] = "finished"
        except Exception as exc:
            provider = {"status": f"unavailable: {exc}"}
    enrichment.to_csv(
        output / "mantel_cosplicex_go_kegg_enrichment.tsv", sep="\t", index=False
    )
    write_json(output / "mantel_cosplicex_enrichment_provider.json", provider)

    cosplice_summary = pd.read_csv(args.cosplice_summary.resolve(), sep="\t")
    cosplice_summary = cosplice_summary.loc[cosplice_summary["method"].eq("CoSplice")].head(1)
    if cosplice_summary.empty:
        raise ValueError("The supplied CoSplice summary does not contain a CoSplice row.")
    cosplice_biology = pd.read_csv(args.cosplice_biology.resolve(), sep="\t")
    cosplice_biology = cosplice_biology.loc[
        cosplice_biology["method"].eq("CoSplice")
        & cosplice_biology["cohort"].astype(str).str.contains("GSE62669")
    ].head(1)
    cosplice = cosplice_summary.copy()
    if not cosplice_biology.empty:
        for column in (
            "functional_annotated_module_fraction",
            "significant_GO_KEGG_term_count",
            "rbp_enriched_module_count_FDR_0.05",
        ):
            cosplice[column] = cosplice_biology.iloc[0].get(column, np.nan)
    mantel = result.method_summary.copy()
    mantel["functional_annotated_module_fraction"] = (
        float(
            enrichment.loc[
                pd.to_numeric(enrichment.get("FDR"), errors="coerce").le(0.05),
                "module",
            ].nunique()
            / max(int(mantel.iloc[0]["n_modules"]), 1)
        )
        if not enrichment.empty
        else 0.0
    )
    mantel["significant_GO_KEGG_term_count"] = (
        int(pd.to_numeric(enrichment.get("FDR"), errors="coerce").le(0.05).sum())
        if not enrichment.empty
        else 0
    )
    mantel["rbp_enriched_module_count_FDR_0.05"] = (
        int(
            pd.to_numeric(
                biology.get("rbp_splicing_factor_enrichment_FDR"), errors="coerce"
            ).le(0.05).sum()
        )
        if not biology.empty
        else 0
    )
    combined = pd.concat([cosplice, mantel], ignore_index=True, sort=False)
    combined.to_csv(output / "strict_comparator_summary.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "method": "CoSplice",
                "biological_unit": "AS/exon-usage event",
                "shared_source": "GSE62669 normalized exon counts",
                "shared_samples": 60,
                "direct_assignment_ARI_allowed": False,
                "reason": "event-level nodes are not interchangeable with gene-level nodes",
            },
            {
                "method": "Mantel-CoSplicEx",
                "biological_unit": "gene-level multi-exon profile",
                "shared_source": "GSE62669 normalized exon counts",
                "shared_samples": 60,
                "direct_assignment_ARI_allowed": False,
                "reason": "gene-level nodes are not interchangeable with event-level nodes",
            },
        ]
    ).to_csv(output / "comparator_scope.tsv", sep="\t", index=False)
    _plot_head_to_head(combined, output / "strict_comparator_scorecard.png")
    _write_report(output, combined, biology, enrichment)
    write_json(
        output / "parameters.json",
        {
            "exon_counts": str(args.exon_counts.resolve()),
            "gtf": str(args.gtf.resolve()),
            "rscript": str(args.rscript.resolve()),
            "max_genes": args.genes,
            "config": config.__dict__,
            "scientific_scope": (
                "same source and samples; native event-level versus gene-level feature universes"
            ),
        },
    )
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
