from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import hypergeom

from katana_cosplice.functional_enrichment import gprofiler_enrichment
from katana_cosplice.module_annotation import BUILTIN_RBP_SYMBOLS
from katana_cosplice.utils import bh_adjust


@dataclass
class BiologicalInterpretabilityResult:
    module_gene_audit: pd.DataFrame
    functional_enrichment: pd.DataFrame
    enrichment_sources: pd.DataFrame
    module_gene_replication: pd.DataFrame
    functional_theme_replication: pd.DataFrame
    method_summary: pd.DataFrame


def gene_from_event_key(event_key: object) -> str:
    """Return the gene component from Katana or TCGA SpliceSeq event keys."""
    text = str(event_key or "").strip()
    if not text:
        return ""
    return text.split("|", 1)[0].split(";", 1)[0].strip()


def _module_gene_audit(
    cohort: str,
    method: str,
    assignments: pd.DataFrame,
    *,
    module_column: str,
) -> pd.DataFrame:
    frame = assignments[["event_key", module_column]].copy()
    frame = frame.rename(columns={module_column: "module"})
    frame["module"] = frame["module"].astype(str)
    frame["gene"] = frame["event_key"].map(gene_from_event_key)
    frame = frame.loc[frame["module"].ne("grey") & frame["gene"].ne("")]
    rbp_symbols = {symbol.casefold() for symbol in BUILTIN_RBP_SYMBOLS}
    universe_genes = set(frame["gene"].str.casefold())
    rbp_universe = universe_genes & rbp_symbols
    rows: list[dict[str, object]] = []
    for module, group in frame.groupby("module", observed=True, sort=True):
        counts = group.groupby("gene", observed=True)["event_key"].nunique().sort_values(
            ascending=False
        )
        event_count = int(counts.sum())
        gene_count = int(len(counts))
        shares = counts.to_numpy(float) / max(event_count, 1)
        positive = shares[shares > 0]
        effective_genes = float(np.exp(-(positive * np.log(positive)).sum()))
        module_genes = set(counts.index.str.casefold())
        rbp_overlap = module_genes & rbp_universe
        pvalue = (
            float(
                hypergeom.sf(
                    len(rbp_overlap) - 1,
                    len(universe_genes),
                    len(rbp_universe),
                    gene_count,
                )
            )
            if universe_genes and rbp_universe and gene_count and rbp_overlap
            else np.nan
        )
        rows.append(
            {
                "cohort": cohort,
                "method": method,
                "module": module,
                "n_events": event_count,
                "n_unique_genes": gene_count,
                "events_per_unique_gene": event_count / max(gene_count, 1),
                "top_gene": str(counts.index[0]) if gene_count else "",
                "top_gene_event_count": int(counts.iloc[0]) if gene_count else 0,
                "top_gene_event_fraction": float(shares[0]) if len(shares) else np.nan,
                "effective_gene_count": effective_genes,
                "effective_gene_fraction": effective_genes / max(gene_count, 1),
                "multi_event_gene_fraction": float(counts.gt(1).mean()) if gene_count else np.nan,
                "rbp_splicing_factor_gene_count": len(rbp_overlap),
                "rbp_splicing_factor_genes": ";".join(
                    sorted(gene for gene in counts.index if gene.casefold() in rbp_overlap)
                ),
                "rbp_splicing_factor_enrichment_pvalue": pvalue,
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["rbp_splicing_factor_enrichment_FDR"] = bh_adjust(
            result["rbp_splicing_factor_enrichment_pvalue"]
        )
    return result


def _functional_enrichment(
    cohort: str,
    method: str,
    assignments: pd.DataFrame,
    *,
    module_column: str,
    organism: str,
    enabled: bool,
) -> tuple[pd.DataFrame, dict[str, object]]:
    frame = assignments[["event_key", module_column]].copy()
    frame = frame.rename(columns={module_column: "module"})
    frame["module"] = frame["module"].astype(str)
    frame["gene"] = frame["event_key"].map(gene_from_event_key)
    frame = frame.loc[frame["module"].ne("grey") & frame["gene"].ne("")]
    module_genes = {
        str(module): sorted(set(group["gene"]))
        for module, group in frame.groupby("module", observed=True, sort=True)
    }
    background = sorted(set(frame["gene"]))
    if not enabled:
        return pd.DataFrame(), {
            "cohort": cohort,
            "method": method,
            "status": "disabled",
            "organism": organism,
        }
    try:
        enrichment, provider = gprofiler_enrichment(
            module_genes,
            organism=organism,
            background=background,
        )
        if "method" in enrichment:
            enrichment = enrichment.rename(columns={"method": "enrichment_method"})
        enrichment.insert(0, "method", method)
        enrichment.insert(0, "cohort", cohort)
        return enrichment, {
            "cohort": cohort,
            "method": method,
            "status": "finished",
            **provider,
        }
    except Exception as exc:
        return pd.DataFrame(), {
            "cohort": cohort,
            "method": method,
            "status": f"unavailable: {exc}",
            "organism": organism,
        }


def _theme_replication(enrichment: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "method",
        "discovery_cohort",
        "discovery_module",
        "validation_cohort",
        "best_validation_module",
        "discovery_significant_terms",
        "validation_significant_terms",
        "shared_significant_terms",
        "best_term_jaccard",
    ]
    if enrichment.empty or "FDR" not in enrichment:
        return pd.DataFrame(columns=columns)
    significant = enrichment.loc[
        pd.to_numeric(enrichment.get("FDR"), errors="coerce").le(0.05)
    ].copy()
    rows: list[dict[str, object]] = []
    if significant.empty:
        return pd.DataFrame(columns=columns)
    methods = sorted(significant["method"].unique())
    for method in methods:
        subset = significant.loc[significant["method"].eq(method)]
        cohorts = sorted(subset["cohort"].unique())
        terms = {
            (str(cohort), str(module)): set(group["term_id"].astype(str))
            for (cohort, module), group in subset.groupby(
                ["cohort", "module"], observed=True
            )
        }
        for discovery, validation in permutations(cohorts, 2):
            discovery_modules = sorted(
                module for cohort, module in terms if cohort == discovery
            )
            validation_modules = sorted(
                module for cohort, module in terms if cohort == validation
            )
            for module in discovery_modules:
                source = terms[(discovery, module)]
                best_module = ""
                best_shared = 0
                best_target_size = 0
                best_jaccard = 0.0
                for target_module in validation_modules:
                    target = terms[(validation, target_module)]
                    shared = len(source & target)
                    union = len(source | target)
                    score = shared / union if union else 0.0
                    if score > best_jaccard:
                        best_module = target_module
                        best_shared = shared
                        best_target_size = len(target)
                        best_jaccard = score
                rows.append(
                    {
                        "method": method,
                        "discovery_cohort": discovery,
                        "discovery_module": module,
                        "validation_cohort": validation,
                        "best_validation_module": best_module,
                        "discovery_significant_terms": len(source),
                        "validation_significant_terms": best_target_size,
                        "shared_significant_terms": best_shared,
                        "best_term_jaccard": best_jaccard,
                    }
                )
    return pd.DataFrame(rows, columns=columns)


def _module_gene_replication(
    assignments: dict[str, dict[str, pd.DataFrame]],
) -> pd.DataFrame:
    columns = [
        "method",
        "discovery_cohort",
        "discovery_module",
        "validation_cohort",
        "best_validation_module",
        "discovery_unique_genes",
        "validation_unique_genes",
        "shared_unique_genes",
        "best_gene_jaccard",
    ]
    gene_sets: dict[tuple[str, str, str], set[str]] = {}
    for cohort, methods in assignments.items():
        for method, frame in methods.items():
            module_column = next(
                name
                for name in ("consensus_module", "wgcna_module", "module")
                if name in frame.columns
            )
            current = frame[["event_key", module_column]].copy()
            current["gene"] = current["event_key"].map(gene_from_event_key)
            current[module_column] = current[module_column].astype(str)
            current = current.loc[
                current[module_column].ne("grey") & current["gene"].ne("")
            ]
            for module, group in current.groupby(module_column, observed=True):
                gene_sets[(cohort, method, str(module))] = set(
                    group["gene"].str.casefold()
                )
    rows: list[dict[str, object]] = []
    methods = sorted({method for _, method, _ in gene_sets})
    for method in methods:
        cohorts = sorted({cohort for cohort, current, _ in gene_sets if current == method})
        for discovery, validation in permutations(cohorts, 2):
            source_modules = sorted(
                module
                for cohort, current, module in gene_sets
                if cohort == discovery and current == method
            )
            target_modules = sorted(
                module
                for cohort, current, module in gene_sets
                if cohort == validation and current == method
            )
            for module in source_modules:
                source = gene_sets[(discovery, method, module)]
                candidates: list[tuple[float, int, str, int]] = []
                for target_module in target_modules:
                    target = gene_sets[(validation, method, target_module)]
                    shared = len(source & target)
                    union = len(source | target)
                    candidates.append(
                        (
                            shared / union if union else 0.0,
                            shared,
                            target_module,
                            len(target),
                        )
                    )
                best = max(candidates, default=(0.0, 0, "", 0))
                rows.append(
                    {
                        "method": method,
                        "discovery_cohort": discovery,
                        "discovery_module": module,
                        "validation_cohort": validation,
                        "best_validation_module": best[2],
                        "discovery_unique_genes": len(source),
                        "validation_unique_genes": best[3],
                        "shared_unique_genes": best[1],
                        "best_gene_jaccard": best[0],
                    }
                )
    return pd.DataFrame(rows, columns=columns)


def evaluate_biological_interpretability(
    assignments: dict[str, dict[str, pd.DataFrame]],
    *,
    organism: str = "hsapiens",
    run_functional_enrichment: bool = True,
) -> BiologicalInterpretabilityResult:
    audits: list[pd.DataFrame] = []
    enrichments: list[pd.DataFrame] = []
    sources: list[dict[str, object]] = []
    for cohort, methods in assignments.items():
        for method, frame in methods.items():
            module_column = next(
                name
                for name in ("consensus_module", "wgcna_module", "module")
                if name in frame.columns
            )
            audits.append(
                _module_gene_audit(
                    cohort,
                    method,
                    frame,
                    module_column=module_column,
                )
            )
            enrichment, source = _functional_enrichment(
                cohort,
                method,
                frame,
                module_column=module_column,
                organism=organism,
                enabled=run_functional_enrichment,
            )
            if not enrichment.empty:
                enrichments.append(enrichment)
            sources.append(source)
    audit = pd.concat(audits, ignore_index=True) if audits else pd.DataFrame()
    enrichment = (
        pd.concat(enrichments, ignore_index=True) if enrichments else pd.DataFrame()
    )
    gene_replication = _module_gene_replication(assignments)
    replication = _theme_replication(enrichment)
    summary_rows: list[dict[str, object]] = []
    if not audit.empty:
        for method, group in audit.groupby("method", observed=True):
            functional = enrichment.loc[enrichment["method"].eq(method)] if not enrichment.empty else pd.DataFrame()
            significant = (
                functional.loc[pd.to_numeric(functional["FDR"], errors="coerce").le(0.05)]
                if not functional.empty
                else pd.DataFrame()
            )
            replicated = replication.loc[replication["method"].eq(method)]
            replicated_genes = gene_replication.loc[
                gene_replication["method"].eq(method)
            ]
            module_keys = set(zip(group["cohort"], group["module"], strict=False))
            annotated_keys = (
                set(zip(significant["cohort"], significant["module"], strict=False))
                if not significant.empty
                else set()
            )
            summary_rows.append(
                {
                    "method": method,
                    "n_modules_audited": len(group),
                    "median_unique_genes_per_module": float(group["n_unique_genes"].median()),
                    "median_effective_gene_fraction": float(group["effective_gene_fraction"].median()),
                    "worst_top_gene_event_fraction": float(group["top_gene_event_fraction"].max()),
                    "fraction_modules_not_gene_dominated": float(
                        group["top_gene_event_fraction"].le(0.20).mean()
                    ),
                    "rbp_enriched_module_count_FDR_0.05": int(
                        pd.to_numeric(
                            group["rbp_splicing_factor_enrichment_FDR"],
                            errors="coerce",
                        ).le(0.05).sum()
                    ),
                    "significant_GO_KEGG_term_count": int(len(significant)),
                    "functional_annotated_module_fraction": (
                        len(annotated_keys & module_keys) / len(module_keys)
                        if module_keys
                        else 0.0
                    ),
                    "median_best_cross_cohort_term_jaccard": (
                        float(replicated["best_term_jaccard"].median())
                        if not replicated.empty
                        else np.nan
                    ),
                    "median_best_cross_cohort_gene_jaccard": (
                        float(replicated_genes["best_gene_jaccard"].median())
                        if not replicated_genes.empty
                        else np.nan
                    ),
                }
            )
    return BiologicalInterpretabilityResult(
        module_gene_audit=audit,
        functional_enrichment=enrichment,
        enrichment_sources=pd.DataFrame(sources),
        module_gene_replication=gene_replication,
        functional_theme_replication=replication,
        method_summary=pd.DataFrame(summary_rows),
    )


def render_biological_interpretability(
    result: BiologicalInterpretabilityResult,
    output: Path,
) -> None:
    if result.method_summary.empty:
        return
    plot = result.method_summary.copy()
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    sns.barplot(
        data=plot,
        x="method",
        y="fraction_modules_not_gene_dominated",
        ax=axes[0],
    )
    axes[0].set_ylim(0, 1)
    axes[0].set_title("Modules not dominated by one gene")
    axes[0].set_ylabel("Fraction with top-gene event share <= 0.20")
    sns.barplot(
        data=plot,
        x="method",
        y="functional_annotated_module_fraction",
        ax=axes[1],
    )
    axes[1].set_ylim(0, 1)
    axes[1].set_title("Modules with GO/KEGG support")
    axes[1].set_ylabel("Fraction with FDR <= 0.05 term")
    replication_metric = (
        "median_best_cross_cohort_term_jaccard"
        if plot["median_best_cross_cohort_term_jaccard"].notna().any()
        else "median_best_cross_cohort_gene_jaccard"
    )
    sns.barplot(
        data=plot,
        x="method",
        y=replication_metric,
        ax=axes[2],
    )
    axes[2].set_ylim(0, 1)
    axes[2].set_title("Cross-cohort functional replication")
    axes[2].set_ylabel(
        "Median best significant-term Jaccard"
        if replication_metric.endswith("term_jaccard")
        else "Median best unique-gene Jaccard"
    )
    for ax in axes:
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=12)
    fig.suptitle(
        "Biological interpretability audit (validation only; not used to discover modules)",
        fontsize=15,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
