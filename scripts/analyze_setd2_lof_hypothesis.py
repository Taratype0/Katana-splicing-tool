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
from scipy import stats

from katana_cosplice.config import FactorContrastConfig
from katana_cosplice.factor_contrasts import analyze_factor_contrasts
from katana_cosplice.functional_enrichment import gprofiler_enrichment
from katana_cosplice.matrix import EVENT_TYPES, _event_frame
from katana_cosplice.module_annotation import BUILTIN_RBP_SYMBOLS
from katana_cosplice.utils import bh_adjust


COLORS = {
    "MUT": "#3167d5",
    "CLK3": "#e68613",
    "Pla": "#138a72",
    "DP": "#8f3fc0",
    "noDP": "#347f9f",
    "ink": "#172033",
    "muted": "#657083",
}


def _read_tsv(path: Path, **kwargs: object) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", low_memory=False, **kwargs)


def _hash_files(paths: list[Path]) -> str:
    digest = hashlib.sha256()
    for path in sorted(paths, key=lambda value: str(value).casefold()):
        digest.update(str(path.resolve()).encode("utf-8"))
        stat = path.stat()
        digest.update(f"{stat.st_size}:{stat.st_mtime_ns}".encode("ascii"))
    return digest.hexdigest()


def _save(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _load_rmats_comparison(rmats_root: Path, comparison: str, mode: str) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for event_type in EVENT_TYPES:
        path = rmats_root / comparison / "rmats_post" / f"{event_type}.MATS.{mode}.txt"
        raw = _read_tsv(path)
        identity = _event_frame(raw, event_type)
        values = identity.merge(
            raw.assign(_row=np.arange(len(raw)))[["_row", "FDR", "IncLevelDifference"]],
            left_index=True,
            right_on="_row",
            how="left",
        ).drop(columns="_row")
        values["FDR"] = pd.to_numeric(values["FDR"], errors="coerce")
        values["dPSI"] = pd.to_numeric(values["IncLevelDifference"], errors="coerce")
        values["comparison_id"] = comparison
        frames.append(values)
    result = pd.concat(frames, ignore_index=True)
    result = result.sort_values(["event_key", "FDR"], na_position="last").drop_duplicates("event_key")
    result["significant"] = result["FDR"].le(0.05) & result["dPSI"].abs().ge(0.10)
    return result


def _event_core_tables(
    comparisons: dict[str, dict[str, str]],
    rmats: dict[str, pd.DataFrame],
    assignments: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summaries: list[dict[str, object]] = []
    cores: list[pd.DataFrame] = []
    triple_rows: list[pd.DataFrame] = []
    for context, mapping in comparisons.items():
        mut = rmats[mapping["MUT"]]
        mut_sig = mut.loc[mut["significant"]].copy()
        context_shared: dict[str, pd.DataFrame] = {}
        for inhibitor in ("CLK3", "Pla"):
            drug = rmats[mapping[inhibitor]]
            drug_sig = drug.loc[drug["significant"]].copy()
            tested = mut[["event_key", "dPSI", "significant"]].merge(
                drug[["event_key", "dPSI", "significant"]],
                on="event_key",
                suffixes=("_MUT", f"_{inhibitor}"),
            )
            shared = mut_sig.merge(
                drug_sig[["event_key", "FDR", "dPSI"]],
                on="event_key",
                suffixes=("_MUT", f"_{inhibitor}"),
            )
            shared["context"] = context
            shared["inhibitor"] = inhibitor
            shared["same_direction"] = shared["dPSI_MUT"].mul(shared[f"dPSI_{inhibitor}"]).gt(0)
            shared["abs_delta_dPSI"] = (shared["dPSI_MUT"] - shared[f"dPSI_{inhibitor}"]).abs()
            shared["core_class"] = np.where(shared["same_direction"], f"{inhibitor}-like_same_direction", f"{inhibitor}-opposite")
            shared = shared.merge(assignments, on="event_key", how="left")
            cores.append(shared)
            context_shared[inhibitor] = shared

            universe = len(tested)
            n_mut_tested = int(tested["significant_MUT"].sum())
            n_drug_tested = int(tested[f"significant_{inhibitor}"].sum())
            n_mut_total = len(mut_sig)
            n_drug_total = len(drug_sig)
            overlap = len(shared)
            pvalue = (
                float(stats.hypergeom.sf(overlap - 1, universe, n_drug_tested, n_mut_tested))
                if universe and n_mut_tested and n_drug_tested
                else np.nan
            )
            same = int(shared["same_direction"].sum())
            broad = tested[["dPSI_MUT", f"dPSI_{inhibitor}"]].dropna()
            summaries.append(
                {
                    "context": context,
                    "inhibitor": inhibitor,
                    "mut_comparison": mapping["MUT"],
                    "inhibitor_comparison": mapping[inhibitor],
                    "n_exact_tested": universe,
                    "n_MUT_significant_total": n_mut_total,
                    "n_inhibitor_significant_total": n_drug_total,
                    "n_MUT_significant_in_exact_tested_universe": n_mut_tested,
                    "n_inhibitor_significant_in_exact_tested_universe": n_drug_tested,
                    "n_shared_significant": overlap,
                    "MUT_contained_pct": 100.0 * overlap / n_mut_total if n_mut_total else np.nan,
                    "overlap_enrichment_pvalue": pvalue,
                    "n_shared_same_direction": same,
                    "shared_same_direction_pct": 100.0 * same / overlap if overlap else np.nan,
                    "shared_dPSI_pearson": shared["dPSI_MUT"].corr(shared[f"dPSI_{inhibitor}"], method="pearson"),
                    "shared_dPSI_spearman": shared["dPSI_MUT"].corr(shared[f"dPSI_{inhibitor}"], method="spearman"),
                    "median_shared_abs_delta_dPSI": shared["abs_delta_dPSI"].median(),
                    "broad_dPSI_spearman": broad["dPSI_MUT"].corr(broad[f"dPSI_{inhibitor}"], method="spearman"),
                }
            )

        clk = context_shared["CLK3"]
        pla = context_shared["Pla"]
        triple = clk.merge(
            pla[["event_key", "FDR_Pla", "dPSI_Pla", "same_direction"]].rename(
                columns={"same_direction": "same_direction_MUT_Pla"}
            ),
            on="event_key",
            how="inner",
        )
        triple = triple.rename(columns={"same_direction": "same_direction_MUT_CLK3"})
        triple["all_three_same_direction"] = (
            triple["dPSI_MUT"].mul(triple["dPSI_CLK3"]).gt(0)
            & triple["dPSI_MUT"].mul(triple["dPSI_Pla"]).gt(0)
        )
        triple["core_class"] = "MUT-CLK3-Pla_triple_core"
        triple_rows.append(triple)

    summary = pd.DataFrame(summaries)
    summary["overlap_enrichment_FDR"] = bh_adjust(summary["overlap_enrichment_pvalue"])
    core = pd.concat(cores, ignore_index=True) if cores else pd.DataFrame()
    triple = pd.concat(triple_rows, ignore_index=True) if triple_rows else pd.DataFrame()
    return summary, core, triple


def _core_module_enrichment(
    core: pd.DataFrame,
    triple: pd.DataFrame,
    assignments: pd.DataFrame,
) -> pd.DataFrame:
    assignment = assignments.loc[
        assignments["module"].notna() & assignments["module"].ne("grey"),
        ["event_key", "module"],
    ].drop_duplicates("event_key")
    universe = set(assignment["event_key"].astype(str))
    module_sets = {
        str(module): set(subset["event_key"].astype(str))
        for module, subset in assignment.groupby("module", observed=True)
    }
    core_sets: list[tuple[str, str, str, set[str]]] = []
    for (context, inhibitor), subset in core.groupby(["context", "inhibitor"], observed=True):
        core_sets.append(("shared_all", str(context), str(inhibitor), set(subset["event_key"].astype(str))))
        same = subset.loc[subset["same_direction"].fillna(False)]
        core_sets.append(
            ("shared_same_direction", str(context), str(inhibitor), set(same["event_key"].astype(str)))
        )
    for context, subset in triple.groupby("context", observed=True):
        core_sets.append(("triple_all", str(context), "CLK3+Pla", set(subset["event_key"].astype(str))))
        same = subset.loc[subset["all_three_same_direction"].fillna(False)]
        core_sets.append(
            ("triple_all_three_same_direction", str(context), "CLK3+Pla", set(same["event_key"].astype(str)))
        )

    rows: list[dict[str, object]] = []
    for core_class, context, inhibitor, raw_events in core_sets:
        events = raw_events & universe
        for module, module_events in module_sets.items():
            overlap = len(events & module_events)
            expected = len(events) * len(module_events) / len(universe) if universe else np.nan
            rows.append(
                {
                    "core_set": core_class,
                    "context": context,
                    "inhibitor": inhibitor,
                    "module": module,
                    "event_universe": len(universe),
                    "module_events": len(module_events),
                    "core_events_total": len(raw_events),
                    "core_events_mapped_to_colored_modules": len(events),
                    "overlap_events": overlap,
                    "expected_overlap": expected,
                    "enrichment_ratio": overlap / expected if expected and expected > 0 else np.nan,
                    "pvalue": (
                        float(stats.hypergeom.sf(overlap - 1, len(universe), len(module_events), len(events)))
                        if overlap and universe
                        else 1.0
                    ),
                }
            )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["FDR"] = result.groupby(["core_set", "context", "inhibitor"], observed=True)["pvalue"].transform(
            lambda value: bh_adjust(value)
        )
        result = result.sort_values(["core_set", "context", "inhibitor", "FDR", "module"])
    return result


def _core_functional_tables(
    comparisons: dict[str, dict[str, str]],
    rmats: dict[str, pd.DataFrame],
    core: pd.DataFrame,
    triple: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
    membership_rows: list[dict[str, object]] = []
    event_type_rows: list[dict[str, object]] = []
    enrichment_frames: list[pd.DataFrame] = []
    providers: list[dict[str, object]] = []
    for context, mapping in comparisons.items():
        context_core = core.loc[core["context"].eq(context)].copy()
        clk = context_core.loc[
            context_core["inhibitor"].eq("CLK3") & context_core["same_direction"].fillna(False)
        ].drop_duplicates("event_key")
        pla = context_core.loc[
            context_core["inhibitor"].eq("Pla") & context_core["same_direction"].fillna(False)
        ].drop_duplicates("event_key")
        triple_same = triple.loc[
            triple["context"].eq(context) & triple["all_three_same_direction"].fillna(False)
        ].drop_duplicates("event_key")
        clk_keys = set(clk["event_key"].astype(str))
        pla_keys = set(pla["event_key"].astype(str))
        sets = {
            "CLK3_like": clk,
            "PlaB_like": pla,
            "MUT_CLK3_PlaB_triple": triple_same,
            "CLK3_like_exclusive": clk.loc[~clk["event_key"].astype(str).isin(pla_keys)],
            "PlaB_like_exclusive": pla.loc[~pla["event_key"].astype(str).isin(clk_keys)],
        }
        gene_sets: dict[str, list[str]] = {}
        for core_set, frame in sets.items():
            frame = frame.loc[frame["gene_symbol"].notna()].copy()
            gene_sets[core_set] = sorted(set(frame["gene_symbol"].astype(str)))
            for (gene_id, gene_symbol), genes in frame.groupby(["gene_id", "gene_symbol"], observed=True):
                membership_rows.append(
                    {
                        "context": context,
                        "core_set": core_set,
                        "gene_id": gene_id,
                        "gene_symbol": gene_symbol,
                        "n_events": genes["event_key"].nunique(),
                        "event_types": ";".join(sorted(set(genes["event_type"].astype(str)))),
                        "event_keys": ";".join(sorted(set(genes["event_key"].astype(str)))),
                    }
                )
            counts = frame.groupby("event_type", observed=True).size()
            for event_type, count in counts.items():
                event_type_rows.append(
                    {
                        "context": context,
                        "core_set": core_set,
                        "event_type": event_type,
                        "n_events": int(count),
                        "n_total_events": len(frame),
                        "fraction": float(count / len(frame)) if len(frame) else np.nan,
                    }
                )

        tested = rmats[mapping["MUT"]][["event_key", "gene_symbol"]].merge(
            rmats[mapping["CLK3"]][["event_key"]], on="event_key"
        ).merge(rmats[mapping["Pla"]][["event_key"]], on="event_key")
        background = sorted(set(tested["gene_symbol"].dropna().astype(str)))
        try:
            enrichment, provider = gprofiler_enrichment(
                gene_sets,
                organism="mmusculus",
                background=background,
                timeout_seconds=90,
            )
            enrichment = enrichment.rename(columns={"module": "core_set"})
            enrichment.insert(0, "context", context)
            enrichment["custom_background_genes"] = len(background)
            enrichment_frames.append(enrichment)
            providers.append({"context": context, **provider})
        except Exception as exc:  # pragma: no cover - remote provider failure is recorded, not hidden
            providers.append(
                {
                    "context": context,
                    "provider": "g:Profiler g:GOSt",
                    "status": "failed",
                    "error": str(exc),
                    "background": "custom",
                    "background_size": len(background),
                }
            )

    membership = pd.DataFrame(membership_rows)
    event_types = pd.DataFrame(event_type_rows)
    enrichment = pd.concat(enrichment_frames, ignore_index=True) if enrichment_frames else pd.DataFrame()
    return membership, event_types, enrichment, providers


def _module_annotation_summary(
    assignments: pd.DataFrame,
    pathways: pd.DataFrame,
    hub_genes: pd.DataFrame,
    rbp: pd.DataFrame,
) -> pd.DataFrame:
    preferred = {
        "M03": "nuclear-speck / RNA-regulatory program",
        "M06": "RNA-splicing program",
    }
    rows: list[dict[str, object]] = []
    for module, subset in assignments.loc[assignments["module"].ne("grey")].groupby("module", observed=True):
        module_pathways = pathways.loc[
            pathways["module"].eq(module) & pathways["FDR"].le(0.05)
        ].sort_values("FDR")
        informative = module_pathways.loc[
            module_pathways["term_name"].astype(str).str.contains(
                "splic|speck|RNA|ribonucleo|nucleoplasm", case=False, regex=True
            )
        ]
        selected_terms = pd.concat([informative.head(5), module_pathways.head(5)]).drop_duplicates("term_id").head(5)
        genes = hub_genes.loc[hub_genes["module"].eq(module)].sort_values("hub_gene_rank").head(10)
        regulators = rbp.loc[
            rbp["module"].eq(module) & rbp["FDR"].le(0.05)
        ].sort_values("FDR").head(8)
        rows.append(
            {
                "module": module,
                "interpretation_label": preferred.get(str(module), "co-splicing program"),
                "n_events": subset["event_key"].nunique(),
                "n_genes": hub_genes.loc[hub_genes["module"].eq(module), "gene_id"].nunique(),
                "top_hub_genes": ";".join(genes["gene_symbol"].dropna().astype(str)),
                "top_pathways": ";".join(selected_terms["term_name"].dropna().astype(str)),
                "top_expression_correlated_RBPs": ";".join(regulators["RBP"].dropna().astype(str)),
                "interpretation_limit": "Program/regulator prioritization only; association is not causal evidence.",
            }
        )
    return pd.DataFrame(rows)


def _write_interpretation_report(
    path: Path,
    summary: pd.DataFrame,
    occlusion: pd.DataFrame,
    enrichment: pd.DataFrame,
    annotation: pd.DataFrame,
    core_functional_enrichment: pd.DataFrame,
) -> None:
    lines = [
        "# Setd2 MUT LOF-like splicing hypothesis closure",
        "",
        "## Executive conclusion",
        "",
        "The RNA-seq data support a **partial inhibitor-like splicing component** of Setd2 MUT, not a complete CLK3/PlaB phenocopy. Exact shared significant events show strong directional and effect-size agreement, but module-level inhibitor responses do not show formal occlusion in M03 or M06. Therefore RNA-seq supports molecular convergence consistent with a possible LOF component, while definitive LOF/GOF classification still requires orthogonal activity or rescue evidence.",
        "",
        "## 1. Exact significant-event evidence",
        "",
        "Events are matched by exact canonical rMATS identity (gene_id, event_type, chromosome, strand and event coordinates). Significance is FDR <= 0.05 and abs(dPSI) >= 0.10.",
        "",
        "| Context | Reference | MUT significant | Shared significant | MUT covered | Same direction | dPSI Spearman | Median abs delta dPSI |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary.itertuples(index=False):
        lines.append(
            f"| {row.context} | {row.inhibitor} | {int(row.n_MUT_significant_total)} | {int(row.n_shared_significant)} | "
            f"{row.MUT_contained_pct:.1f}% | {row.shared_same_direction_pct:.1f}% | "
            f"{row.shared_dPSI_spearman:.2f} | {row.median_shared_abs_delta_dPSI:.3f} |"
        )
    lines.extend(
        [
            "",
            "Interpretation: the shared subset is strongly directionally concordant, but it represents only part of the MUT signature. CLK3-like and PlaB-like cores must therefore be reported separately rather than collapsed into one inhibitor category.",
            "",
            "## 2. MUT + inhibitor combination tests",
            "",
            "Occlusion is tested as `[MUT inhibitor - MUT baseline] - [WT inhibitor - WT baseline]`. Convergence asks whether the absolute MUT-WT module-score gap becomes smaller after inhibitor treatment.",
            "",
            "| Module | Context | Inhibitor | WT effect | MUT effect | Interaction FDR | Occlusion supported | Convergence (descriptive) |",
            "|---|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    focus = occlusion.loc[occlusion["module"].isin(["M03", "M06"])].copy()
    for row in focus.itertuples(index=False):
        lines.append(
            f"| {row.module} | {row.DP_status} | {row.inhibitor} | {row.WT_inhibitor_effect:+.3f} | "
            f"{row.MUT_inhibitor_effect:+.3f} | {row.occlusion_interaction_FDR:.3g} | "
            f"{'yes' if row.formal_occlusion_supported else 'no'} | {'yes' if row.descriptive_convergence else 'no'} |"
        )
    lines.extend(
        [
            "",
            "M03/M06 show robust inhibitor movement, but neither module passes the formal occlusion criterion. M03 under DP + PlaB shows significant enhancement rather than attenuation. This argues against presenting these modules as proof of a simple shared-pathway LOF mechanism.",
            "",
            "## 3. Co-splicing programs",
            "",
            "Modules are exact-WGCNA co-splicing programs calculated from replicate PSI. They are supporting program-level evidence and are not equivalent to the exact significant-event overlap above.",
            "",
            "| Module | Interpretation | Events | Genes | Top hub genes | Top pathways | Prioritized RBPs |",
            "|---|---|---:|---:|---|---|---|",
        ]
    )
    for row in annotation.loc[annotation["module"].isin(["M03", "M06"])].itertuples(index=False):
        lines.append(
            f"| {row.module} | {row.interpretation_label} | {int(row.n_events)} | {int(row.n_genes)} | "
            f"{row.top_hub_genes} | {row.top_pathways} | {row.top_expression_correlated_RBPs} |"
        )
    lines.extend(
        [
            "",
            "## 4. Do the exact inhibitor-like cores concentrate in M03/M06?",
            "",
            "The enrichment table formally tests this question. A small raw overlap must not be called module support unless its within-core enrichment FDR is significant.",
            "",
            "| Core | Context | Reference | Module | Core mapped | Overlap | Enrichment ratio | FDR |",
            "|---|---|---|---|---:|---:|---:|---:|",
        ]
    )
    shown = enrichment.loc[
        enrichment["module"].isin(["M03", "M06"])
        & enrichment["core_set"].isin(["shared_same_direction", "triple_all_three_same_direction"])
    ]
    for row in shown.itertuples(index=False):
        lines.append(
            f"| {row.core_set} | {row.context} | {row.inhibitor} | {row.module} | "
            f"{int(row.core_events_mapped_to_colored_modules)} | {int(row.overlap_events)} | "
            f"{row.enrichment_ratio:.2f} | {row.FDR:.3g} |"
        )
    lines.extend(
        [
            "",
            "## Reporting boundary",
            "",
            "- Strong claim supported: a finite MUT AS subset converges with CLK3/PlaB inhibition in event identity, direction and dPSI magnitude.",
            "- Claim not supported: MUT is a complete inhibitor phenocopy or M03/M06 demonstrate formal inhibitor occlusion.",
            "- LOF classification: RNA-seq is supportive molecular evidence only; add SETD2 activity/mark measurement, rescue, or another orthogonal functional assay for a definitive conclusion.",
        ]
    )
    if not core_functional_enrichment.empty:
        lines.extend(
            [
                "",
                "## Functional interpretation of inhibitor-like event cores",
                "",
                "GO/KEGG over-representation uses unique genes and a custom background restricted to genes represented by exact events testable in MUT, CLK3 and PlaB contrasts in the same DP context. Results are functional prioritization, not proof that the pathway caused the AS event.",
                "",
                "| Context | Core | Source | Term | FDR | Overlap genes |",
                "|---|---|---|---|---:|---:|",
            ]
        )
        functional = core_functional_enrichment.sort_values(["context", "core_set", "FDR"])
        functional = functional.groupby(["context", "core_set"], observed=True).head(5)
        for row in functional.itertuples(index=False):
            lines.append(
                f"| {row.context} | {row.core_set} | {row.source} | {row.term_name} | "
                f"{row.FDR:.3g} | {int(row.intersection_size)} |"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _program_scores_long(scores: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    long = scores.melt(id_vars="sample_id", var_name="module", value_name="module_score")
    return long.merge(samples, on="sample_id", how="left", validate="many_to_one")


def _condition_and_occlusion_tables(
    scores: pd.DataFrame,
    samples: pd.DataFrame,
    module_contrasts: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    long = _program_scores_long(scores, samples)
    condition = (
        long.groupby(["module", "DP_status", "genotype", "treatment"], observed=True)["module_score"]
        .agg(mean="mean", sd="std", n="count")
        .reset_index()
    )
    interactions = module_contrasts.loc[
        module_contrasts["contrast_type"].eq("perturbation_occlusion_interaction")
    ].copy()
    interactions = interactions.rename(
        columns={
            "context": "context_label",
            "effect": "occlusion_interaction_effect",
            "FDR": "occlusion_interaction_FDR",
            "pvalue": "occlusion_interaction_pvalue",
        }
    )
    rows: list[dict[str, object]] = []
    for module in scores.columns.drop("sample_id"):
        for dp in ("noDP", "DP"):
            for inhibitor in ("CLK3", "Pla"):
                subset = condition.loc[condition["module"].eq(module) & condition["DP_status"].eq(dp)]
                lookup = {
                    (str(row.genotype), str(row.treatment)): float(row.mean)
                    for row in subset.itertuples(index=False)
                }
                required = [("WT", "Baseline"), ("MUT", "Baseline"), ("WT", inhibitor), ("MUT", inhibitor)]
                if any(key not in lookup for key in required):
                    continue
                wt_base, mut_base = lookup[("WT", "Baseline")], lookup[("MUT", "Baseline")]
                wt_drug, mut_drug = lookup[("WT", inhibitor)], lookup[("MUT", inhibitor)]
                wt_effect = wt_drug - wt_base
                mut_effect = mut_drug - mut_base
                gap_base = mut_base - wt_base
                gap_drug = mut_drug - wt_drug
                formal = interactions.loc[
                    interactions["module"].eq(module)
                    & interactions["context_label"].eq(f"DP_status={dp}")
                    & interactions["perturbation"].eq(inhibitor)
                ]
                formal_row = formal.iloc[0] if not formal.empty else pd.Series(dtype=object)
                same = wt_effect * mut_effect > 0
                attenuation = 1.0 - abs(mut_effect) / abs(wt_effect) if same and wt_effect else np.nan
                rows.append(
                    {
                        "module": module,
                        "DP_status": dp,
                        "inhibitor": inhibitor,
                        "WT_baseline_score": wt_base,
                        "MUT_baseline_score": mut_base,
                        "WT_inhibitor_score": wt_drug,
                        "MUT_inhibitor_score": mut_drug,
                        "baseline_MUT_minus_WT": gap_base,
                        "inhibitor_MUT_minus_WT": gap_drug,
                        "WT_inhibitor_effect": wt_effect,
                        "MUT_inhibitor_effect": mut_effect,
                        "same_inhibitor_response_direction": same,
                        "descriptive_occlusion_fraction": attenuation,
                        "absolute_genotype_gap_reduction": abs(gap_base) - abs(gap_drug),
                        "descriptive_convergence": abs(gap_drug) < abs(gap_base),
                        "occlusion_interaction_effect": formal_row.get("occlusion_interaction_effect", np.nan),
                        "occlusion_interaction_pvalue": formal_row.get("occlusion_interaction_pvalue", np.nan),
                        "occlusion_interaction_FDR": formal_row.get("occlusion_interaction_FDR", np.nan),
                        "formal_occlusion_supported": bool(
                            pd.notna(formal_row.get("occlusion_interaction_FDR", np.nan))
                            and float(formal_row["occlusion_interaction_FDR"]) <= 0.05
                            and same
                            and pd.notna(attenuation)
                            and attenuation > 0
                        ),
                    }
                )
    return condition, pd.DataFrame(rows)


def _hub_tables(
    psi: pd.DataFrame,
    scores: pd.DataFrame,
    assignments: pd.DataFrame,
    metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    score_matrix = scores.set_index("sample_id")
    rows: list[pd.DataFrame] = []
    for module, subset in assignments.groupby("module", observed=True):
        keys = psi.index.intersection(subset["event_key"].astype(str))
        if module not in score_matrix or len(keys) == 0:
            continue
        module_score = pd.to_numeric(score_matrix[module], errors="coerce").reindex(psi.columns)
        correlations = psi.loc[keys].T.corrwith(module_score).rename("kME").reset_index()
        correlations = correlations.rename(columns={correlations.columns[0]: "event_key"})
        correlations["module"] = module
        rows.append(correlations)
    events = pd.concat(rows, ignore_index=True)
    events["abs_kME"] = events["kME"].abs()
    events = events.merge(metadata, on="event_key", how="left")
    events["hub_event_rank"] = events.groupby("module")["abs_kME"].rank(method="min", ascending=False)
    genes = (
        events.groupby(["module", "gene_id", "gene_symbol"], as_index=False)
        .agg(
            n_module_events=("event_key", "nunique"),
            event_types=("event_type", lambda value: ";".join(sorted(set(value.astype(str))))),
            max_abs_kME=("abs_kME", "max"),
            mean_abs_kME=("abs_kME", "mean"),
        )
    )
    genes["hub_gene_score"] = genes["max_abs_kME"] * np.log1p(genes["n_module_events"])
    genes["hub_gene_rank"] = genes.groupby("module")["hub_gene_score"].rank(method="min", ascending=False)
    return events.sort_values(["module", "hub_event_rank"]), genes.sort_values(["module", "hub_gene_rank"])


def _rbp_tables(
    scores: pd.DataFrame,
    hub_genes: pd.DataFrame,
    expression_path: Path,
    motif_hits_path: Path,
    assignments: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    score_matrix = scores.set_index("sample_id")
    expression = _read_tsv(expression_path, index_col=0)
    expression.index = expression.index.astype(str).str.casefold()
    expression = expression.apply(pd.to_numeric, errors="coerce")
    rbp_keys = {symbol.casefold(): symbol for symbol in BUILTIN_RBP_SYMBOLS}
    corr_rows: list[dict[str, object]] = []
    for key in sorted(set(expression.index) & set(rbp_keys)):
        values = expression.loc[key].reindex(score_matrix.index)
        for module in score_matrix.columns:
            aligned = pd.concat([values, score_matrix[module]], axis=1).dropna()
            if len(aligned) < 4 or aligned.iloc[:, 0].nunique() < 2:
                continue
            correlation, pvalue = stats.pearsonr(aligned.iloc[:, 0], aligned.iloc[:, 1])
            corr_rows.append(
                {
                    "module": module,
                    "RBP": rbp_keys[key],
                    "expression_module_correlation": correlation,
                    "pvalue": pvalue,
                    "n_samples": len(aligned),
                }
            )
    correlations = pd.DataFrame(corr_rows)
    correlations["FDR"] = bh_adjust(correlations["pvalue"])
    member = hub_genes.loc[
        hub_genes["gene_symbol"].astype(str).str.casefold().isin(rbp_keys),
        ["module", "gene_symbol", "n_module_events", "max_abs_kME", "hub_gene_rank"],
    ].rename(columns={"gene_symbol": "RBP"})
    correlations = correlations.merge(member, on=["module", "RBP"], how="outer")
    correlations["RBP_is_module_gene"] = correlations["n_module_events"].notna()

    hits = _read_tsv(motif_hits_path)[["event_key", "RBP"]].drop_duplicates()
    universe = set(assignments["event_key"].astype(str))
    motif_rows: list[dict[str, object]] = []
    for rbp, rbp_hits in hits.groupby("RBP", observed=True):
        hit_set = set(rbp_hits["event_key"].astype(str)) & universe
        for module, members in assignments.groupby("module", observed=True):
            member_set = set(members["event_key"].astype(str)) & universe
            overlap = member_set & hit_set
            pvalue = (
                float(stats.hypergeom.sf(len(overlap) - 1, len(universe), len(hit_set), len(member_set)))
                if overlap
                else np.nan
            )
            motif_rows.append(
                {
                    "module": module,
                    "RBP_or_motif": rbp,
                    "event_universe": len(universe),
                    "motif_events": len(hit_set),
                    "module_events": len(member_set),
                    "overlap_events": len(overlap),
                    "pvalue": pvalue,
                }
            )
    motif = pd.DataFrame(motif_rows)
    motif["FDR"] = bh_adjust(motif["pvalue"])
    return correlations.sort_values(["module", "FDR"], na_position="last"), motif.sort_values(["module", "FDR"], na_position="last")


def _plot_event_evidence(summary: pd.DataFrame, path: Path) -> None:
    labels = (summary["context"] + "\n" + summary["inhibitor"]).tolist()
    x = np.arange(len(labels))
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    axes[0].bar(x, summary["MUT_contained_pct"], color=[COLORS[value] for value in summary["inhibitor"]])
    axes[0].set_title("MUT significant-event containment")
    axes[0].set_ylabel("MUT events also significant in inhibitor (%)")
    axes[1].bar(x, summary["shared_same_direction_pct"], color=[COLORS[value] for value in summary["inhibitor"]])
    axes[1].axhline(50, color=COLORS["muted"], linestyle="--", linewidth=1)
    axes[1].set_title("Direction agreement in exact shared events")
    axes[1].set_ylabel("Same standardized dPSI direction (%)")
    axes[1].set_ylim(0, 105)
    axes[2].bar(x, summary["shared_dPSI_pearson"], color=[COLORS[value] for value in summary["inhibitor"]])
    axes[2].axhline(0, color=COLORS["muted"], linewidth=1)
    axes[2].set_title("Effect-size agreement in exact shared events")
    axes[2].set_ylabel("Pearson r of standardized dPSI")
    axes[2].set_ylim(-1, 1)
    for ax in axes:
        ax.set_xticks(x, labels, rotation=35, ha="right")
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle("Setd2 MUT versus pharmacological LOF-reference AS signatures", fontsize=17, fontweight="bold")
    fig.text(0.5, 0.01, "Exact rMATS event identity; FDR <= 0.05 and abs(dPSI) >= 0.10 in both contrasts.", ha="center", color=COLORS["muted"])
    fig.tight_layout(rect=[0.02, 0.05, 0.98, 0.93])
    _save(fig, path)


def _plot_module_scores(condition: pd.DataFrame, path: Path, focus: list[str]) -> None:
    selected = condition.loc[condition["module"].isin(focus)].copy()
    modules = [module for module in focus if module in set(selected["module"])]
    fig, axes = plt.subplots(len(modules), 2, figsize=(13, 4.3 * len(modules)), squeeze=False, sharex=False)
    order = ["Baseline", "CLK3", "Pla"]
    for row_index, module in enumerate(modules):
        for col_index, dp in enumerate(["noDP", "DP"]):
            ax = axes[row_index, col_index]
            part = selected.loc[selected["module"].eq(module) & selected["DP_status"].eq(dp)]
            for genotype, color, marker in [("WT", "#7b8797", "o"), ("MUT", COLORS["MUT"], "s")]:
                line = part.loc[part["genotype"].eq(genotype)].set_index("treatment").reindex(order)
                ax.errorbar(order, line["mean"], yerr=line["sd"], label=genotype, color=color, marker=marker, linewidth=2, capsize=3)
            ax.axhline(0, color="#9aa3af", linewidth=0.8)
            ax.set_title(f"{module} | {dp}")
            ax.set_ylabel("Module eigengene / program score")
            ax.grid(axis="y", alpha=0.2)
            ax.legend(frameon=False)
    fig.suptitle("Exact-WGCNA module projection across WT/MUT and inhibitor conditions", fontsize=16, fontweight="bold")
    fig.text(0.5, 0.01, "Module-score sign is orientation-dependent; compare within-module condition shifts, not absolute sign between modules.", ha="center", color=COLORS["muted"])
    fig.tight_layout(rect=[0.02, 0.05, 0.98, 0.94])
    _save(fig, path)


def _plot_module_occlusion(table: pd.DataFrame, path: Path, focus: list[str]) -> None:
    part = table.loc[table["module"].isin(focus)].copy()
    part["label"] = part["module"] + " | " + part["DP_status"] + " | " + part["inhibitor"]
    x = np.arange(len(part))
    fig, axes = plt.subplots(2, 1, figsize=(14, 9), sharex=True)
    width = 0.36
    axes[0].bar(x - width / 2, part["WT_inhibitor_effect"], width, label="WT inhibitor effect", color="#7b8797")
    axes[0].bar(x + width / 2, part["MUT_inhibitor_effect"], width, label="MUT inhibitor effect", color=COLORS["MUT"])
    axes[0].axhline(0, color=COLORS["ink"], linewidth=0.8)
    axes[0].set_ylabel("Inhibitor minus baseline module score")
    axes[0].set_title("Occlusion: is the added inhibitor response attenuated in MUT?")
    axes[0].legend(frameon=False)
    axes[1].bar(x - width / 2, part["baseline_MUT_minus_WT"].abs(), width, label="abs(MUT-WT), baseline", color="#b7c1ce")
    axes[1].bar(x + width / 2, part["inhibitor_MUT_minus_WT"].abs(), width, label="abs(MUT-WT), inhibitor", color=COLORS["Pla"])
    axes[1].set_ylabel("Absolute genotype gap in module score")
    axes[1].set_title("Convergence: does inhibitor reduce the MUT-WT module-score gap?")
    axes[1].legend(frameon=False)
    axes[1].set_xticks(x, part["label"], rotation=55, ha="right")
    for ax in axes:
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle("Exact-WGCNA module occlusion and convergence diagnostics", fontsize=16, fontweight="bold")
    fig.text(0.5, 0.01, "Formal evidence requires the difference-in-differences interaction FDR; bar similarity alone is descriptive.", ha="center", color=COLORS["muted"])
    fig.tight_layout(rect=[0.02, 0.06, 0.98, 0.94])
    _save(fig, path)


def _plot_module_core_mapping(core: pd.DataFrame, triple: pd.DataFrame, path: Path) -> None:
    mapped = core.copy()
    mapped["module"] = mapped["module"].fillna("not_in_exact_network")
    counts = mapped.groupby(["context", "inhibitor", "module"], observed=True).size().rename("events").reset_index()
    top_modules = counts.groupby("module")["events"].sum().nlargest(8).index
    counts["module_plot"] = np.where(counts["module"].isin(top_modules), counts["module"], "other_modules")
    pivot = counts.pivot_table(index=["context", "inhibitor"], columns="module_plot", values="events", aggfunc="sum", fill_value=0)
    fig, ax = plt.subplots(figsize=(12, 6))
    pivot.plot(kind="bar", stacked=True, ax=ax, colormap="tab20")
    ax.set_title("Where do exact MUT-inhibitor shared significant events map in the exact-WGCNA network?", fontweight="bold")
    ax.set_ylabel("Shared significant events")
    ax.set_xlabel("Context and LOF-reference inhibitor")
    ax.legend(title="Exact-WGCNA module", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    ax.grid(axis="y", alpha=0.2)
    fig.text(0.5, 0.01, f"Triple-core rows available: {len(triple):,}. Network mapping is support context, not proof that the entire module phenocopies inhibition.", ha="center", color=COLORS["muted"])
    fig.tight_layout(rect=[0.02, 0.05, 0.84, 0.96])
    _save(fig, path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Close the Setd2 MUT LOF-like splicing hypothesis with event and exact-WGCNA evidence.")
    parser.add_argument("--project-root", type=Path, default=Path(r"E:\Project\RNASEQ\Result\Setd2_results"))
    parser.add_argument("--rmats-mode", choices=["JC", "JCEC"], default="JC")
    args = parser.parse_args()

    project = args.project_root
    network_root = project / "katana_output" / "09_cosplice_network"
    exact_root = network_root / "exact_wgcna_full_stable_mut_effect"
    output = network_root / "06i_lof_hypothesis_closure"
    plots = output / "plots"
    output.mkdir(parents=True, exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)

    comparisons = {
        "DP": {"MUT": "MUT_DP-WT_DP", "CLK3": "WT_DP_CLK3-WT_DP", "Pla": "WT_DP_Pla-WT_DP"},
        "noDP": {"MUT": "MUT_NO_DP-WT_NO_DP", "CLK3": "WT_NO_DP_CLK3-WT_NO_DP", "Pla": "WT_NO_DP_Pla-WT_NO_DP"},
    }
    comparison_ids = sorted({value for mapping in comparisons.values() for value in mapping.values()})
    rmats_root = project / "Setd2_Splicing_genomebam" / "star_salmon" / "rmats"
    rmats = {comparison: _load_rmats_comparison(rmats_root, comparison, args.rmats_mode) for comparison in comparison_ids}

    assignment_path = exact_root / "event_module_assignments.tsv"
    raw_assignments = _read_tsv(assignment_path).rename(columns={"wgcna_module": "module"})
    raw_assignments = raw_assignments[["event_key", "module"]].drop_duplicates("event_key")
    event_metadata = _read_tsv(network_root / "02_matrix" / "event_metadata.tsv")
    summary, core, triple = _event_core_tables(comparisons, rmats, raw_assignments)
    core_module_enrichment = _core_module_enrichment(core, triple, raw_assignments)
    core_genes, core_event_types, core_functional_enrichment, core_enrichment_providers = _core_functional_tables(
        comparisons,
        rmats,
        core,
        triple,
    )
    summary.to_csv(output / "exact_shared_event_lof_summary.tsv", sep="\t", index=False)
    core.to_csv(output / "exact_shared_significant_event_cores.tsv", sep="\t", index=False)
    triple.to_csv(output / "exact_triple_significant_event_core.tsv", sep="\t", index=False)
    core_module_enrichment.to_csv(output / "exact_core_module_enrichment.tsv", sep="\t", index=False)
    core_genes.to_csv(output / "exact_core_gene_membership.tsv", sep="\t", index=False)
    core_event_types.to_csv(output / "exact_core_event_type_summary.tsv", sep="\t", index=False)
    core_functional_enrichment.to_csv(output / "exact_core_go_kegg_enrichment.tsv", sep="\t", index=False)
    (output / "exact_core_go_kegg_provider.json").write_text(
        json.dumps(core_enrichment_providers, indent=2), encoding="utf-8"
    )

    psi = _read_tsv(network_root / "03_stable_as_qc" / "stable_psi_matrix.tsv", index_col=0)
    psi.index = psi.index.astype(str)
    samples = _read_tsv(network_root / "00_config" / "factorial_sample_metadata.tsv")
    scores = _read_tsv(exact_root / "module_program_scores.tsv")
    factor_settings = FactorContrastConfig(
        enabled=True,
        primary_factor="genotype",
        target_level="MUT",
        reference_level="WT",
        context_factors=["DP_status"],
        perturbation_factor="treatment",
        baseline_level="Baseline",
        perturbation_levels=["CLK3", "Pla"],
        network_source="exact_wgcna_full_stable_mut_effect",
        fdr_cutoff=0.05,
        min_abs_delta_psi=0.10,
        min_samples_per_group=3,
    )
    factor = analyze_factor_contrasts(
        psi,
        scores.set_index("sample_id"),
        samples,
        event_metadata,
        factor_settings,
    )
    factor.event_contrasts.to_csv(output / "replicate_event_factor_contrasts.tsv", sep="\t", index=False)
    factor.module_contrasts.to_csv(output / "exact_module_factor_contrasts.tsv", sep="\t", index=False)
    condition, occlusion = _condition_and_occlusion_tables(scores, samples, factor.module_contrasts)
    condition.to_csv(output / "exact_module_condition_scores.tsv", sep="\t", index=False)
    occlusion.to_csv(output / "exact_module_occlusion_convergence.tsv", sep="\t", index=False)

    hub_events, hub_genes = _hub_tables(psi, scores, raw_assignments, event_metadata)
    hub_events.to_csv(output / "exact_module_hub_events.tsv", sep="\t", index=False)
    hub_genes.to_csv(output / "exact_module_hub_genes.tsv", sep="\t", index=False)
    rbp, motif = _rbp_tables(
        scores,
        hub_genes,
        network_root / "06e_module_interpretation" / "rbp_gene_expression_matrix.tsv",
        network_root / "06e_module_interpretation" / "event_rbp_motif_hits.tsv",
        raw_assignments,
    )
    rbp.to_csv(output / "exact_module_rbp_regulators.tsv", sep="\t", index=False)
    motif.to_csv(output / "exact_module_rbp_motif_enrichment.tsv", sep="\t", index=False)
    pathways = _read_tsv(exact_root / "module_go_kegg_enrichment.tsv")
    pathways.to_csv(output / "exact_module_pathway_enrichment.tsv", sep="\t", index=False)
    module_annotation = _module_annotation_summary(raw_assignments, pathways, hub_genes, rbp)
    module_annotation.to_csv(output / "exact_module_annotation_summary.tsv", sep="\t", index=False)
    _write_interpretation_report(
        output / "LOF_HYPOTHESIS_INTERPRETATION.md",
        summary,
        occlusion,
        core_module_enrichment,
        module_annotation,
        core_functional_enrichment,
    )

    _plot_event_evidence(summary, plots / "01_exact_event_lof_evidence.png")
    _plot_module_scores(condition, plots / "02_M03_M06_condition_scores.png", ["M03", "M06"])
    _plot_module_occlusion(occlusion, plots / "03_M03_M06_occlusion_convergence.png", ["M03", "M06"])
    _plot_module_core_mapping(core, triple, plots / "04_shared_core_module_mapping.png")

    input_paths = [assignment_path, exact_root / "module_program_scores.tsv", network_root / "03_stable_as_qc" / "stable_psi_matrix.tsv"]
    input_paths.extend(rmats_root / comparison / "rmats_post" / f"{event_type}.MATS.{args.rmats_mode}.txt" for comparison in comparison_ids for event_type in EVENT_TYPES)
    parameters = {
        "analysis": "Setd2 MUT LOF-like splicing hypothesis closure",
        "scientific_question": "Does MUT resemble WT plus CLK3 inhibitor or WT plus PlaB inhibitor at event and co-splicing-module levels?",
        "rmats_mode": args.rmats_mode,
        "significant_event_rule": "FDR <= 0.05 and abs(IncLevelDifference) >= 0.10",
        "event_matching": "exact canonical event_key: gene_id + event_type + chromosome + strand + rMATS coordinates",
        "module_branch": "exact_wgcna_full_stable_mut_effect",
        "module_branch_fingerprint": _hash_files(input_paths),
        "comparisons": comparisons,
        "occlusion_definition": "[MUT inhibitor - MUT baseline] - [WT inhibitor - WT baseline]",
        "convergence_definition": "descriptive reduction in abs(MUT-WT module-score gap) under inhibitor; formal inference uses the interaction FDR",
        "generated_tables": sorted(path.name for path in output.glob("*.tsv")),
        "generated_plots": sorted(path.name for path in plots.glob("*.png")),
        "generated_reports": ["LOF_HYPOTHESIS_INTERPRETATION.md"],
        "core_functional_enrichment": {
            "method": "g:Profiler g:GOSt GO/KEGG over-representation with a custom exact-AS-testable gene background",
            "gene_unit": "unique gene_symbol; multiple AS events from one gene are counted once",
            "providers": core_enrichment_providers,
        },
        "limits": [
            "Module eigengene sign is arbitrary; only within-module condition differences are interpreted.",
            "Shared significant-event direction/effect agreement supports partial molecular convergence, not complete phenocopy by itself.",
            "LOF/GOF functional classification still requires orthogonal activity/rescue evidence.",
        ],
    }
    (output / "parameters.json").write_text(json.dumps(parameters, indent=2), encoding="utf-8")
    print(f"Finished: {output}")
    print(summary.to_string(index=False))
    print("\nM03/M06 formal occlusion rows:")
    print(occlusion.loc[occlusion["module"].isin(["M03", "M06"])].to_string(index=False))


if __name__ == "__main__":
    main()
