from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import linear_sum_assignment
from scipy.stats import hypergeom, norm

from katana_cosplice.config import ExternalCohortConfig, MetaAnalysisConfig
from katana_cosplice.event_harmonization import annotation_aware_event_alignment, exact_event_alignment
from katana_cosplice.utils import bh_adjust


@dataclass
class CrossCohortMetaResult:
    event_alignment: pd.DataFrame
    module_matching: pd.DataFrame
    module_state: pd.DataFrame
    cohort_effects: pd.DataFrame
    random_effects_meta: pd.DataFrame
    cohort_audit: pd.DataFrame
    multi_cohort_consensus: pd.DataFrame


def _module_column(frame: pd.DataFrame) -> str:
    for column in ("consensus_module", "module"):
        if column in frame:
            return column
    raise ValueError("Module assignments require consensus_module or module.")


def _load_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="," if path.suffix.lower() == ".csv" else "\t", low_memory=False)


def _jaccard_matching(
    reference: pd.DataFrame,
    target: pd.DataFrame,
    cohort: str,
    event_universe: set[str],
    fdr_cutoff: float,
) -> pd.DataFrame:
    ref_col = _module_column(reference)
    target_col = _module_column(target)
    reference_sets = {
        str(module): set(part["event_key"].astype(str)) & event_universe
        for module, part in reference.loc[reference[ref_col].astype(str).ne("grey")].groupby(ref_col, observed=True)
    }
    target_sets = {
        str(module): set(part["event_key"].astype(str)) & event_universe
        for module, part in target.loc[target[target_col].astype(str).ne("grey")].groupby(target_col, observed=True)
    }
    if not reference_sets or not target_sets:
        return pd.DataFrame()
    ref_names = sorted(reference_sets)
    target_names = sorted(target_sets)
    similarity = np.zeros((len(ref_names), len(target_names)), dtype=float)
    overlap = np.zeros_like(similarity)
    for i, left in enumerate(ref_names):
        for j, right in enumerate(target_names):
            shared = reference_sets[left] & target_sets[right]
            union = reference_sets[left] | target_sets[right]
            similarity[i, j] = len(shared) / len(union) if union else 0.0
            overlap[i, j] = len(shared)
    all_pairs: list[dict[str, object]] = []
    for i, left in enumerate(ref_names):
        for j, right in enumerate(target_names):
            shared_count = int(overlap[i, j])
            left_size = len(reference_sets[left])
            right_size = len(target_sets[right])
            pvalue = float(
                hypergeom.sf(shared_count - 1, len(event_universe), left_size, right_size)
            ) if shared_count and event_universe else 1.0
            all_pairs.append({
                "cohort": cohort,
                "reference_module": left,
                "target_module": right,
                "shared_events": shared_count,
                "jaccard": float(similarity[i, j]),
                "reference_module_size": left_size,
                "target_module_size": right_size,
                "reference_event_coverage": shared_count / left_size if left_size else 0.0,
                "target_event_coverage": shared_count / right_size if right_size else 0.0,
                "overlap_pvalue": pvalue,
            })
    pair_table = pd.DataFrame(all_pairs)
    pair_table["overlap_FDR"] = bh_adjust(pair_table["overlap_pvalue"])
    rows, columns = linear_sum_assignment(1.0 - similarity)
    selected = pair_table.loc[
        [
            (pair_table["reference_module"].eq(ref_names[i]) & pair_table["target_module"].eq(target_names[j])).idxmax()
            for i, j in zip(rows, columns, strict=False)
        ]
    ].copy()
    selected["significant_exact_event_overlap"] = selected["overlap_FDR"].le(fdr_cutoff)
    selected["matching_method"] = (
        "one-to-one Hungarian matching on exact-event Jaccard; overlap evidence from hypergeometric test with BH-FDR"
    )
    return selected.reset_index(drop=True)


def _fit_module_projection(psi: pd.DataFrame, assignments: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, dict[str, object]]]:
    module_col = _module_column(assignments)
    scores: dict[str, pd.Series] = {}
    models: dict[str, dict[str, object]] = {}
    for module, part in assignments.loc[assignments[module_col].astype(str).ne("grey")].groupby(module_col, observed=True):
        events = psi.index.intersection(part["event_key"].astype(str))
        if len(events) < 3:
            continue
        values = psi.loc[events].T.apply(pd.to_numeric, errors="coerce")
        means = values.mean(axis=0)
        scales = values.std(axis=0).replace(0, np.nan)
        standardized = ((values - means) / scales).fillna(0.0)
        _, _, right = np.linalg.svd(standardized.to_numpy(dtype=float), full_matrices=False)
        loadings = right[0]
        raw_score = standardized.to_numpy(dtype=float) @ loadings
        mean_profile = standardized.mean(axis=1).to_numpy(dtype=float)
        if np.corrcoef(raw_score, mean_profile)[0, 1] < 0:
            loadings = -loadings
            raw_score = -raw_score
        scores[str(module)] = pd.Series(raw_score, index=values.index)
        models[str(module)] = {
            "events": events.tolist(),
            "loadings": loadings,
        }
    return pd.DataFrame(scores), models


def _project_module_scores(psi: pd.DataFrame, models: dict[str, dict[str, object]]) -> pd.DataFrame:
    scores: dict[str, pd.Series] = {}
    for module, model in models.items():
        events = pd.Index(model["events"]).intersection(psi.index)
        if len(events) < 3:
            continue
        positions = [model["events"].index(event) for event in events]
        loadings = np.asarray(model["loadings"], dtype=float)[positions]
        norm = np.linalg.norm(loadings)
        if norm <= 0:
            continue
        loadings = loadings / norm
        values = psi.loc[events].T.apply(pd.to_numeric, errors="coerce")
        standardized = ((values - values.mean(axis=0)) / values.std(axis=0).replace(0, np.nan)).fillna(0.0)
        scores[module] = pd.Series(standardized.to_numpy(dtype=float) @ loadings, index=values.index)
    return pd.DataFrame(scores)


def _hedges_effects(scores: pd.DataFrame, metadata: pd.DataFrame, cohort: str) -> pd.DataFrame:
    if "sample_id" not in metadata or "condition" not in metadata:
        return pd.DataFrame()
    metadata = metadata.set_index("sample_id").reindex(scores.index)
    levels = [str(value) for value in metadata["condition"].dropna().unique()]
    if len(levels) < 2:
        return pd.DataFrame()
    baseline = levels[0]
    rows: list[dict[str, object]] = []
    for level in levels[1:]:
        for module in scores.columns:
            left = pd.to_numeric(scores.loc[metadata["condition"].astype(str).eq(level), module], errors="coerce").dropna()
            right = pd.to_numeric(scores.loc[metadata["condition"].astype(str).eq(baseline), module], errors="coerce").dropna()
            if len(left) < 2 or len(right) < 2:
                continue
            pooled_df = len(left) + len(right) - 2
            pooled_sd = np.sqrt(((len(left) - 1) * left.var(ddof=1) + (len(right) - 1) * right.var(ddof=1)) / pooled_df)
            if not np.isfinite(pooled_sd) or pooled_sd <= 0:
                continue
            d = float((left.mean() - right.mean()) / pooled_sd)
            correction = 1.0 - 3.0 / (4.0 * pooled_df - 1.0) if pooled_df > 1 else 1.0
            g = correction * d
            variance = (len(left) + len(right)) / (len(left) * len(right)) + (g * g) / (2.0 * pooled_df)
            rows.append({
                "cohort": cohort, "module": module, "contrast": f"{level}_vs_{baseline}",
                "effect_size_hedges_g": g, "standard_error": float(np.sqrt(variance)),
                "n_experiment": len(left), "n_reference": len(right),
            })
    return pd.DataFrame(rows)


def _random_effects(frame: pd.DataFrame, settings: MetaAnalysisConfig) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if frame.empty:
        return pd.DataFrame()
    for (module, contrast), part in frame.groupby(["module", "contrast"], observed=True):
        part = part.dropna(subset=["effect_size_hedges_g", "standard_error"])
        if len(part) < settings.min_cohorts:
            continue
        effect = part["effect_size_hedges_g"].to_numpy(dtype=float)
        variance = np.square(part["standard_error"].to_numpy(dtype=float))
        fixed_weights = 1.0 / variance
        fixed = float(np.sum(fixed_weights * effect) / np.sum(fixed_weights))
        q = float(np.sum(fixed_weights * np.square(effect - fixed)))
        df = len(effect) - 1
        c = float(np.sum(fixed_weights) - np.sum(np.square(fixed_weights)) / np.sum(fixed_weights))
        tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
        weights = 1.0 / (variance + tau2)
        pooled = float(np.sum(weights * effect) / np.sum(weights))
        se = float(np.sqrt(1.0 / np.sum(weights)))
        z = pooled / se if se > 0 else np.nan
        pvalue = float(2.0 * norm.sf(abs(z))) if np.isfinite(z) else np.nan
        i2 = max(0.0, (q - df) / q) if q > 0 else 0.0
        rows.append({
            "module": module, "contrast": contrast, "n_cohorts": len(part),
            "random_effect_hedges_g": pooled, "standard_error": se,
            "ci_low": pooled - 1.96 * se, "ci_high": pooled + 1.96 * se,
            "z": z, "pvalue": pvalue, "Q": q, "Q_df": df, "I2": i2,
            "tau2_DerSimonian_Laird": tau2,
            "method": "DerSimonian-Laird random-effects meta-analysis",
        })
    result = pd.DataFrame(rows)
    if not result.empty:
        result["FDR"] = bh_adjust(result["pvalue"])
        result["significant_at_configured_FDR"] = result["FDR"].le(settings.fdr_cutoff)
        result = result.sort_values(["FDR", "pvalue"], na_position="last")
    return result


def meta_analyze_modules(
    reference_psi: pd.DataFrame,
    reference_assignments: pd.DataFrame,
    reference_metadata: pd.DataFrame,
    external_cohorts: list[ExternalCohortConfig],
    settings: MetaAnalysisConfig,
    *,
    reference_name: str = "discovery",
    reference_genome_build: str = "",
) -> CrossCohortMetaResult:
    alignments: list[pd.DataFrame] = []
    matching: list[pd.DataFrame] = []
    states: list[dict[str, object]] = []
    effects: list[pd.DataFrame] = []
    audits: list[dict[str, object]] = []
    reference_psi.index = reference_psi.index.astype(str)
    reference_scores, projection_models = _fit_module_projection(reference_psi, reference_assignments)
    effects.append(_hedges_effects(reference_scores, reference_metadata, reference_name))
    reference_modules = set(reference_assignments[_module_column(reference_assignments)].astype(str)) - {"grey"}

    if not external_cohorts:
        audits.append({
            "cohort": "",
            "status": "not_configured",
            "message": "No independent external cohort configured; cross-cohort preservation and meta-analysis were not assessed.",
        })

    for cohort in external_cohorts:
        if not cohort.psi_matrix.exists() or not cohort.sample_metadata.exists():
            audits.append({"cohort": cohort.name, "status": "missing_input", "message": "PSI matrix or sample metadata missing"})
            continue
        target_psi = pd.read_csv(
            cohort.psi_matrix,
            sep="," if cohort.psi_matrix.suffix.lower() == ".csv" else "\t",
            index_col=0,
        )
        target_psi.index = target_psi.index.astype(str)
        target_metadata = _load_table(cohort.sample_metadata)
        explicit_mapping = (
            _load_table(cohort.event_mapping_table)
            if cohort.event_mapping_table is not None and cohort.event_mapping_table.exists()
            else None
        )
        alignment = (
            exact_event_alignment(reference_psi.index, target_psi.index, cohort.name)
            if settings.event_matching == "exact"
            else annotation_aware_event_alignment(
                reference_psi.index, target_psi.index, cohort.name,
                reference_genome_build=reference_genome_build,
                target_genome_build=cohort.genome_build,
                mapping_table=explicit_mapping,
            )
        )
        alignments.append(alignment.mapping)
        shared = alignment.shared_events
        if len(shared) < settings.min_shared_events:
            audits.append({"cohort": cohort.name, "status": "insufficient_shared_events", "message": f"{len(shared)} exact shared events"})
            continue
        target_aligned = target_psi.loc[[alignment.reference_to_target[event] for event in shared]].copy()
        target_aligned.index = shared
        projected_scores = _project_module_scores(target_aligned, projection_models)
        effects.append(_hedges_effects(projected_scores, target_metadata, cohort.name))
        if cohort.effect_table is not None and cohort.effect_table.exists():
            external_effects = _load_table(cohort.effect_table)
            required = {"module", "contrast", "effect_size_hedges_g", "standard_error"}
            if required.issubset(external_effects):
                external_effects["cohort"] = cohort.name
                effects.append(external_effects[["cohort", *sorted(required)]])
        if cohort.module_assignments is None or not cohort.module_assignments.exists():
            for module in sorted(reference_modules):
                states.append({
                    "cohort": cohort.name, "reference_module": module, "target_module": "",
                    "module_state": "not_assessed", "reason": "independent module assignments not configured",
                })
            audits.append({"cohort": cohort.name, "status": "projected_effects_only", "message": "Independent modules unavailable; gained/lost/rewired not inferred"})
            continue
        target_assignments = _load_table(cohort.module_assignments)
        target_assignments = target_assignments.copy()
        target_assignments["event_key"] = target_assignments["event_key"].astype(str).map(
            {target: reference for reference, target in alignment.reference_to_target.items()}
        )
        target_assignments = target_assignments.dropna(subset=["event_key"])
        matched = _jaccard_matching(
            reference_assignments,
            target_assignments,
            cohort.name,
            set(map(str, shared)),
            settings.fdr_cutoff,
        )
        matching.append(matched)
        significant_matches = matched.loc[matched["significant_exact_event_overlap"]] if not matched.empty else pd.DataFrame()
        matched_target = set(significant_matches["target_module"]) if not significant_matches.empty else set()
        target_modules = set(target_assignments[_module_column(target_assignments)].astype(str)) - {"grey"}
        preservation = (
            _load_table(cohort.preservation_table)
            if cohort.preservation_table is not None and cohort.preservation_table.exists()
            else pd.DataFrame()
        )
        for module in sorted(reference_modules):
            row = matched.loc[matched["reference_module"].eq(module)].head(1) if not matched.empty else pd.DataFrame()
            jaccard = float(row.iloc[0]["jaccard"]) if not row.empty else 0.0
            significant_overlap = bool(row.iloc[0]["significant_exact_event_overlap"]) if not row.empty else False
            preservation_row = (
                preservation.loc[preservation["consensus_module"].astype(str).eq(module)].head(1)
                if not preservation.empty and "consensus_module" in preservation
                else pd.DataFrame()
            )
            connectivity = float(preservation_row.iloc[0].get("connectivity_correlation", np.nan)) if not preservation_row.empty else np.nan
            connectivity_p = float(preservation_row.iloc[0].get("connectivity_empirical_pvalue", np.nan)) if not preservation_row.empty else np.nan
            connectivity_fdr = float(preservation_row.iloc[0].get("connectivity_empirical_FDR", np.nan)) if not preservation_row.empty else np.nan
            connectivity_significance = connectivity_fdr if np.isfinite(connectivity_fdr) else connectivity_p
            rewired = bool(
                np.isfinite(connectivity)
                and connectivity < 0
                and np.isfinite(connectivity_significance)
                and connectivity_significance <= settings.fdr_cutoff
            )
            if rewired:
                module_state = "rewired_supported_by_connectivity"
            elif significant_overlap:
                module_state = "consensus"
            else:
                module_state = "lost_or_not_recovered"
            states.append({
                "cohort": cohort.name, "reference_module": module,
                "target_module": str(row.iloc[0]["target_module"]) if not row.empty else "",
                "module_state": module_state,
                "jaccard": jaccard,
                "overlap_FDR": float(row.iloc[0]["overlap_FDR"]) if not row.empty else np.nan,
                "reason": (
                    "negative connectivity with permutation-adjusted support"
                    if rewired and np.isfinite(connectivity_fdr)
                    else "negative connectivity with empirical permutation support"
                    if rewired
                    else "significant exact-event module overlap after BH-FDR"
                    if significant_overlap
                    else "no statistically supported exact-event module match"
                ),
                "rewiring_scope": "rewiring is not inferred from overlap alone; use preservation/connectivity diagnostics",
                "connectivity_correlation": connectivity,
                "connectivity_empirical_pvalue": connectivity_p,
                "connectivity_empirical_FDR": connectivity_fdr,
                "connectivity_significance_used": connectivity_significance,
            })
        for module in sorted(target_modules - matched_target):
            states.append({
                "cohort": cohort.name,
                "reference_module": "",
                "target_module": module,
                "module_state": "gained_or_not_recovered",
                "jaccard": 0.0,
                "reason": "target module lacks a statistically supported exact-event reference match",
                "rewiring_scope": "descriptive structural follow-up; not causal rewiring",
            })
        audits.append({"cohort": cohort.name, "status": "finished", "message": f"{len(shared)} exact shared events; independent modules matched"})

    alignment_columns = [
        "cohort", "event_key", "in_reference", "in_target", "exact_match",
        "reference_event_key", "target_event_key", "matched", "matching_rule",
        "reference_genome_build", "target_genome_build",
    ]
    matching_columns = [
        "cohort", "reference_module", "target_module", "shared_events", "jaccard",
        "reference_module_size", "target_module_size", "reference_event_coverage",
        "target_event_coverage", "overlap_pvalue", "overlap_FDR",
        "significant_exact_event_overlap", "matching_method",
    ]
    state_columns = [
        "cohort", "reference_module", "target_module", "module_state", "jaccard",
        "overlap_FDR", "reason", "rewiring_scope",
        "connectivity_correlation", "connectivity_empirical_pvalue",
        "connectivity_empirical_FDR", "connectivity_significance_used",
    ]
    effect_columns = [
        "cohort", "module", "contrast", "effect_size_hedges_g", "standard_error",
        "n_experiment", "n_reference",
    ]
    meta_columns = [
        "module", "contrast", "n_cohorts", "random_effect_hedges_g", "standard_error",
        "ci_low", "ci_high", "z", "pvalue", "Q", "Q_df", "I2",
        "tau2_DerSimonian_Laird", "method", "FDR", "significant_at_configured_FDR",
    ]
    audit_columns = ["cohort", "status", "message"]
    cohort_effects = (
        pd.concat([frame for frame in effects if not frame.empty], ignore_index=True)
        if any(not frame.empty for frame in effects)
        else pd.DataFrame(columns=effect_columns)
    )
    meta = _random_effects(cohort_effects, settings)
    state_frame = pd.DataFrame(states).reindex(columns=state_columns)
    consensus_rows: list[dict[str, object]] = []
    finished_cohorts = {row["cohort"] for row in audits if row.get("status") == "finished"}
    for module in sorted(reference_modules):
        part = state_frame.loc[state_frame["reference_module"].astype(str).eq(module)]
        supported = int(part["module_state"].eq("consensus").sum())
        rewired = int(part["module_state"].eq("rewired_supported_by_connectivity").sum())
        denominator = len(finished_cohorts)
        fraction = supported / denominator if denominator else np.nan
        rewired_fraction = rewired / denominator if denominator else np.nan
        if denominator == 0:
            consensus_state = "not_assessed"
        elif fraction >= settings.min_supported_cohort_fraction:
            consensus_state = "preserved_across_cohorts"
        elif rewired_fraction >= settings.min_supported_cohort_fraction:
            consensus_state = "rewired_across_cohorts"
        elif supported > 0:
            consensus_state = "partially_preserved"
        elif rewired > 0:
            consensus_state = "partially_rewired"
        else:
            consensus_state = "not_preserved"
        consensus_rows.append({
            "reference_module": module, "n_external_cohorts_assessed": denominator,
            "n_supported_cohorts": supported, "supported_cohort_fraction": fraction,
            "n_connectivity_rewired_cohorts": rewired,
            "rewired_cohort_fraction": rewired_fraction,
            "multi_cohort_state": consensus_state,
            "state_rule": "exact/annotation-aware event matching plus FDR-supported module overlap; rewiring requires negative permutation-supported connectivity",
        })
    return CrossCohortMetaResult(
        event_alignment=(pd.concat(alignments, ignore_index=True) if alignments else pd.DataFrame()).reindex(columns=alignment_columns),
        module_matching=(pd.concat(matching, ignore_index=True) if matching else pd.DataFrame()).reindex(columns=matching_columns),
        module_state=state_frame,
        cohort_effects=cohort_effects.reindex(columns=effect_columns),
        random_effects_meta=meta.reindex(columns=meta_columns),
        cohort_audit=pd.DataFrame(audits).reindex(columns=audit_columns),
        multi_cohort_consensus=pd.DataFrame(consensus_rows),
    )


def render_meta_analysis_plots(result: CrossCohortMetaResult, output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    if not result.module_matching.empty:
        matrix = result.module_matching.pivot_table(index="reference_module", columns=["cohort", "target_module"], values="jaccard", aggfunc="max")
        fig, ax = plt.subplots(figsize=(max(9, matrix.shape[1] * 0.55), max(5, matrix.shape[0] * 0.45)))
        sns.heatmap(matrix, cmap="crest", vmin=0, vmax=1, annot=True, fmt=".2f", ax=ax, cbar_kws={"label": "Exact-event Jaccard"})
        ax.set_title("Cross-cohort CoSplice module matching")
        fig.tight_layout()
        path = output_dir / "cross_cohort_module_matching.png"
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
    if not result.random_effects_meta.empty:
        frame = result.random_effects_meta.nsmallest(40, "FDR").copy()
        frame["label"] = frame["module"].astype(str) + " | " + frame["contrast"].astype(str)
        fig, ax = plt.subplots(figsize=(9, max(5, len(frame) * 0.3)))
        ax.errorbar(frame["random_effect_hedges_g"], np.arange(len(frame)), xerr=1.96 * frame["standard_error"], fmt="o", color="#0f766e", ecolor="#94a3b8")
        ax.axvline(0, color="#1f2937", linewidth=1)
        ax.set_yticks(np.arange(len(frame)), frame["label"])
        ax.set_xlabel("Random-effects Hedges g (95% CI)")
        ax.set_title("Cross-cohort module program-score meta-analysis")
        fig.tight_layout()
        path = output_dir / "cross_cohort_module_effect_meta.png"
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
    return paths
