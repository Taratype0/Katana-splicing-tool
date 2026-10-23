"""Versioned execution adapter around the frozen MotifTrace kernels.

Only execution budgeting differs from frozen analyze_consequences: exact
pre-scoring dimensions and explicitly bounded allocation. Models, matching,
statistics, permutations and the complete correction family are unchanged.
"""
from __future__ import annotations
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from katana_next.science.motiftrace_consequences import (
    POLICY, build_paired_datasets, paired_cluster_family, EVENT_TYPES,
    extract_event_rna, match_event_controls, paired_form_scores, trace_instances,
    validate_models,
)
from katana_next.reports.motif_native_index import render_budgeted_consequences as render_consequences
from katana_next.science.motiftrace_budget import (
    DEFAULT_MAX_TESTS, DEFAULT_MAX_MATRIX_CELLS, execution_limits,
    estimate_paired_family, enforce_paired_budget,
)

EXECUTION_ADAPTER_VERSION = "1.0.1"


def analyze_consequences(complete, selected_ids, models, fasta, destination, *,
                         permutations=999, seed=20261015, check_cancelled=lambda: None,
                         provenance=None, max_tests=DEFAULT_MAX_TESTS,
                         max_matrix_cells=DEFAULT_MAX_MATRIX_CELLS):
    """Write a new immutable result directory. Never rewrite a source result."""
    validate_models(models)
    limits = execution_limits(max_tests, max_matrix_cells)
    needed = {"event_id", "gene_id", "comparison_id", "event_type", "dpsi", "fdr",
              "stable_significant", "measurement_eligible", "numerator_mean_reads", "denominator_mean_reads", "IncLevel1", "IncLevel2"}
    if not needed <= set(complete):
        raise ValueError("Missing full measured-event contract: " + ", ".join(sorted(needed-set(complete))))
    selected_ids = set(map(str, selected_ids))
    if complete.event_id.duplicated().any() or complete.comparison_id.nunique() != 1:
        raise ValueError("One full comparison with unique event IDs is required")
    selected = complete[complete.event_id.isin(selected_ids)]
    if not selected_ids or set(selected.event_id) != selected_ids or not selected.event_type.isin(EVENT_TYPES).all():
        raise ValueError("Select exact supported DAS events from the full measured catalog")
    if not selected.stable_significant.fillna(False).all() or not selected.measurement_eligible.fillna(False).all() or not (np.isfinite(selected.dpsi) & selected.dpsi.ne(0)).all():
        raise ValueError("Foreground must retain existing supported nonzero DAS qualification")
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    plan = {**POLICY, "UTC_before_sequence_scoring": datetime.now(timezone.utc).isoformat(),
            "comparison_id": str(complete.comparison_id.iloc[0]), "permutations": permutations,
            "seed": seed, "minimum_independent_gene_pairs": 20, "window_nt": [10, 20, 50], "step": 1,
            "selected_event_ids": sorted(selected_ids), "models": models, "provenance": provenance or {},
            "execution_budget": limits, "execution_adapter_version": EXECUTION_ADAPTER_VERSION, "statistical_family_split": False,
            "external_exon_anchor_nt": 100, "regulatory_intron_nt": 200,
            "input_event_table_sha256": hashlib.sha256(complete.to_json(orient="split", double_precision=15).encode()).hexdigest()}
    (destination/"ANALYSIS_PLAN_BEFORE_SCORING.json").write_text(json.dumps(plan, indent=2, allow_nan=False), encoding="utf-8")
    selected.to_csv(destination/"FROZEN_SELECTED_EVENTS.tsv", sep="\t", index=False)
    selected_genes = set(complete.loc[complete.stable_significant.fillna(False) | complete.fdr.le(.05) | complete.event_id.isin(selected_ids), "gene_id"])
    candidates = complete[complete.event_id.isin(selected_ids) | (
        complete.event_type.isin(set(selected.event_type)) & complete.measurement_eligible.fillna(False)
        & complete.fdr.gt(.1) & complete.dpsi.abs().lt(.05) & ~complete.gene_id.isin(selected_genes))].copy()
    features = extract_event_rna(candidates, fasta, feature_only=True, check_cancelled=check_cancelled)
    matching = match_event_controls(candidates, selected_ids, features)
    matching.to_csv(destination/"MATCHING_ALL_SELECTED_WITH_LOSSES.tsv", sep="\t", index=False)
    matched_ids = set(matching.loc[matching.status.eq("matched"), "background_event"])
    used = candidates[candidates.event_id.isin(selected_ids | matched_ids)].copy()
    used["role"] = np.where(used.event_id.isin(selected_ids), "foreground", "background")
    rna = extract_event_rna(used, fasta, check_cancelled=check_cancelled)
    estimate = estimate_paired_family(rna, matching, models)
    budget = {"limits": limits, "estimate": estimate, "checked_before_motif_scoring": True}
    (destination/"EXECUTION_BUDGET_BEFORE_SCORING.json").write_text(json.dumps(budget, indent=2, allow_nan=False), encoding="utf-8")
    enforce_paired_budget(estimate, limits, permutations=permutations)
    hits, changes, weighted = trace_instances(rna[rna.role.eq("foreground")], selected, models, check_cancelled=check_cancelled)
    scores = paired_form_scores(rna[rna.role.eq("foreground")], selected, models, check_cancelled=check_cancelled)
    for name, frame in (("ACTUAL_RNA_FORMS", rna), ("ACTUAL_LOCUS_HITS", hits),
                        ("PAIRED_FORM_INSTANCE_CHANGES", changes), ("PER_REPLICATE_SEQUENCE_OPPORTUNITY", weighted),
                        ("PAIRED_FORM_NORMALIZED_SCORES", scores)):
        if not len(frame.columns):
            frame = pd.DataFrame(columns=["event_id", "model_id", "status"])
        frame.to_csv(destination/(name+".tsv.gz"), sep="\t", index=False, compression="gzip")
    with (destination/"ACTUAL_LOCAL_RNA.fa").open("x", encoding="ascii") as stream:
        for i, row in rna.iterrows():
            stream.write(f">rna_{i:06d} event={row.event_id} region={row.region} strand={row.strand} SHA256_DNA={row.sequence_sha256}\n{row.sequence_RNA}\n")
    datasets, pairs = build_paired_datasets(rna, matching, models, **limits, check_cancelled=check_cancelled)
    tests, clusters, null_max = paired_cluster_family(datasets, pairs, permutations=permutations,
        seed=seed, minimum_pairs=20, check_cancelled=check_cancelled)
    if tests.empty:
        tests = pd.DataFrame(columns=["event_type", "direction", "model_id", "method", "region", "window_nt", "raw_P_new_paired_swaps", "q_new_window_family", "status"])
    if clusters.empty:
        clusters = pd.DataFrame(columns=["event_type", "direction", "model_id", "method", "region", "window_nt", "cluster_P_FWER_new_family"])
    tests.to_csv(destination/"NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz", sep="\t", index=False, compression="gzip")
    clusters.to_csv(destination/"NEW_CLUSTER_FWER.tsv", sep="\t", index=False)
    np.save(destination/"MAX_CLUSTER_NULL.npy", null_max, allow_pickle=False)
    summary = {"schema": POLICY["schema"], "comparison_id": plan["comparison_id"], "selected_events": len(selected),
        "selected_event_types": selected.event_type.value_counts().to_dict(), "matched_independent_gene_pairs": len(pairs),
        "foreground_events_including_unmatched": int(rna.loc[rna.role.eq("foreground"), "event_id"].nunique()),
        "eligible_new_window_tests": int(tests.raw_P_new_paired_swaps.notna().sum()),
        "new_window_q_lt_0_05": int(tests.q_new_window_family.lt(.05).sum()),
        "cluster_P_FWER_lt_0_05": int(clusters.cluster_P_FWER_new_family.lt(.05).sum()),
        "cluster_rows": len(clusters), "instance_rows": len(changes),
        "matching_loss_counts": matching.status.value_counts().to_dict(),
        "inference_status": "computed_conditional_on_matching_exchangeability" if tests.raw_P_new_paired_swaps.notna().any() else "not_estimable_insufficient_independent_matched_pairs",
        "historical_results_modified": False, "quantification_rerun": False,
        "CLIP_or_functional_validation": "not_computed_by_this_workflow", "policy": plan}
    (destination/"RESULTS.json").write_text(json.dumps(summary, indent=2, allow_nan=False), encoding="utf-8")
    render_consequences(rna, changes, tests, clusters, summary, destination/"MOTIFTRACE_EVENT_CONSEQUENCES.html", paired_scores=scores)
    members = {}
    for path in sorted(destination.iterdir()):
        if path.is_file():
            with path.open("rb") as stream:
                sha = hashlib.file_digest(stream, "sha256").hexdigest()
            members[path.name] = {"bytes": path.stat().st_size, "sha256": sha}
    (destination/"COMPLETE.json").write_text(json.dumps({"complete": True, "members": members, "quantification_rerun": False}, indent=2), encoding="utf-8")
    return summary

