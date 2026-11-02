"""Production integration of the frozen, publicly validated EventGraph kernels."""
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
from katana_motiftrace.eventgraph.adapter import match_controls
from katana_motiftrace.eventgraph.scanning import scan_family
from katana_next.reports.motif_eventgraph import render_eventgraph_report
from katana_next.science.motiftrace_budget import (
    DEFAULT_MAX_TESTS, DEFAULT_MAX_MATRIX_CELLS, execution_limits,
    estimate_paired_family,
)

EXECUTION_ADAPTER_VERSION = "2.0.0-eventgraph.20261101"


def analyze_consequences(complete, selected_ids, models, fasta, destination, *,
                         permutations=999, seed=20261015, check_cancelled=lambda: None,
                         provenance=None, max_tests=DEFAULT_MAX_TESTS,
                         max_matrix_cells=DEFAULT_MAX_MATRIX_CELLS):
    """Write a new result directory using the frozen EventGraph algorithms."""
    validate_models(models)
    if any(m["kind"] != "PWM" for m in models.values()):
        raise ValueError("EventGraph requires PWM models. Select analysis_engine=legacy for consensus libraries; no model is silently omitted.")
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
            "seed": seed, "minimum_independent_gene_blocks": 20, "analysis_engine": "eventgraph_20261030", "window_nt": [10, 20, 50], "step": 1,
            "selected_event_ids": sorted(selected_ids), "models": models, "provenance": provenance or {},
            "execution_budget": limits, "execution_limits_enforced": False,
            "execution_adapter_version": EXECUTION_ADAPTER_VERSION, "statistical_family_split": False,
            "external_exon_anchor_nt": 100, "regulatory_intron_nt": 200,
            "input_event_table_sha256": hashlib.sha256(complete.to_json(orient="split", double_precision=15).encode()).hexdigest()}
    plan.update(schema="katana.motiftrace_eventgraph/v1", matching="frozen global maximum-cardinality/minimum-cost event matching; original foreground matches protected",
        independence="foreground gene blocks; all matched events retained; conditional joint central symmetry",
        window_family="BH across all estimable scalar/window gene-block P values in one comparison",
        statistic="event signed ranks summed within foreground gene blocks; sum(block sums)/sqrt(sum squared block sums)",
        exchangeability="conditional joint central symmetry within independent foreground gene blocks",
        methods="PWM log2(1+mean odds), no hit cutoff; consensus libraries use the explicit legacy engine",
        short_sequences="all estimable scalar/window positions; unavailable groups recorded in TEST_FAMILY_ESTIMABILITY.tsv",
        historical_results_modified=False)
    (destination/"ANALYSIS_PLAN_BEFORE_SCORING.json").write_text(json.dumps(plan, indent=2, allow_nan=False), encoding="utf-8")
    selected.to_csv(destination/"FROZEN_SELECTED_EVENTS.tsv", sep="\t", index=False)
    selected_genes = set(complete.loc[complete.stable_significant.fillna(False) | complete.fdr.le(.05) | complete.event_id.isin(selected_ids), "gene_id"])
    candidates = complete[complete.event_id.isin(selected_ids) | (
        complete.event_type.isin(set(selected.event_type)) & complete.measurement_eligible.fillna(False)
        & complete.fdr.gt(.1) & complete.dpsi.abs().lt(.05) & ~complete.gene_id.isin(selected_genes))].copy()
    features = extract_event_rna(candidates, fasta, feature_only=True, check_cancelled=check_cancelled)
    if candidates.event_id.isin(selected_ids).all():
        # No assignment or hypothesis exists without a background pool. Keep
        # every selected event explicit; avoid the legacy fast matcher's empty
        # boolean-array limitation without changing either matching algorithm.
        if selected.gene_id.isna().any() or selected.gene_id.astype(str).str.strip().eq('').any():
            raise ValueError('Gene identity required to define dependence')
        legacy = match_event_controls(candidates, selected_ids, features)
        matching = legacy.copy()
        matching['status'] = 'unmatched'
        matching['reason'] = 'no_compatible_background'
        matching['candidate_background_events'] = 0
        matching['match_cost'] = np.nan
        certificate = {'maximum_cardinality_certificate': 0, 'candidate_gene_edges': 0,
                       'protected_original_events': 0, 'total_matched_cost': 0.0,
                       'status': 'empty_background_pool_no_assignment'}
    else:
        matching, certificate, legacy = match_controls(candidates, selected_ids, features)
    (destination/"MATCHING_CERTIFICATE.json").write_text(json.dumps(certificate, indent=2), encoding="utf-8")
    legacy.to_csv(destination/"LEGACY_MATCHING_REFERENCE.tsv", sep="\t", index=False)
    matching.to_csv(destination/"MATCHING_ALL_SELECTED_WITH_LOSSES.tsv", sep="\t", index=False)
    matched_ids = set(matching.loc[matching.status.eq("matched"), "background_event"])
    used = candidates[candidates.event_id.isin(selected_ids | matched_ids)].copy()
    used["role"] = np.where(used.event_id.isin(selected_ids), "foreground", "background")
    rna = extract_event_rna(used, fasta, check_cancelled=check_cancelled)
    estimate = {"selected_events": len(selected), "matched_events": int(matching.status.eq("matched").sum()),
                "models": len(models), "pre_scoring_only": True}
    budget = {"limits": limits, "estimate": estimate, "checked_before_motif_scoring": True,
              "enforced": False, "capacity_policy": "no_fixed_computation_limits"}
    (destination/"EXECUTION_BUDGET_BEFORE_SCORING.json").write_text(json.dumps(budget, indent=2, allow_nan=False), encoding="utf-8")
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
    statistics, clusters = scan_family(rna, matching, models, destination,
        permutations=permutations, seed=seed, check_cancelled=check_cancelled)
    summary = {"schema": plan["schema"], "comparison_id": plan["comparison_id"], "selected_events": len(selected),
        "selected_event_types": selected.event_type.value_counts().to_dict(), "matched_independent_gene_pairs": None, "analysis_engine": "eventgraph_20261030", **statistics,
        "foreground_events_including_unmatched": int(rna.loc[rna.role.eq("foreground"), "event_id"].nunique()),
        "cluster_rows": len(clusters), "instance_rows": len(changes),
        "matching_loss_counts": matching.status.value_counts().to_dict(),
        "inference_status": "computed_conditional_on_gene_block_symmetry" if statistics["eligible_new_window_tests"] else "not_estimable_insufficient_independent_gene_blocks",
        "historical_results_modified": False, "quantification_rerun": False,
        "CLIP_or_functional_validation": "not_computed_by_this_workflow", "policy": plan}
    (destination/"RESULTS.json").write_text(json.dumps(summary, indent=2, allow_nan=False), encoding="utf-8")
    render_eventgraph_report(summary, clusters, destination/"MOTIFTRACE_EVENT_CONSEQUENCES.html")
    members = {}
    for path in sorted(destination.iterdir()):
        if path.is_file():
            with path.open("rb") as stream:
                sha = hashlib.file_digest(stream, "sha256").hexdigest()
            members[path.name] = {"bytes": path.stat().st_size, "sha256": sha}
    (destination/"COMPLETE.json").write_text(json.dumps({"complete": True, "members": members, "quantification_rerun": False, "analysis_engine": "eventgraph_20261030", "statistical_family_split": False}, indent=2), encoding="utf-8")
    return summary

