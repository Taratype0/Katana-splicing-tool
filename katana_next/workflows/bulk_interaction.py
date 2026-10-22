from itertools import combinations
from dataclasses import replace
import numpy as np
import pandas as pd

from katana_next.science.bulk_interaction import interaction_observations, fit_read_interaction
from katana_next.science.statistics import bh_adjust
from katana_next.workflows.base import WorkflowBase, WorkflowDefinition, ParameterDefinition, ArtifactWriter
from katana_next.workflows.bulk import _latest_table, _plot
from katana_next.workflows.sources import selection_rows, selected_source


class BulkInteractionWorkflow(WorkflowBase):
    version = "2.3.0"
    definition = WorkflowDefinition(
        workflow_id="differential.replicate_interaction", workspace_id="differential",
        label="Exploratory beta-binomial contrast interaction (legacy)", supports_contrast_selection=True,
        question="Does the exact event's splicing effect differ between two declared contrasts?",
        method="Beta-binomial likelihood-ratio test of (A numerator - A denominator) - (B numerator - B denominator); unique sample observations and rMATS effective-length offset.",
        output_summary="One filterable union table: A-specific, B-specific, shared, not significant and failed/not-testable statuses; model diagnostics, raw support and per-pair BH FDR.",
        analysis_unit="exact event, contrast pair, independent biological sample",
        parameters=(ParameterDefinition("min_reads", "Minimum IJC + SJC per sample", "number", 10, "Applied before model fitting in every unique arm.", minimum=1, step=1),
                    ParameterDefinition("min_replicates", "Supported replicates per arm", "number", 2, "Independent biological samples, not repeated comparison observations.", minimum=2, step=1),
                    ParameterDefinition("event_set", "Caller classification", "select", "stable_significant", "Classification only; models test all eligible exact common events, not only significant ones.", ("stable_significant", "significant")),
                    ParameterDefinition("selection_id", "Optional exact event selection", "selection", "", "Empty tests every available event. Selected-event results are exploratory: within-selection BH does not correct for choosing events using the same DAS data.")),
        limitations="Asymptotic one-df LRT, not a donor GLMM. Shared identical arms are deduplicated; partial overlaps, conflicting counts and repeated biological units are not silently treated as independent. An exclusive caller result is not proof of zero effect in the other contrast. Event selection makes this an exploratory restricted-family analysis, not a selection-adjusted or genome-wide FDR claim.",
    )

    def parent_artifacts(self, context, parameters):
        if parameters.get("selection_id"):
            parent = selected_source(context, parameters)
            if parent is None or parent.analysis_id != "differential.as":
                raise ValueError("Interaction requires an exact-event selection from Differential AS")
            return [parent]
        parent = _latest_table(context, "differential.as", "Differential AS events", parameters)
        return [parent] if parent else []

    def input_paths(self, context, parameters):
        return [a.path for a in self.parent_artifacts(context, parameters)]

    def run(self, context, job, run_fingerprint, parameters):
        parents = self.parent_artifacts(context, parameters)
        if not parents:
            raise ValueError("Run Differential AS first; interaction never starts upstream workflows.")
        events = pd.read_parquet(parents[0].path)
        if parameters.get("selection_id"):
            events = selection_rows(context, parameters, parents[0], events)
        ids = parameters.get("contrast_ids") or list(events.comparison_id.unique())
        contrasts = {c.contrast_id: c for c in context.manifest.contrasts}
        if len(ids) != 2:
            raise ValueError("Select exactly two biologically justified contrasts; automatic all-pair fitting is disabled")
        frames = {c: events.loc[events.comparison_id.eq(c)].set_index("event_id", verify_integrity=True) for c in ids}
        tasks = [(a, b, sorted(set(frames[a].index) | set(frames[b].index))) for a, b in combinations(ids, 2)]
        rows, observations = [], []
        total = sum(len(keys) for _, _, keys in tasks)
        for a, b, keys in tasks:
            for event in keys:
                left = frames[a].loc[event] if event in frames[a].index else None
                right = frames[b].loc[event] if event in frames[b].index else None
                key = f"{a}::{b}::{event}"
                sa = left is not None and pd.notna(left.get(parameters["event_set"])) and bool(left.get(parameters["event_set"], False))
                sb = right is not None and pd.notna(right.get(parameters["event_set"])) and bool(right.get(parameters["event_set"], False))
                row = {"row_id": key, "event_id": event, "comparison_A": a, "comparison_B": b,
                       "contrast_A": contrasts[a].label, "contrast_B": contrasts[b].label,
                       "caller_class": "shared_significant" if sa and sb else "A_specific_significant" if sa else "B_specific_significant" if sb else "neither_significant",
                       "status": "counterpart_not_found", "pvalue": np.nan}
                for side, entry in (("A", left), ("B", right)):
                    for column in ("gene_id", "gene_symbol", "event_type", "dpsi", "fdr"):
                        row[f"{column}_{side}"] = entry.get(column) if entry is not None else None
                if left is not None and right is not None:
                    row["delta_dpsi_A_minus_B"] = left.dpsi - right.dpsi
                    row["direction_class"] = "same" if left.dpsi*right.dpsi > 0 else "opposite" if left.dpsi*right.dpsi < 0 else "zero_or_unknown"
                    try:
                        obs, vector, offset = interaction_observations(left, right, contrasts[a], contrasts[b], context.manifest.samples)
                        row.update(fit_read_interaction(obs, vector, offset, min_reads=parameters["min_reads"], min_replicates=parameters["min_replicates"]))
                        obs["interaction_id"] = key
                        obs["row_id"] = key + "::" + obs.sample_id
                        observations.extend(obs.to_dict("records"))
                    except ValueError as exc:
                        row["status"] = str(exc)
                rows.append(row)
                if len(rows) % 20 == 0:
                    job.check_cancelled()
                    job.report(.05 + .85*len(rows)/max(total, 1), f"Beta-binomial event-pair {len(rows)}/{total}; shared samples counted once")
        result = pd.DataFrame(rows)
        if result.empty:
            raise ValueError("No events available in the selected contrast family")
        result["interaction_fdr"] = result.groupby(["comparison_A", "comparison_B"])["pvalue"].transform(lambda s: bh_adjust(s))
        scope = "restricted_selection_exploratory" if parameters.get("selection_id") else "all_available_eligible_events_per_pair"
        result["inference_scope"] = scope
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
                                common_metadata={"method": self.definition.method, "limitations": self.definition.limitations,
                                                 "inference_scope": scope,
                                                 "fdr_family": "all successfully modeled events within each declared pair; not just caller-significant events"})
        parent_ids = tuple(a.artifact_id for a in parents)
        table = writer.table("All replicate-aware event interactions", "replicate_interactions.parquet", result, entity_type="splice_event", parents=parent_ids)
        if observations:
            writer.table("Deduplicated sample read support", "interaction_observations.parquet", pd.DataFrame(observations), entity_type="sample_event", parents=(table.artifact_id,))
        summary = result.groupby(["contrast_A", "contrast_B", "caller_class", "status"], dropna=False).size().reset_index(name="events")
        summary["row_id"] = summary.index.astype(str)
        writer.table("Interaction execution and caller classes", "interaction_summary.parquet", summary, entity_type="contrast_pair", parents=(table.artifact_id,))
        writer.plot("A-specific and B-specific event support", "interaction_classes.vl.json", _plot(summary.to_dict("records"), "bar",
                    {"x": {"field": "caller_class", "type": "nominal", "axis": {"labelAngle": -30}},
                     "y": {"field": "events", "type": "quantitative"}, "color": {"field": "status", "type": "nominal"},
                     "row": {"field": "contrast_A", "type": "nominal"}, "column": {"field": "contrast_B", "type": "nominal"}},
                    title="Exact events by caller status and model eligibility", description="Each panel compares the two explicitly named contrasts. A-specific and B-specific categories are separate; caller exclusivity is not an interaction test.", height=240), entity_type="interaction_plot", parents=(table.artifact_id,))
        return writer.finish()


class ExploratoryInteractionWorkflow(BulkInteractionWorkflow):
    definition = replace(BulkInteractionWorkflow.definition,
        workflow_id='differential.exploratory_event_interaction',
        label='Exploratory event interaction (beta-binomial)',
        question='For two explicitly selected contrasts, is the event-level response different?')
