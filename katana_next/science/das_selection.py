"""Selection of completed DAS evidence without changing caller statistics.

The two flags remain separate in every published event table. Legacy kernels
which consume ``stable_significant`` receive a private compatibility copy only.
This contract changes foreground membership, never a caller P/FDR, PSI matrix,
network measurement population, or a completed upstream result.
"""
from __future__ import annotations

from typing import Any

import pandas as pd

from katana_next.science.bulk import DEFAULT_MIN_READS, READ_SUPPORT_POLICY_VERSION


DAS_SELECTION_POLICY_VERSION = "das-event-set-v1"
DEFAULT_EVENT_SET = "stable_significant"
EVENT_SET_LABELS = {"stable_significant": "RS-DAS", "significant": "DAS"}

SWITCHABLE_WORKFLOWS = frozenset({
    "differential.as", "differential.cross_comparison", "prioritize.evidence",
    "prioritize.enrichment", "prioritize.rbp_support",
    "interpretation.joint_response", "interpretation.joint_pathways",
    "differential.exploratory_event_interaction", "reports.bundle", "reports.study",
})
RS_REQUIRED_WORKFLOWS = frozenset({
    "prioritize.rbp_mechanism", "prioritize.event_motifs",
    "prioritize.binding_enrichment", "prioritize.se_rbp_impact",
    "prioritize.motiftrace_consequences",
})
MEASUREMENT_WORKFLOWS = frozenset({
    "quality.stable_as", "landscape.bulk", "programs.cosplice",
    "programs.psi_network", "programs.djexpress", "programs.psi_projection",
    "programs.module_functions", "programs.hub_events",
    "programs.module_explanation", "programs.module_diagnostics",
    "programs.explain_expression", "programs.diagnostics",
    "programs.external_preservation", "programs.external_validation",
    "single_cell.programs", "prioritize.selected_gene_views",
})
# These are registrations/views of completed result collections, not a request
# to compute a new ordinary-DAS foreground from those historical collections.
FIXED_RESULT_WORKFLOWS = frozenset({
    "prioritize.rbp_sequence_results", "prioritize.motiftrace_parallel",
    "prioritize.motiftrace_usage",
})


def normalize_event_set(value: Any = None) -> str:
    token = DEFAULT_EVENT_SET if value is None else str(value).strip().lower()
    aliases = {"rs-das": "stable_significant", "rs_das": "stable_significant",
               "das": "significant", "ordinary": "significant"}
    token = aliases.get(token, token)
    if token not in EVENT_SET_LABELS:
        raise ValueError("event_set must be stable_significant (RS-DAS) or significant (DAS).")
    return token


def workflow_das_capability(workflow_id: str) -> dict[str, Any]:
    mode = ("switchable" if workflow_id in SWITCHABLE_WORKFLOWS else
            "rs_required" if workflow_id in RS_REQUIRED_WORKFLOWS else
            "measurement_population" if workflow_id in MEASUREMENT_WORKFLOWS else
            "not_applicable")
    description = {
        "switchable": "Select RS-DAS (default) or ordinary DAS from unchanged caller flags. Supporting downstream artifacts must use the same selection.",
        "rs_required": "RS-DAS is required for new selected-event motif/CLIP analyses. Ordinary-DAS selections are rejected; completed imported collections retain their recorded source scope.",
        "measurement_population": "Uses its declared measured PSI/read-support population, not the significant DAS/RS-DAS foreground. The DAS switch does not change or refit this population.",
        "not_applicable": "This module does not select a differential-event foreground using the bulk DAS switch.",
    }[mode]
    if workflow_id in FIXED_RESULT_WORKFLOWS:
        description = "Preserves its completed result collection or transcript-usage scope and original P/FDR; the DAS switch does not refilter or relabel it."
    if workflow_id == "prioritize.selected_gene_views":
        description = "Uses the saved gene selection and full measured PSI for heatmaps. Existing gene-card RS-event views keep their original flags; this module does not switch the event foreground."
    if workflow_id == "differential.exploratory_event_interaction":
        description += " The switch changes caller-status classification only; the interaction model continues to test all eligible events in its declared family."
    return {"mode": mode, "default_event_set": DEFAULT_EVENT_SET if mode in {"switchable", "rs_required"} else None,
            "label": "RS-DAS required" if mode == "rs_required" else "DAS event selection",
            "description": description, "policy_version": DAS_SELECTION_POLICY_VERSION}


def normalize_workflow_selection(workflow_id: str, parameters: dict[str, Any]) -> dict[str, Any]:
    values = dict(parameters)
    capability = workflow_das_capability(workflow_id)
    if capability["mode"] in {"switchable", "rs_required"}:
        values["event_set"] = normalize_event_set(values.get("event_set"))
        recorded = values.get("das_selection_policy", DAS_SELECTION_POLICY_VERSION)
        if recorded != DAS_SELECTION_POLICY_VERSION:
            raise ValueError("The requested DAS selection policy is historical or unknown; preserve the old result and start a new explicit downstream selection.")
        values["das_selection_policy"] = DAS_SELECTION_POLICY_VERSION
        if capability["mode"] == "rs_required" and values["event_set"] != DEFAULT_EVENT_SET:
            raise ValueError(f"{workflow_id} requires RS-DAS; event_set=significant is not supported. No ordinary-to-RS fallback is performed.")
    elif "event_set" in values:
        raise ValueError(f"{workflow_id} uses its own measurement/result scope; the DAS event_set switch does not apply.")
    return values


def selection_metadata(event_set: Any = None) -> dict[str, Any]:
    key = normalize_event_set(event_set)
    return {"policy_version": DAS_SELECTION_POLICY_VERSION, "event_set": key,
            "label": EVENT_SET_LABELS[key], "default_event_set": DEFAULT_EVENT_SET,
            "caller_statistics": "preserved; no P/FDR or effect refitting",
            "read_support_applied_to_foreground": key == DEFAULT_EVENT_SET,
            "as_background_eligibility": "caller-testable and read-supported genes" if key == DEFAULT_EVENT_SET else "caller-testable genes without the RS foreground gate",
            "read_support_policy": READ_SUPPORT_POLICY_VERSION,
            "minimum_counts": DEFAULT_MIN_READS,
            "read_support_scope": "IJC + SJC >= 10 in every declared sample in both comparison arms"}


def _flag(values: pd.Series, name: str) -> pd.Series:
    """Reject strings such as 'False' instead of Python's truthy-string cast."""
    observed = values.dropna()
    if not observed.map(lambda item: isinstance(item, (bool, int, float)) and item in (True, False, 0, 1)).all():
        # numpy boolean scalars are accepted by numeric equality, explicitly
        # excluding text and arbitrary objects.
        if not observed.map(lambda item: not isinstance(item, str) and item in (True, False, 0, 1)).all():
            raise ValueError(f"{name} must contain explicit boolean/0/1 flags, not text or arbitrary values.")
    return values.fillna(False).astype(bool)


def selected_event_mask(events: pd.DataFrame, event_set: Any = None) -> pd.Series:
    key = normalize_event_set(event_set)
    if key not in events:
        raise ValueError(f"Completed event source lacks {key}; missing RS evidence cannot fall back to ordinary DAS. Create current downstream read-support evidence explicitly.")
    chosen = _flag(events[key], key)
    if key == DEFAULT_EVENT_SET:
        if "significant" not in events:
            raise ValueError("RS-DAS requires the separate original significant flag.")
        if (chosen & ~_flag(events["significant"], "significant")).any():
            raise ValueError("Invalid source: RS-DAS is not a subset of ordinary DAS.")
        if "measurement_eligible" not in events or "read_support_policy" not in events:
            raise ValueError("RS-DAS requires measurement_eligible and the recorded every-sample read-support policy; a historical stable flag alone is insufficient.")
        if not events["read_support_policy"].eq(READ_SUPPORT_POLICY_VERSION).all():
            raise ValueError("RS-DAS source uses a missing or historical read-support policy; preserve it and create current evidence explicitly.")
        if (chosen & ~_flag(events["measurement_eligible"], "measurement_eligible")).any():
            raise ValueError("Invalid source: RS-DAS contains events failing every-sample read support.")
        if "criterion_min_reads" in events and not events["criterion_min_reads"].eq(DEFAULT_MIN_READS).all():
            raise ValueError("RS-DAS source does not implement the agreed ten-count threshold.")
    return chosen


def selected_event_view(events: pd.DataFrame, event_set: Any = None) -> pd.DataFrame:
    """Full immutable-source audit view with a separate chosen-membership flag."""
    result = events.copy()
    result["selected_das"] = selected_event_mask(events, event_set)
    result["das_event_set"] = normalize_event_set(event_set)
    return result


def legacy_kernel_event_view(events: pd.DataFrame, event_set: Any = None) -> pd.DataFrame:
    """Private input adapter for frozen kernels. Never publish as raw events."""
    result = selected_event_view(events, event_set)
    result["source_stable_significant"] = events.get("stable_significant", pd.Series(pd.NA, index=events.index))
    result["stable_significant"] = result["selected_das"]
    return result


def policy_matches_artifact(artifact: Any, event_set: Any) -> bool:
    """Full DAS tables are dual-policy sources; derived foregrounds are not."""
    if artifact.analysis_id == "differential.as" and artifact.name == "Differential AS events":
        return True
    if workflow_das_capability(artifact.analysis_id)["mode"] != "switchable":
        return True
    recorded = artifact.metadata.get("das_selection", {})
    return (recorded.get("policy_version") == DAS_SELECTION_POLICY_VERSION
            and recorded.get("event_set") == normalize_event_set(event_set))
