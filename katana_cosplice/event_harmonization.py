from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass
class EventAlignmentResult:
    mapping: pd.DataFrame
    shared_events: list[str]
    reference_to_target: dict[str, str]


def splice_site_signature(event_key: str) -> str:
    fields = str(event_key).split("|")
    if len(fields) < 5:
        return ""
    return "|".join(fields[1:])


def exact_event_alignment(reference_events: pd.Index, target_events: pd.Index, cohort: str) -> EventAlignmentResult:
    """Align canonical rMATS event keys without coordinate tolerance or gene-only collapsing."""
    reference = pd.Index(reference_events.astype(str)).drop_duplicates()
    target = pd.Index(target_events.astype(str)).drop_duplicates()
    shared = reference.intersection(target).tolist()
    union = reference.union(target)
    mapping = pd.DataFrame(
        {
            "cohort": cohort,
            "event_key": union,
            "in_reference": union.isin(reference),
            "in_target": union.isin(target),
            "exact_match": union.isin(shared),
            "matching_rule": "exact canonical event_key (gene_id + event_type + strand + rMATS coordinates)",
        }
    )
    return EventAlignmentResult(mapping=mapping, shared_events=shared, reference_to_target={key: key for key in shared})


def annotation_aware_event_alignment(
    reference_events: pd.Index,
    target_events: pd.Index,
    cohort: str,
    *,
    reference_genome_build: str = "",
    target_genome_build: str = "",
    mapping_table: pd.DataFrame | None = None,
) -> EventAlignmentResult:
    reference = pd.Index(reference_events.astype(str)).drop_duplicates()
    target = pd.Index(target_events.astype(str)).drop_duplicates()
    pairs: dict[str, str] = {key: key for key in reference.intersection(target)}
    methods: dict[str, str] = {key: "exact canonical event_key" for key in pairs}
    if mapping_table is not None and not mapping_table.empty:
        required = {"reference_event_key", "target_event_key"}
        if not required.issubset(mapping_table):
            raise ValueError("Event mapping table requires reference_event_key and target_event_key.")
        for _, row in mapping_table.iterrows():
            left, right = str(row["reference_event_key"]), str(row["target_event_key"])
            if left in reference and right in target:
                pairs[left] = right
                methods[left] = "explicit validated event/liftover mapping"
    same_build = bool(reference_genome_build and target_genome_build and reference_genome_build == target_genome_build)
    if same_build or (not reference_genome_build and not target_genome_build):
        target_by_signature: dict[str, list[str]] = {}
        for event in target:
            target_by_signature.setdefault(splice_site_signature(event), []).append(event)
        for event in reference:
            if event in pairs:
                continue
            candidates = target_by_signature.get(splice_site_signature(event), [])
            if len(candidates) == 1:
                pairs[event] = candidates[0]
                methods[event] = "same-build exact splice-site signature; annotation label may differ"
    rows = []
    for event in reference:
        target_event = pairs.get(event, "")
        matching_rule = methods.get(event, "not matched")
        rows.append({
            "cohort": cohort,
            "event_key": event,
            "in_reference": True,
            "in_target": bool(target_event),
            "exact_match": matching_rule == "exact canonical event_key",
            "reference_event_key": event,
            "target_event_key": target_event,
            "matched": bool(target_event),
            "matching_rule": matching_rule,
            "reference_genome_build": reference_genome_build,
            "target_genome_build": target_genome_build,
        })
    return EventAlignmentResult(
        mapping=pd.DataFrame(rows), shared_events=list(pairs), reference_to_target=pairs,
    )
