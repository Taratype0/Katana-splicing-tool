"""Shared study design policy for native Nextflow, CLI and browser execution."""
from __future__ import annotations

RECOMMENDED_UNITS = 12


def validate_study_design(samples, *, paired=False):
    """Count biological units, without imposing a minimum study-size gate.

    Pairing is explicitly enabled; repeated observations never increase the
    independent unit count. Individual statistical models still report NA when
    their design is not estimable.
    """
    if not isinstance(paired, bool):
        raise ValueError("paired must be true or false")
    required = ["sample_id", "biological_sample_id", "condition"]
    if samples.empty or not set(required).issubset(samples.columns):
        raise ValueError("Provide sample_id, biological_sample_id and condition")
    values = samples[required]
    if values.isna().any().any() or values.astype(str).apply(lambda x: x.str.strip().eq("").any()).any():
        raise ValueError("Sample, biological-unit and condition identities cannot be empty")
    if not samples.sample_id.is_unique:
        raise ValueError("Each sample_id must occur exactly once")
    if samples.duplicated(["biological_sample_id", "condition"]).any():
        raise ValueError("Merge technical replicates within biological unit and condition before network analysis")
    repeated = samples.biological_sample_id.duplicated().any()
    if repeated and not paired:
        raise ValueError("Repeated biological units require explicit paired mode; observations cannot be treated as independent")
    if paired and not repeated:
        raise ValueError("Paired mode requires shared biological_sample_id values across conditions")
    units = int(samples.biological_sample_id.nunique())
    return {"paired": paired, "samples": len(samples), "biological_units": units,
            "recommended_biological_units": RECOMMENDED_UNITS,
            "below_recommended_units": units < RECOMMENDED_UNITS,
            "minimum_units_gate": False,
            "sample_size_note": "At least 12 independent biological units are recommended (project default, not a significance or eligibility threshold). Repeated measurements count once."}


def project_study_design(manifest, *, paired=False, association_contrast_ids=None):
    """Translate explicit per-comparison pairing, preserving numerator direction."""
    import pandas as pd
    samples = pd.DataFrame([{"sample_id": s.sample_id,
                             "biological_sample_id": s.biological_sample_id,
                             "condition": s.condition} for s in manifest.samples])
    designs = manifest.metadata.get("comparison_designs", {})
    enabled = paired or any(d.get("paired") is True for d in designs.values())
    policy = validate_study_design(samples, paired=enabled)
    rows = []
    selected = None
    if association_contrast_ids is not None:
        if not isinstance(association_contrast_ids, (list, tuple)) or any(not isinstance(c, str) for c in association_contrast_ids):
            raise ValueError("association_contrast_ids must be a list of comparison IDs")
        selected = set(association_contrast_ids)
        if len(selected) != len(association_contrast_ids) or not selected <= {c.contrast_id for c in manifest.contrasts}:
            raise ValueError("Association comparisons must be unique declared comparison IDs")
    by_id = samples.set_index("sample_id")
    for c in manifest.contrasts:
        if selected is not None and c.contrast_id not in selected:
            continue
        arms = [by_id.loc[list(arm.sample_ids)] for arm in (c.numerator, c.denominator)]
        if any(arm.condition.nunique() != 1 for arm in arms):
            raise ValueError("Network association requires one declared condition in each contrast arm")
        a, b = [str(arm.condition.iloc[0]) for arm in arms]
        if a == b:
            raise ValueError("Network contrast arms must have different condition labels")
        declared = designs.get(c.contrast_id, {}).get("paired", paired)
        if not isinstance(declared, bool):
            raise ValueError("Comparison paired mode must be true or false")
        units_a, units_b = [set(arm.biological_sample_id) for arm in arms]
        if declared:
            if units_a != units_b or any(not arm.biological_sample_id.is_unique for arm in arms):
                raise ValueError(f"{c.contrast_id}: paired analysis needs one observation per unit in each arm")
        elif units_a & units_b:
            raise ValueError(f"{c.contrast_id}: repeated units require paired=true for this comparison")
        row = {"contrast_id": c.contrast_id, "factor": "condition", "a": a, "b": b, "paired": declared}
        selected_samples = list(c.numerator.sample_ids) + list(c.denominator.sample_ids)
        if set(selected_samples) != set(samples.loc[samples.condition.isin([a, b]), "sample_id"]):
            row["sample_ids"] = selected_samples
        rows.append(row)
    return policy, rows
