"""Design-preserving biological-unit bootstrap shared by discovery and benchmarks."""
from dataclasses import dataclass
import hashlib
import json

import numpy as np
import pandas as pd


@dataclass
class BootstrapPlan:
    blocks: tuple
    strata: tuple
    columns: tuple
    unit_field: str

    def draw(self, rng):
        chosen = [rng.choice(group, size=len(group), replace=True) for group in self.strata]
        return np.concatenate([self.blocks[int(i)] for group in chosen for i in group])

    def audit(self):
        movable = sum(len(group) for group in self.strata if len(group) > 1)
        return {"resampling_unit_field": self.unit_field,
            "resampling_units": len(self.blocks), "resampling_strata": len(self.strata),
            "resampleable_units": movable, "fixed_singleton_strata": sum(len(g) == 1 for g in self.strata),
            "resampling_scope": "declared biological blocks" if self.unit_field != "sample_id" else "sample IDs assumed independent; donor identity unavailable",
            "resampling_status": ("no_between_unit_resampling" if not movable else
                "partially_resampleable_fixed_strata" if movable < len(self.blocks) else "available")}

    def draw_hash(self, indices):
        return hashlib.sha256(json.dumps([self.columns[int(i)] for i in indices]).encode()).hexdigest()

    def require_resampleable(self):
        if not any(len(group) > 1 for group in self.strata):
            raise ValueError("No biological units can be resampled within matched condition profiles; bootstrap confidence is not estimable")


def bootstrap_plan(samples, columns):
    columns = pd.Index(columns)
    if not len(columns) or not columns.is_unique:
        raise ValueError("Bootstrap requires nonempty, unique sample columns")
    if samples is None:
        samples = pd.DataFrame({"sample_id": columns})
    if "sample_id" not in samples or samples.sample_id.isna().any() or samples.sample_id.duplicated().any():
        raise ValueError("Bootstrap sample IDs must be unique and nonmissing")
    if not set(columns) <= set(samples.sample_id):
        raise ValueError("Bootstrap columns are absent from sample metadata")
    ordered = samples.set_index("sample_id").loc[columns]
    field = next((key for key in ("donor_id", "biological_sample_id") if key in ordered), "sample_id")
    units = pd.Series(columns, index=columns) if field == "sample_id" else ordered[field]
    if units.isna().any() or units.astype(str).str.strip().eq("").any():
        raise ValueError("Declared bootstrap biological IDs must be complete; no guessed pairing")
    conditions = ordered.get("condition", pd.Series("all", index=columns))
    if conditions.isna().any() or conditions.astype(str).str.strip().eq("").any():
        raise ValueError("Bootstrap condition metadata must be complete")
    units, conditions = units.astype(str).to_numpy(), conditions.astype(str).to_numpy()
    blocks = tuple(np.flatnonzero(units == unit) for unit in pd.unique(units))
    strata = {}
    for i, block in enumerate(blocks):
        # Keep every paired measurement together. Matching condition multisets
        # also preserves condition counts when units have different observation patterns.
        profile = tuple(sorted(conditions[block]))
        strata.setdefault(profile, []).append(i)
    return BootstrapPlan(blocks, tuple(np.asarray(g) for g in strata.values()),
                         tuple(map(str, columns)), field)


def bootstrap_indices(samples, columns, rng):
    return bootstrap_plan(samples, columns).draw(rng)
