"""Execution-only bounds for the complete, unchanged MotifTrace test family.

No sequence scores or significance results are computed here. Counts follow
the existing paired-dataset dimensions, including explicit unavailable windows.
"""
from __future__ import annotations

DEFAULT_MAX_TESTS = 1_000_000
DEFAULT_MAX_MATRIX_CELLS = 8_000_000


def execution_limits(max_tests=DEFAULT_MAX_TESTS, max_matrix_cells=DEFAULT_MAX_MATRIX_CELLS):
    values = {"max_tests": max_tests, "max_matrix_cells": max_matrix_cells}
    for name, value in values.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)) or value < 1 or not float(value).is_integer():
            raise ValueError(f"{name} must be a finite positive integer execution budget")
    return {key: int(value) for key, value in values.items()}


def estimate_paired_family(rna, matching, models, *, widths=(10, 20, 50)):
    """Exact allocated tests/cells, before motif scanning or scoring."""
    from katana_next.science.motiftrace_events import validate_models
    validate_models(models)
    if not widths or any(not isinstance(w, int) or w < 1 for w in widths):
        raise ValueError("Explicit positive integer window widths required")
    pairs = matching[matching.status.eq("matched")]
    if (pairs.foreground_gene.duplicated().any() or pairs.background_gene.duplicated().any() or
        set(pairs.foreground_gene) & set(pairs.background_gene)):
        raise ValueError("Dependent/reused foreground/background genes are not independent pairs")
    regions = rna[rna.layer.isin(["affected_sequence", "regulatory_pre_mRNA"])]
    lookup = {(r.event_id, r.region): len(r.sequence) for r in regions.itertuples()}
    tests = cells = datasets = 0
    for (kind, _direction), group in pairs.groupby(["event_type", "direction"], sort=True):
        for region in sorted(set(regions.loc[regions.event_type.eq(kind), "region"])):
            lengths = [min(lookup.get((r.foreground_event, region), 0),
                           lookup.get((r.background_event, region), 0)) for r in group.itertuples()]
            for model in models.values():
                methods = 1 if model["kind"] == "PWM" else 2
                for width in (0, *sorted(set(widths))):
                    span = 1 if not width else max(1, max((max(0, n-width+1) for n in lengths), default=0))
                    tests += methods * span
                    cells += methods * len(group) * span
                    datasets += methods
    return {"models": len(models), "matched_gene_pairs": len(pairs), "datasets": datasets,
            "scalar_and_window_tests": tests, "paired_matrix_cells": cells,
            "paired_matrix_float64_bytes": cells * 8,
            "memory_note": "Matrix storage only; Python tables, permutation arrays and reports require additional memory.",
            "statistical_family_split": False, "model_or_window_truncation": False}


def enforce_paired_budget(estimate, limits, *, permutations):
    for count, limit in (("scalar_and_window_tests", "max_tests"), ("paired_matrix_cells", "max_matrix_cells")):
        if estimate[count] > limits[limit]:
            raise ValueError(f"MotifTrace pre-scoring execution budget exceeded: {count}={estimate[count]} > "
                             f"{limit}={limits[limit]}. Configure a sufficient execution budget; "
                             "no models, windows, events or statistical families were truncated or scored.")
    if permutations * estimate["matched_gene_pairs"] > 8_000_000:
        raise ValueError("MotifTrace pre-scoring whole-pair permutation memory budget exceeded")
