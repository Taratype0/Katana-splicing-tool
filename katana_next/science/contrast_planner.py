"""Cell-mean contrasts and row-space estimability, without fitting outcomes."""
from itertools import combinations, product
import json
import math

import numpy as np
import pandas as pd


def plan_contrasts(manifest, factors=None, *, include_batch=False, pairing="auto", reference_levels=None):
    factors = list(factors or ["condition"])
    if len(set(factors)) != len(factors) or not 1 <= len(factors) <= 3:
        raise ValueError("Choose one to three distinct categorical factors explicitly")
    if pairing not in {"auto", "paired", "unpaired"}:
        raise ValueError("Pairing must be auto, paired or unpaired")
    rows = []
    for sample in manifest.samples:
        values = tuple(sample.condition if f == "condition" else sample.attributes.get(f, "") for f in factors)
        if any(not v.strip() for v in values):
            raise ValueError(f"Incomplete declared factors for {sample.sample_id}; no label guessing or imputation")
        rows.append({"sample_id": sample.sample_id, "unit": sample.biological_sample_id,
                     "batch": sample.batch or "", "cell": values})
    if not rows:
        raise ValueError("Declare samples before generating contrasts")
    frame = pd.DataFrame(rows)
    levels = [sorted({v[i] for v in frame.cell}) for i in range(len(factors))]
    reference_levels = {} if reference_levels is None else reference_levels
    if not isinstance(reference_levels, dict) or set(reference_levels)-set(factors):
        raise ValueError("Reference levels must name selected factors only")
    for i,factor in enumerate(factors):
        if factor in reference_levels:
            reference = reference_levels[factor]
            if reference not in levels[i]:
                raise ValueError(f"Unknown reference level for {factor}: {reference}")
            levels[i] = [reference, *[level for level in levels[i] if level != reference]]
    if math.prod(map(len, levels)) > 128:
        raise ValueError("More than 128 factorial cells; explicitly narrow the factors rather than truncating contrasts")
    cells = list(product(*levels))
    names = [json.dumps(dict(zip(factors, c)), sort_keys=True) for c in cells]
    if include_batch and (frame.batch.eq("").any() or frame.groupby(["unit", "cell"]).batch.nunique().gt(1).any()):
        raise ValueError("Batch must be complete and unique within biological-unit/factor-cell; aggregate technical partitions explicitly")
    observations = frame.drop_duplicates(["unit", "cell"]).reset_index(drop=True)
    repeated = observations.groupby("unit").cell.nunique().gt(1).any()
    paired = pairing == "paired" or (pairing == "auto" and bool(repeated))
    if pairing == "paired" and not repeated:
        raise ValueError("Paired model requested but no biological unit is observed in more than one cell")
    if pairing == "unpaired" and repeated:
        raise ValueError("Repeated biological units cannot be treated as independent; use auto or paired")
    x = np.column_stack([observations.cell.map(lambda v: v == c).astype(float) for c in cells])
    columns = list(names)
    for key, enabled in (("unit", paired), ("batch", include_batch)):
        if enabled:
            dummy = pd.get_dummies(observations[key], drop_first=True, dtype=float)
            x = np.column_stack([x, dummy.to_numpy()])
            columns.extend(f"{key}:{v}" for v in dummy.columns)
    rank = int(np.linalg.matrix_rank(x))
    residual_df = len(observations)-rank
    projection = np.linalg.pinv(x) @ x
    results = []

    def emit(kind, label, weights):
        c = np.zeros(x.shape[1])
        for cell, weight in weights.items():
            c[cells.index(cell)] = weight
        estimable = bool(np.allclose(c, c @ projection, atol=1e-9, rtol=1e-9))
        samples_by_cell = {names[cells.index(cell)]: sorted(frame.loc[frame.cell.map(lambda v: v == cell), "sample_id"])
                           for cell in weights}
        units = {cell: set(frame.loc[frame.cell.map(lambda v: v == cell), "unit"]) for cell in weights}
        minimum = min(map(len, units.values()))
        complete_units = len(set.intersection(*units.values()))
        matched = []
        if kind == "simple_effect":
            positive = set(next(v for cell, v in samples_by_cell.items() if weights[cells[names.index(cell)]] > 0))
            negative = set(next(v for cell, v in samples_by_cell.items() if weights[cells[names.index(cell)]] < 0))
            for contrast in manifest.contrasts:
                a, b = set(contrast.numerator.sample_ids), set(contrast.denominator.sample_ids)
                if (a, b) == (positive, negative):
                    matched.append({"contrast_id": contrast.contrast_id, "direction_multiplier": 1})
                elif (a, b) == (negative, positive):
                    matched.append({"contrast_id": contrast.contrast_id, "direction_multiplier": -1})
        results.append({"kind": kind, "label": label, "estimable": estimable,
                        "residual_df": residual_df, "minimum_cell_biological_units": minimum,
                        "complete_paired_units": complete_units if paired else None,
                        "status": "not_estimable_missing_or_confounded_cells" if not estimable else
                                  "no_residual_degrees_of_freedom" if residual_df <= 0 else
                                  "limited_replication" if minimum < 2 else "design_estimable",
                        "coefficients": {names[cells.index(cell)]: w for cell, w in weights.items()},
                        "sample_arms": samples_by_cell, "existing_imported_contrasts": matched})

    for i, factor in enumerate(factors):
        others = [j for j in range(len(factors)) if j != i]
        for baseline, target in combinations(levels[i], 2):
            marginal = {}
            contexts = list(product(*(levels[j] for j in others)))
            for context in contexts:
                a = [None]*len(factors)
                for j, v in zip(others, context):
                    a[j] = v
                a[i] = target; positive = tuple(a)
                a[i] = baseline; negative = tuple(a)
                label = f"{factor}: {target} minus {baseline}"
                if others:
                    label += " | " + ", ".join(f"{factors[j]}={v}" for j, v in zip(others, context))
                emit("simple_effect", label, {positive: 1., negative: -1.})
                marginal[positive], marginal[negative] = 1./len(contexts), -1./len(contexts)
            if others:
                emit("marginal_effect", f"{factor}: {target} minus {baseline} | equal-weight average over other factor cells", marginal)
    for i, j in combinations(range(len(factors)), 2):
        remaining = [k for k in range(len(factors)) if k not in (i, j)]
        for (a0, a1), (b0, b1), rest in product(combinations(levels[i], 2), combinations(levels[j], 2), product(*(levels[k] for k in remaining))):
            weights = {}
            for a, b, w in ((a1,b1,1.), (a0,b1,-1.), (a1,b0,-1.), (a0,b0,1.)):
                cell = [None]*len(factors); cell[i], cell[j] = a,b
                for k,v in zip(remaining,rest):
                    cell[k] = v
                weights[tuple(cell)] = w
            label = f"({factors[i]}: {a1}-{a0}) at {factors[j]}={b1} minus at {factors[j]}={b0}"
            if remaining:
                label += " | " + ", ".join(f"{factors[k]}={v}" for k,v in zip(remaining,rest))
            emit("interaction", label, weights)
    return {"factors": factors, "levels": dict(zip(factors,levels)), "paired_block": paired,
            "include_batch": include_batch, "observations": len(observations), "input_samples": len(frame),
            "technical_partitions_collapsed": len(frame)-len(observations), "rank": rank, "columns": columns,
            "residual_df": residual_df, "design_matrix": x.tolist(), "contrasts": results,
            "reference_levels": reference_levels,
            "scope": "Design estimability only, before event-specific missingness. Declared references come first; otherwise alphabetical level order is shown, not an inferred biological baseline. No statistical outputs or imported contrasts are changed. Donor blocking is fixed-effect, not a mixed model."}
