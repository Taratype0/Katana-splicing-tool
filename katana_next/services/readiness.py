"""Read-only run planning. Never infer experimental factors from filenames."""
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


def experimental_design(manifest):
    samples = [{"sample_id": s.sample_id, "biological_sample_id": s.biological_sample_id,
                "condition": s.condition, "batch": s.batch or "", **{
                    f"factor:{k}": v for k, v in s.attributes.items()}} for s in manifest.samples]
    by_id = {s.sample_id: s for s in manifest.samples}
    contrasts = []
    for c in manifest.contrasts:
        left = {by_id[s].biological_sample_id for s in c.numerator.sample_ids}
        right = {by_id[s].biological_sample_id for s in c.denominator.sample_ids}
        common = left & right
        contrasts.append({"contrast_id": c.contrast_id, "label": c.label,
            "numerator": c.numerator.label, "denominator": c.denominator.label,
            "effect_definition": f"{c.numerator.label} minus {c.denominator.label}",
            "numerator_samples": list(c.numerator.sample_ids), "denominator_samples": list(c.denominator.sample_ids),
            "numerator_biological_units": len(left), "denominator_biological_units": len(right),
            "shared_biological_units": len(common),
            "pairing": "complete" if common and common == left == right else "partial" if common else "unpaired",
            "replicate_warning": min(len(left), len(right)) < 3})
    pairs = []
    for a, b in combinations(manifest.contrasts, 2):
        ids_a = set(a.numerator.sample_ids) | set(a.denominator.sample_ids)
        ids_b = set(b.numerator.sample_ids) | set(b.denominator.sample_ids)
        shared = ids_a & ids_b
        shared_bio = {by_id[s].biological_sample_id for s in ids_a} & {by_id[s].biological_sample_id for s in ids_b}
        pairs.append({"contrast_A": a.contrast_id, "contrast_B": b.contrast_id,
            "label_A": a.label, "label_B": b.label, "shared_samples": sorted(shared),
            "shared_biological_units": sorted(shared_bio),
            "matched_denominator": set(a.denominator.sample_ids) == set(b.denominator.sample_ids),
            "independent_contrasts": not bool(shared_bio),
            "interpretation": "Shared biological units: overlap enrichment is not an independence-adjusted test."
                if shared_bio else "No shared biological units declared; biological independence still requires design review."})
    frame = pd.DataFrame(samples)
    batch = {"included": False, "estimable": None, "reason": "No batch supplied"}
    if len(frame) and frame.batch.ne("").any():
        if frame.batch.eq("").any():
            batch = {"included": False, "estimable": False, "reason": "Partial batch metadata; no implicit imputation"}
        else:
            x = pd.get_dummies(frame[["condition", "batch"]], drop_first=True, dtype=float)
            design = np.column_stack([np.ones(len(frame)), x.to_numpy()])
            rank = int(np.linalg.matrix_rank(design))
            batch = {"included": False, "estimable": rank == design.shape[1] and len(frame) > rank,
                     "design_columns": design.shape[1], "rank": rank, "residual_df": len(frame)-rank,
                     "reason": "Condition + optional batch rank check; batch is not automatically enabled"}
    return {"samples": samples, "contrasts": contrasts, "comparison_pairs": pairs, "batch_design": batch,
            "factor_policy": "Only declared attributes are used; genotype, DP and inhibitor are not guessed from names."}


def preflight(service, project_id, workflow_id, parameters=None):
    context = service.projects.get(project_id)
    workflow = service.registry.get(workflow_id)
    errors, warnings, inputs, parents = [], [], [], []
    normalized = {}
    try:
        normalized = workflow.normalize_parameters(parameters)
        from katana_next.services.selection_routes import validate_route
        validate_route(context, workflow_id, normalized)
        normalized = service._normalize_contrast_scope(context, workflow, normalized)
        workflow.validate_assets(context)
        parents = workflow.parent_artifacts(context, normalized)
        paths = workflow.input_paths(context, normalized)
        for path in dict.fromkeys(Path(p).resolve() for p in paths):
            if not path.is_file():
                errors.append(f"Input missing: {path}")
            else:
                inputs.append({"path": str(path), "bytes": path.stat().st_size})
    except (ValueError, KeyError, OSError, TypeError) as exc:
        errors.append(str(exc))
    if workflow_id.startswith("programs."):
        n = len({s.biological_sample_id for s in context.manifest.samples})
        if n < 15:
            warnings.append(f"Only {n} declared biological units. WGCNA recommends >=15, preferably >=20; this is not a calibrated cutoff for splicing.")
    if workflow_id.startswith("single_cell."):
        warnings.append("Assay coverage and donor metadata determine eligibility. Cell counts do not increase biological replicate counts; absent coverage is not zero usage.")
    design = experimental_design(context.manifest)
    if design["batch_design"]["estimable"] is False:
        warnings.append("Condition + batch is incomplete, confounded or has no residual degrees of freedom.")
    return {"workflow_id": workflow_id, "label": workflow.definition.label,
        "ready": not errors, "errors": errors, "warnings": warnings,
        "inputs": inputs, "parent_artifacts": [{"artifact_id": p.artifact_id, "name": p.name,
            "workflow": p.analysis_id, "run_fingerprint": p.run_fingerprint} for p in parents],
        "parameters": normalized, "input_bytes": sum(p["bytes"] for p in inputs),
        "memory_estimate": "Input bytes are not a peak RAM estimate; dense network memory grows quadratically with events.",
        "scope": "Read-only prerequisite check; no analysis, tool installation or download. It does not validate every input record or executable dependency."}
