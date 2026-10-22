"""Diagnostics for a fixed partition, never a new module-selection objective."""
from itertools import combinations
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.stats import beta as beta_distribution

from katana_cosplice.uncertainty import estimate_beta_binomial_uncertainty
from katana_cosplice.resampling import bootstrap_plan
from katana_cosplice.partition_metrics import UNASSIGNED


def calibration_grid(*, events=200, seed=1701, replicates=3):
    """Known latent *read-fraction* truth, not length-normalized rMATS PSI.

    Coverage targets per-sample latent proportions. This tests the existing
    effective-count approximation, not a new fully fitted beta-binomial model.
    """
    if events < 10 or replicates < 2:
        raise ValueError("Use >=10 independent simulated events and >=2 replicates")
    rng = np.random.default_rng(seed)
    metadata = pd.DataFrame({"sample_id": [f"s{i}" for i in range(2*replicates)],
                             "condition": ["A"]*replicates + ["B"]*replicates})
    records = []
    for depth, rho in ((n, r) for n in (5, 20, 100, 500) for r in (0., .05, .2)):
        mean = rng.uniform(.05, .95, (events, 2))
        mean = np.repeat(mean, replicates, axis=1)
        truth = mean if rho == 0 else rng.beta(mean*(1/rho-1), (1-mean)*(1/rho-1))
        ijc = pd.DataFrame(rng.binomial(depth, truth), columns=metadata.sample_id)
        sjc = depth-ijc
        fit = estimate_beta_binomial_uncertainty(ijc, sjc, metadata, include_sample_metrics=False)
        for method, a, b in (("Jeffreys binomial", ijc.to_numpy()+.5, sjc.to_numpy()+.5),
                             ("CoSplice effective-count approximation", fit.posterior_alpha.to_numpy(), fit.posterior_beta.to_numpy())):
            lower, upper = beta_distribution.ppf(.025, a, b), beta_distribution.ppf(.975, a, b)
            covered = (lower <= truth) & (truth <= upper)
            per_event = covered.mean(axis=1)
            records.append({"method": method, "depth": depth, "true_rho": rho,
                "replicates_per_condition": replicates, "simulated_events": events,
                "nominal_coverage": .95, "observed_coverage": float(covered.mean()),
                "coverage_mc_se": float(per_event.std(ddof=1)/np.sqrt(events)),
                "mean_interval_width": float((upper-lower).mean()),
                "mean_absolute_error": float(np.abs(a/(a+b)-truth).mean()),
                "mean_reliability": float(fit.event_metrics.network_reliability.mean()),
                "target": "per-sample latent read fraction; not differential significance or normalized PSI"})
    return pd.DataFrame(records)


def thinning_diagnostic(ijc, sjc, samples, *, fractions=(.25, .5), seed=1701, iterations=10):
    """Conditional split-read reproducibility, not independent biological truth."""
    if not ijc.index.equals(sjc.index) or not ijc.columns.equals(sjc.columns):
        raise ValueError("Count axes must match exactly")
    i, s = ijc.to_numpy(float), sjc.to_numpy(float)
    if iterations < 1:
        raise ValueError("At least one thinning iteration is required")
    if np.any(np.isinf(i) | np.isinf(s) | (i < 0) | (s < 0)):
        raise ValueError("Raw counts must be nonnegative and finite or missing")
    valid = np.isfinite(i) & np.isfinite(s) & (i >= 0) & (s >= 0)
    if np.any(valid & ((i != np.floor(i)) | (s != np.floor(s)))):
        raise ValueError("Read thinning requires integer raw counts, not effective counts")
    rng = np.random.default_rng(seed)
    total = np.where(valid, i+s, 0).astype(np.int64)
    i = np.where(valid, i, 0).astype(np.int64)
    records = []
    for fraction in fractions:
        if not 0 < fraction < 1:
            raise ValueError("Read fraction must be strictly between zero and one")
        for iteration in range(iterations):
            n = np.floor(total*fraction).astype(np.int64)
            chosen = rng.hypergeometric(i, total-i, n)
            held = total-n
            eligible = valid & (n > 0) & (held > 0)
            train = np.divide(chosen, n, out=np.full(n.shape, np.nan), where=n>0)
            test = np.divide(i-chosen, held, out=np.full(n.shape, np.nan), where=held>0)
            fit = estimate_beta_binomial_uncertainty(pd.DataFrame(chosen, index=ijc.index, columns=ijc.columns),
                pd.DataFrame(n-chosen, index=ijc.index, columns=ijc.columns), samples, include_sample_metrics=False)
            estimate = fit.posterior_psi.to_numpy()
            records.append({"fraction": fraction, "iteration": iteration,
                "total_event_sample_cells": int(valid.sum()), "evaluable_split_cells": int(eligible.sum()),
                "evaluable_fraction": float(eligible.sum()/max(1, valid.sum())),
                "raw_split_mae": float(np.abs(train-test)[eligible].mean()) if eligible.any() else np.nan,
                "approximation_split_mae": float(np.abs(estimate-test)[eligible].mean()) if eligible.any() else np.nan,
                "reference": "disjoint remaining reads; noisy conditional reference, not independent validation"})
    return pd.DataFrame(records)


def event_junctions(row):
    chrom, strand = str(row.get("chr", "")), str(row.get("strand", ""))
    if chrom.lower() in {"", "nan", "none", "<na>"} or strand not in {"+", "-"}:
        return set()
    def site(a, b):
        start, end = float(row[a]), float(row[b])
        if not np.isfinite([start, end]).all() or start < 0 or start != int(start) or end != int(end) or start >= end:
            raise ValueError("Invalid junction interval")
        return (chrom, int(start), int(end), strand)
    try:
        kind = row.get("event_type")
        if kind == "SE":
            return {site("upstreamEE", "exonStart_0base"), site("exonEnd", "downstreamES"), site("upstreamEE", "downstreamES")}
        if kind == "MXE":
            return {site("upstreamEE", f"{exon}ExonStart_0base") for exon in ("1st", "2nd")} | {
                site(f"{exon}ExonEnd", "downstreamES") for exon in ("1st", "2nd")}
        if kind == "RI":
            return {site("upstreamEE", "downstreamES")}
        if kind in {"A3SS", "A5SS"}:
            if int(row["flankingES"]) < int(row["longExonStart_0base"]):
                return {site("flankingEE", "longExonStart_0base"), site("flankingEE", "shortES")}
            return {site("longExonEnd", "flankingES"), site("shortEE", "flankingES")}
    except (KeyError, ValueError, TypeError, OverflowError):
        return set()
    return set()


def _pairs(n, maximum, rng):
    count = n*(n-1)//2
    if count <= maximum:
        return list(combinations(range(n), 2))
    found = set()
    while len(found) < maximum:
        for a, b in rng.integers(0, n, size=(min(4096, maximum-len(found)+100), 2)):
            if a != b:
                found.add((int(min(a, b)), int(max(a, b))))
            if len(found) == maximum:
                break
    return sorted(found)


def _corr(a, b):
    valid = np.isfinite(a) & np.isfinite(b)
    if valid.sum() < 4 or np.std(a[valid]) < 1e-10 or np.std(b[valid]) < 1e-10:
        return np.nan
    return float(np.corrcoef(a[valid], b[valid])[0, 1])


def _partial_corr(a, b, design):
    """Refit both events on exactly the same measured samples."""
    shared = np.isfinite(a) & np.isfinite(b)
    x = design[shared]
    rank = int(np.linalg.matrix_rank(x)) if shared.any() else 0
    df = int(shared.sum()) - rank
    if rank < design.shape[1]:
        return np.nan, df, "incomplete_pair_design"
    if df < 3:
        return np.nan, df, "insufficient_residual_df"
    y = np.column_stack([a[shared], b[shared]])
    residual = y - x @ np.linalg.lstsq(x, y, rcond=None)[0]
    corr = _corr(residual[:, 0], residual[:, 1])
    return corr, df, "available" if np.isfinite(corr) else "constant_residual"


def _denominator_key(row):
    chrom, strand, cluster = row.get("chr"), row.get("strand"), row.get("cluster_id")
    if (pd.isna(chrom) or not str(chrom).strip() or strand not in {"+", "-"}
            or pd.isna(cluster) or not str(cluster).strip()):
        return None
    return (str(chrom), str(strand), str(cluster))


def structural_components(event_ids, structural, records):
    """All-event inverted index; connected structure does not imply independence."""
    parents = list(range(len(event_ids)))
    def root(i):
        while parents[i] != i:
            parents[i] = parents[parents[i]]
            i = parents[i]
        return i
    owners = {}
    assessable = []
    for i, (junctions, record) in enumerate(zip(structural, records)):
        denominator = _denominator_key(record)
        keys = [("junction", *j) for j in junctions]
        if denominator is not None:
            keys.append(("denominator", *denominator))
        assessable.append(bool(keys))
        for key in keys:
            if key in owners:
                parents[root(i)] = root(owners[key])
            else:
                owners[key] = i
    groups = {}
    for i in range(len(event_ids)):
        groups.setdefault(root(i), []).append(i)
    rows = []
    for positions in groups.values():
        members = sorted(str(event_ids[i]) for i in positions)
        identity = "structure_" + hashlib.sha256(json.dumps(members).encode()).hexdigest()[:16]
        for i in positions:
            rows.append({"event_key": event_ids[i], "structural_component": identity if assessable[i] else None,
                "structural_component_events": len(positions) if assessable[i] else 0,
                "structural_status": "linked" if len(positions) > 1 else "isolated_known" if assessable[i] else "not_assessable"})
    return pd.DataFrame(rows, columns=["event_key", "structural_component", "structural_component_events", "structural_status"])


def fixed_module_audit(psi, samples, events, assignments, *, max_pairs=5000, seed=1701, include_batch=False):
    if max_pairs < 1:
        raise ValueError("max_pairs must be positive")
    if not psi.index.is_unique or not psi.columns.is_unique:
        raise ValueError("PSI axes must be unique")
    if not {"event_key", "module"} <= set(assignments):
        raise ValueError("Assignments require event_key and module")
    for frame, key in ((samples, "sample_id"), (events, "event_key"), (assignments, "event_key")):
        if key not in frame or frame[key].duplicated().any():
            raise ValueError(f"Missing or ambiguous {key} mapping")
    if not set(psi.columns) <= set(samples.sample_id):
        raise ValueError("PSI samples are absent from metadata")
    if not set(psi.index) <= set(events.event_key):
        raise ValueError("PSI events are absent from structural metadata")
    if not set(assignments.event_key) <= set(psi.index):
        raise ValueError("Module assignments contain events outside the supplied PSI universe")
    ordered = samples.set_index("sample_id").loc[psi.columns]
    factors = ["condition"] + (["batch"] if include_batch else [])
    if any(c not in ordered or ordered[c].isna().any() or ordered[c].astype(str).eq("").any() for c in factors):
        raise ValueError("Explicit complete condition/batch metadata is required")
    x = np.column_stack([np.ones(len(ordered)), pd.get_dummies(ordered[factors].astype(str), drop_first=True, dtype=float)])
    rank = int(np.linalg.matrix_rank(x))
    if rank < x.shape[1]:
        raise ValueError("Condition/batch design is confounded; residual correlations are not identifiable")
    values = psi.to_numpy(float)
    if np.isinf(values).any() or np.any(np.isfinite(values) & ((values < 0) | (values > 1))):
        raise ValueError("PSI must lie in [0, 1]; missing values must remain missing")
    r2 = np.full(len(values), np.nan)
    for k, y in enumerate(values):
        valid = np.isfinite(y)
        local_rank = np.linalg.matrix_rank(x[valid]) if valid.any() else 0
        if valid.sum()-local_rank < 3 or local_rank < x.shape[1]:
            continue
        fit = x[valid] @ np.linalg.lstsq(x[valid], y[valid], rcond=None)[0]
        total = np.sum((y[valid]-y[valid].mean())**2)
        if total > 1e-15:
            r2[k] = 1-np.sum((y[valid]-fit)**2)/total
    annotations = events.set_index("event_key").reindex(psi.index).assign(_audit_event_key=psi.index)
    structural = [event_junctions(row) for row in annotations.to_dict("records")]
    annotation_records = annotations.to_dict("records")
    components = structural_components(psi.index, structural, annotation_records).set_index("event_key").reindex(psi.index).rename_axis("event_key")
    labels = assignments.set_index("event_key").module.reindex(psi.index).fillna("grey").astype(str)
    labels = labels.mask(labels.isin(UNASSIGNED), "grey")
    rng = np.random.default_rng(seed)
    pairs, modules = [], []
    for module in sorted(set(labels)-{"grey", "gray", "0", "unassigned"}):
        positions = np.flatnonzero(labels.to_numpy() == module)
        chosen = _pairs(len(positions), max_pairs, rng)
        rows = []
        for a, b in chosen:
            a, b = positions[a], positions[b]
            ma, mb = annotation_records[a], annotation_records[b]
            cluster_a, cluster_b = _denominator_key(ma), _denominator_key(mb)
            same_cluster = cluster_a is not None and cluster_a == cluster_b
            known = bool(structural[a] and structural[b])
            common = structural[a] & structural[b]
            identical = known and structural[a] == structural[b] and ma.get("event_type") == mb.get("event_type")
            ca, cb = ma.get("event_coordinates"), mb.get("event_coordinates")
            # A junction set alone does not encode inclusion/skipping roles.
            identical = bool(identical and (ma.get("event_type") == "SE"
                or (pd.notna(ca) and pd.notna(cb) and ca == cb)))
            relation = "duplicate_structure" if identical else "shared_denominator_cluster" if same_cluster else "shared_junction" if common else "no_known_shared_junction" if known else "not_assessable"
            shared = np.isfinite(values[a]) & np.isfinite(values[b])
            complementary = bool(shared.sum() >= 4 and np.std(values[a, shared]) > 1e-10 and np.max(np.abs(values[a, shared]+values[b, shared]-1)) < 1e-8)
            partial_r, partial_df, partial_status = _partial_corr(values[a], values[b], x)
            row = {"module": module, "event_A": str(psi.index[a]), "event_B": str(psi.index[b]),
                   "relation": relation, "shared_junction_count": len(common),
                   "numerically_complementary": complementary,
                   "observed_r": _corr(values[a], values[b]), "condition_adjusted_r": partial_r,
                   "condition_adjusted_status": partial_status, "pair_residual_df": partial_df,
                   "shared_measured_samples": int(shared.sum()),
                   "shared_residual_samples": int(shared.sum()) if partial_status == "available" else 0,
                   "complete_design_residual_df": len(ordered)-rank}
            rows.append(row)
        pairs.extend(rows)
        genes = sorted({str(value) for i in positions
                        for value in [next((annotation_records[i].get(c) for c in ("gene_symbol", "gene_id")
                            if pd.notna(annotation_records[i].get(c)) and str(annotation_records[i].get(c)).strip()), "")]
                        if value})
        same = sum(r["relation"] in {"duplicate_structure", "shared_denominator_cluster", "shared_junction"} for r in rows)
        unrelated = [r for r in rows if r["relation"] == "no_known_shared_junction" and not r["numerically_complementary"]]
        local_components = components.iloc[positions]
        modules.append({"module": module, "n_events": len(positions), "genes": ";".join(genes), "n_genes": len(genes),
            "all_possible_pairs": len(positions)*(len(positions)-1)//2, "evaluated_pairs": len(rows),
            "known_structural_dependency_pairs": same,
            "unassessable_structure_pairs": sum(r["relation"] == "not_assessable" for r in rows),
            "all_event_structurally_linked_events": int(local_components.structural_status.eq("linked").sum()),
            "all_event_known_structure_components": local_components.structural_component.nunique(),
            "all_event_structure_unknown_events": int(local_components.structural_status.eq("not_assessable").sum()),
            "pairs_without_known_dependency": len(unrelated),
            "adjusted_evaluable_pairs": sum(np.isfinite(r["condition_adjusted_r"]) for r in rows),
            "adjusted_evaluable_pairs_without_known_dependency": sum(np.isfinite(r["condition_adjusted_r"]) for r in unrelated),
            "mean_condition_r_squared": float(np.nanmean(r2[positions])) if np.isfinite(r2[positions]).any() else np.nan,
            "mean_observed_r": _finite_mean([r["observed_r"] for r in rows]),
            "mean_observed_r_on_adjusted_pairs": _finite_mean([r["observed_r"] for r in rows if np.isfinite(r["condition_adjusted_r"])]),
            "mean_condition_adjusted_r": _finite_mean([r["condition_adjusted_r"] for r in rows]),
            "mean_adjusted_r_without_known_dependency": _finite_mean([r["condition_adjusted_r"] for r in unrelated]),
            "interpretation": "Fixed-module descriptive audit; not a rewiring test, module refit or causal claim"})
    event_table = pd.DataFrame({"event_key": psi.index, "module": labels.to_numpy(),
        "condition_r_squared": r2, "measured_samples": np.isfinite(values).sum(axis=1),
        "structural_junctions": [json.dumps(sorted(j)) for j in structural]})
    event_table = event_table.merge(components.reset_index(), on="event_key", validate="one_to_one")
    resampling = bootstrap_plan(samples, psi.columns).audit()
    coverage = pd.DataFrame([{"input_events": len(psi), "assigned_events": int((~labels.isin(["grey", "gray", "0", "unassigned"])).sum()),
        "retention_fraction": float((~labels.isin(["grey", "gray", "0", "unassigned"])).mean()),
        "independent_units_declared": resampling["resampling_units"], **resampling,
        "condition_count": ordered.condition.nunique(), "design_rank": rank, "residual_df": len(ordered)-rank,
        "max_pairs_per_module": max_pairs, "seed": seed,
        "limitations": "All events checked for connected structure; component count is not an effective independent-event count. Correlation pairs sampled; partial correlation uses pair-complete refits. No automatic module deletion or preservation claim"}])
    return {"modules": pd.DataFrame(modules), "pairs": pd.DataFrame(pairs), "events": event_table, "coverage": coverage}


def _finite_mean(values):
    a = np.asarray(values, float)
    return float(a[np.isfinite(a)].mean()) if np.isfinite(a).any() else np.nan
