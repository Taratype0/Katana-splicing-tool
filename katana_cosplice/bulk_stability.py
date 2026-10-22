"""Label-invariant diagnostics for frozen bulk module-refit experiments.

These are audit utilities, not a module-selection objective. No phenotype or
functional annotation is used to choose members.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


def checked_labels(labels):
    values = np.asarray(labels, dtype=str)
    if values.ndim != 1 or len(values) == 0:
        raise ValueError("Expected a nonempty one-dimensional partition")
    return values


def align_and_classify(reference, candidate):
    """Maximize overlap without matching grey; unmatched splits stay distinct."""
    a, b = checked_labels(reference), checked_labels(candidate)
    if len(a) != len(b):
        raise ValueError("Partitions must have identical event order and length")
    am, bm = sorted(set(a) - {"grey"}), sorted(set(b) - {"grey"})
    mapping = {}
    if am and bm:
        overlaps = np.array([[np.sum((a == x) & (b == y)) for y in bm] for x in am])
        ii, jj = linear_sum_assignment(-overlaps)
        mapping = {bm[j]: am[i] for i, j in zip(ii, jj) if overlaps[i, j] > 0}
    aligned = np.array([mapping.get(x, f"unmatched:{x}") if x != "grey" else "grey" for x in b])
    status = np.full(len(a), "switched", dtype=object)
    status[(a == "grey") & (b == "grey")] = "shared_grey"
    status[(a == "grey") & (b != "grey")] = "newly_assigned"
    status[(a != "grey") & (b == "grey")] = "dropped"
    status[(a != "grey") & (b != "grey") & (a == aligned)] = "retained_same"
    return aligned, status


def module_topology(reference, candidate):
    """Partition each reference module's event pairs into retained/lost/split.

    Pair fractions avoid arbitrary split-size thresholds and label matching.
    Merge pairs are reported relative to candidate modules that contain members
    of this reference module; they must not be summed across reference modules.
    """
    a, b = checked_labels(reference), checked_labels(candidate)
    if len(a) != len(b):
        raise ValueError("Partitions must have identical length")
    choose2 = lambda n: int(n) * (int(n) - 1) // 2
    rows = []
    for module in sorted(set(a) - {"grey"}):
        mask = a == module
        n, retained = int(mask.sum()), int(np.sum(mask & (b != "grey")))
        fragments = pd.Series(b[mask & (b != "grey")]).value_counts()
        total_pairs = choose2(n)
        preserved_pairs = sum(choose2(k) for k in fragments)
        assigned_pairs = choose2(retained)
        rows.append({
            "reference_module": module, "reference_events": n,
            "retained_events": retained, "dropped_events": n - retained,
            "candidate_fragments": len(fragments),
            "largest_fragment_fraction": float(fragments.max() / n) if len(fragments) else 0.,
            "pair_preserved_fraction": preserved_pairs / total_pairs if total_pairs else np.nan,
            "pair_split_fraction": (assigned_pairs - preserved_pairs) / total_pairs if total_pairs else np.nan,
            "pair_loss_to_grey_fraction": (total_pairs - assigned_pairs) / total_pairs if total_pairs else np.nan,
            "merging_cross_reference_pairs": sum(
                int(count) * int(np.sum((b == label) & (a != module) & (a != "grey")))
                for label, count in fragments.items()),
            "module_completely_lost": retained == 0,
        })
    return pd.DataFrame(rows)


def boundary_stages(membership, settings):
    """Reconstruct saved refinement stages, validating against final labels.

    A cache is accepted only when both the retained mask and the partition are
    exactly reproduced. This also guards unresolvable best-profile ties.
    """
    from sklearn.metrics import adjusted_rand_score

    original = membership.observed_module.astype(str).to_numpy()
    own = membership.consensus_membership_confidence.to_numpy(float)
    other = membership.next_best_module_coassignment.to_numpy(float)
    margin = membership.consensus_confidence_margin.to_numpy(float)
    assigned = original != "grey"
    core = assigned & (own >= other) & (own >= settings["core_membership_min"])
    peripheral = assigned & (own >= other) & ~core & (own >= settings["peripheral_membership_min"]) & (margin >= settings["peripheral_membership_margin"])
    pre = np.where(core | peripheral, original, "grey")
    final = pre.copy()
    counts = pd.Series(pre[pre != "grey"]).value_counts()
    small = counts.index[counts < settings["min_module_size"]]
    final[np.isin(final, small)] = "grey"
    saved = membership.consensus_module.astype(str).to_numpy()
    if not np.array_equal(final == "grey", saved == "grey") or adjusted_rand_score(final, saved) != 1:
        raise ValueError("Saved consensus cannot be reconstructed exactly; refit required")
    reason = np.full(len(original), "retained", dtype=object)
    reason[~assigned] = "backbone_grey"
    reason[assigned & (pre == "grey")] = "coassignment_threshold"
    reason[(pre != "grey") & (final == "grey")] = "small_module_after_filter"
    return pre, final, reason


def retain_ranked(labels, scores, *, fraction=None, target=None):
    """Within-backbone ranked retention; ties use the fixed input event order."""
    labels = checked_labels(labels)
    scores = np.nan_to_num(np.asarray(scores, dtype=float), nan=-np.inf)
    if len(scores) != len(labels) or (fraction is None) == (target is None):
        raise ValueError("Supply matching scores and exactly one retention rule")
    if fraction is not None and not 0 <= fraction <= 1:
        raise ValueError("fraction must be in [0, 1]")
    if target is not None and len(target) != len(labels):
        raise ValueError("target must match event order")
    result = np.full(len(labels), "grey", dtype=object)
    for module in sorted(set(labels) - {"grey"}):
        pos = np.flatnonzero(labels == module)
        count = int(np.floor(len(pos) * fraction)) if target is None else int(np.sum(np.asarray(target)[pos] != "grey"))
        keep = pos[np.argsort(-scores[pos], kind="stable")[:count]]
        result[keep] = module
    return result


def core_kme_rescue(membership, settings):
    """Development-only alternative using existing core-evidence thresholds.

    Rescue may retain backbone members, never create or move modules. Missing
    evidence and nonpositive separation cannot rescue an event. This is not a
    calibrated probability or a production default.
    """
    pre, _, _ = boundary_stages(membership, settings)
    original = membership.observed_module.astype(str).to_numpy()
    allowed = (
        (original != "grey")
        & (membership.consensus_confidence_margin.to_numpy(float) > 0)
        & (membership.bootstrap_mean_abs_kME.to_numpy(float) >= settings["core_kme_mean_min"])
        & (membership.bootstrap_q10_abs_kME.to_numpy(float) >= settings["core_kme_lower_quantile_min"])
        & (membership.bootstrap_kME_sign_consistency.to_numpy(float) >= settings["core_kme_sign_consistency_min"])
    )
    result = np.where((pre != "grey") | allowed, original, "grey")
    counts = pd.Series(result[result != "grey"]).value_counts()
    result[np.isin(result, counts.index[counts < settings["min_module_size"]])] = "grey"
    return result
