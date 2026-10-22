from __future__ import annotations

from dataclasses import dataclass
from itertools import product, permutations

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import chi2, t

from katana_cosplice.module_models import _fit_model


@dataclass(frozen=True)
class SharedControlPermutationResult:
    observed_rank_overlap: int
    observed_rank_same_direction: int
    null_distribution: pd.DataFrame
    empirical_overlap_pvalue: float
    empirical_same_direction_pvalue: float


def _top_indices(values: np.ndarray, size: int) -> np.ndarray:
    finite = np.flatnonzero(np.isfinite(values))
    if size <= 0 or finite.size == 0:
        return np.asarray([], dtype=int)
    size = min(size, finite.size)
    scores = values[finite]
    selected = np.argpartition(scores, -size)[-size:]
    return finite[selected]


def _rank_overlap(
    values: np.ndarray,
    wt_indices: np.ndarray,
    mut_indices: np.ndarray,
    drug_indices: np.ndarray,
    n_mut: int,
    n_drug: int,
) -> tuple[int, int]:
    wt_mean = np.nanmean(values[:, wt_indices], axis=1)
    mut_effect = np.nanmean(values[:, mut_indices], axis=1) - wt_mean
    drug_effect = np.nanmean(values[:, drug_indices], axis=1) - wt_mean
    mut_top = _top_indices(np.abs(mut_effect), n_mut)
    drug_top = _top_indices(np.abs(drug_effect), n_drug)
    shared = np.intersect1d(mut_top, drug_top, assume_unique=False)
    same_direction = int(np.sum(mut_effect[shared] * drug_effect[shared] > 0))
    return int(shared.size), same_direction


def blocked_shared_control_permutation(
    psi: pd.DataFrame,
    *,
    n_mut: int,
    n_drug: int,
) -> SharedControlPermutationResult:
    """Enumerate condition-label permutations within three replicate blocks.

    Columns must be ordered as WT, MUT, drug for replicate 1, followed by the
    same order for replicates 2 and 3. Fixing the selected-set sizes isolates
    the extra overlap caused by shared controls and event-level effect
    propensity instead of conflating it with different marginal set sizes.
    """

    if psi.shape[1] != 9:
        raise ValueError("Shared-control permutation requires three groups x three replicates.")
    values = psi.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    actual_wt = np.asarray([0, 3, 6], dtype=int)
    actual_mut = np.asarray([1, 4, 7], dtype=int)
    actual_drug = np.asarray([2, 5, 8], dtype=int)
    observed_overlap, observed_same = _rank_overlap(
        values,
        actual_wt,
        actual_mut,
        actual_drug,
        n_mut,
        n_drug,
    )

    block_permutations = tuple(permutations((0, 1, 2)))
    rows: list[dict[str, int]] = []
    for permutation_id, block_assignment in enumerate(product(block_permutations, repeat=3)):
        wt_indices = np.asarray(
            [block * 3 + block_assignment[block][0] for block in range(3)],
            dtype=int,
        )
        mut_indices = np.asarray(
            [block * 3 + block_assignment[block][1] for block in range(3)],
            dtype=int,
        )
        drug_indices = np.asarray(
            [block * 3 + block_assignment[block][2] for block in range(3)],
            dtype=int,
        )
        overlap, same_direction = _rank_overlap(
            values,
            wt_indices,
            mut_indices,
            drug_indices,
            n_mut,
            n_drug,
        )
        rows.append(
            {
                "permutation_id": permutation_id,
                "rank_matched_overlap": overlap,
                "rank_matched_same_direction": same_direction,
            }
        )
    null = pd.DataFrame(rows)
    overlap_pvalue = float(
        (1 + null["rank_matched_overlap"].ge(observed_overlap).sum()) / (len(null) + 1)
    )
    same_direction_pvalue = float(
        (1 + null["rank_matched_same_direction"].ge(observed_same).sum()) / (len(null) + 1)
    )
    return SharedControlPermutationResult(
        observed_rank_overlap=observed_overlap,
        observed_rank_same_direction=observed_same,
        null_distribution=null,
        empirical_overlap_pvalue=overlap_pvalue,
        empirical_same_direction_pvalue=same_direction_pvalue,
    )


def _likelihood_ratio_pvalue(full_nll: float, null_nll: float, df: int = 1) -> float:
    statistic = max(0.0, 2.0 * (null_nll - full_nll))
    return float(chi2.sf(statistic, df))


def welch_equivalence_test(
    left: np.ndarray,
    right: np.ndarray,
    *,
    margin: float,
) -> dict[str, float | bool]:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    left = left[np.isfinite(left)]
    right = right[np.isfinite(right)]
    if left.size < 2 or right.size < 2:
        return {
            "direct_delta_psi": np.nan,
            "equivalence_ci90_low": np.nan,
            "equivalence_ci90_high": np.nan,
            "equivalence_pvalue": np.nan,
            "practical_similarity": False,
        }
    delta = float(np.mean(left) - np.mean(right))
    left_term = float(np.var(left, ddof=1) / left.size)
    right_term = float(np.var(right, ddof=1) / right.size)
    variance = left_term + right_term
    if variance <= 0:
        equivalent = abs(delta) < margin
        return {
            "direct_delta_psi": delta,
            "equivalence_ci90_low": delta,
            "equivalence_ci90_high": delta,
            "equivalence_pvalue": 0.0 if equivalent else 1.0,
            "practical_similarity": bool(abs(delta) <= margin),
        }
    standard_error = float(np.sqrt(variance))
    denominator = 0.0
    if left.size > 1:
        denominator += left_term**2 / (left.size - 1)
    if right.size > 1:
        denominator += right_term**2 / (right.size - 1)
    degrees_freedom = variance**2 / denominator if denominator > 0 else left.size + right.size - 2
    critical = float(t.ppf(0.95, degrees_freedom))
    lower = delta - critical * standard_error
    upper = delta + critical * standard_error
    lower_test = (delta + margin) / standard_error
    upper_test = (delta - margin) / standard_error
    pvalue = max(
        float(t.sf(lower_test, degrees_freedom)),
        float(t.cdf(upper_test, degrees_freedom)),
    )
    return {
        "direct_delta_psi": delta,
        "equivalence_ci90_low": lower,
        "equivalence_ci90_high": upper,
        "equivalence_pvalue": pvalue,
        "practical_similarity": bool(abs(delta) <= margin),
    }


def fit_shared_control_beta_binomial(
    ijc: np.ndarray,
    sjc: np.ndarray,
    psi: np.ndarray,
    *,
    equivalence_margin: float,
    max_iterations: int = 300,
) -> dict[str, object]:
    """Fit WT, MUT and inhibitor replicate counts in one beta-binomial model."""

    ijc = np.asarray(ijc, dtype=float)
    sjc = np.asarray(sjc, dtype=float)
    psi = np.asarray(psi, dtype=float)
    totals = ijc + sjc
    valid = np.isfinite(ijc) & np.isfinite(totals) & (totals >= 1)
    group = np.repeat(np.arange(3), 3)
    valid_counts = pd.Series(group[valid]).value_counts()
    if len(valid_counts) < 3 or (valid_counts < 2).any():
        return {
            "joint_model_status": "insufficient_read_supported_replicates",
            "joint_model_warning": str(valid_counts.to_dict()),
        }

    y = ijc[valid]
    n = totals[valid]
    group_valid = group[valid]
    full = np.column_stack(
        [
            np.ones(valid.sum(), dtype=float),
            (group_valid == 1).astype(float),
            (group_valid == 2).astype(float),
        ]
    )
    null_mut_wt = np.column_stack(
        [np.ones(valid.sum(), dtype=float), (group_valid == 2).astype(float)]
    )
    null_drug_wt = np.column_stack(
        [np.ones(valid.sum(), dtype=float), (group_valid == 1).astype(float)]
    )
    null_mut_drug = np.column_stack(
        [np.ones(valid.sum(), dtype=float), (group_valid != 0).astype(float)]
    )
    try:
        full_parameters, full_nll, full_ok, full_message = _fit_model(
            full,
            y,
            n,
            max_iterations,
        )
        _, mut_null_nll, mut_ok, mut_message = _fit_model(
            null_mut_wt,
            y,
            n,
            max_iterations,
        )
        _, drug_null_nll, drug_ok, drug_message = _fit_model(
            null_drug_wt,
            y,
            n,
            max_iterations,
        )
        _, direct_null_nll, direct_ok, direct_message = _fit_model(
            null_mut_drug,
            y,
            n,
            max_iterations,
        )
        coefficients = full_parameters[:-1]
        model_means = expit(
            np.asarray(
                [
                    coefficients[0],
                    coefficients[0] + coefficients[1],
                    coefficients[0] + coefficients[2],
                ]
            )
        )
        rho = 1e-7 + (0.95 - 1e-7) * expit(full_parameters[-1])
        equivalence = welch_equivalence_test(
            psi[3:6],
            psi[6:9],
            margin=equivalence_margin,
        )
        return {
            "joint_model_status": "converged"
            if full_ok and mut_ok and drug_ok and direct_ok
            else "optimizer_warning",
            "joint_model_warning": "; ".join(
                [
                    f"full={full_message}",
                    f"mut_vs_wt={mut_message}",
                    f"drug_vs_wt={drug_message}",
                    f"mut_vs_drug={direct_message}",
                ]
            ),
            "beta_binomial_rho": float(rho),
            "model_mean_WT": float(model_means[0]),
            "model_mean_MUT": float(model_means[1]),
            "model_mean_inhibitor": float(model_means[2]),
            "joint_MUT_vs_WT_pvalue": _likelihood_ratio_pvalue(full_nll, mut_null_nll),
            "joint_inhibitor_vs_WT_pvalue": _likelihood_ratio_pvalue(
                full_nll,
                drug_null_nll,
            ),
            "direct_MUT_vs_inhibitor_pvalue": _likelihood_ratio_pvalue(
                full_nll,
                direct_null_nll,
            ),
            "joint_same_direction": bool(
                (model_means[1] - model_means[0])
                * (model_means[2] - model_means[0])
                > 0
            ),
            **equivalence,
        }
    except Exception as exc:
        return {
            "joint_model_status": "failed",
            "joint_model_warning": str(exc),
        }
