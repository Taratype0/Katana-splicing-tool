"""Read-level beta-binomial contrasts without duplicating shared controls."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.linalg import null_space
from scipy.optimize import minimize
from scipy.special import betaln, expit, gammaln, logit
from scipy.stats import chi2

from katana_next.science.bulk import parse_number_list


def interaction_observations(row_a, row_b, contrast_a, contrast_b, samples):
    records, arm_sets, offsets = {}, [], []
    biological = {s.sample_id: s.biological_sample_id for s in samples}
    for row, contrast in ((row_a, contrast_a), (row_b, contrast_b)):
        li, ls = float(row.get("IncFormLen", np.nan)), float(row.get("SkipFormLen", np.nan))
        if not np.isfinite([li, ls]).all() or min(li, ls) <= 0:
            raise ValueError("missing_effective_lengths")
        offsets.append(float(np.log(li / ls)))
        for arm in ("numerator", "denominator"):
            ids = tuple(getattr(contrast, arm).sample_ids)
            arm_sets.append(frozenset(ids))
            index = 1 if (arm == "numerator") != bool(row.get("reverse_direction", False)) else 2
            ijc = parse_number_list(row.get(f"IJC_SAMPLE_{index}", ""))
            sjc = parse_number_list(row.get(f"SJC_SAMPLE_{index}", ""))
            if len(ijc) != len(ids) or len(sjc) != len(ids):
                raise ValueError("replicate_id_count_mismatch")
            for sample, inc, skip in zip(ids, ijc, sjc, strict=True):
                if not np.isfinite([inc, skip]).all():
                    continue
                if min(inc, skip) < 0 or inc != int(inc) or skip != int(skip):
                    raise ValueError("invalid_integer_counts")
                value = (inc, inc + skip)
                if sample in records and records[sample] != value:
                    raise ValueError("conflicting_shared_sample_counts")
                records[sample] = value
    if not np.isclose(offsets[0], offsets[1], atol=1e-10, rtol=0):
        raise ValueError("inconsistent_effective_lengths")
    groups = list(dict.fromkeys(arm_sets))
    if any(a & b for i, a in enumerate(groups) for b in groups[i+1:]):
        raise ValueError("partially_overlapping_arms_require_explicit_design")
    observed_ids = list(records)
    bio = [biological.get(s, s) for s in observed_ids]
    if len(bio) != len(set(bio)):
        raise ValueError("repeated_biological_units_require_blocked_model")
    weights = np.zeros(len(groups))
    for arm, sign in zip(arm_sets, (1, -1, -1, 1), strict=True):
        weights[groups.index(arm)] += sign
    frame = pd.DataFrame([{"sample_id": s, "IJC": value[0], "total": value[1],
                           "group": next(i for i, group in enumerate(groups) if s in group)}
                          for s, value in records.items()], columns=["sample_id", "IJC", "total", "group"])
    return frame, weights, offsets[0]


def fit_read_interaction(observations, contrast, offset=0.0, *, min_reads=10, min_replicates=2):
    """H0: (A_num-A_den)-(B_num-B_den)=0 on length-adjusted logit PSI.

    The one-df chi-square LRT is asymptotic. Low replicate numbers still limit
    calibration; this is not a permutation-calibrated or causal test.
    """
    frame = observations.loc[observations.total.ge(min_reads)].copy()
    g = len(contrast)
    counts = frame.groupby("group").size().reindex(range(g), fill_value=0)
    result = {"status": "insufficient_read_supported_replicates", "pvalue": np.nan,
              "effect_log_odds_A_minus_B": np.nan, "rho": np.nan,
              "n_unique_samples": len(frame), "min_group_replicates": int(counts.min()),
              "likelihood_ratio": np.nan, "full_loglik": np.nan, "null_loglik": np.nan}
    if counts.min() < min_replicates:
        return result
    c = np.asarray(contrast, float)
    if not np.any(c):
        result["status"] = "identical_contrasts"
        return result
    y, n = frame.IJC.to_numpy(float), frame.total.to_numpy(float)
    design = np.eye(g)[frame.group.to_numpy(int)]
    basis = null_space(c[None, :])
    group_rates = frame.groupby("group").apply(lambda x: (x.IJC.sum()+.5)/(x.total.sum()+1), include_groups=False).reindex(range(g)).to_numpy()
    initial = logit(np.clip(group_rates, 1e-5, 1-1e-5)) - offset
    combinations = gammaln(n+1)-gammaln(y+1)-gammaln(n-y+1)

    def objective(p, x):
        mu = np.clip(expit(x @ p[:-1] + offset), 1e-8, 1-1e-8)
        concentration = np.exp(p[-1])
        a, b = mu*concentration, (1-mu)*concentration
        return -float((combinations + betaln(y+a, n-y+b)-betaln(a, b)).sum())

    def fit(x, init):
        candidates = [minimize(objective, np.r_[init, concentration], args=(x,), method="L-BFGS-B",
                               bounds=[(-25, 25)]*x.shape[1]+[(-4, 12)],
                               options={"maxiter": 350, "ftol": 1e-10}) for concentration in (3., 7.)]
        valid = [fit for fit in candidates if fit.success and np.isfinite(fit.fun)]
        return min(valid, key=lambda fit: fit.fun) if valid else None

    null = fit(design @ basis, basis.T @ initial)
    if null is None:
        result["status"] = "null_optimizer_failed"
        return result
    full = fit(design, initial)
    if full is None or full.fun > null.fun + 1e-5:
        result["status"] = "full_optimizer_failed_or_non_nested_likelihood"
        return result
    statistic = max(0., 2*(null.fun-full.fun))
    result.update(status="fit", pvalue=float(chi2.sf(statistic, 1)),
                  effect_log_odds_A_minus_B=float(c @ full.x[:-1]), rho=float(1/(1+np.exp(full.x[-1]))),
                  likelihood_ratio=statistic, full_loglik=-float(full.fun), null_loglik=-float(null.fun),
                  rho_at_boundary=bool(abs(full.x[-1]+4) < .01 or abs(full.x[-1]-12) < .01))
    return result
