"""Count-conditional uncertainty for effective-length-normalized rMATS PSI.

Intervals describe an individual sample's latent PSI, not a condition mean,
biological dispersion, or a differential-splicing test.
"""
import numpy as np
import pandas as pd
from scipy.stats import beta


def length_correct(q, inclusion_length, skipping_length):
    q, li, ls = np.broadcast_arrays(np.asarray(q, float), inclusion_length, skipping_length)
    if np.any(~np.isfinite(li) | ~np.isfinite(ls) | (li <= 0) | (ls <= 0)):
        raise ValueError("Effective inclusion/skipping lengths must be finite and positive")
    if np.any(np.isfinite(q) & ((q < 0) | (q > 1))):
        raise ValueError("Read fraction must be between zero and one")
    return q * ls / (q * ls + (1 - q) * li)


def psi_intervals(ijc, sjc, inclusion_length, skipping_length, level=.95):
    i, s, li, ls = np.broadcast_arrays(np.asarray(ijc, float), np.asarray(sjc, float),
                                      inclusion_length, skipping_length)
    if not 0 < level < 1:
        raise ValueError("Interval level must be between zero and one")
    for counts in (i, s):
        if np.any(np.isinf(counts) | (counts < 0) | (np.isfinite(counts) & (counts != np.floor(counts)))):
            raise ValueError("IJC/SJC must be nonnegative integer counts or missing")
    length_correct(np.zeros(i.shape), li, ls)
    valid = np.isfinite(i) & np.isfinite(s) & ((i + s) > 0)
    n = i + s
    q = np.divide(i, n, out=np.full(i.shape, np.nan), where=valid)
    # Monotone transformation preserves posterior interval probability.
    lower = beta.ppf((1-level)/2, i+.5, s+.5)
    upper = beta.ppf(1-(1-level)/2, i+.5, s+.5)
    return {"psi": length_correct(q, li, ls),
            "psi_posterior_median": np.where(valid, length_correct(beta.ppf(.5, i+.5, s+.5), li, ls), np.nan),
            "psi_lower": np.where(valid, length_correct(lower, li, ls), np.nan),
            "psi_upper": np.where(valid, length_correct(upper, li, ls), np.nan),
            "measurable": valid}


def length_calibration_grid(events=1000, replicates=3, seed=1701):
    """Predeclared latent sample-PSI truth, with varying lengths and dispersion."""
    rng = np.random.default_rng(seed)
    rows = []
    for ratio in (.25, .5, 1., 2., 4.):
        for depth in (5, 20, 100, 500):
            for rho in (0., .05, .2):
                center = rng.uniform(.05, .95, (events, 1))
                truth = (np.broadcast_to(center, (events, replicates)).copy() if rho == 0
                         else rng.beta(center*(1/rho-1), (1-center)*(1/rho-1), (events, replicates)))
                # Inverse length transform: sample reads, not normalized pseudo-counts.
                q = truth*ratio / (truth*ratio + 1-truth)
                i = rng.binomial(depth, q)
                result = psi_intervals(i, depth-i, ratio, 1.)
                hit = ((result["psi_lower"] <= truth) & (truth <= result["psi_upper"])).mean(axis=1)
                rows.append({"length_ratio": ratio, "depth": depth, "biological_psi_rho": rho,
                    "events": events, "replicates": replicates, "nominal_coverage": .95,
                    "observed_coverage": float(hit.mean()), "coverage_MCSE": float(hit.std(ddof=1)/np.sqrt(events)),
                    "mean_interval_width": float(np.mean(result["psi_upper"]-result["psi_lower"])),
                    "raw_PSI_MAE": float(np.mean(abs(result["psi"]-truth))),
                    "target": "latent individual-sample length-normalized PSI; not a condition mean"})
    return pd.DataFrame(rows)


def audit_rmats_lengths(frame):
    """Verify caller parity and expose measured intervals without changing calls."""
    required = {"IncFormLen", "SkipFormLen", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1",
                "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel2"}
    if not required <= set(frame):
        raise ValueError(f"Missing rMATS columns: {sorted(required-set(frame))}")
    rows = []
    def vector(value):
        return [float(x) if x not in {"NA", "nan", ""} else np.nan for x in str(value).split(",")]
    for index, row in frame.iterrows():
        for arm in (1, 2):
            i, s, caller = [vector(row[col]) for col in (f"IJC_SAMPLE_{arm}", f"SJC_SAMPLE_{arm}", f"IncLevel{arm}")]
            if len({len(i), len(s), len(caller)}) != 1:
                raise ValueError(f"Unmatched replicate vectors at row {index}, arm {arm}")
            for j in range(len(i)):
                rows.append({"source_row": str(index), "arm": arm, "replicate_index": j,
                    "IJC": i[j], "SJC": s[j], "IncFormLen": row.IncFormLen, "SkipFormLen": row.SkipFormLen,
                    "caller_PSI": caller[j]})
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    metrics = psi_intervals(result.IJC, result.SJC, result.IncFormLen, result.SkipFormLen)
    for key, value in metrics.items():
        result[key] = value
    result["caller_parity"] = (result.psi.isna() & result.caller_PSI.isna()) | (abs(result.psi-result.caller_PSI) <= .000500001)
    return result
