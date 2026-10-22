"""Held-out biological-repeat prediction; never a latent-PSI truth claim."""
import hashlib

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.special import betaln

from katana_next.science.bulk import parse_number_list
from katana_next.science.psi_calibration import length_correct
from katana_next.science.hierarchical_prediction import hierarchical_draws, paired_moment_rho


def gene_fold(gene, folds):
    return int(hashlib.sha256(str(gene).encode()).hexdigest()[:16], 16) % folds


def train_dispersion(i, s):
    """Shared rho likelihood conditional on training-only empirical event means."""
    n = i+s
    valid = np.isfinite(n) & (n > 0)
    usable = valid.sum(axis=1) >= 2
    i, n, valid = i[usable], n[usable], valid[usable]
    if len(i) < 10:
        return None
    mu = (np.where(valid, i, 0).sum(axis=1)+.5)/(np.where(valid, n, 0).sum(axis=1)+1)
    def objective(log_concentration):
        k = np.exp(log_concentration)
        a, b = mu[:, None]*k, (1-mu[:, None])*k
        terms = betaln(i+a, n-i+b)-betaln(a, b)
        return -float(np.where(valid, terms, 0).sum())
    fit = minimize_scalar(objective, bounds=(-2, 14), method="bounded")
    if not fit.success or not np.isfinite(fit.fun):
        return None
    return {"rho": float(1/(1+np.exp(fit.x))), "training_events": int(len(i)),
            "rho_at_search_boundary": bool(fit.x < -1.99 or fit.x > 13.99)}


def predictive_draws(training_i, training_s, held_n, *, rho, rng, draws):
    a = np.nansum(training_i, axis=1)+.5
    b = np.nansum(training_s, axis=1)+.5
    q = rng.beta(a[:, None], b[:, None], (len(a), draws))
    if rho > 0:
        k = 1/rho-1
        q = rng.beta(np.clip(q, 1e-10, 1-1e-10)*k, np.clip(1-q, 1e-10, 1-1e-10)*k)
    return rng.binomial(held_n[:, None].astype(np.int64), q)


def heldout_arm(frame, sample_ids, biological_ids, *, index, folds=5, draws=2000, seed=1701, report=None, model="legacy"):
    if model not in {"legacy", "hierarchical", "all"}:
        raise ValueError("Unknown held-repeat prediction model")
    if len(set(biological_ids)) != len(biological_ids):
        raise ValueError("Aggregate technical partitions before biological-repeat prediction")
    if len(sample_ids) < 3:
        raise ValueError("At least three independent biological replicates per arm are required")
    if folds < 2 or draws < 100:
        raise ValueError("At least two gene folds and 100 predictive draws are required")
    frame = frame.sort_values("event_id").copy()
    if frame.event_id.duplicated().any():
        raise ValueError("Duplicate exact events in one comparison arm")
    vectors = [frame[f"{kind}_SAMPLE_{index}"].map(parse_number_list).tolist() for kind in ("IJC", "SJC")]
    if any(len(v) != len(sample_ids) for values in vectors for v in values):
        raise ValueError("Declared sample order and count vectors differ")
    i, s = map(lambda x: np.array(x, dtype=float), vectors)
    for counts in (i,s):
        if np.any(np.isinf(counts) | (counts < 0) | (np.isfinite(counts) & (counts != np.floor(counts)))):
            raise ValueError("Read counts must be integer, nonnegative or missing")
    li, ls = frame.IncFormLen.to_numpy(float), frame.SkipFormLen.to_numpy(float)
    length_correct(np.zeros(len(frame)), li, ls)
    known_gene = frame.gene_id.notna() & frame.gene_id.astype(str).str.strip().ne("")
    assignments = frame.gene_id.map(lambda g: gene_fold(g, folds)).to_numpy()
    records, fits = [], []
    for held, sample in enumerate(sample_ids):
        train = [j for j in range(len(sample_ids)) if j != held]
        ti, ts = i[:, train].copy(), s[:, train].copy()
        valid_train = np.isfinite(ti+ts) & ((ti+ts) > 0)
        ti[~valid_train], ts[~valid_train] = np.nan, np.nan
        valid = known_gene.to_numpy() & (valid_train.sum(axis=1) >= 2) & np.isfinite(i[:, held]+s[:, held]) & ((i[:, held]+s[:, held]) > 0)
        for f in range(folds):
            training_mask = known_gene.to_numpy() & (assignments != f)
            fit = train_dispersion(ti[training_mask], ts[training_mask]) if model != "hierarchical" else None
            moment = paired_moment_rho(ti[training_mask], ts[training_mask]) if model != "legacy" else None
            fit_methods = ([("conditional_plugin", fit)] if model != "hierarchical" else [])
            if model != "legacy":
                fit_methods.append(("paired_moment_hierarchical", moment))
            for estimator, fitted in fit_methods:
                fits.append({"held_sample_id":sample, "test_gene_fold":f, "estimator":estimator,
                    "status":"fit" if fitted else "insufficient_training_or_optimizer_failure",
                    "excluded_from_fit": "held sample in all genes and every test-fold gene", **(fitted or {})})
            methods = [("pooled_count_only", 0., predictive_draws)] if model != "hierarchical" else []
            if model != "hierarchical":
                methods.append(("training_gene_crossfit_dispersion", fit["rho"] if fit else None, predictive_draws))
            if model != "legacy":
                methods.append(("training_gene_crossfit_hierarchical", moment["rho"] if moment else None, hierarchical_draws))
            chosen = np.flatnonzero((assignments == f) & valid)
            rng_seed = seed+held*1009+f*97
            for begin in range(0, len(chosen), 128):
                if report:
                    report(sample, f, begin, len(chosen))
                ix = chosen[begin:begin+128]
                n = (i+s)[ix,held]
                truth = length_correct(i[ix,held]/n, li[ix],ls[ix])
                prediction = length_correct((np.nansum(ti[ix],axis=1)+.5)/(np.nansum(ti[ix]+ts[ix],axis=1)+1),li[ix],ls[ix])
                for method,rho,sampler in methods:
                    if rho is None:
                        for event_index in ix:
                            records.append({"event_id":frame.event_id.iloc[event_index],
                                "gene_id":frame.gene_id.iloc[event_index], "held_sample_id":sample,
                                "test_gene_fold":f, "method":method,
                                "reason":"insufficient_training_or_dispersion_failure", "covered":None})
                        continue
                    samples = sampler(ti[ix],ts[ix],n,rho=rho,
                        rng=np.random.default_rng(rng_seed+begin),draws=draws)
                    lo,hi = np.quantile(samples,[.025,.975],axis=1,method="inverted_cdf")
                    lo,hi = length_correct(lo/n,li[ix],ls[ix]),length_correct(hi/n,li[ix],ls[ix])
                    for pos,event_index in enumerate(ix):
                        records.append({"event_id":frame.event_id.iloc[event_index],"gene_id":frame.gene_id.iloc[event_index],
                            "held_sample_id":sample,"test_gene_fold":f,"method":method,"rho":rho,
                            "training_sample_ids":";".join(sample_ids[j] for j in train),
                            "held_reads":n[pos],"observed_held_PSI":truth[pos],"predicted_PSI":prediction[pos],
                            "predictive_lower":lo[pos],"predictive_upper":hi[pos],
                            "covered":bool(lo[pos] <= truth[pos] <= hi[pos]),"width":hi[pos]-lo[pos],
                            "absolute_error":abs(prediction[pos]-truth[pos]),"draws":draws,
                            "point_estimator":"pooled-count Jeffreys mean, unchanged across interval methods"})
        for pos in np.flatnonzero(~valid):
            records.append({"event_id":frame.event_id.iloc[pos],"gene_id":frame.gene_id.iloc[pos],"held_sample_id":sample,
                "method":"not_evaluable", "reason":"missing_gene" if not known_gene.iloc[pos] else "insufficient_training_or_held_reads"})
    return pd.DataFrame(records), pd.DataFrame(fits)


def prediction_summary(frame):
    rows=[]
    for (comparison, arm, method), part in frame.groupby(["comparison_id","arm","method"],dropna=False):
        eligible = part[part.covered.notna()] if "covered" in part else part.iloc[:0]
        rows.append({"comparison_id":comparison,"arm":arm,"method":method,"evaluated_predictions":len(eligible),
            "events":eligible.event_id.nunique(),"held_biological_replicates":eligible.held_sample_id.nunique(),
            "observed_predictive_coverage":eligible.covered.astype(float).mean() if len(eligible) else np.nan,
            "mean_interval_width":eligible.width.mean() if len(eligible) else np.nan,
            "MAE":eligible.absolute_error.mean() if len(eligible) else np.nan,
            "not_evaluable_predictions":len(part)-len(eligible),
            "scope":"Observed held-repeat PSI prediction, not latent PSI truth; overlapping folds and shared controls are not independent tests"})
    return pd.DataFrame(rows)
