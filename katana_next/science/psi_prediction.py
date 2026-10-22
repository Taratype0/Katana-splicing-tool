"""Discovery-only multi-module prediction with explicit inner-CV regularization."""
from __future__ import annotations
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold

ALPHAS = (0.1, 1., 10., 100., 1000.)


def _standardize_fit(x):
    center = x.mean(axis=0)
    scale = x.std(axis=0, ddof=0)
    scale = np.where(scale < 1e-12, 1., scale)
    return center, scale


def fit_multimodule_predictor(scores, traits):
    if not scores.index.is_unique or not traits.index.is_unique or set(scores.index) != set(traits.index):
        raise ValueError("Training scores and traits require identical unique sample IDs")
    x = scores.to_numpy(float)
    if not np.isfinite(x).all():
        raise ValueError("Discovery module scores must be finite")
    entries = []
    for trait in traits:
        y = pd.to_numeric(traits.loc[scores.index, trait], errors="raise").to_numpy(float)
        if np.isinf(y).any():
            raise ValueError("Traits cannot contain infinity")
        good = np.isfinite(y)
        xx, yy = x[good], y[good]
        if len(yy) < 15 or np.std(yy) < 1e-12:
            continue
        options = (None,) + ALPHAS  # None is the intercept-only candidate.
        errors = np.zeros((5, len(options)))
        for fold, (tr, va) in enumerate(KFold(5, shuffle=True, random_state=20261017).split(yy)):
            center, scale = _standardize_fit(xx[tr])
            for col, alpha in enumerate(options):
                pred = np.full(len(va), yy[tr].mean())
                if alpha is not None and xx.shape[1]:
                    fit = Ridge(alpha=alpha).fit((xx[tr]-center)/scale, yy[tr])
                    pred = fit.predict((xx[va]-center)/scale)
                errors[fold, col] = np.mean((yy[va]-pred)**2)
        avg = errors.mean(axis=0)
        best = int(np.argmin(avg))
        bound = avg[best] + errors[:, best].std(ddof=1)/np.sqrt(5)
        eligible = np.flatnonzero(avg <= bound + 1e-15)
        # Prefer intercept-only, otherwise the strongest regularization within 1 SE.
        chosen = 0 if 0 in eligible else int(eligible[-1])
        alpha = options[chosen]
        center, scale = _standardize_fit(xx)
        coefficients = np.zeros(xx.shape[1])
        intercept = float(yy.mean())
        if alpha is not None and xx.shape[1]:
            fit = Ridge(alpha=alpha).fit((xx-center)/scale, yy)
            coefficients = fit.coef_
            intercept = float(fit.intercept_)
        entries.append(dict(trait=str(trait), modules=list(scores.columns),
            centers=center.tolist(), scales=scale.tolist(), coefficients=coefficients.tolist(),
            intercept=intercept, training_mean=float(yy.mean()), training_n=len(yy), alpha=alpha,
            selection="five_inner_folds_one_SE_conditional_on_frozen_discovery_modules",
            inner_candidates=["intercept"]+list(ALPHAS), inner_mean_MSE=avg.tolist(),
            inner_selected_MSE=float(avg[chosen]), inner_one_SE_bound=float(bound)))
    return entries


def predict_multimodule(entry, scores):
    coef = np.array(entry["coefficients"])
    used = np.abs(coef) > 0
    if not used.any():
        return np.full(len(scores), entry["training_mean"])
    x = scores.reindex(columns=entry["modules"]).to_numpy(float)
    z = (x[:, used]-np.array(entry["centers"])[used])/np.array(entry["scales"])[used]
    pred = entry["intercept"] + z @ coef[used]
    pred[~np.isfinite(z).all(axis=1)] = np.nan
    return pred
