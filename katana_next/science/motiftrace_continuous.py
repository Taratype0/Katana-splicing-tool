"""Continuous usage association with expression weights and gene-block nulls.

This is a distinct, versioned alternative to matched foreground enrichment.
Its resampling units are genes, not extra biological replicates. P values are
obtained from a predeclared maximum-statistic calibration, not a normal tail.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def continuous_inputs(events, inclusion, *, source_group1_is_numerator):
    if not source_group1_is_numerator:
        raise ValueError("Explicit, source-supported RNA group identity required")
    if events.event_id.duplicated().any() or events.gene_id.isna().any():
        raise ValueError("Unique event identity and nonmissing gene blocks required")
    if inclusion.index.tolist() != events.event_id.tolist():
        raise ValueError("Continuous measurements must match event IDs, never position-guess")
    if inclusion.shape[1] != 8:
        raise ValueError("This validation adapter expects four included and four total TPM columns")
    included, total = inclusion.iloc[:, :4].to_numpy(float), inclusion.iloc[:, 4:].to_numpy(float)
    if not np.isfinite(included).all() or not np.isfinite(total).all() or (total <= 0).any():
        raise ValueError("Invalid or missing compatible-transcript TPM")
    if ((included < 0) | (included > total + 1e-8)).any():
        raise ValueError("Included TPM must be a subset of compatible transcript TPM")
    fractions = included / total
    tau = float(np.median(total))
    sample_weights = total / (total + tau)
    # Symmetric over samples: null swaps cannot change the weight or its scale.
    expression_weight = 4. / np.sum(1. / sample_weights, axis=1)

    def slots(column):
        result = []
        for value in column:
            tokens = str(value).split(",")
            if len(tokens) != 2:
                raise ValueError("Require exactly two released PSI slots per source arm")
            values = np.array([float(t) if t.lower() not in {"na", "nan"} else np.nan for t in tokens])
            if ((values[np.isfinite(values)] < 0) | (values[np.isfinite(values)] > 1)).any():
                raise ValueError("PSI must be a fraction or explicitly missing")
            result.append(values)
        return np.asarray(result)

    kd, ctrl = slots(events.IncLevel1), slots(events.IncLevel2)
    frame = events[["event_id", "gene_id", "role", "dpsi", "fdr"]].copy()
    for j, label in enumerate(("Ctrl1", "Ctrl2", "KD1", "KD2")):
        frame[f"transcript_inclusion_fraction_{label}"] = fractions[:, j]
        frame[f"compatible_transcript_TPM_{label}"] = total[:, j]
    for j in range(2):
        # Source PSI slot order is not silently equated to RSEM replicate IDs.
        frame[f"PSI_KD_source_slot{j+1}"] = kd[:, j]
        frame[f"PSI_Ctrl_source_slot{j+1}"] = ctrl[:, j]
    frame["delta_PSI_KD_minus_Ctrl"] = kd.mean(axis=1) - ctrl.mean(axis=1)
    frame["delta_transcript_fraction_KD_minus_Ctrl"] = fractions[:, 2:].mean(axis=1) - fractions[:, :2].mean(axis=1)
    frame["expression_weight"] = expression_weight
    frame["expression_saturation_TPM"] = tau
    frame["PSI_eligible"] = np.isfinite(kd).all(axis=1) & np.isfinite(ctrl).all(axis=1)
    return frame


def gene_contributions(scores, response, expression_weight, genes):
    scores = np.asarray(scores, float)
    response = np.asarray(response, float)
    weight = np.asarray(expression_weight, float)
    genes = np.asarray(genes)
    if scores.ndim != 2 or len(scores) != len(response) or len(weight) != len(response):
        raise ValueError("Event rows and measurements must align")
    if not np.isfinite(scores).all():
        raise ValueError("Sequence opportunity missingness requires an explicit shared mask")
    levels = sorted(set(genes))
    lookup = {g: j for j, g in enumerate(levels)}
    codes = np.array([lookup[g] for g in genes])
    valid = np.isfinite(response) & np.isfinite(weight) & (weight > 0)
    if len(set(genes[valid])) < 20:
        return np.full((len(levels), scores.shape[1]), np.nan), levels, valid
    x, w, y = scores[valid], weight[valid], response[valid]
    centered_x = x - np.average(x, axis=0, weights=w)
    contribution = centered_x * (w*y)[:, None]
    result = np.zeros((len(levels), scores.shape[1]), float)
    np.add.at(result, codes[valid], contribution)
    return result, levels, valid


def score_from_contributions(contribution, flips=None):
    contribution = np.asarray(contribution, float)
    scale = np.sqrt(np.sum(contribution**2, axis=0))
    numerator = contribution.sum(axis=0) if flips is None else np.asarray(flips, float) @ contribution
    result = np.divide(numerator, scale, out=np.full_like(numerator, np.nan), where=np.isfinite(scale) & (scale > 0))
    result[..., scale == 0] = 0
    return result


def calibrated_tail(observed, calibration_maxima):
    maxima = np.asarray(calibration_maxima, float)
    if maxima.ndim != 1 or not len(maxima) or not np.isfinite(maxima).all():
        raise ValueError("A complete independent finite calibration family is required")
    values = np.asarray(observed, float)
    greater_equal = len(maxima) - np.searchsorted(np.sort(maxima), values, side="left")
    result = (1. + greater_equal) / (1. + len(maxima))
    return np.where(np.isfinite(values), result, np.nan)
