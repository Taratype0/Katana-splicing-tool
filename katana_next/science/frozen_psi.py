"""Portable discovery-frozen PSI module projection, without validation-set fitting."""
from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

SCHEMA = "katana-frozen-psi-1.0"
EXT_SCHEMA = "katana-frozen-psi-1.1"
ROBUST_SCHEMA = "katana-frozen-psi-1.2"
SCHEMAS = {SCHEMA,EXT_SCHEMA,ROBUST_SCHEMA}


def _digest(model):
    payload = {k: v for k, v in model.items() if k != "sha256"}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, allow_nan=False,
                                    separators=(",", ":")).encode()).hexdigest()


def _samples(matrix, samples, *, allow_empty_events=False):
    if not {"sample_id", "biological_sample_id"}.issubset(samples):
        raise ValueError("Declare sample_id and biological_sample_id")
    if not matrix.index.is_unique or not matrix.columns.is_unique or len(matrix.columns)==0 or (matrix.empty and not allow_empty_events):
        raise ValueError("PSI requires unique nonempty event and sample axes")
    for key in ("sample_id", "biological_sample_id"):
        if samples[key].isna().any() or samples[key].astype(str).str.strip().eq("").any():
            raise ValueError(f"Missing {key}")
        if not samples[key].is_unique:
            raise ValueError("Aggregate technical replicates before frozen module validation")
    if set(samples.sample_id) != set(matrix.columns):
        raise ValueError("PSI columns must match sample metadata exactly")
    x = matrix.to_numpy(float)
    if np.isinf(x).any() or np.any(np.isfinite(x) & ((x < 0) | (x > 1))):
        raise ValueError("PSI must lie in [0,1]; missing values remain missing")
    return samples.set_index("sample_id").loc[matrix.columns]


def _mean_pairwise(measured):
    if measured.shape[1] < 15 or len(measured) < 2:
        return None, 0
    cor = measured.T.corr(min_periods=15).to_numpy()
    vals = cor[np.triu_indices(len(cor), 1)]
    vals = vals[np.isfinite(vals)]
    return (float(vals.mean()) if len(vals) else None), len(vals)


def freeze_psi_modules(measured, means, baseline, samples, *, feature_space_id,
                       network, filters, traits=None, trait_prediction_method="single_module",
                       reference_scores=None,reference_space_id=None,filter_profile='legacy'):
    """Freeze exactly the feature order and SVD convention of program_scores."""
    if trait_prediction_method not in {"single_module", "multimodule_ridge",'validated_hierarchy'}:
        raise ValueError("Unknown trait_prediction_method")
    ordered = _samples(measured, samples,allow_empty_events=True)
    if not str(feature_space_id).strip():
        raise ValueError("Declare feature_space_id (build, event definition and quantification)")
    selected = baseline.selected_matrix
    labels = baseline.assignments.set_index("event_key").reindex(selected.index).module
    if labels.isna().any() or not set(selected.index).issubset(measured.index):
        raise ValueError("Discovery membership and measured PSI do not align")
    model = dict(schema=SCHEMA, feature_space_id=str(feature_space_id),
                 training_samples=list(selected.columns),
                 training_biological_ids=ordered.biological_sample_id.astype(str).tolist(),
                 network=asdict(network), filters=asdict(filters), modules=[], trait_models=[],
                 score_definition="training-standardized PSI @ oriented_PC1_loading / training_singular_value",
                 filter_profile=filter_profile,
                 network_status='modules_available' if baseline.program_scores.shape[1] else 'no_estimable_modules',
                 limitations="Projection and measured-event coherence are not proof of biological regulation or WGCNA/NetRep preservation.")
    for module in sorted(set(labels).difference({"grey",'gray','0',0,'0.0'})):
        keys = selected.index[labels.eq(module)]
        complete = selected.loc[keys]
        center = complete.mean(axis=1).to_numpy(float)
        scale = complete.std(axis=1, ddof=1).to_numpy(float)
        scale = np.where(scale < 1e-12, 1., scale)
        z = (complete.to_numpy(float)-center[:, None])/scale[:, None]
        u, singular, vt = np.linalg.svd(z, full_matrices=False)
        if singular[0] <= 1e-12:
            raise ValueError(f"Cannot freeze degenerate module {module}")
        loading = u[:, 0]
        if np.corrcoef(vt[0], z.mean(axis=0))[0, 1] < 0:
            loading = -loading
        scores = z.T @ loading / singular[0]
        error = float(np.max(np.abs(scores-baseline.program_scores[module].to_numpy())))
        if error > 1e-8:
            raise ValueError("Frozen loading failed to reproduce discovery scores")
        coherence, pairs = _mean_pairwise(measured.loc[keys])
        model["modules"].append(dict(module=str(module), events=keys.tolist(),
            imputation_means=means.reindex(keys).tolist(), centers=center.tolist(),
            scales=scale.tolist(), loading=loading.tolist(), singular_value=float(singular[0]),
            training_reconstruction_max_error=error, training_mean_pairwise_r=coherence,
            training_valid_pairs=pairs))
        from katana_next.science.psi_robust import calibrate_guarded_projection
        model['modules'][-1]['guarded_projection']=calibrate_guarded_projection(model['modules'][-1],measured)
    if traits is not None:
        if not traits.index.is_unique or set(traits.index) != set(selected.columns):
            raise ValueError("Training trait rows must match discovery samples exactly")
        for trait in traits:
            y = pd.to_numeric(traits.loc[selected.columns, trait], errors="raise").to_numpy(float)
            if np.isinf(y).any():
                raise ValueError("Traits cannot contain infinity")
            valid = np.isfinite(y)
            if valid.sum() < 15 or np.std(y[valid]) < 1e-12:
                continue
            associations = {}
            for module in baseline.program_scores:
                x = baseline.program_scores[module].to_numpy(float)[valid]
                if np.std(x) > 1e-12:
                    associations[module] = abs(float(np.corrcoef(x, y[valid])[0, 1]))
            best = max(associations, key=associations.get) if associations else None
            mean = float(y[valid].mean())
            coefficients = [mean, 0.]
            if best is not None:
                x = baseline.program_scores[best].to_numpy(float)[valid]
                coefficients = np.linalg.lstsq(np.column_stack([np.ones(len(x)), x]), y[valid], rcond=None)[0].tolist()
            model["trait_models"].append(dict(trait=str(trait), module=best,
                intercept=float(coefficients[0]), slope=float(coefficients[1]), training_mean=mean,
                training_n=int(valid.sum()), selection="max_abs_Pearson_on_discovery_only"))
    if trait_prediction_method == "multimodule_ridge":
        if traits is None:
            raise ValueError("Supply training traits for multimodule_ridge")
        from katana_next.science.psi_prediction import fit_multimodule_predictor
        model["schema"] = EXT_SCHEMA
        model["trait_prediction_method"] = trait_prediction_method
        model["multimodule_trait_models"] = fit_multimodule_predictor(baseline.program_scores, traits)
    if trait_prediction_method=='validated_hierarchy':
        if traits is None:raise ValueError('Supply training traits for validated_hierarchy')
        if reference_scores is not None and not str(reference_space_id or '').strip():
            raise ValueError('Declare reference_space_id for supplied reference scores')
        from katana_next.science.psi_robust import fit_reliable_predictors
        model['schema']=ROBUST_SCHEMA;model['trait_prediction_method']=trait_prediction_method
        model['reference_space_id']=reference_space_id
        model['reliable_trait_models']=fit_reliable_predictors(baseline.program_scores,traits,reference_scores)
    model["sha256"] = _digest(model)
    return model


def save_frozen_model(model, path):
    if model.get("sha256") != _digest(model):
        raise ValueError("Frozen model checksum mismatch")
    Path(path).write_text(json.dumps(model, indent=2, allow_nan=False), encoding="utf-8")


def load_frozen_model(path):
    model = json.loads(Path(path).read_text(encoding="utf-8"))
    if model.get("schema") not in SCHEMAS or model.get("sha256") != _digest(model):
        raise ValueError("Unsupported or modified frozen PSI model")
    return model


def project_psi_modules(model, matrix, samples, *, feature_space_id,
                        min_observed_fraction=.8, traits=None,
                        projection_method="mean_impute", min_loading_fraction=.8,
                        reference_scores=None,reference_space_id=None):
    """Project independent samples; missing PSI uses discovery means only."""
    if model.get("schema") not in SCHEMAS or model.get("sha256") != _digest(model):
        raise ValueError("Unsupported or modified frozen PSI model")
    if projection_method not in {"mean_impute", "observed_loading_ls",'guarded_loading'}:
        raise ValueError("Unknown projection_method")
    if not np.isfinite(min_loading_fraction) or not 0 < min_loading_fraction <= 1:
        raise ValueError("min_loading_fraction must be within (0,1]")
    ordered = _samples(matrix, samples)
    if feature_space_id != model["feature_space_id"]:
        raise ValueError("feature_space_id differs: no coordinate/quantification conversion is performed")
    if set(matrix.columns) & set(model["training_samples"]) or set(ordered.biological_sample_id.astype(str)) & set(model["training_biological_ids"]):
        raise ValueError("Validation overlaps discovery sample/biological IDs")
    if not np.isfinite(min_observed_fraction) or not 0 < min_observed_fraction <= 1:
        raise ValueError("min_observed_fraction must be within (0,1]")
    scores = pd.DataFrame(index=matrix.columns)
    scores.index.name = "sample_id"
    audits, coverage, summaries = [], [], []
    for entry in model["modules"]:
        module, keys = entry["module"], entry["events"]
        target = matrix.reindex(keys)
        observed = target.notna().mean(axis=0)
        filled = target.T.fillna(pd.Series(entry["imputation_means"], index=keys)).T
        z = (filled.to_numpy(float)-np.array(entry["centers"])[:, None])/np.array(entry["scales"])[:, None]
        loading = np.array(entry["loading"])
        information = target.notna().to_numpy().T @ (loading**2) / np.sum(loading**2)
        eligible = observed.ge(min_observed_fraction)
        score = z.T @ loading/entry["singular_value"]
        route=np.full(len(matrix.columns),'mean_impute',dtype=object)
        if projection_method == "observed_loading_ls":
            # Missing standardized values have no contribution; solve over measured entries.
            z = np.where(target.notna().to_numpy(), z, 0.)
            score = np.divide(z.T @ loading, entry["singular_value"] * information,
                              out=np.full(len(matrix.columns), np.nan), where=information > 1e-12)
            eligible &= information >= min_loading_fraction - 1e-12
            route[:] = 'observed_loading_ls'
        elif projection_method=='guarded_loading':
            mean_score=score.copy()
            z=np.where(target.notna().to_numpy(),z,0.)
            ls=np.divide(z.T@loading,entry['singular_value']*information,
                         out=mean_score.copy(),where=information>1e-12)
            weight=entry.get('guarded_projection',{}).get('weight',0.)
            use=eligible.to_numpy()&(information>=min_loading_fraction-1e-12)
            score=np.where(use,mean_score+weight*(ls-mean_score),mean_score)
            route=np.where(use& (weight>0),'training_calibrated_blend',
                           np.where(use,'training_selected_mean','mean_fallback_loading_information'))
        scores[module] = np.where(eligible, score, np.nan)
        for key in keys:
            present = key in matrix.index
            count = int(target.loc[key].notna().sum())
            audits.append(dict(module=module, event_key=key, present=present, observed_samples=count,
                projection_status="absent_event_imputed_with_discovery_mean" if not present else
                                  "all_missing_imputed_with_discovery_mean" if count == 0 else "measured"))
        for i, sample in enumerate(matrix.columns):
            coverage.append(dict(module=module, sample_id=sample, observed_fraction=float(observed[sample]),
                required_fraction=min_observed_fraction, score_status="scored" if eligible[sample] else
                    "insufficient_measured_events" if observed[sample] < min_observed_fraction else "insufficient_loading_information",
                observed_loading_fraction=float(information[i]), projection_method=projection_method,
                required_loading_fraction=min_loading_fraction if projection_method in {'observed_loading_ls','guarded_loading'} else None,
                scoring_route=str(route[i])))
        # Coherence uses measured PSI only, and only samples eligible for scoring.
        coherence, pairs = _mean_pairwise(target.loc[:, eligible])
        summaries.append(dict(module=module, training_events=len(keys),
            present_events=int(target.index.isin(matrix.index).sum()),
            scored_samples=int(eligible.sum()), validation_samples=len(matrix.columns),
            mean_observed_fraction=float(observed.mean()),
            training_mean_pairwise_r=entry["training_mean_pairwise_r"],
            validation_mean_pairwise_r=coherence, validation_valid_pairs=pairs,
            coherence_status="measured_Pearson_diagnostic" if pairs else "insufficient_measurements",
            regulation_status="not_tested", preservation_significance="not_computed"))
    predictions, performance = [], []
    if traits is not None:
        if not traits.index.is_unique or set(traits.index) != set(matrix.columns):
            raise ValueError("Validation trait rows must match validation samples exactly")
    multi = model.get("trait_prediction_method") == "multimodule_ridge"
    robust=model.get('trait_prediction_method')=='validated_hierarchy'
    if robust and reference_scores is not None and reference_space_id!=model.get('reference_space_id'):
        raise ValueError('reference_space_id differs from the frozen predictor reference')
    entries = model['reliable_trait_models'] if robust else model["multimodule_trait_models"] if multi else model["trait_models"]
    for entry in entries:
        audit=None
        if robust:
            from katana_next.science.psi_robust import predict_reliable
            pred,audit=predict_reliable(entry,scores,reference_scores)
            module=';'.join(entry['module_columns']) or None
        elif multi:
            from katana_next.science.psi_prediction import predict_multimodule
            module = ";".join(entry["modules"]) if entry["alpha"] is not None else None
            pred = predict_multimodule(entry, scores)
        else:
            module = entry["module"]
            pred = (entry["intercept"]+entry["slope"]*scores[module].to_numpy()) if module else np.full(len(matrix.columns), entry["training_mean"])
        y = np.full(len(pred), np.nan)
        if traits is not None and entry["trait"] in traits:
            y = pd.to_numeric(traits.loc[matrix.columns, entry["trait"]], errors="raise").to_numpy(float)
            if np.isinf(y).any():
                raise ValueError("Traits cannot contain infinity")
        for i, sample in enumerate(matrix.columns):
            record=dict(sample_id=sample, trait=entry["trait"], selected_module=module,
                prediction=pred[i], observation=y[i], training_mean=entry["training_mean"])
            if audit is not None:record.update(audit.loc[sample].to_dict())
            predictions.append(record)
        mask = np.isfinite(pred) & np.isfinite(y)
        base = float(np.sum((y[mask]-entry["training_mean"])**2))
        err = float(np.sum((y[mask]-pred[mask])**2))
        performance.append(dict(trait=entry["trait"], matched_samples=int(mask.sum()),
            total_samples=len(pred), prediction_coverage=float(np.isfinite(pred).mean()),
            relative_MSE_improvement=1-err/base if base > 0 else np.nan,
            sum_squared_error=err if mask.any() else np.nan,
            baseline_sum_squared_error=base if mask.any() else np.nan,
            interpretation="heldout_prediction_not_regulatory_accuracy"))
    return dict(scores=scores, event_audit=pd.DataFrame(audits, columns=["module", "event_key", "present", "observed_samples", "projection_status"]),
                sample_coverage=pd.DataFrame(coverage, columns=["module", "sample_id", "observed_fraction", "required_fraction", "score_status", "observed_loading_fraction", "projection_method", "required_loading_fraction","scoring_route"]),
                module_validation=pd.DataFrame(summaries, columns=["module", "training_events", "present_events", "scored_samples", "validation_samples", "mean_observed_fraction", "training_mean_pairwise_r", "validation_mean_pairwise_r", "validation_valid_pairs", "coherence_status", "regulation_status", "preservation_significance"]),
                predictions=pd.DataFrame(predictions, columns=["sample_id", "trait", "selected_module", "prediction", "observation", "training_mean"]+(['prediction_route','uses_module','uses_reference','primary_available','available_module_count'] if robust else [])),
                trait_validation=pd.DataFrame(performance, columns=["trait", "matched_samples", "total_samples", "prediction_coverage", "relative_MSE_improvement", "sum_squared_error", "baseline_sum_squared_error", "interpretation"]))
