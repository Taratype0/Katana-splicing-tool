"""Full original-module validation inputs, with no core trimming or filling."""
from dataclasses import dataclass

import numpy as np
import pandas as pd

from katana_next.science.frozen_event_transport import full_frozen_member_gate


@dataclass
class FullFamilyInput:
    reference: pd.DataFrame
    target: pd.DataFrame
    labels: pd.DataFrame
    members: list[str]


def strict_frozen_scores(matrix, loadings):
    """Project fixed A weights only when every original member is measured."""
    if not {"event_key", "center", "scale", "loading"} <= set(loadings):
        raise ValueError("Frozen event centers, scales and weights required")
    if loadings.event_key.duplicated().any() or loadings.event_key.isna().any():
        raise ValueError("Frozen score membership must be unique and explicit")
    if matrix.index.has_duplicates or matrix.columns.has_duplicates:
        raise ValueError("Duplicate event/sample score identity")
    model = loadings[["center", "scale", "loading"]].to_numpy(float)
    if not np.isfinite(model).all() or (model[:, 1] <= 0).any():
        raise ValueError("Frozen model must be finite with positive scales")
    values = matrix.reindex(loadings.event_key).to_numpy(float)
    if np.any(np.isinf(values)) or np.any(np.isfinite(values) & ((values < 0) | (values > 1))):
        raise ValueError("Measured score PSI outside [0, 1] or infinite")
    complete = np.isfinite(values).all(axis=0)
    scores = np.full(values.shape[1], np.nan)
    scores[complete] = model[:, 2] @ ((values[:, complete] - model[:, 0, None]) / model[:, 1, None])
    return pd.Series(scores, index=matrix.columns, name="frozen_score")


def prepare_full_frozen_families(reference, target, membership, source_samples,
                                target_samples, exact_identity, *,
                                minimum_patients=15, minimum_members=20, variance_floor=1e-8):
    """Freeze the estimable family on measurements, before any P values.

    The primary R invocation has exactly one full target module. Other usable
    original-universe events are background, not smaller modules to be tested.
    """
    if membership.index.has_duplicates or membership.empty or membership.isna().any().any():
        raise ValueError("Original unique explicit memberships required")
    for matrix in (reference, target):
        if matrix.index.has_duplicates or matrix.columns.has_duplicates:
            raise ValueError("Duplicate event/sample matrix identity")
        if set(matrix.index) - set(membership.index):
            raise ValueError("Matrix contains events outside original frozen universe")
        values = matrix.to_numpy(float)
        if np.isinf(values).any() or np.any(np.isfinite(values) & ((values < 0) | (values > 1))):
            raise ValueError("Invalid measured PSI")
    for samples in (source_samples, target_samples):
        if not {"sample_id", "patient_id", "condition"} <= set(samples):
            raise ValueError("Explicit biological identity and scope required")
        if samples[["sample_id", "patient_id", "condition"]].isna().any().any():
            raise ValueError("Missing biological identity or scope")
        if samples.sample_id.duplicated().any() or samples.duplicated(["patient_id", "condition"]).any():
            raise ValueError("Technical replicates must not become independent patients")
        if samples[["sample_id", "patient_id", "condition"]].astype(str).apply(lambda s: s.str.strip().eq("")).any().any():
            raise ValueError("Empty biological identity or scope")
    if set(source_samples.patient_id.astype(str)) & set(target_samples.patient_id.astype(str)):
        raise ValueError("Biological IDs overlap between discovery and external studies")
    scopes = sorted(set(source_samples.condition) | set(target_samples.condition))
    exact = pd.Series({key: isinstance(exact_identity.get(key), (bool, np.bool_))
                       and bool(exact_identity[key]) for key in membership.index})
    attempts, inputs, losses = [], {}, []
    for scope in scopes:
        a_samples = source_samples.loc[source_samples.condition.eq(scope)]
        b_samples = target_samples.loc[target_samples.condition.eq(scope)]
        for method in membership:
            for module in sorted(set(membership[method]) - {"grey", "0"}):
                keys = membership.index[membership[method].eq(module)].tolist()
                gates = [full_frozen_member_gate(matrix, keys, samples, identity,
                    minimum_patients, minimum_members, variance_floor)
                    for matrix, samples, identity in (
                        (reference, a_samples, dict.fromkeys(keys, True)),
                        (target, b_samples, exact_identity))]
                reasons = [f"{study}:{reason}" for study, gate in zip(("A", "B"), gates)
                           for reason in gate["reasons"]]
                row = dict(scope=scope, method=method, module=module, original_members=len(keys),
                    source_complete_patients=gates[0]["complete_case_patients"],
                    target_complete_patients=gates[1]["complete_case_patients"],
                    source_variable_members=gates[0]["variable_members_in_complete_cases"],
                    target_variable_members=gates[1]["variable_members_in_complete_cases"],
                    target_exact_members=gates[1]["exact_identity_members"],
                    full_module_input_eligible=False, preservation_pvalue=None,
                    preservation_fdr=None, imputed_cells=0, dropped_original_members=0)
                if not reasons:
                    chosen = []
                    for matrix, samples, gate in zip((reference, target), (a_samples, b_samples), gates):
                        selected = samples.loc[samples.patient_id.astype(str).isin(gate["complete_case_patient_ids"]), "sample_id"]
                        chosen.append(matrix.reindex(membership.index).loc[:, selected])
                    a, b = chosen
                    complete_a = np.isfinite(a.to_numpy(float)).all(axis=1)
                    complete_b = np.isfinite(b.to_numpy(float)).all(axis=1)
                    variable_a = a.var(axis=1, ddof=1).to_numpy() > variance_floor
                    variable_b = b.var(axis=1, ddof=1).to_numpy() > variance_floor
                    usable = exact.reindex(membership.index).to_numpy() & complete_a & complete_b & variable_a & variable_b
                    background = membership.index[usable]
                    for i, key in enumerate(membership.index):
                        losses.append(dict(scope=scope, method=method, module=module, event_key=key,
                            primary_original_member=key in keys, exact_identity=bool(exact[key]),
                            complete_A=bool(complete_a[i]), complete_B=bool(complete_b[i]),
                            variable_A=bool(variable_a[i]), variable_B=bool(variable_b[i]),
                            retained_background=bool(usable[i])))
                    row["background_events"] = len(background)
                    if not set(keys) <= set(background):
                        raise RuntimeError("Full original target members vanished from background")
                    if len(background) <= len(keys):
                        reasons.append("no_external_permutation_background")
                    else:
                        labels = pd.DataFrame({"event_key": background,
                            "wgcna_module": [module if key in keys else "grey" for key in background],
                            "netrep_module": [module if key in keys else "0" for key in background]})
                        inputs[(scope, method, module)] = FullFamilyInput(a.loc[background], b.loc[background], labels, keys)
                        row["full_module_input_eligible"] = True
                row["reasons"] = ";".join(reasons)
                attempts.append(row)
    return pd.DataFrame(attempts), inputs, pd.DataFrame(losses)
