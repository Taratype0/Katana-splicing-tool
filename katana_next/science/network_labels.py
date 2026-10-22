"""Descriptive network labels are opt-in views, never membership filters.

The fitted partition and continuous evidence remain unchanged. Historical
composite confidence scores include threshold-coded split-half roles and must
not be presented as threshold-independent measurements or calibrated probability.
"""
from __future__ import annotations

import pandas as pd


NETWORK_LABEL_POLICY_ID = 'network-descriptive-labels-opt-in-v1'


def network_label_policy(apply_threshold_labels=False, *, label_thresholds=None):
    if not isinstance(apply_threshold_labels, bool):
        raise ValueError('apply_threshold_labels must be true or false')
    policy = {
        'policy_id': NETWORK_LABEL_POLICY_ID,
        'apply_threshold_labels': apply_threshold_labels,
        'hard_filter': False,
        'default_reporting': 'continuous_metrics_and_full_ranks',
        'fdr_default': .05,
        'label_thresholds': dict(label_thresholds or {}),
        'label_threshold_role': 'descriptive_only' if apply_threshold_labels else 'configured_not_applied',
        'partition': 'full_data_wgcna_backbone',
        'inferential_fdr_default': .05,
        'ranking': 'Complete raw co-assignment, kME and margin rankings; no descriptive-label eligibility filter.',
        'legacy_composite': 'Threshold-coded split-half roles contribute to the historical composite score; '
                            'retained with an explicit legacy name, never a calibrated probability.',
        'literature_basis': {
            'kME_and_connectivity': 'WGCNA defines module eigengene membership and intramodular connectivity: '
                                   'https://doi.org/10.1186/1471-2105-9-559',
            'exact_member_labels': 'The .70/.45/.05 co-assignment thresholds, bootstrap kME cutoffs and '
                                   'split-half confirmation are Katana descriptive rules, not independently calibrated literature standards.',
            'stability_labels': 'Bootstrap/Jaccard recovery is supported by resampling methodology; the Katana '
                                '.75/.50 labels use a composite recovery measure, not a statistical significance test.',
            'inference': 'Native P values, BH FDR families and the .05 default are retained separately from descriptive labels.',
        },
    }
    policy['basis'] = policy['literature_basis']
    return policy


def configured_member_thresholds(consensus, validation=None):
    names = ('core_membership_min', 'peripheral_membership_min', 'peripheral_membership_margin',
             'core_kme_mean_min', 'core_kme_lower_quantile_min', 'core_kme_sign_consistency_min',
             'peripheral_kme_mean_min', 'peripheral_kme_lower_quantile_min', 'peripheral_kme_sign_consistency_min')
    values = {name: getattr(consensus, name) for name in names}
    if validation is not None:
        values.update(robust_recovery_min=validation.robust_jaccard,
                      moderate_recovery_min=validation.moderate_jaccard)
    return values


_LEGACY_COMPOSITES = {
    'continuous_membership_confidence': 'legacy_threshold_dependent_membership_score',
    'continuous_confidence_rank_within_module': 'legacy_threshold_dependent_rank_within_module',
    'continuous_confidence_percentile_within_module': 'legacy_threshold_dependent_percentile_within_module',
    'confidence_split_half_component': 'legacy_threshold_coded_split_half_component',
    'mean_continuous_confidence': 'legacy_threshold_dependent_mean_confidence',
    'median_continuous_confidence': 'legacy_threshold_dependent_median_confidence',
    'min_continuous_confidence': 'legacy_threshold_dependent_min_confidence',
    'mean_continuous_membership_confidence': 'legacy_threshold_dependent_mean_membership_score',
    'median_continuous_membership_confidence': 'legacy_threshold_dependent_median_membership_score',
    'min_continuous_membership_confidence': 'legacy_threshold_dependent_min_membership_score',
}
_MEMBER_LABELS = (
    'membership_role', 'membership_evidence', 'split_half_membership_role_a',
    'split_half_membership_role_b', 'stability_class',
    'observed_residual_program_class',
)
_LABEL_COUNTS = (
    'n_core_events', 'n_peripheral_events', 'n_uncertain_events',
    'n_uncertain_candidate_events', 'n_low_confidence_assigned_events',
    'robust_event_fraction',
)


def member_label_view(frame, apply_threshold_labels=False):
    """Return a nonmutating annotation view with all rows and numeric evidence.

    Opting in exposes the engine's existing descriptive rules. It never refines
    the partition, removes uncertain events, or changes a native test result.
    """
    network_label_policy(apply_threshold_labels)
    result = frame.copy(deep=True)
    recognized = set(_MEMBER_LABELS) | set(_LEGACY_COMPOSITES) | set(_LABEL_COUNTS) | {
        'split_half_core_confirmed', 'consensus_membership_confidence',
    }
    if not recognized.intersection(result.columns):
        return result
    for old, new in _LEGACY_COMPOSITES.items():
        if old in result:
            if new not in result:
                result[new] = result[old]
            result = result.drop(columns=old)
    if apply_threshold_labels:
        result['member_label_status'] = ('applied_descriptive_only'
            if any(column in result for column in _MEMBER_LABELS) else 'not_assessed')
    else:
        for column in _MEMBER_LABELS:
            if column in result:
                result[column] = 'not_applied'
        for column in _LABEL_COUNTS:
            if column in result:
                result[column] = pd.Series(pd.NA, index=result.index, dtype='Float64')
        if 'split_half_core_confirmed' in result:
            result['split_half_core_confirmed'] = pd.Series(pd.NA, index=result.index, dtype='boolean')
        result['member_label_status'] = 'not_applied'
        if 'analysis_assignment_scope' in result:
            module = result.get('consensus_module', result.get('analysis_module'))
            if module is not None:
                result['analysis_assignment_scope'] = module.astype(str).map(
                    lambda name: 'unassigned' if name in {'grey', 'gray', '0'} else 'wgcna_backbone_member')
    if 'consensus_module' in result:
        assigned = ~result.consensus_module.astype(str).isin(['grey', 'gray', '0'])
        for column, rank in (
            ('consensus_membership_confidence', 'coassignment_rank_within_module'),
            ('bootstrap_mean_abs_kME', 'kme_rank_within_module'),
            ('consensus_confidence_margin', 'margin_rank_within_module'),
        ):
            if column in result:
                values = pd.to_numeric(result[column], errors='coerce').where(assigned)
                result[rank] = values.groupby(result.consensus_module, observed=True).rank(
                    ascending=False, method='min')
    result['network_label_policy_id'] = NETWORK_LABEL_POLICY_ID
    return result
