"""Display completed, explicitly scoped native ORA without refitting statistics."""
from __future__ import annotations


def native_enrichment_scope(artifact):
    """Do not substitute unbound historical plots for a current DAS selection."""
    if artifact.get('analysis_id') != 'prioritize.enrichment':
        return None
    meta = artifact.get('metadata', {})
    comparisons = artifact.get('contrast_ids') or []
    foreground = meta.get('foreground_source')
    if (len(comparisons) != 1 or foreground not in {'deg', 'das', 'dtu', 'dexseq'}
            or not meta.get('native_background_provenance')
            or meta.get('report_completion_evidence') != 'completed_workflow_manifest'):
        return None
    selection = meta.get('das_selection') or {}
    policy = None
    if foreground == 'das':
        policy = selection.get('event_set')
        if policy not in {'significant', 'stable_significant'}:
            return None
        if policy == 'stable_significant':
            from katana_next.science.bulk import READ_SUPPORT_POLICY_VERSION
            if selection.get('read_support_policy') != READ_SUPPORT_POLICY_VERSION:
                return None
    return {'contrast_id': comparisons[0], 'foreground': foreground,
            'policy': policy, 'background': meta.get('background_source', foreground)}
