"""Reconcile deferred report statuses with complete, explicitly scoped evidence.

This runs after all views are collected. Numerical values and producer files
are never changed; a non-significant completed result is still completed.
"""
from copy import deepcopy
from pathlib import Path, PurePosixPath, PureWindowsPath
import re

from .sources import digest


# Identity alone is insufficient: a producer must also declare the component
# or a completion scope appropriate to that particular evidence layer.
_RULES = {
    'clip': ({'prioritize.binding_enrichment', 'prioritize.rbp_mechanism'},
             {'clip', 'reference_clip', 'binding_enrichment'},
             {'clip', 'reference_clip', 'binding_enrichment', 'full_binding_enrichment'}),
    'rbp_identity': ({'prioritize.rbp_support'}, {'rbp_identity'}, {'rbp_identity'}),
    'motif_enrichment': ({'prioritize.event_motifs', 'prioritize.se_rbp_impact',
                          'prioritize.rbp_sequence_results', 'prioritize.rbp_mechanism',
                          'prioritize.motiftrace_consequences'},
                         {'motif_enrichment'}, {'motif_enrichment', 'full_motif_enrichment', 'full_boundary_association'}),
    'motiftrace_sequence': ({'prioritize.rbp_sequence_results', 'prioritize.motiftrace_consequences',
                            'prioritize.motiftrace_usage'},
                           {'motiftrace_sequence'}, {'motiftrace_sequence', 'sequence_consequences_only'}),
    'mapp_parallel': ({'prioritize.motiftrace_parallel'}, {'mapp_parallel', 'calibrated_parallel'},
                      {'calibrated_parallel', 'completed_calibrated_parallel'}),
    'network_bootstrap': ({'programs.cosplice', 'flow.network.bootstrap'}, {'network_bootstrap'},
                          {'bootstrap_confidence', 'full_network_with_bootstrap'}),
}


def _expected_contrasts(builder, state):
    cid = state.get('contrast_id', 'all')
    if cid and cid != 'all':
        return [str(cid)]
    project = builder.source.get('project', {})
    values = builder.source.get('report_selection', {}).get('contrast_ids')
    if values is None:
        values = project.get('report_contrast_ids')
    if values is None:
        values = [row['contrast_id'] for row in project.get('contrasts', [])]
    return list(dict.fromkeys(map(str, values or [])))


def _is_component(panel, family):
    meta = panel.get('metadata', {})
    workflows, markers, scopes = _RULES[family]
    if meta.get('analysis_id') not in workflows:
        return False
    if meta.get('completion_scope') is not None and meta['completion_scope'] not in scopes:
        return False
    if family == 'motif_enrichment' and meta.get('completion_scope') == 'sequence_consequences_only':
        return False
    if family == 'mapp_parallel' and meta.get('calibrated_motiftrace_parallel') is not True:
        return False
    if family == 'network_bootstrap' and (meta.get('bootstrap_enabled') is False
                                          or meta.get('completion_scope') not in scopes):
        return False
    return (meta.get('component') in markers or meta.get('report_module') in markers
            or meta.get('completion_scope') in scopes)


def _completed_metadata(builder, panel):
    meta = panel.get('metadata', {})
    if meta.get('completion_verified') is False:
        return False
    if (meta.get('completion_verified') is not True and
            meta.get('completion_evidence') not in {'completed_workflow_manifest', 'verified_collection_receipt',
                                                     'native_terminal_completion'}):
        return False
    if panel.get('status') not in (None, 'available', 'complete', 'completed'):
        return False
    binding = meta.get('report_source_binding')
    if binding is not None:
        project = builder.source.get('project', {})
        if (binding.get('schema') != 'katana.report_source_binding/v1'
                or binding.get('biological_design_identical') is not True
                or binding.get('source_files_modified') is not False
                or binding.get('scientific_computation_performed') is not False
                or binding.get('target_project_id') != project.get('project_id')):
            return False
        target = builder.source.get('native_project_manifest_sha256')
        if target and binding.get('target_manifest_sha256') != target:
            return False
    if meta.get('project_id') and meta['project_id'] != builder.source.get('project', {}).get('project_id'):
        return False
    return True


def _native_links(builder, panel, family, cache):
    """Return exact bound data identities, or None without clearing a status."""
    required = {'full_results', 'methods'}
    if family == 'clip':
        required |= {'foreground', 'background'}
    cid, cell = str(panel.get('contrast_id', 'all')), str(panel.get('cell_type', 'all'))
    records, covered = [], set()
    root = Path(builder.output).resolve()
    for link in panel.get('data_links', []):
        role = link.get('role')
        if role not in required:
            continue
        if str(link.get('contrast_id', 'all')) != cid or str(link.get('cell_type', cell)) != cell:
            return None
        relative = link.get('path')
        identity = link.get('portable_sha256', link.get('source_sha256'))
        if (not isinstance(relative, str) or not relative or '\\' in relative
                or PurePosixPath(relative).is_absolute() or PureWindowsPath(relative).drive
                or '..' in PurePosixPath(relative).parts
                or not isinstance(identity, str) or not re.fullmatch('[0-9a-f]{64}', identity)):
            return None
        path = (root / relative).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            return None
        stat = path.stat()
        key = (str(path), stat.st_size, stat.st_mtime_ns)
        if key not in cache:
            cache[key] = digest(path)
        if cache[key] != identity:
            return None
        covered.add(role)
        records.append({'role': role, 'path': relative, 'sha256': identity})
    return records if required <= covered else None


def _input_supplied(builder, family):
    workflows = _RULES[family][0]
    records = [*builder.panels, *builder.tables, *builder.source.get('artifacts', []),
               *builder.source.get('extra_assets', [])]
    return any(row.get('metadata', {}).get('analysis_id', row.get('analysis_id')) in workflows for row in records)


def reconcile_component_statuses(builder):
    """Resolve only fully proven obsolete statuses; retain their exact history.

    Mutates only the new report's ``exclusions`` and ``derived_audits``. Call
    from presentation preparation after every native and supplemental view.
    Aggregate panels are never assumed to cover individual comparison/cell
    scopes, and partial scopes never erase a project-wide required gap.
    """
    kept, decisions, cache = [], [], {}
    for original in builder.exclusions:
        state = deepcopy(original)
        meta = state.get('metadata', {})
        family = meta.get('component')
        if state.get('status') != 'not_computed' or family not in _RULES:
            kept.append(original)
            continue
        expected = _expected_contrasts(builder, state)
        cell = str(state.get('cell_type', 'all'))
        evidence = []
        for panel in builder.panels:
            if (str(panel.get('contrast_id', 'all')) not in expected
                    or str(panel.get('cell_type', 'all')) != cell
                    or not _is_component(panel, family) or not _completed_metadata(builder, panel)):
                continue
            links = _native_links(builder, panel, family, cache)
            if links is not None:
                evidence.append({'panel_id': panel.get('id'), 'contrast_id': panel['contrast_id'],
                                 'cell_type': cell, 'analysis_id': panel['metadata']['analysis_id'],
                                 'completion_scope': panel['metadata'].get('completion_scope'),
                                 'bound_native_files': links})
        complete = {row['contrast_id'] for row in evidence}
        missing = [cid for cid in expected if cid not in complete]
        decision = {'original_status': state, 'component': family, 'required_contrast_ids': expected,
                    'cell_type': cell, 'verified_evidence': evidence, 'missing_contrast_ids': missing,
                    'original_native_statistics_preserved': True}
        if expected and not missing:
            decision.update(action='resolved_completed', status='completed',
                            significance_used_as_completion_gate=False)
        elif (meta.get('requirement') == 'conditional_import'
                and meta.get('optional_source_supplied') is False and not _input_supplied(builder, family)):
            # Absence stays a fact, rather than being promoted to completion.
            decision.update(action='optional_absence_retained_in_provenance', status='not_computed',
                            optional_source_supplied=False, main_unavailable_display=False)
        else:
            kept.append(original)
            decision.update(action='retained_unresolved_required_scope', status='not_computed')
        decisions.append(decision)
    audit = {'schema': 'katana.completed_component_status_reconciliation/v1', 'decisions': decisions,
             'scientific_computation_performed': False, 'source_files_modified': False,
             'native_statistics_modified': False}
    builder.exclusions = kept
    if decisions:
        builder.derived_audits.append(audit)
    return audit
