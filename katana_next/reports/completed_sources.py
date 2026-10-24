"""Assemble explicitly declared completed result sources without scientific jobs."""
from copy import deepcopy
from pathlib import Path
import json

SCHEMA = 'katana.completed_result_set/v1'


def assemble_result_set(value, *, base, source_files, mapper):
    from .sources import _load_source, digest, _native_paths, report_contrast_selection
    from .path_mapping import identity_projection
    location = mapper.resolve(value['project_manifest'], base)
    raw = json.loads(location.read_text('utf-8-sig'))
    parent_sha = digest(location)
    parent = _native_paths(raw, location.parent, mapper=mapper)
    selection = report_contrast_selection(parent, value)
    wanted = set(selection['contrast_ids'])
    rows = value.get('sources', [])
    if not rows or len({s['source_id'] for s in rows}) != len(rows):
        raise ValueError('Completed result set needs uniquely named explicit sources')
    all_artifacts, extras, views, audits, scopes, seen = [], [], [], [], {}, {}
    files = [*source_files, location]
    for row in rows:
        source_path = mapper.resolve(row['path'], base)
        kind = row.get('kind', 'descriptor')
        if kind == 'collection':
            from katana_next.flow.collect import materialize_report_descriptor
            child_descriptor = materialize_report_descriptor(source_path)
            child = _load_source(child_descriptor, inline_base=source_path, mapper=mapper)
            receipt_path = source_path / 'collection.receipt.json'
            if source_path.is_file(): receipt_path = source_path.parent / 'collection.receipt.json'
            receipt = json.loads(receipt_path.read_text('utf-8-sig'))
            files.append(receipt_path)
            completion = {'receipt_sha256':receipt['receipt_sha256'], 'task_count':len(receipt['tasks']),
                          'status':'completed', 'receipt_file_sha256':digest(receipt_path)}
        elif kind == 'descriptor':
            child = _load_source(source_path, inline_base=source_path.parent, mapper=mapper)
            completion = {'status':'completed_artifact_descriptor'}
        elif kind == 'views':
            # Authenticate a saved view against its original project revision
            # before the shared design/target-SHA gates below permit rebinding.
            # Never replace its native identity with the target just to load it.
            view_manifest = mapper.resolve(row['project_manifest'], base) if row.get('project_manifest') else location
            child = _load_source({'project_manifest':str(view_manifest),'artifacts':[],
                                  'extra_views':[str(source_path)]},inline_base=base,mapper=mapper)
            completion = {'status':'completed_view_descriptor'}
        else:
            raise ValueError('Supported result sources are descriptor, views and verified collection')
        bound_evidence=[]
        for evidence in row.get('completion_evidence', []):
            path=mapper.resolve(evidence['path'],base)
            if digest(path)!=evidence['sha256']: raise ValueError('Completed-source evidence SHA differs')
            proof=json.loads(path.read_text('utf-8-sig'))
            if not proof.get('schema') or not (proof.get('passed') is True or proof.get('status') in ('completed','verified','passed')):
                raise ValueError('Completion evidence does not declare verified/completed results')
            expected=proof.get('target_report_project_manifest_sha256',proof.get('target_project_manifest_sha256',proof.get('report_manifest_sha256')))
            if expected and expected!=parent_sha: raise ValueError('Completion evidence binds another target study')
            bound_evidence.append({'path':str(path),'sha256':evidence['sha256'],'schema':proof['schema'],
                                   'evidence_kind':'caller_bound_completed_outputs'})
            files.append(path)
        if bound_evidence: completion['bound_completion_evidence']=bound_evidence
        original = child.get('native_project', child['project'])
        if identity_projection(original) != identity_projection(parent):
            raise ValueError('Completed result source has different biological samples, contrast arms or reference: '+row['source_id'])
        source_sha = child['native_project_manifest_sha256']
        if source_sha != parent_sha and row.get('target_project_manifest_sha256') != parent_sha:
            raise ValueError('A different project revision requires an explicit exact target manifest SHA binding')
        if original.get('project_id') != parent.get('project_id') and not row.get('allow_project_rebind', False):
            raise ValueError('Different project IDs require explicit presentation-only rebinding')
        binding = {'schema':'katana.report_source_binding/v1', 'source_id':row['source_id'],
                   'native_project_id':original.get('project_id'), 'native_manifest_sha256':source_sha,
                   'target_project_id':parent['project_id'], 'target_manifest_sha256':parent_sha,
                   'biological_design_identical':True, 'source_files_modified':False,
                   'scientific_computation_performed':False, **completion}
        excluded = set(row.get('exclude_analysis_ids', []))
        ids = set()
        def include(record):
            return record.get('analysis_id', record.get('metadata', {}).get('analysis_id')) not in excluded
        for extra, records in ((False,child['artifacts']), (True,child.get('extra_assets', []))):
            for original_record in records:
                if not include(original_record): continue
                record = deepcopy(original_record)
                ids.update(record.get('contrast_ids') or [])
                record['project_id'] = parent['project_id']
                meta = record.setdefault('metadata', {})
                meta['report_source_binding'] = binding
                if extra: meta['source_project_manifest_sha256'] = parent_sha
                key = record.get('artifact_id')
                if key:
                    signature = (record.get('analysis_id'),str(record['path']))
                    if key in seen:
                        if seen[key] != signature: raise ValueError('Conflicting completed artifact identities: '+str(key))
                        continue
                    seen[key] = signature
                (extras if extra else all_artifacts).append(record)
        for native_view in child.get('extra_views', []):
            view = deepcopy(native_view)
            view['panels'] = [p for p in view.get('panels', []) if include(p)]
            view['tables'] = [t for t in view.get('tables', []) if include(t)]
            ids.update(str(p['contrast_id']) for p in view['panels'] if p.get('contrast_id') not in (None,'all'))
            view.update(project_id=parent['project_id'],source_project_manifest_sha256=parent_sha)
            view['_report_source_id'] = row['source_id']
            for panel in view['panels']: panel.setdefault('metadata', {})['report_source_binding'] = binding
            if row.get('source_collection_root'):
                root=str(mapper.resolve(row['source_collection_root'],base))
                for panel in view['panels']:
                    if panel.get('metadata',{}).get('analysis_id')=='prioritize.mapp_native':
                        for link in panel.get('data_links',[]):
                            if link.get('role')=='full_results': link.setdefault('source_collection_root',root)
            views.append(view)
        scopes[row['source_id']] = {'verified_collection':kind=='collection' or bool(bound_evidence), 'contrasts':ids}
        audits.append(binding)
        files.extend(child.get('source_files', []))
    superseded = []
    for rule in value.get('superseded_statuses', []):
        replacements = rule.get('replacement_source_ids', [])
        if not replacements or any(s not in scopes or not scopes[s]['verified_collection'] for s in replacements):
            raise ValueError('A superseded status requires verified completed replacement collections')
        covered = set().union(*(scopes[s]['contrasts'] for s in replacements))
        required = set(rule.get('contrast_ids', selection['contrast_ids']))
        if not required <= covered: raise ValueError('Replacement collections do not cover the superseded status scope')
        matches = 0
        for view in views:
            if view['_report_source_id'] != rule['source_id']: continue
            retained = []
            for state in view.get('unavailable', []):
                if state.get('title') == rule['title']:
                    superseded.append({'original_status':deepcopy(state),'replacement_source_ids':replacements,
                                       'contrast_ids':sorted(required),'source_history_preserved':True})
                    matches += 1
                else: retained.append(state)
            view['unavailable'] = retained
        if matches != 1: raise ValueError('Superseded status must identify exactly one recorded old status')
    for view in views: view.pop('_report_source_id', None)
    display = deepcopy(parent)
    by_id = {c['contrast_id']:c for c in parent.get('contrasts', [])}
    display['contrasts'] = [by_id[c] for c in selection['contrast_ids']]
    display['report_contrast_ids'] = selection['contrast_ids']
    display['report_contrast_groups'] = selection['groups']
    return {'project':display,'native_project':parent,'native_project_manifest_sha256':parent_sha,
            'report_selection':selection,'artifacts':all_artifacts,'extra_assets':extras,'extra_views':views,
            'supplemental':[],'source_files':list(dict.fromkeys(files)), 'search_roots':[],
            'dummy':bool(value.get('dummy',False)), 'mode':value.get('mode','bulk'),
            'report_native_figures_in_main':bool(value.get('report_native_figures_in_main',False)),
            'report_require_plot_data':bool(value.get('report_require_plot_data',True)),
            'assembly_audit':{'schema':'katana.completed_report_source_assembly/v1','sources':audits,
                              'superseded_status_history':superseded,'scientific_jobs_submitted':0}}
