"""Recover completed ORA display inputs from a portable report, never refit ORA."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path, PurePosixPath

from .sources import digest, read_table, _artifact_comparisons


def portable_enrichment_artifacts(manifest, provenance, root):
    root = Path(root).resolve()
    recorded = {r.get('artifact_id'): r for source in provenance.get('source_projects', [])
                for r in source.get('parameters', []) if r.get('analysis_id') == 'prioritize.enrichment'}
    hashes = {r['report_path']: r['source_sha256'] for r in manifest.get('source_files', [])
              if r.get('report_path') and r.get('source_sha256')}
    names = {'Gene-set enrichment': 'gene_set', 'Enrichment gene inputs by comparison': 'gene_input_audit',
             'Enrichment background eligibility': 'background_audit'}
    artifacts = []
    checked = {}
    for table in manifest.get('complete_results', []):
        if table.get('analysis_id') != 'prioritize.enrichment':
            continue
        saved = deepcopy(table.get('completed_enrichment') or {})
        entity = saved.get('entity_type') or names.get(table.get('name'))
        if entity not in names.values():
            continue
        previous = recorded.get(table.get('artifact_id'), {})
        meta = saved.get('metadata') or deepcopy(previous.get('method_metadata') or {})
        if meta.get('report_completion_evidence') != 'completed_workflow_manifest':
            continue
        relative = table['source']
        path = (root / relative).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            raise ValueError('Completed enrichment input is outside the portable report or missing')
        if relative not in checked: checked[relative] = digest(path)
        if not hashes.get(relative) or checked[relative] != hashes[relative]:
            raise ValueError('Portable enrichment input differs from its recorded result')
        parameters = saved.get('parameters') or deepcopy(previous.get('parameters') or {})
        # Older manifests omitted run IDs on tables. Their identical recorded
        # parameters and lineage provide a display grouping key, not a new fit.
        export_roots = sorted({str(PurePosixPath(str(p).replace('\\', '/')).parent)
                               for p in (meta.get('exports') or {}).values() if p})
        run = saved.get('run_fingerprint') or hashlib.sha256(json.dumps(
            [parameters, meta.get('native_background_provenance'), export_roots],
            sort_keys=True).encode()).hexdigest()
        artifact = dict(analysis_id='prioritize.enrichment', kind='table', entity_type=entity,
                        artifact_id=table.get('artifact_id'), name=table['name'], path=str(path),
                        metadata=meta, parameters=parameters, run_fingerprint=run)
        frame = read_table(path)
        ids = sorted(set(frame.comparison_id.dropna().astype(str))) if 'comparison_id' in frame else []
        declared = saved.get('contrast_ids') or ([table['contrast_id']] if table.get('contrast_id') else [])
        if declared and set(ids) - set(declared):
            raise ValueError('Portable enrichment comparison differs from recorded identity')
        artifact['contrast_ids'] = ids or declared or list(_artifact_comparisons(artifact))
        artifacts.append(artifact)
    return artifacts
