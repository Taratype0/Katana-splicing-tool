"""Import completed portable evidence into a separate browsing workspace."""
from pathlib import Path, PureWindowsPath
import hashlib
import json
import pandas as pd

from katana_next.domain.models import ProjectManifest, AnalysisArtifact
from katana_next.services.project_service import ProjectService


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def member(root, name):
    p = Path(name)
    if p.is_absolute() or PureWindowsPath(name).drive or '..' in p.parts or '\\' in name:
        raise ValueError('Report member must be a bounded relative path')
    current = root
    for part in p.parts:
        current = current / part
        if current.is_symlink():
            raise ValueError('Report symlinks are not supported')
    if not current.resolve().is_relative_to(root) or not current.is_file():
        raise ValueError(f'Missing report member: {name}')
    return current


def import_report(report_path, output_root, projects=None):
    source = Path(report_path).expanduser().resolve()
    if source.is_file():
        source = source.parent
    destination = Path(output_root).expanduser().resolve()
    if destination == source or destination.is_relative_to(source) or source.is_relative_to(destination):
        raise ValueError('Choose a separate empty workspace outside the report')
    if destination.exists() and any(destination.iterdir()):
        raise ValueError('Workspace must be new or empty')
    manifest_path = source / 'report-manifest.json'
    manifest = json.loads(manifest_path.read_text('utf-8-sig'))
    if manifest.get('schema') != 'katana.approved_report/v1':
        raise ValueError('Unsupported report schema')
    display = manifest.get('display_project')
    if not display or display.get('project_id') != manifest.get('project_id'):
        raise ValueError('Report needs recorded sample/design metadata with matching project identity')
    ledger = {}
    for row in manifest.get('source_files', []):
        if row.get('report_path'):
            ledger[row['report_path']] = row.get('portable_sha256') or row.get('source_sha256')
    for panel in manifest.get('panels', []):
        for link in panel.get('data_links', []):
            if link.get('path'):
                ledger[link['path']] = link.get('portable_sha256') or link.get('source_sha256')
    tables, seen = [], set()
    # Complete tables only; figure previews and Top 10 never define the universe.
    for entry in manifest.get('complete_results', []):
        relative = entry.get('source', '')
        if relative in seen:
            continue
        seen.add(relative)
        if not relative.lower().endswith(('.tsv', '.tsv.gz', '.csv', '.csv.gz', '.parquet')):
            continue
        path = member(source, relative)
        digest = sha(path)
        expected = entry.get('metadata', {}).get('portable_sha256') or ledger.get(relative)
        if expected and expected != digest:
            raise ValueError(f'Report checksum mismatch: {relative}')
        tables.append((entry, path, digest, bool(expected)))
    if not tables:
        raise ValueError('No complete machine-readable result tables in this report')
    contrasts = [{**c, 'label': c.get('label', c['contrast_id'])} for c in display.get('contrasts', [])]
    samples = [{k: v for k, v in s.items() if k in {'sample_id','biological_sample_id','condition','batch','attributes'}}
               for s in display.get('samples', [])]
    origin = sha(manifest_path)
    project = ProjectManifest(project_id=display['project_id'], name=display.get('name', display['project_id']),
        input_root=destination, output_root=destination / 'results', samples=samples, contrasts=contrasts,
        annotation={'genome_build': display.get('annotation', {}).get('genome_build', 'unknown')},
        metadata={'report_workspace': True, 'source_report_sha256': origin,
                  'source_report_directory': str(source), 'assay_mode': manifest.get('mode', 'bulk'),
                  'coordinate_system': display.get('metadata', {}).get('coordinate_system', '0-based-half-open'),
                  'coordinate_contract': 'Canonical rMATS *_0base/ES/EE fields only; other geometry remains unsupported'})
    projects = projects or ProjectService()
    projects._validate_output_root(destination)
    destination.mkdir(parents=True, exist_ok=True)
    local_manifest = destination / 'katana.project.json'
    local_manifest.write_text(project.model_dump_json(indent=2), encoding='utf-8')
    context = projects.open_project(local_manifest)
    audit = []
    for entry, path, digest, verified in tables:
        if path.suffix == '.parquet':
            frame = pd.read_parquet(path)
        else:
            # Preserve identifier strings, including leading zeroes and literal NA.
            sep = ',' if '.csv' in path.suffixes else '\t'
            header = pd.read_csv(path, sep=sep, nrows=0).columns
            identifiers = {c: str for c in header if c in {'gene_id','gene_symbol','GeneID','geneSymbol','event_id','event_uid','row_id','ID','comparison_id','contrast_id','sample_id'}}
            frame = pd.read_csv(path, sep=sep, dtype=identifiers, keep_default_na=False)
        for old, new in [('GeneID','gene_id'), ('geneSymbol','gene_symbol')]:
            if old in frame and new not in frame:
                frame[new] = frame[old].astype(str).str.strip('"')
        cid = entry.get('contrast_id')
        if cid and cid != 'all' and 'comparison_id' not in frame:
            frame['comparison_id'] = cid
        if 'row_id' in frame:
            frame['source_row_id'] = frame['row_id']
        frame['row_id'] = [f'{digest[:16]}:{i}' for i in range(len(frame))]
        target = destination / 'results' / 'imported' / f'{digest}.parquet'
        context.results.write_dataframe(frame, target)
        source_metadata = entry.get('metadata', {})
        metadata = {'source_metadata': source_metadata,
            **{k:v for k,v in source_metadata.items() if k in {'count_mode','rmats_mode','selected_rmats_mode','das_selection'}},
            'project_snapshot_fingerprint': context.snapshot.manifest_fingerprint,
            'source_analysis_id': entry.get('analysis_id'), 'source_report_sha256': origin,
            'source_relative_path': path.relative_to(source).as_posix(), 'source_table_sha256': digest,
            'imported_table_sha256': sha(target), 'report_hash_verified': verified,
            'completion_verification': 'report_recorded_result', 'chapter': entry.get('chapter'),
            'statistical_scope': 'Original completed result; browsing never changes its test family'}
        artifact = AnalysisArtifact(project_id=project.project_id, analysis_id='imported.report',
            name=entry.get('name', path.name), kind='table', entity_type='report_result',
            fingerprint=digest, run_fingerprint=origin, path=target, row_count=len(frame),
            contrast_ids=(cid,) if cid and cid != 'all' else (),
            parameters=entry.get('parameters', {}), metadata=metadata)
        context.catalog.register_artifact(artifact)
        audit.append({'artifact_id': artifact.artifact_id, 'rows': len(frame), **metadata})
    (destination / 'report-import.json').write_text(json.dumps({'source_manifest_sha256': origin,
        'tables': audit, 'upstream_runs': 0}, indent=2), encoding='utf-8')
    return context
