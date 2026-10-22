"""Read project manifests and artifact exports without opening their database."""
from pathlib import Path
import copy
import hashlib, json


def digest(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def load_source(project, catalog=None):
    source_files = []
    if isinstance(project, dict):
        value = dict(project)
        base = Path.cwd()
    else:
        location = Path(project).expanduser().resolve()
        value = json.loads(location.read_text('utf-8-sig'))
        base = location.parent
        source_files.append(location)
    if value.get('project_manifest'):
        manifest = Path(value['project_manifest'])
        if not manifest.is_absolute():
            manifest = base / manifest
        native = json.loads(manifest.read_text('utf-8-sig'))
        source_files.append(manifest)
        base = manifest.parent
    else:
        native = value.get('project', value)
    records = value.get('artifacts')
    catalog = catalog or value.get('artifact_catalog') or value.get('catalog_export') or value.get('native_artifact_catalog')
    if isinstance(records, str):
        catalog, records = records, None
    if records is None:
        candidates = [Path(catalog)] if catalog else [base / 'native-artifact-catalog.json', Path(native.get('output_root', base)) / '.katana/artifact_manifest.json']
        for path in candidates:
            if not path.is_absolute():
                path = base / path
            if path.is_file():
                records = json.loads(path.read_text('utf-8-sig'))
                records = records.get('artifacts', records) if isinstance(records, dict) else records
                source_files.append(path)
                break
    if records is None:
        db = Path(native.get('output_root', base)) / '.katana/catalog.duckdb'
        if not db.is_file():
            raise ValueError('No completed artifact export found; supply catalog= with a JSON artifact export.')
        import duckdb
        with duckdb.connect(str(db), read_only=True) as connection:
            rows = connection.execute('SELECT record_json FROM artifacts WHERE project_id = ? ORDER BY created_at DESC', [native['project_id']]).fetchall()
            records = [json.loads(row[0]) for row in rows]
        source_files.append(db)
    if not isinstance(records, list):
        raise ValueError('Artifact export must contain a list of completed artifact records.')
    artifacts, seen = [], set()
    for raw in records:
        record = copy.deepcopy(raw)
        if record.get('status', 'complete') not in ('complete', 'completed'):
            continue
        if record.get('project_id') and record['project_id'] != native.get('project_id'):
            raise ValueError('Artifact belongs to a different project; use an explicitly verified supplemental source.')
        path = Path(record.get('path', ''))
        if not path.is_absolute():
            path = base / path
        if not path.is_file():
            record['missing_reason'] = 'Recorded completed artifact file is unavailable.'
        record['path'] = str(path.resolve())
        if path.is_file():
            from .rich_results import _completion
            complete, evidence = _completion(record, path)
            if not complete:
                continue
            record.setdefault('metadata', {})['report_completion_evidence'] = evidence
        key = (record.get('analysis_id'), record.get('name'), tuple(record.get('contrast_ids', [])), record.get('parameters', {}).get('sequence_scope'), record.get('parameters', {}).get('event_set'))
        if key in seen:
            continue
        seen.add(key)
        artifacts.append(record)
    native_location = manifest if value.get('project_manifest') else (location if not isinstance(project, dict) and 'project' not in value else None)
    manifest_identity = digest(native_location) if native_location else hashlib.sha256(json.dumps(native, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()
    extras = copy.deepcopy(value.get('extra_assets', []))
    views = []
    for asset in extras:
        meta = asset.get('metadata', {})
        if asset.get('project_id') != native.get('project_id') or meta.get('source_project_manifest_sha256', meta.get('source_project_manifest_SHA256')) != manifest_identity:
            raise ValueError('Extra report assets require matching project identity and parent manifest hash.')
        if asset.get('status', 'complete') not in ('complete', 'completed'):
            raise ValueError('Extra report assets must be explicitly completed.')
    for raw in value.get('extra_views', []):
        if isinstance(raw, str):
            path = Path(raw)
            if not path.is_absolute(): path = base/path
            raw = json.loads(path.read_text('utf-8-sig'))
            source_files.append(path)
        if raw.get('project_id') != native.get('project_id') or raw.get('source_project_manifest_sha256') != manifest_identity:
            raise ValueError('Extra report views require matching project identity and parent manifest hash.')
        if raw.get('status', 'complete') not in ('complete', 'completed'):
            raise ValueError('Extra report views must be completed; running or failed collections cannot be presented as completed results.')
        views.append(copy.deepcopy(raw))
    supplemental = []
    for item in value.get('supplemental_sources', []):
        if isinstance(item, str):
            item = str((base / item).resolve()) if not Path(item).is_absolute() else item
        child = load_source(item)
        parent = child['project'].get('metadata', {})
        parent_id = parent.get('source_project_id')
        expected = parent.get('source_project_manifest_SHA256')
        actual = manifest_identity
        same_samples = [(s['sample_id'], s.get('biological_sample_id'), s.get('condition')) for s in native.get('samples', [])] == [(s['sample_id'], s.get('biological_sample_id'), s.get('condition')) for s in child['project'].get('samples', [])]
        same_contrasts = native.get('contrasts') == child['project'].get('contrasts')
        if not (same_samples and same_contrasts and expected and expected == actual):
            raise ValueError('Supplemental results require the exact parent manifest hash, samples and comparisons.')
        if parent_id and parent_id != native.get('project_id'):
            raise ValueError('Supplemental result parent identity does not match the source study.')
        supplemental.append(child)
    synthetic = bool(value.get('dummy', value.get('synthetic', native.get('metadata', {}).get('dummy', native.get('metadata', {}).get('synthetic', False)))))
    return {'project': native, 'artifacts': artifacts, 'supplemental': supplemental, 'extra_assets': extras, 'extra_views': views, 'dummy': synthetic,
            'search_roots': value.get('search_roots', []), 'source_files': source_files,
            'mode': value.get('mode', 'single-cell' if str(native.get('metadata', {}).get('assay', '')).startswith(('sc', 'sn', 'single')) else 'bulk')}


def read_table(path, *, columns=None, filters=None):
    path = Path(path)
    import pandas as pd
    if path.suffix == '.parquet':
        return pd.read_parquet(path, columns=columns, filters=filters)
    return pd.read_csv(path, sep='\t' if any(s.lower() in ('.tsv', '.txt') for s in path.suffixes) else ',', usecols=columns)


def table_schema(path):
    path = Path(path)
    if path.suffix == '.parquet':
        import pyarrow.parquet as pq
        meta = pq.ParquetFile(path)
        return {'columns': meta.schema_arrow.names, 'rows': meta.metadata.num_rows}
    import gzip
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8-sig') as handle:
        header = handle.readline().rstrip('\r\n').split('\t' if any(s.lower() in ('.tsv', '.txt') for s in path.suffixes) else ',')
        count = sum(1 for _ in handle)
    return {'columns': header, 'rows': count}
