"""Read project manifests and artifact exports without opening their database."""
from pathlib import Path
import copy
import hashlib, json


def digest(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def _absolute_path(value, base, mapper=None):
    if mapper is not None:
        return mapper.resolve(value, base)
    path = Path(value).expanduser()
    return (path if path.is_absolute() else base / path).resolve()


def _is_remote_uri(value):
    return '://' in str(value)


def _native_paths(raw, base, mapper=None):
    """Resolve local manifest path fields on a copy, never arbitrary labels."""
    native = mapper.project(raw) if mapper else copy.deepcopy(raw)

    def normalize(record, keys):
        for key in keys:
            if record.get(key) and not _is_remote_uri(record[key]):
                record[key] = str(_absolute_path(record[key], base, mapper))

    normalize(native, ('input_root', 'output_root'))
    normalize(native.get('annotation') or {}, ('gtf', 'genome_fasta', 'transcript_fasta'))
    for asset in native.get('assets', []):
        normalize(asset, ('path',))
    for sample in native.get('samples', []):
        normalize(sample.get('assets') or {}, ('quant_sf',))
    return native


def _view_paths(raw, base, mapper=None):
    """Resolve only declared view file references, without changing metadata."""
    views = mapper.project(raw) if mapper else copy.deepcopy(raw)
    for section, keys in (('panels', ('path', 'image', 'svg')),
                          ('tables', ('path', 'source', 'view')),
                          ('top10', ('full_path', 'full_table'))):
        for record in views.get(section, []):
            for key in keys:
                if record.get(key):
                    record[key] = str(_absolute_path(record[key], base, mapper))
            for link in record.get('data_links', []):
                if link.get('path'):
                    link['path'] = str(_absolute_path(link['path'], base, mapper))
                if link.get('source_collection_root'):
                    link['source_collection_root'] = str(_absolute_path(link['source_collection_root'], base, mapper))
    return views


def report_contrast_selection(native, descriptor):
    """Validate a report-only projection, retaining original contrast identity.

    No group names or biological meanings are inferred. In the absence of
    explicit groups, stable ordered groups of at most four axes are used.
    """
    contrasts = native.get('contrasts', [])
    native_ids = [row.get('contrast_id') for row in contrasts]
    if any(not isinstance(cid, str) or not cid.strip() for cid in native_ids) or len(native_ids) != len(set(native_ids)):
        raise ValueError('Native contrast IDs must be nonempty and unique before report selection.')
    chosen = descriptor.get('report_contrast_ids', native.get('report_contrast_ids', native_ids))
    if not isinstance(chosen, list) or any(not isinstance(cid, str) for cid in chosen) or len(chosen) != len(set(chosen)):
        raise ValueError('report_contrast_ids must be a unique list of native contrast IDs.')
    if set(chosen) - set(native_ids):
        raise ValueError('Report selection contains an unknown native contrast ID.')
    if native_ids and not chosen:
        raise ValueError('A nonempty source study requires at least one selected report contrast.')
    supplied = descriptor.get('report_contrast_groups', native.get('report_contrast_groups'))
    groups = []
    if supplied is None:
        for start in range(0, len(chosen), 4):
            ids = chosen[start:start + 4]
            groups.append({'group_id': 'comparisons-' + str(start // 4 + 1),
                           'label': 'Comparisons ' + str(start // 4 + 1), 'contrast_ids': ids})
    else:
        if not isinstance(supplied, list):
            raise ValueError('report_contrast_groups must be a list of explicit group records.')
        for raw in supplied:
            if not isinstance(raw, dict) or not isinstance(raw.get('group_id'), str) or not raw['group_id'].strip():
                raise ValueError('Each report contrast group requires a nonempty group_id.')
            ids = raw.get('contrast_ids')
            if not isinstance(ids, list) or not 1 <= len(ids) <= 4 or any(not isinstance(cid, str) for cid in ids) or len(ids) != len(set(ids)):
                raise ValueError('Each report contrast group needs one to four distinct contrast IDs.')
            if not set(ids) <= set(chosen):
                raise ValueError('Report group axes must belong to the selected native contrasts.')
            groups.append({'group_id': raw['group_id'], 'label': str(raw.get('label', raw['group_id'])), 'contrast_ids': list(ids)})
    group_ids = [group['group_id'] for group in groups]
    axes = [cid for group in groups for cid in group['contrast_ids']]
    if len(group_ids) != len(set(group_ids)) or len(axes) != len(set(axes)) or set(axes) != set(chosen):
        raise ValueError('Report groups must have unique IDs and partition the selected comparisons exactly once.')
    return {'contrast_ids': list(chosen), 'excluded_contrast_ids': [cid for cid in native_ids if cid not in chosen],
            'groups': groups, 'applied': chosen != native_ids or supplied is not None,
            'native_contrast_ids': native_ids, 'statistical_results_changed': False}


def load_source(project, catalog=None, *, path_mappings=None):
    """Load a descriptor using its own directory, or cwd for inline input.

    An explicit catalog argument follows command-line cwd semantics. Paths in
    a descriptor follow that descriptor, while native artifact paths follow the
    native project manifest. A file of extra views owns its own relative paths.
    """
    from .path_mapping import PathMappings
    mapper = PathMappings(path_mappings)
    return _load_source(project, catalog, inline_base=Path.cwd(), mapper=mapper)


def _artifact_comparisons(record):
    """Use recorded comparison identity, including native enrichment lineage."""
    comparisons = record.get('contrast_ids') or record.get('comparison_ids') or []
    if not comparisons:
        comparisons = [scope.get(key) for scope in (record, record.get('parameters') or {}, record.get('metadata') or {})
                       for key in ('contrast_id', 'comparison_id') if scope.get(key)]
    if not comparisons:
        provenance = (record.get('metadata') or {}).get('native_background_provenance') or {}
        comparisons = [item['contrast_id'] for item in provenance.get('native_inputs', []) if item.get('contrast_id')]
    return tuple(sorted(set(comparisons)))


def _load_source(project, catalog=None, *, inline_base, mapper=None):
    from .path_mapping import PathMappings
    mapper = mapper or PathMappings()
    source_files = []
    if isinstance(project, dict):
        value = copy.deepcopy(project)
        descriptor_base = inline_base
    else:
        location = Path(project).expanduser().resolve()
        value = json.loads(location.read_text('utf-8-sig'))
        descriptor_base = location.parent
        source_files.append(location)
    if value.get('schema') == 'katana.completed_result_set/v1':
        if catalog is not None: raise ValueError('A completed result set declares each source catalog explicitly')
        from .completed_sources import assemble_result_set
        return assemble_result_set(value, base=descriptor_base, source_files=source_files, mapper=mapper)
    value = mapper.project(value)
    native_base = descriptor_base
    if value.get('project_manifest'):
        manifest = _absolute_path(value['project_manifest'], descriptor_base, mapper)
        native = json.loads(manifest.read_text('utf-8-sig'))
        source_files.append(manifest)
        native_base = manifest.parent
    else:
        native = value.get('project', value)
    # Identity belongs to the exact source JSON, before paths are normalized in
    # a report-only copy. Existing supplementary evidence binds that identity.
    native_location = manifest if value.get('project_manifest') else (location if not isinstance(project, dict) and 'project' not in value else None)
    manifest_identity = digest(native_location) if native_location else hashlib.sha256(json.dumps(native, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()
    native = _native_paths(native, native_base, mapper=mapper)
    records = value.get('artifacts')
    explicit_catalog = catalog is not None
    catalog_base = Path.cwd() if catalog is not None else descriptor_base
    catalog = catalog if catalog is not None else (value.get('artifact_catalog') or value.get('catalog_export') or value.get('native_artifact_catalog'))
    if isinstance(records, str):
        if not explicit_catalog:
            catalog, catalog_base = records, descriptor_base
        records = None
    output_value = native.get('output_root', native_base)
    output_root = None if _is_remote_uri(output_value) else _absolute_path(output_value, native_base, mapper)
    if records is None:
        candidates = [_absolute_path(catalog, catalog_base, mapper)] if catalog else [native_base / 'native-artifact-catalog.json']
        if not catalog and output_root is not None:
            candidates.append(output_root / '.katana/artifact_manifest.json')
        for path in candidates:
            if path.is_file():
                records = json.loads(path.read_text('utf-8-sig'))
                records = records.get('artifacts', records) if isinstance(records, dict) else records
                source_files.append(path)
                break
    if records is None:
        if output_root is None:
            raise ValueError('Remote result URIs are not read by the report loader; supply a local completed artifact export.')
        db = output_root / '.katana/catalog.duckdb'
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
        record = mapper.project(raw)
        if record.get('status', 'complete') not in ('complete', 'completed'):
            continue
        if record.get('project_id') and record['project_id'] != native.get('project_id'):
            raise ValueError('Artifact belongs to a different project; use an explicitly verified supplemental source.')
        path = _absolute_path(record.get('path', ''), native_base, mapper)
        if not path.is_file():
            record['missing_reason'] = 'Recorded completed artifact file is unavailable.'
        record['path'] = str(path.resolve())
        if path.is_file():
            from .rich_results import _completion
            complete, evidence = _completion(record, path)
            if not complete:
                continue
            record.setdefault('metadata', {})['report_completion_evidence'] = evidence
        parameters = record.get('parameters') or {}
        key = (record.get('analysis_id'), record.get('name'), record.get('kind'), record.get('entity_type'),
               _artifact_comparisons(record), parameters.get('sequence_scope'), parameters.get('event_set'),
               parameters.get('foreground_source', (record.get('metadata') or {}).get('foreground_source')))
        if key in seen:
            continue
        seen.add(key)
        artifacts.append(record)
    extras = copy.deepcopy(value.get('extra_assets', []))
    views = []
    for asset in extras:
        meta = asset.get('metadata', {})
        if asset.get('project_id') != native.get('project_id') or meta.get('source_project_manifest_sha256', meta.get('source_project_manifest_SHA256')) != manifest_identity:
            raise ValueError('Extra report assets require matching project identity and parent manifest hash.')
        if asset.get('status', 'complete') not in ('complete', 'completed'):
            raise ValueError('Extra report assets must be explicitly completed.')
        if asset.get('path'):
            asset['path'] = str(_absolute_path(asset['path'], descriptor_base, mapper))
        for key, path in meta.get('exports', {}).items():
            if key in ('png', 'svg') and path:
                meta['exports'][key] = str(_absolute_path(path, descriptor_base, mapper))
        for link in meta.get('data_links', []):
            if link.get('path'):
                link['path'] = str(_absolute_path(link['path'], descriptor_base, mapper))
    for raw in value.get('extra_views', []):
        view_base = descriptor_base
        if isinstance(raw, str):
            path = _absolute_path(raw, descriptor_base, mapper)
            raw = json.loads(path.read_text('utf-8-sig'))
            source_files.append(path)
            view_base = path.parent
        if raw.get('project_id') != native.get('project_id') or raw.get('source_project_manifest_sha256') != manifest_identity:
            raise ValueError('Extra report views require matching project identity and parent manifest hash.')
        if raw.get('status', 'complete') not in ('complete', 'completed'):
            raise ValueError('Extra report views must be completed; running or failed collections cannot be presented as completed results.')
        views.append(_view_paths(raw, view_base, mapper=mapper))
    supplemental = []
    for item in value.get('supplemental_sources', []):
        if isinstance(item, str):
            item = str(_absolute_path(item, descriptor_base, mapper))
        child = _load_source(item, inline_base=descriptor_base, mapper=mapper)
        parent = child['project'].get('metadata', {})
        parent_id = parent.get('source_project_id')
        expected = parent.get('source_project_manifest_SHA256')
        actual = manifest_identity
        child_native = child.get('native_project', child['project'])
        same_samples = [(s['sample_id'], s.get('biological_sample_id'), s.get('condition')) for s in native.get('samples', [])] == [(s['sample_id'], s.get('biological_sample_id'), s.get('condition')) for s in child_native.get('samples', [])]
        same_contrasts = native.get('contrasts') == child_native.get('contrasts')
        if not (same_samples and same_contrasts and expected and expected == actual):
            raise ValueError('Supplemental results require the exact parent manifest hash, samples and comparisons.')
        if parent_id and parent_id != native.get('project_id'):
            raise ValueError('Supplemental result parent identity does not match the source study.')
        supplemental.append(child)
    synthetic = bool(value.get('dummy', value.get('synthetic', native.get('metadata', {}).get('dummy', native.get('metadata', {}).get('synthetic', False)))))
    selection = report_contrast_selection(native, value)
    display_project = copy.deepcopy(native)
    if native.get('contrasts') or 'report_contrast_ids' in value or 'report_contrast_groups' in value:
        indexed = {row['contrast_id']: row for row in native.get('contrasts', [])}
        display_project['contrasts'] = [copy.deepcopy(indexed[cid]) for cid in selection['contrast_ids']]
        display_project['report_contrast_ids'] = selection['contrast_ids']
        display_project['report_contrast_groups'] = selection['groups']
        display_project['report_native_contrast_ids'] = selection['native_contrast_ids']
    return {'project': display_project, 'native_project': native, 'native_project_manifest_sha256': manifest_identity,
            'report_selection': selection, 'artifacts': artifacts, 'supplemental': supplemental, 'extra_assets': extras, 'extra_views': views, 'dummy': synthetic,
            'report_native_figures_in_main': bool(value.get('report_native_figures_in_main', True)),
            'report_require_plot_data': bool(value.get('report_require_plot_data', False)),
            'search_roots': value.get('search_roots', []), 'source_files': source_files,
            'mode': value.get('mode', 'single-cell' if str(native.get('metadata', {}).get('assay', '')).startswith(('sc', 'sn', 'single')) else 'bulk')}


def read_table(path, *, columns=None, filters=None):
    path = Path(path)
    import pandas as pd
    if path.suffix == '.parquet':
        import pyarrow.dataset as ds
        import pyarrow.parquet as pq
        # Large completed tables often put several contrasts in one row group.
        # Parallel row-group prefetch can expand the whole study before the
        # requested contrast is returned. Bound scan buffers without changing
        # predicates, row order, types or native statistical values.
        scanner = ds.dataset(path, format='parquet').scanner(
            columns=columns, filter=pq.filters_to_expression(filters) if filters else None,
            batch_size=16384, batch_readahead=1, fragment_readahead=1, use_threads=False)
        return scanner.to_table().to_pandas(use_threads=False)
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
