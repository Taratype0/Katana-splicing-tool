"""Exact data bindings for simple, self-contained saved Vega-Lite figures.

This adapter copies the producer's recorded rows. It does not infer a parent
table, refilter a result, evaluate a transform, or run statistical analysis.
Unsupported specifications keep the normal strict missing-data failure.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _simple_values(spec: dict) -> list[dict] | None:
    if not isinstance(spec, dict) or 'vega-lite' not in str(spec.get('$schema', '')):
        return None
    # These constructs need their own explicit post-transform producer binding.
    if any(key in spec for key in ('transform', 'layer', 'concat', 'hconcat', 'vconcat', 'facet', 'repeat', 'datasets')):
        return None
    if not isinstance(spec.get('mark'), (str, dict)):
        return None
    data = spec.get('data')
    if not isinstance(data, dict) or set(data) != {'values'}:
        return None
    rows = data['values']
    if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
        return None
    def transforms(value):
        if isinstance(value, dict):
            return any(key in value for key in ('aggregate', 'bin', 'timeUnit', 'stack')) or any(transforms(v) for v in value.values())
        return isinstance(value, list) and any(transforms(v) for v in value)
    if transforms(spec.get('encoding', {})):
        return None
    return rows


def recorded_plot_data_links(artifact: dict, panel: dict, output: Path) -> list[dict]:
    """Materialize all inline rows only for the figure's exact saved source.

    An image export is accepted only with a binding explicitly propagated from
    its own producer artifact. Filenames or neighbouring tables are not used.
    JSON keeps original nulls/types; TSV provides convenient complete browsing.
    """
    path = Path(artifact['path']).resolve()
    metadata = artifact.get('metadata', {})
    source = path if path.name.endswith('.vl.json') else None
    binding = metadata.get('saved_plot_specification')
    if source is None and isinstance(binding, dict):
        source = Path(binding['path']).resolve()
        if binding.get('export_sha256') != _hash(path):
            raise ValueError('Saved plot export SHA does not match its producer binding.')
    if source is None or not source.is_file():
        return []
    before = source.stat()
    original = source.read_bytes()
    sha = hashlib.sha256(original).hexdigest()
    if binding and binding.get('source_sha256') != sha:
        raise ValueError('Saved plot specification SHA does not match its producer binding.')
    missing_values = []
    def invalid_constant(value):
        if value == 'NaN':
            missing_values.append(value)
            return None
        raise ValueError('Saved plot specification contains a non-finite number: ' + value)
    spec = json.loads(original.decode('utf-8-sig'), parse_constant=invalid_constant)
    rows = _simple_values(spec)
    if rows is None:
        return []
    contrast = str(panel.get('contrast_id') or 'all')
    policy = panel.get('policy')
    scope = panel.get('metadata', {}).get('das_scope', {})
    if scope.get('mode') in {'rs_required', 'fixed_completed_result'} and scope.get('recorded_event_set') in {'significant', 'stable_significant'}:
        policy = scope['recorded_event_set']
    for row in rows:
        for name in ('contrast_id', 'comparison_id'):
            if name in row and contrast != 'all' and str(row[name]) != contrast:
                raise ValueError('Saved plot rows contain a different comparison.')
        for name in ('report_policy', 'report_event_set', 'das_event_set'):
            if name in row and row[name] != policy:
                raise ValueError('Saved plot rows contain a different DAS policy.')
    target = Path(output).resolve() / 'outputs' / 'saved-plot-data' / sha
    target.mkdir(parents=True, exist_ok=True)
    json_path = target / 'values.json'
    values_bytes = (json.dumps(rows, ensure_ascii=False, allow_nan=False, indent=2) + '\n').encode('utf-8')
    json_path.write_bytes(values_bytes)
    columns = list(dict.fromkeys(key for row in rows for key in row))
    links = []
    common = {'contrast_id': contrast, 'policy': policy, 'unit': 'saved figure data row; field units recorded in specification',
              'producer_specification_sha256': sha, 'producer_artifact_id': metadata.get('source_artifact_id', artifact.get('artifact_id'))}
    links.append({**common, 'role': 'plot_data', 'path': str(json_path), 'label': 'Exact figure data (JSON)',
                  'source_sha256': _hash(json_path)})
    if columns:
        tsv_path = target / 'values.tsv'
        with tsv_path.open('w', encoding='utf-8', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            for row in rows:
                writer.writerow({key: json.dumps(value, ensure_ascii=False, allow_nan=False) if isinstance(value, (list, dict)) else value for key, value in row.items()})
        links.append({**common, 'role': 'plot_data', 'path': str(tsv_path), 'label': 'Exact figure data (TSV)',
                      'source_sha256': _hash(tsv_path)})
    links.append({**common, 'role': 'methods', 'path': str(source), 'label': 'Original figure specification', 'source_sha256': sha})
    receipt = {'schema': 'katana.saved_plot_data/v1', 'source_specification_sha256': sha,
               'producer_artifact_id': common['producer_artifact_id'], 'contrast_id': contrast, 'policy': policy,
               'rows': len(rows), 'columns': columns, 'data_json_sha256': _hash(json_path),
               'new_scientific_computation': False, 'all_recorded_plot_rows_retained': True,
               'statistical_values_preserved': True,
               'missing_value_encoding': {'source_NaN_to_JSON_null': len(missing_values), 'imputed_values': 0}}
    (target / 'binding.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8', newline='\n')
    after = source.stat()
    if (after.st_size, after.st_mtime_ns) != (before.st_size, before.st_mtime_ns) or _hash(source) != sha:
        raise RuntimeError('Saved plot specification changed during data collection.')
    return links
