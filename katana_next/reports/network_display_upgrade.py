"""Upgrade copied legacy network views using their unchanged completed proof."""
from __future__ import annotations

from pathlib import Path
import json
import re


def _completed_inputs(metadata, *, require_display_parameters=True):
    from .network_views import _completed_network_scope, _fingerprint
    from katana_next.flow.runtime import canonical_hash

    completion = Path(metadata['completion_manifest']).resolve()
    proof = json.loads(completion.read_text(encoding='utf-8-sig'))
    if proof.get('schema') != 'katana.flow.network_view_inputs/v1':
        raise ValueError('Automatic network display upgrade needs the exact completed flow proof')
    root = completion.parent.parent
    receipt_path = root / 'collection.receipt.json'
    receipt = json.loads(receipt_path.read_text(encoding='utf-8-sig'))
    if (receipt.get('schema') != 'katana.flow.collection_receipt/v1'
            or receipt.get('status') != 'completed'
            or receipt.get('receipt_sha256') != canonical_hash({k: v for k, v in receipt.items() if k != 'receipt_sha256'})
            or receipt.get('plan_id') != proof.get('plan_id')):
        raise ValueError('Network display upgrade needs the unchanged completed collection receipt')
    identity = proof.get('network_identity')
    if metadata.get('network_identity') != identity:
        raise ValueError('Network display metadata and completed sample/contrast scope differ')
    names = ('event_module_assignments.tsv', 'module_program_scores.tsv', 'selected_event_matrix.tsv',
             'event_metadata.tsv', 'sample_metadata_resolved.tsv', 'network_edges.tsv')
    inputs = []
    for name in names:
        rows = [row for row in proof.get('files', []) if Path(row['path']).name == name]
        if len(rows) != 1:
            raise ValueError('Network display requires one exact completed input: ' + name)
        path = (root / rows[0]['path']).resolve()
        if not path.is_relative_to(root):
            raise ValueError('Completed network display input escapes its collection')
        inputs.append(path)
    scope, _, precomputed, authenticated = _completed_network_scope(completion, inputs, proof['contrasts'])
    if scope != 'flow_verified_network_complete':
        raise ValueError('Network display input scope is not a verified completed flow')
    # Legacy view proofs omit settings/diagnostics; the optional bridge above
    # authenticates them against this exact fit's native output ledger.
    parameters = inputs[0].with_name('parameters.json')
    required = {parameters, *[parameters.with_name(name) for name in (
        'network_diagnostics.tsv', 'event_module_assignments.tsv', 'selected_event_matrix.tsv', 'network_edges.tsv')]}
    if require_display_parameters and not required <= set(authenticated):
        raise ValueError('Complete saved settings and fitted numerical inputs are required for the approved module graph')
    binding = metadata.get('report_source_binding', {})
    if binding.get('native_manifest_sha256') and _fingerprint(root / 'project.json')['sha256'] != binding['native_manifest_sha256']:
        raise ValueError('Network display source binding differs from its actual native project')
    return root, receipt, completion, proof, parameters


def _verify_registered_browser(source, metadata, root, receipt):
    """Bind even an already portable HTML copy to its actual native browser."""
    from .network_views import _fingerprint
    native = (Path(metadata['portable_root']) / metadata['relative_path']).resolve()
    if not native.is_relative_to(root) or not native.is_file():
        raise ValueError('The native completed network browser is missing or unbounded')
    relative = native.relative_to(root).as_posix()
    rows = [r for r in receipt['outputs'] if r.get('relative_path', r['path']) == relative]
    actual = _fingerprint(native)
    if len(rows) != 1 or (rows[0]['sha256'], rows[0]['size_bytes']) != (actual['sha256'], actual['bytes']):
        raise ValueError('Native network browser differs from its completed collection ledger')
    import re
    def payload(path):
        match = re.search(r'<script\s+id=["\']data["\']\s+type=["\']application/json["\']>(.*?)</script>',
                          path.read_text(encoding='utf-8-sig'), re.S)
        if match is None:
            raise ValueError('Completed browser has no full native network payload')
        return json.loads(match.group(1))
    if payload(source) != payload(native):
        raise ValueError('Portable network browser differs from its authenticated native numerical payload')


def upgrade_copied_network_browser(asset, target):
    """Refresh a destination copy only; a saved collection is never rewritten."""
    from .network_views import refresh_completed_network_browser
    metadata = asset.get('metadata', {})
    if metadata.get('network_view_role') != 'module_browser':
        return None
    source, target = Path(asset['path']).resolve(), Path(target).resolve()
    if source == target:
        raise ValueError('A network display upgrade must write an independent copy')
    # A fresh current collection already contains the approved browser. Keep it
    # byte-identical instead of repeating its drawing export during packaging.
    text = source.read_text(encoding='utf-8-sig')
    graph_match = re.search(r'<script\s+id=["\']module-display-data["\'][^>]*>(.*?)</script>', text, re.S)
    graph = json.loads(graph_match.group(1)) if graph_match else None
    current_controls = 'id="gene-body"' in text and 'id="network-display-size"' in text
    current_graph = isinstance(graph, dict) and isinstance(graph.get('modules'), dict)
    # Zero assigned modules is a valid completed graph, not missing analysis.
    # Historical bootstrap reports can already have the current complete-gene
    # UI while honestly retaining only their original exported native edges.
    # They have no no-bootstrap |kME| graph to re-export during packaging.
    current_native_edges = graph_match is not None and graph is None and metadata.get('bootstrap_enabled') is True
    if current_controls and (current_graph or current_native_edges):
        return None
    if metadata.get('completion_scope') != 'flow_verified_network_complete':
        # Older externally imported networks can gain complete gene browsing;
        # an unproved native graph is not invented for them.
        return refresh_completed_network_browser(source, target)
    root, receipt, completion, proof, parameters = _completed_inputs(metadata)
    _verify_registered_browser(source, metadata, root, receipt)
    result = refresh_completed_network_browser(source, target, parameters_path=parameters)
    result['authentication'] = {'completion_manifest': str(completion),
        'plan_id': proof['plan_id'], 'collection_receipt_sha256': receipt['receipt_sha256'],
        'network_identity': proof.get('network_identity'), 'original_collection_modified': False}
    return result


def bind_legacy_network_table_scopes(manifest):
    """Recover missing table scope from authenticated native record identities.

    No contrast is inferred from the text of a title or a filename. Artifact ID
    is preferred; exact registered name plus content SHA supports old attachment
    entries that omitted their artifact ID. Ambiguous records stay unmodified.
    """
    from katana_next.flow.collect import materialize_report_descriptor
    from .network_views import _fingerprint

    by_id, by_record, visited = {}, {}, set()
    for browser in manifest.get('network_browsers', []):
        metadata = browser.get('metadata', {})
        if metadata.get('completion_scope') != 'flow_verified_network_complete':
            continue
        completion = Path(metadata.get('completion_manifest', ''))
        if not completion.is_file():
            continue  # Portable-only refresh cannot manufacture absent proof.
        root, _, _, _, _ = _completed_inputs(metadata, require_display_parameters=False)
        if root in visited:
            continue
        visited.add(root)
        descriptor = materialize_report_descriptor(root)
        for artifact in [*descriptor.get('artifacts', []), *descriptor.get('extra_assets', [])]:
            ids = artifact.get('contrast_ids') or []
            if len(ids) != 1:
                continue
            cid = ids[0]
            if artifact.get('artifact_id'):
                by_id.setdefault(artifact['artifact_id'], set()).add(cid)
            path = Path(artifact['path'])
            if path.is_file():
                key = (artifact.get('analysis_id', ''), artifact.get('name', path.stem), _fingerprint(path)['sha256'])
                by_record.setdefault(key, set()).add(cid)
    hashes = {row['report_path']: row.get('portable_sha256', row['source_sha256'])
              for row in manifest.get('source_files', [])}
    bindings = []
    for table in manifest.get('complete_results', []):
        if table.get('contrast_id') not in (None, '') or table.get('analysis_id') == 'rich_browser':
            continue
        ids = by_id.get(table.get('artifact_id'), set())
        if not ids and table.get('source') in hashes:
            ids = by_record.get((table.get('analysis_id', ''), table.get('name'), hashes[table['source']]), set())
        if len(ids) == 1:
            cid = next(iter(ids))
            table['contrast_id'] = cid
            table.setdefault('metadata', {})['legacy_scope_binding'] = {
                'source': 'verified native artifact identity and exact content SHA',
                'contrast_id': cid, 'title_parsed_for_contrast': False, 'statistics_recomputed': False}
            bindings.append({'artifact_id': table.get('artifact_id'), 'source': table.get('source'), 'contrast_id': cid})
    return {'schema': 'katana.legacy_network_table_scope_binding/v1', 'bound_tables': bindings,
            'verified_collections': len(visited), 'statistics_recomputed': False}
