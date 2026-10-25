"""Refresh a portable completed report without invoking scientific workflows."""
from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import zipfile

from .engine import ASSET_ROOT, Builder
from .sources import digest


def write_archive(output, archive_path=None):
    output = Path(output).resolve()
    archive_path = Path(archive_path).resolve() if archive_path else output / 'REPORT.zip'
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    if archive_path.exists():
        raise FileExistsError('Preserve the existing archive; choose a new report directory.')
    with zipfile.ZipFile(archive_path, 'w', zipfile.ZIP_DEFLATED,
                         compresslevel=6, strict_timestamps=False) as archive:
        for path in sorted(output.rglob('*')):
            if path.is_file() and path != archive_path:
                archive.write(path, path.relative_to(output).as_posix())
    return archive_path


def _project(manifest, project):
    if project is None:
        value = deepcopy(manifest.get('display_project'))
        if not value:
            raise ValueError('This older report needs --project for its display-only sample/design metadata.')
    elif isinstance(project, dict):
        value = deepcopy(project)
    else:
        location = Path(project).expanduser().resolve()
        value = json.loads(location.read_text('utf-8-sig'))
        if value.get('project_manifest'):
            parent = Path(value['project_manifest'])
            if not parent.is_absolute(): parent = location.parent / parent
            value = json.loads(parent.read_text('utf-8-sig'))
        else:
            value = value.get('project', value)
    if value.get('project_id') != manifest.get('project_id'):
        raise ValueError('Report and display project have different identities.')
    contrasts = value.get('contrasts', [])
    by_id = {c['contrast_id']: c for c in contrasts}
    selected = manifest.get('contrasts', [])
    if len(by_id) != len(contrasts) or not set(selected) <= set(by_id):
        raise ValueError('Display project must retain every exact selected comparison.')
    value['contrasts'] = [by_id[cid] for cid in selected]
    recorded = manifest.get('display_project')
    if recorded:
        from .presentation import portable_display_project
        current_display = portable_display_project(value)
        if (current_display['contrasts'] != recorded.get('contrasts')
                or current_display['samples'] != recorded.get('samples')):
            raise ValueError('Presentation refresh cannot change the recorded contrast arms or biological sample design.')
    return value


def refresh_report(report_dir, output_dir, *, project=None, create_archive=True,
                   compact=False, network_parameters=None):
    """Refresh existing figures/data in a new directory; never refit a result.

    ``network_parameters`` may declare verified original network parameters for
    a report-only graph export. No network discovery or statistical model is run.
    """
    source = Path(report_dir).expanduser().resolve()
    if source.is_file(): source = source.parent
    output = Path(output_dir).expanduser().resolve()
    if output == source or output.is_relative_to(source) or source.is_relative_to(output):
        raise ValueError('Presentation refresh must use a separate report directory.')
    if output.exists():
        raise FileExistsError('Completed reports are preserved; select a new directory.')
    manifest_path, provenance_path = source / 'report-manifest.json', source / 'provenance.json'
    manifest = json.loads(manifest_path.read_text('utf-8-sig'))
    if manifest.get('schema') != 'katana.approved_report/v1':
        raise ValueError('A complete portable Katana report manifest is required.')
    project_value = _project(manifest, project)
    control_paths = [manifest_path, provenance_path, source / 'REPORT.html']
    previous_receipt = source / 'presentation-refresh.json'
    if previous_receipt.is_file():
        control_paths.append(previous_receipt)
    before = {str(p): digest(p) for p in control_paths}
    before_stats = {str(p): (p.stat().st_size, p.stat().st_mtime_ns) for p in control_paths}

    def original_unchanged():
        return all(p.is_file() and digest(p) == before[str(p)]
                   and (p.stat().st_size, p.stat().st_mtime_ns) == before_stats[str(p)]
                   for p in control_paths)

    provenance = json.loads(provenance_path.read_text('utf-8-sig'))
    output.mkdir(parents=True)
    excluded = {'REPORT.html', 'REPORT.zip', 'report-manifest.json', 'provenance.json'}
    for path in source.rglob('*'):
        if not path.is_file(): continue
        relative = path.relative_to(source)
        if (len(relative.parts) == 1 and path.name in excluded) or relative.parts[0] == 'assets' or relative.as_posix() == 'outputs/index.html':
            continue
        destination = output / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        try: os.link(path, destination)
        except OSError: shutil.copy2(path, destination)
    shutil.copytree(ASSET_ROOT, output / 'assets')

    from .network_views import refresh_completed_network_browser
    from .network_display_upgrade import bind_legacy_network_table_scopes, upgrade_copied_network_browser
    scope_binding = bind_legacy_network_table_scopes(manifest)
    network_receipts = []
    for table in manifest.get('complete_results', []):
        if (table.get('chapter') != 6 or table.get('analysis_id') != 'rich_browser'
                or (table.get('metadata', {}).get('network_view_role') != 'module_browser'
                    and not table.get('name', '').startswith('All network modules:'))):
            continue
        relative = Path(table['source'])
        if relative.is_absolute() or not (source / relative).resolve().is_relative_to(source):
            raise ValueError('Network browser escapes the portable report.')
        options = {}
        if network_parameters and table.get('contrast_id') in network_parameters:
            options['parameters_path'] = network_parameters[table['contrast_id']]
        if options:
            receipt = refresh_completed_network_browser(source / relative, output / relative, **options)
        else:
            asset = {'path': str(source / relative), 'metadata': table.get('metadata', {})}
            receipt = upgrade_copied_network_browser(asset, output / relative)
            if receipt is None:
                receipt = {'schema': 'katana.completed_network_browser_refresh/v1', 'passed': True,
                           'already_current': True, 'network_discovery_runs': 0,
                           'association_tests': 0, 'enrichment_tests': 0}
        receipt['contrast_id'] = table.get('contrast_id')
        network_receipts.append(receipt)
        table.setdefault('metadata', {})['presentation_refresh'] = receipt
        table['metadata']['portable_sha256'] = digest(output / relative)

    from .motif_display import apply_motif_significance_display
    manifest, motif_audit = apply_motif_significance_display(manifest, output)
    builder = Builder({'project': project_value, 'artifacts': [], 'mode': manifest.get('mode', 'bulk'),
                       'report_selection': manifest.get('report_selection', {}),
                       'report_require_plot_data': bool(manifest.get('figure_data_links', {}).get('strict_plot_data_required'))},
                      output, manifest.get('event_set', 'stable_significant'))
    for attribute, key in [('panels', 'panels'), ('tables', 'complete_results'), ('top10', 'top10'),
                           ('policy_data', 'policy_data'), ('audit', 'source_files'),
                           ('enrichment_summaries', 'rbp_enrichment_summaries'),
                           ('motif_summaries', 'motif_summaries'), ('exclusions', 'unavailable')]:
        if key in manifest: setattr(builder, attribute, deepcopy(manifest[key]))
    builder.data_link_audit = deepcopy(manifest.get('figure_data_links', {}).get('bindings', []))
    builder.derived_audits = deepcopy(provenance.get('derived_views', []))
    builder.preserved_source_projects = deepcopy(provenance.get('source_projects', []))
    builder.derived_audits.append({'schema': 'katana.report_presentation_refresh/v1',
        'source_manifest_sha256': before[str(manifest_path)], 'network_browsers': network_receipts,
        'motif_significance': motif_audit, 'network_table_scope_binding': scope_binding, 'completed_analyses_rerun': 0})
    if manifest.get('release_status'): builder.release_status = deepcopy(manifest['release_status'])
    builder.create_archive = False
    result = builder.write()
    if compact:
        from .compact import compact_report
        compact_report(output)
        result = json.loads((output / 'report-manifest.json').read_text('utf-8'))
    receipt = {'schema': 'katana.completed_report_refresh/v1', 'status': 'completed',
        'utc': datetime.now(timezone.utc).isoformat(), 'source_manifest_sha256': before[str(manifest_path)],
        'network_browsers': network_receipts, 'motif_significance': motif_audit,
        'completed_analyses_rerun': 0, 'original_report_unchanged': original_unchanged(),
        'compact': compact}
    if not receipt['original_report_unchanged']:
        raise RuntimeError('Original report changed during refresh.')
    # The copied lifecycle receipt can share an inode with the source. Replace
    # its directory entry instead of writing through that retained hard link.
    from .compact import _atomic_write
    _atomic_write(output / 'presentation-refresh.json',
                  (json.dumps(receipt, indent=2, default=str) + '\n').encode('utf-8'))
    if create_archive: write_archive(output)
    if not original_unchanged():
        raise RuntimeError('Original report changed during refresh.')
    return result
