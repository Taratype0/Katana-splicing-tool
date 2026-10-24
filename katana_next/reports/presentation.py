"""Concise report presentation from completed, comparison-scoped evidence."""


def prepare_presentation(builder):
    """Keep native outputs while exposing networks and significant motifs."""
    for panel in builder.panels:
        if panel.get('chapter') == 6 and 'dendrogram' in panel.get('title', '').lower():
            panel['output_only'] = True
            panel.setdefault('metadata', {})['main_display_note'] = 'Dendrogram retained in outputs by request.'
    browsers, seen = [], set()
    for record in builder.tables:
        if (record.get('chapter') != 6 or record.get('analysis_id') != 'rich_browser'
                or (record.get('metadata', {}).get('network_view_role') != 'module_browser'
                    and not record.get('name', '').startswith('All network modules:'))):
            continue
        cid = record.get('contrast_id', 'all')
        if cid in seen:
            raise ValueError('One completed module browser is required per independent network scope.')
        path = (builder.output / record['source']).resolve()
        if not path.is_relative_to(builder.output) or not path.is_file():
            raise ValueError('Main network browser must be a complete bounded report member.')
        seen.add(cid)
        browsers.append({'id': 'network-browser-' + str(cid), 'chapter': 6, 'contrast_id': cid,
                         'source': record['source'], 'title': 'Module network and complete gene lists',
                         'metadata': record.get('metadata', {})})
    builder.network_browsers = browsers


def portable_display_project(project):
    """Save only display metadata, without workstation references or BAM paths."""
    return {'project_id': project.get('project_id'), 'name': project.get('name', 'Katana report'),
            'samples': [{k: s[k] for k in ('sample_id', 'biological_sample_id', 'condition', 'attributes') if k in s}
                        for s in project.get('samples', [])],
            'contrasts': project.get('contrasts', []),
            'annotation': {k: project.get('annotation', {})[k]
                           for k in ('genome_build', 'annotation_release', 'species') if k in project.get('annotation', {})},
            'metadata': {k: project.get('metadata', {})[k]
                         for k in ('annotation_release', 'gene_filter_policy') if k in project.get('metadata', {})}}
