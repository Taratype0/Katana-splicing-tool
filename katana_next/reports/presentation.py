"""Concise report presentation from completed, comparison-scoped evidence."""


def prepare_presentation(builder):
    """Keep native outputs while exposing networks and significant motifs."""
    selected = {row['contrast_id'] for row in builder.source['project'].get('contrasts', [])}
    # Completed sources can contain more comparisons than the requested report.
    # Their native tables remain in outputs; only selected summaries belong in
    # the approved one-comparison presentation.
    builder.enrichment_summaries = [
        row for row in builder.enrichment_summaries
        if not selected or row.get('contrast_id') in selected or row.get('contrast_id') == 'all'
    ]
    rbp_scopes = {(row['contrast_id'], row.get('cell_type', 'all'))
                  for row in builder.enrichment_summaries
                  if str(row.get('id', '')).startswith('rbp-enrichment-')}
    for panel in builder.panels:
        metadata = panel.get('metadata', {})
        # Legacy completed ORA views did not carry analysis_id. Their explicit
        # RBP-identity title and exact comparison/cell scope identify the saved
        # redundant charts. Keep every chart and its data in outputs, and show
        # the compact significant-set summary with complete contributors.
        rbp_identity = (metadata.get('enrichment_category') == 'RBP_identity'
                        or panel.get('title', '').startswith('RBP identity |'))
        if (panel.get('chapter') == 7 and rbp_identity
                and (panel.get('contrast_id'), panel.get('cell_type', 'all')) in rbp_scopes):
            panel['output_only'] = True
            panel.setdefault('metadata', {})['main_display_note'] = (
                'Compact RBP enrichment and complete hit genes are shown in the main view; '
                'the recorded identity chart and all tests remain in outputs.')
        recorded_tree = (metadata.get('network_view_role') == 'dendrogram' or
                         (metadata.get('flow_stage') == 'network.fit' and
                          metadata.get('flow_artifact_key') == 'dendrogram'))
        if (panel.get('chapter') == 6 and
                (recorded_tree or 'dendrogram' in panel.get('title', '').lower())):
            panel['output_only'] = True
            panel.setdefault('metadata', {})['main_display_note'] = 'Dendrogram retained in outputs by request.'
        if (panel.get('chapter') == 7 and panel.get('contrast_id', 'all') == 'all'
                and 'reference clip' in panel.get('title', '').lower()
                and 'status' in panel.get('title', '').lower()):
            panel['output_only'] = True
            panel.setdefault('metadata', {})['main_display_note'] = 'Whole-study completion status remains in outputs; the main view is comparison-specific.'
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
