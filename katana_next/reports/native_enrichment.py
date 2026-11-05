"""Display completed, explicitly scoped native ORA without refitting statistics."""
from __future__ import annotations


def native_enrichment_scope(artifact):
    """Do not substitute unbound historical plots for a current DAS selection."""
    if artifact.get('analysis_id') != 'prioritize.enrichment':
        return None
    meta = artifact.get('metadata', {})
    from .sources import _artifact_comparisons
    comparisons = _artifact_comparisons(artifact)
    foreground = meta.get('foreground_source')
    if (len(comparisons) != 1 or foreground not in {'deg', 'das', 'dtu', 'dexseq'}
            or not meta.get('native_background_provenance')
            or meta.get('report_completion_evidence') != 'completed_workflow_manifest'):
        return None
    selection = meta.get('das_selection') or {}
    policy = None
    if foreground == 'das':
        policy = selection.get('event_set')
        if policy not in {'significant', 'stable_significant'}:
            return None
        if policy == 'stable_significant':
            from katana_next.science.bulk import READ_SUPPORT_POLICY_VERSION
            if selection.get('read_support_policy') != READ_SUPPORT_POLICY_VERSION:
                return None
    return {'contrast_id': comparisons[0], 'foreground': foreground,
            'policy': policy, 'background': meta.get('background_source', foreground)}


def build_native_enrichment_facets(artifacts, destination):
    """Render full completed native ORA tables, not the old truncated plot data.

    Comparison, cell type, foreground and original run remain separate. GO and
    KEGG retain their original correction families; this never runs ORA or BH.
    """
    import hashlib
    from pathlib import Path
    import pandas as pd
    from .sources import read_table
    from .enrichment_facets import render_enrichment_facets, FACETS

    panels, audit, replaced = [], [], []
    for artifact in artifacts:
        scope = native_enrichment_scope(artifact)
        if not scope or artifact.get('entity_type') != 'gene_set' or artifact.get('kind') != 'table':
            continue
        source = Path(artifact['path'])
        frame = read_table(source)
        required = {'term', 'ontology_category', 'overlap_genes', 'selected_genes', 'fdr'}
        if not required <= set(frame):
            audit.append({'source': str(source), 'status': 'not_rendered', 'reason': 'Missing recorded ontology/count/FDR fields'})
            continue
        if 'comparison_id' in frame and set(frame.comparison_id.dropna().astype(str)) - {scope['contrast_id']}:
            raise ValueError('Native enrichment table contains a different comparison')
        meta = artifact.get('metadata', {})
        run = artifact.get('run_fingerprint') or meta.get('run_fingerprint', '')
        cell_type = meta.get('cell_type', 'all')
        key = (run, scope['contrast_id'], scope['foreground'], scope['policy'], cell_type)
        identity = 'ora-facets-' + hashlib.sha256((str(source.resolve()) + repr(key)).encode()).hexdigest()[:20]
        folder = Path(destination).resolve() / 'outputs/enrichment' / identity
        folder.mkdir(parents=True, exist_ok=True)
        data = frame.copy(deep=True)
        data['category'] = data.ontology_category.replace({'GO BP': 'BP', 'GO CC': 'CC', 'GO MF': 'MF'})
        classified = data.category.isin(FACETS)
        library_counts = []
        if 'library' in frame:
            for library, rows in frame.groupby('library', sort=True):
                q = pd.to_numeric(rows.fdr, errors='coerce')
                library_counts.append({'library': str(library), 'tests': len(rows),
                                       'significant': int(q.le(.05).sum()),
                                       'mapped_input_genes': int(pd.to_numeric(rows.selected_genes).max()),
                                       'mapped_background_genes': int(pd.to_numeric(rows.background_genes).max()) if 'background_genes' in rows else None})
        empty_messages = {}
        unknown_go = ~classified & data.get('library', pd.Series('', index=data.index)).astype(str).str.contains('GO', case=False)
        if unknown_go.any():
            for category in ('BP', 'CC', 'MF'):
                if not (data.category == category).any():
                    empty_messages[category] = 'GO category metadata unavailable\nSee complete GO tests'
        for category in FACETS:
            rows = data.loc[data.category.eq(category)]
            if len(rows) and pd.to_numeric(rows.selected_genes).eq(0).all():
                empty_messages[category] = 'No input genes mapped to this library'
        data = data.loc[classified].copy()
        data['Count'] = pd.to_numeric(data.overlap_genes)
        denominator = pd.to_numeric(data.selected_genes)
        # A zero-size foreground has no meaningful ratio; retain it in the full
        # table, rather than inventing a plotted zero or a significance result.
        data = data.loc[denominator.gt(0) & pd.to_numeric(data.fdr, errors='coerce').notna()].copy()
        data['GeneRatio'] = data.Count / pd.to_numeric(data.selected_genes)
        data['p.adjust'] = pd.to_numeric(data.fdr)
        label = ('RS-DAS' if scope['policy'] == 'stable_significant' else 'DAS') if scope['foreground'] == 'das' else scope['foreground'].upper()
        title = f"GO / KEGG enrichment | {label} | {scope['contrast_id']}"
        image = folder / 'facets.png'
        plotted = render_enrichment_facets(data, image, title=title,
                                          available_categories=set(data.category), empty_messages=empty_messages)
        plot_path = folder / 'plotted.tsv'
        plotted.to_csv(plot_path, sep='\t', index=False)
        links = []
        def link(role, path, label):
            links.append(dict(role=role, path=str(path), label=label,
                              contrast_id=scope['contrast_id'], policy=scope['policy'], unit='gene_set' if role in {'plot_data','full_results'} else 'gene'))
        link('plot_data', plot_path, 'Plotted terms, counts, adjusted P and hit genes')
        link('full_results', source, 'All native tests and original correction families')
        for sibling in artifacts:
            if native_enrichment_scope(sibling) != scope or sibling.get('kind') != 'table':
                continue
            sm = sibling.get('metadata', {})
            if (sibling.get('run_fingerprint') or sm.get('run_fingerprint', ''), sm.get('cell_type', 'all')) != (run, cell_type):
                continue
            role = {'gene_input_audit': 'foreground', 'background_audit': 'background'}.get(sibling.get('entity_type'))
            if role:
                link(role, sibling['path'], 'Input genes and selection flags' if role == 'foreground' else 'Detected/testable background and eligibility')
        panels.append(dict(id=identity, chapter=5, title=title, path=str(image),
                           contrast_id=scope['contrast_id'], policy=scope['policy'], cell_type=cell_type,
                           data_links=links, source_run=run,
                           metadata={'renderer':'go-kegg-facets-v1','foreground_source':scope['foreground'],
                                     'library_counts': library_counts, 'unclassified_terms': int((~classified).sum()),
                                     'facet_order':list(FACETS),'adjusted_p_recomputed':False},
                           caption='BP, CC, MF and KEGG; Gene Ratio = hits / input genes mapped to the library. Bubble area = Count; colour = original adjusted P. Up to five significant terms per facet; complete native tests, hit genes and correction families retained. No enrichment or statistical model was refitted. '
                           + ' '.join(f"{r['library']}: {r['significant']} significant / {r['tests']} recorded tests; {r['mapped_input_genes']} mapped input genes / {r['mapped_background_genes']} mapped background genes." for r in library_counts)
                           + (f" {int((~classified).sum())} tests lack recorded BP/CC/MF category metadata; their complete statistics remain linked. An empty category is not an unperformed GO analysis." if (~classified).any() else '')))
        replaced.append(key)
        audit.append(dict(source=str(source), status='rendered_completed_tests', eligible_tests=len(frame),
                          displayed_terms=len(plotted), unclassified_terms=int((~classified).sum()), new_inference=False))
    return dict(panels=panels, audit=audit, replaced=replaced)
