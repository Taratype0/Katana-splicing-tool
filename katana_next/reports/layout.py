"""Approved report layout; presentation only, no statistical interpretation."""
from pathlib import Path
from datetime import datetime, timezone
import json
import re
import zipfile
from .engine import CHAPTERS, esc, safe, chapter_for


def _scope(record):
    return ' data-contrast="'+esc(record.get('contrast_id', record.get('comparison_id', 'all')))+'" data-celltype="'+esc(record.get('cell_type', 'all'))+'" data-policy="'+esc(record.get('policy') or 'common')+'" data-comparisons="'+esc(json.dumps(record.get('comparison_ids') or []))+'"'


def _policy_control(identity):
    return '<label>DAS selection <select class="das-select" id="'+identity+'"><option value="stable_significant">RS-DAS (default)</option><option value="significant">Ordinary DAS</option></select></label>'


def _link(path, label):
    return '<a href="'+esc(path)+'">'+esc(label)+'</a>'


def _data_links(record, prefix=''):
    links = record.get('data_links', [])
    if not links:
        return '<span class="selection-note">Exact plot data not linked for this imported figure.</span>'
    return ''.join('<a data-data-role="'+esc(link['role'])+'" href="'+esc(prefix+link['path'])+'" title="'+esc('Unit: '+str(link['unit'])+'; comparison: '+str(link['contrast_id'])+'; policy: '+str(link['policy'] or 'fixed/common'))+'">'+esc(link['label'])+'</a>' for link in links)


def _measurement_note(record):
    scope = record.get('metadata', {}).get('das_scope', {})
    if scope.get('mode') == 'rs_required':
        return '<span class="selection-note">RS-DAS required; fixed recorded measurement scope.</span>'
    if scope.get('mode') == 'fixed_completed_result':
        name = {'stable_significant': 'RS-DAS', 'significant': 'ordinary DAS'}.get(scope.get('recorded_event_set'), 'original collection')
        return '<span class="selection-note">Fixed completed '+esc(name)+' scope; unchanged by the DAS selector.</span>'
    return ''


def _cell(value):
    if value is None: return 'NA'
    if isinstance(value, float): return f'{value:.5g}'
    return esc(value)


def _analysis(record):
    title = record.get('title', '')
    for name in ('DEXSeq', 'DTU', 'DEG'):
        if re.search(r'\b'+name+r'\b', title, re.I): return name
    if re.search(r'\b(?:RS-DAS|DAS|AS volcano)\b', title, re.I): return 'DAS'
    return 'integration'


def _ordered_panels(builder, chapter):
    records = [r for r in builder.panels if r['chapter'] == chapter and not r.get('output_only')]
    if chapter == 1:
        return sorted(records, key=lambda r: 0 if r['id'] == 'biological-design' else 1)
    if chapter == 2:
        def rank(record):
            title = record['title'].lower()
            if any(x in title for x in ('raw', 'read support', 'read-support', 'event composition', 'event types')): return 0
            if 'pca' in title: return 1
            if 'heatmap' in title: return 2
            return 3
        return sorted(records, key=rank)
    if chapter == 3:
        contrasts = {c['contrast_id']: i for i,c in enumerate(builder.source['project'].get('contrasts', []))}
        levels = {'DEG': 0, 'DAS': 1, 'DTU': 2, 'DEXSeq': 3, 'integration': 4}
        return sorted(records, key=lambda r: (contrasts.get(r.get('contrast_id'), 10**6), levels[_analysis(r)]))
    if chapter == 7:
        def evidence_layer(record):
            title = record['title'].lower()
            if 'motiftrace' in title:
                return 3 if any(x in title for x in ('window', 'cluster', 'statistics')) else 2
            return 1 if any(x in title for x in ('motif', 'clip', 'mapp')) else 0
        return sorted(records, key=evidence_layer)
    return records


def _top10(record, output):
    columns = record.get('columns', [])
    raw_rows = record.get('rows', [])
    rows = [r if isinstance(r, dict) else dict(zip(columns, r)) for r in raw_rows]
    if not columns and rows: columns = list(rows[0])
    visible = [c for c in columns if c in ('gene_symbol', 'gene_id', 'event_id', 'transcript_id', 'feature_id', 'featureID', 'event_type', 'log2fc', 'dpsi', 'fdr', 'deg_fdr', 'dtu_fdr', 'effect', 'native_effect', 'padj')]
    visible = visible[:7] or columns[:7]
    body = ''
    for rank, row in enumerate(rows[:10], 1):
        details = '<dl>'+''.join('<dt>'+esc(k)+'</dt><dd>'+_cell(row.get(k))+'</dd>' for k in columns)+'</dl>'
        body += '<tr><td>'+str(rank)+'</td>'+''.join('<td title="'+esc(row.get(k))+'"><span class="top10-cell-value">'+_cell(row.get(k))+'</span></td>' for k in visible)+'<td><details><summary>Details</summary>'+details+'</details></td></tr>'
    if not rows: body = '<tr><td colspan="20">No qualifying result for this selection. No nonsignificant result has been promoted to fill ten rows.</td></tr>'
    full = record.get('full_path', record.get('full_table', ''))
    if full:
        full = Path(full)
        full = full.relative_to(output).as_posix() if full.is_absolute() and full.is_relative_to(output) else str(full)
    return '<article class="exhibit top10 wide" id="'+esc(record['id'])+'"'+_scope(record)+'><div class="exhibit-head"><span class="exhibit-no">TOP 10</span><h3>'+esc(record.get('kind', record.get('analysis', 'Result')))+' | '+esc(record.get('contrast_id', record.get('comparison_id', '')) )+'</h3></div><p class="fig-notes">'+esc(record.get('caption', record.get('ranking', 'Native FDR ascending; absolute native effect descending when available; stable object identity.')))+' Unit: '+esc(record.get('unit', 'native result'))+'.</p>'+(_policy_control(safe(record['id'])+'-selection') if record.get('policy') else '')+'<div class="table-scroll"><table class="metadata"><tr><th>Rank</th>'+''.join('<th>'+esc(k)+'</th>' for k in visible)+'<th>All fields</th></tr>'+body+'</table></div>'+(_link(full, 'Complete result table') if full else '')+'</article>'


def _rbp_summary(record):
    rows = ''
    for row in record['rows']:
        genes = ', '.join(row['hit_genes'])
        rows += '<tr><td>'+esc(row['term'])+'</td><td>'+str(row['hits'])+'</td><td>'+_cell(row['fold_enrichment'])+'</td><td>'+_cell(row['pvalue'])+'</td><td>'+_cell(row['fdr'])+'</td><td><details><summary>All '+str(row['hits'])+' hit genes</summary><p>'+esc(genes)+'</p></details></td></tr>'
    if not rows:
        rows = '<tr><td colspan="6">No recorded RBP set passes FDR ≤ 0.05. All tested sets remain in the complete results.</td></tr>'
    families = '; '.join(dict.fromkeys(row['fdr_family'] for row in record['rows']))
    return '<article class="exhibit wide rbp-enrichment" id="'+esc(record['id'])+'"'+_scope(record)+'><div class="exhibit-head"><h3>'+esc(record.get('title', 'RBP identity / declared-set enrichment'))+' | '+esc(record['contrast_id'])+'</h3></div><p class="fig-notes">'+esc(str(len(record['rows']))+' significant / '+str(record['tested_sets'])+' tested sets. Foreground: '+record['foreground_source']+'. Background: '+record['background_source']+'. '+record['resource_semantics'])+'</p><div class="table-scroll"><table class="metadata"><tr><th>RBP set</th><th>Hit genes</th><th>'+esc(record.get('effect_label', 'Fold enrichment'))+'</th><th>P</th><th>FDR</th><th>Recorded contributors</th></tr>'+rows+'</table></div><p class="fig-notes">'+esc(families)+'</p>'+_link(record['full_path'], 'Complete native enrichment statistics')+_link(record.get('all_hits_path', record['full_path']), 'All tested sets and full hit genes')+'</article>'


def _motif_summary(record):
    status = record.get('status')
    heading = 'Significant motif sequences' if status == 'significant' else 'No significant motif' if status == 'no_significant_motif' else 'Motif tests not estimable'
    minimum = record.get('min_native_fdr')
    note = (' Minimum recorded test FDR: '+_cell(minimum)+'.') if minimum is not None else ''
    full = _link(record['ranking_path'], 'All models and original P / FDR') if record.get('ranking_path') else ''
    return '<article class="exhibit wide motif-significance" id="motif-status-'+esc(record['contrast_id'])+'"'+_scope(record)+'><h3>MotifTrace: '+heading+'</h3><p>'+esc(record['notice'])+note+'</p><p class="fig-notes">'+esc(record['native_test_family'])+'</p>'+full+'</article>'


def _network_browser(record):
    return '<article class="exhibit wide network-browser" id="'+esc(record['id'])+'"'+_scope(record)+'><div class="exhibit-head"><h3>Module network and complete gene lists | '+esc(record['contrast_id'])+'</h3></div><p class="fig-notes">Select a module below. The graph uses saved module membership and native adjacency; display-node and edge limits are stated in the browser. Gene symbols, IDs and all event memberships are searchable and downloadable.</p>'+_link(record['source'], 'Open network and module genes in a full window')+'<iframe loading="lazy" data-src="'+esc(record['source'])+'" title="Module network and gene lists" class="network-frame"></iframe></article>'


def render(builder):
    source, output = builder.source, builder.output
    project = source['project']
    sections = []
    for number, (title, note) in CHAPTERS.items():
        blocks = []
        if number == 6:
            blocks.extend(_network_browser(r) for r in builder.network_browsers)
        if number == 4 and len(project.get('contrasts', [])) == 1:
            blocks.append('<div class="status-panel"><strong>Cross-comparison: not applicable</strong>'
                          'One comparison was selected. Overlap and between-comparison direction analyses require '
                          'at least two comparisons. Any separately declared interaction model is reported independently.</div>')
        if number == 7:
            blocks.extend(_rbp_summary(r) for r in builder.enrichment_summaries)
            blocks.extend(_motif_summary(r) for r in builder.motif_summaries)
        chapter_items = _ordered_panels(builder, number)
        inserted = set()
        for index, item in enumerate(chapter_items, 1):
            links = _link(item['image'], 'Figure')
            if item.get('svg'): links += _link(item['svg'], 'SVG')
            links += _data_links(item)
            links += _measurement_note(item)
            if item.get('selection_scope'): links += '<span class="selection-note">Recorded result scope; see provenance</span>'
            blocks.append('<figure class="exhibit wide" id="'+esc(item['id'])+'"'+_scope(item)+'><div class="exhibit-head"><span class="exhibit-no">'+str(number)+'.'+str(index)+'</span><h3>'+esc(item['title'])+'</h3></div><button class="plot-button"><img loading="lazy" src="'+esc(item['image'])+'" alt="'+esc(item['title'])+'"></button><figcaption><p class="fig-notes">'+esc(item.get('caption', ''))+'</p><div class="fig-links">'+links+'</div></figcaption></figure>')
            if number == 3:
                group = (item.get('contrast_id'), _analysis(item))
                following = chapter_items[index] if index < len(chapter_items) else None
                if following is None or (following.get('contrast_id'), _analysis(following)) != group:
                    for record in builder.top10:
                        if (record.get('contrast_id', record.get('comparison_id')), record.get('analysis', record.get('kind'))) == group:
                            blocks.append(_top10(record, output)); inserted.add(record['id'])
        if number == 3:
            blocks.extend(_top10(r, output) for r in builder.top10 if r['id'] not in inserted)
        missing = [r for r in builder.exclusions if r.get('chapter') == number]
        # A source's absent component is not absent if a verified supplemental
        # result supplies that exact component; unrelated optional analyses stay missing.
        available_analyses = {p.get('metadata', {}).get('analysis_id') for p in builder.panels}
        if number == 6 and any(str(x).startswith('programs.') for x in available_analyses):
            missing = [r for r in missing if 'network' not in str(r.get('id', '')).lower()]
        if not blocks:
            blocks.append('<div class="status-panel"><strong>No completed result available</strong>This component is not computed. No synthetic data or other study has been substituted.</div>')
        seen_missing = set()
        for entry in missing:
            text = entry.get('reason') or entry.get('caption') or 'No completed result supplied.'
            key = (entry.get('title'), text)
            if key in seen_missing: continue
            seen_missing.add(key)
            blocks.append('<div class="status-panel"><strong>'+esc(entry.get('title', 'Not computed'))+'</strong>'+esc(text)+'</div>')
        controls = '<div class="local-controls">'+_policy_control('chapter3-das-selection')+'</div>' if number == 3 else ''
        sections.append('<section class="chapter" id="chapter-'+str(number)+'"><div class="chapter-head"><span class="chapter-number">'+str(number).zfill(2)+'</span><h2>'+esc(title)+'</h2></div><p class="chapter-note">'+esc(note)+'</p>'+controls+'<p>'+_link('outputs/index.html#chapter-'+str(number), 'Complete data, figures and interactive browsers for this section')+'</p><div class="figure-grid">'+''.join(blocks)+'</div><p class="empty-scope">No result for this selected view.</p></section>')
    contrasts = project.get('contrasts', [])
    sample_rows = ''.join('<tr><td>'+esc(s['sample_id'])+'</td><td>'+esc(s.get('biological_sample_id', s['sample_id']))+'</td><td>'+esc(s.get('condition', ''))+'</td></tr>' for s in project.get('samples', []))
    contrast_rows = ''.join('<tr><td>'+esc(c.get('label', c['contrast_id']))+'</td><td>'+esc(c['numerator'].get('label', ''))+'</td><td>'+esc(c['denominator'].get('label', ''))+'</td><td>'+str(len(c['numerator'].get('sample_ids', [])))+' / '+str(len(c['denominator'].get('sample_ids', [])))+'</td></tr>' for c in contrasts)
    annotation = project.get('annotation', {})
    design = '<h3>Samples and declared comparisons</h3><p>Effects follow the explicitly recorded numerator / denominator. Independent biological units are recorded per sample; no pairing is inferred from replicate numbers.</p><details><summary>All biological samples ('+str(len(project.get('samples', [])))+')</summary><table class="metadata"><tr><th>Sample</th><th>Biological unit</th><th>Condition</th></tr>'+sample_rows+'</table></details><details open><summary>Declared comparisons ('+str(len(contrasts))+')</summary><div class="table-scroll"><table class="metadata"><tr><th>Comparison</th><th>Numerator</th><th>Denominator</th><th>n / n</th></tr>'+contrast_rows+'</table></div></details><p>Reference: '+esc(annotation.get('genome_build', 'Not recorded'))+' · '+esc(project.get('metadata', {}).get('annotation_release', 'Annotation release not recorded'))+'. Reference file identities and actual method parameters are in provenance.</p>'
    sections[0] = sections[0].replace('<div class="figure-grid">', design+'<div class="figure-grid">', 1)
    toc = ''.join('<a href="#chapter-'+str(n)+'"><b>'+str(n).zfill(2)+'</b>'+esc(value[0])+'</a>' for n, value in CHAPTERS.items())
    options = ''.join('<option value="'+esc(c['contrast_id'])+'">'+esc(c.get('label', c['contrast_id']))+'</option>' for c in contrasts)
    footnote = project.get('metadata', {}).get('gene_filter_policy', {}).get('footnote_en', 'Actual source filtering and fitting scope are retained in provenance. Missing filtering records are not inferred.')
    appendix = '<section class="chapter" id="appendix"><div class="chapter-head"><span class="chapter-number">A</span><h2>Methods, parameters and complete results</h2></div><p>'+esc(footnote)+'</p><p>RS-DAS: native caller FDR/effect criteria plus IJC + SJC ≥ 10 in every declared biological sample in both contrast arms. Ordinary DAS: the same native caller criteria without the read-support gate. Original caller P/FDR and effects are preserved. Historical analyses with another measurement policy are retained as fixed completed outputs and never silently relabelled.</p><p>Top 10: native FDR ascending, absolute native effect descending when available, stable object identity. Missing effect sizes remain unavailable; transcript, exon and gene units are separate.</p><p>Networks require their declared measurable PSI population, not just significant DAS. Motif occurrences, group enrichment, window statistics and CLIP binding evidence retain separate meanings.</p><p>'+_link('outputs/index.html', 'Complete results')+' · '+_link('provenance.json', 'Actual parameters and provenance')+' · '+_link('REPORT.zip', 'Download report package')+'</p><p>Katana Splicing Tool · '+_link('https://github.com/Taratype0/Katana-splicing-tool', 'Software and citation information')+'. No publication DOI is claimed.</p></section>'
    page = '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+esc(project.get('name', 'Katana report'))+'</title><link rel="stylesheet" href="assets/report.css"><header class="topbar"><a class="wordmark" href="REPORT.html">KATANA</a><small>RNA SPLICING ANALYSIS</small><div class="top-actions">'+_link('outputs/index.html', 'Output files')+_link('REPORT.zip', 'Download package')+'<button id="printReport">Print / PDF</button></div></header><section class="cover"><p class="eyebrow">COMPLETED RESULTS · '+esc(source['mode'].upper())+'</p><h1>'+esc(project.get('name', 'Katana report'))+'</h1><p>'+str(len(project.get('samples', [])))+' biological samples · '+str(len(contrasts))+' declared comparisons · '+esc(annotation.get('genome_build', 'Reference recorded in source'))+'</p></section><div class="shell"><nav class="toc"><p class="toc-title">REPORT CONTENTS</p>'+toc+'<a href="#appendix">Methods and outputs</a></nav><main><div class="scopebar"><label>Comparison <select id="contrastSelect">'+options+'</select></label>'+_policy_control('report-das-selection')+'<button id="resetScope">Reset view</button><span id="scopeSummary" class="scope-summary"></span></div><p class="policy-explanation">Switching DAS selection changes its volcano, Top 10 and matched downstream views together. Fixed completed imports retain their original input population.</p>'+''.join(sections)+appendix+'</main></div><dialog id="figureModal"><div class="modal-top"><span id="modalTitle"></span><div class="modal-actions"><button id="fitFigure">Fit</button><button id="actualFigure">100%</button><a id="downloadFigure" download>Figure</a><button id="closeFigure">Close</button></div></div><div class="modal-scroll" id="modalScroll"><img id="modalImage" alt=""></div></dialog><script id="report-settings" type="application/json">'+json.dumps({'event_set': builder.event_set, 'initial_comparison': contrasts[0]['contrast_id'] if contrasts else 'all'})+'</script><script src="assets/report.js"></script></html>'
    (output/'REPORT.html').write_text(page, encoding='utf-8')
    _outputs(builder)
    payload = {'schema': 'katana.approved_report/v1', 'specification_version': '1.1', 'project_id': project.get('project_id'),
               'project_name': project.get('name'), 'generated_utc': datetime.now(timezone.utc).isoformat(), 'event_set': builder.event_set,
               'samples': len(project.get('samples', [])), 'contrasts': [c['contrast_id'] for c in contrasts], 'panels': builder.panels,
               'top10': builder.top10, 'policy_data': builder.policy_data, 'complete_results': builder.tables,
               'rbp_enrichment_summaries': builder.enrichment_summaries,
               'source_files': builder.audit, 'unavailable': builder.exclusions, 'no_upstream_rerun': True,
               'no_source_statistics_refit': True, 'dummy': bool(source.get('dummy', False))}
    from .presentation import portable_display_project
    payload['display_project'] = portable_display_project(project)
    payload['mode'] = source['mode']
    payload['motif_summaries'] = builder.motif_summaries
    payload['network_browsers'] = builder.network_browsers
    payload['presentation'] = {'default_comparison': contrasts[0]['contrast_id'] if contrasts else 'all',
                               'view': 'one_comparison', 'network_dendrogram_in_main': False,
                               'motif_sequence_rule': 'native selected-test FDR <= 0.05; at most 10 distinct models'}
    if hasattr(builder, 'release_status'): payload['release_status'] = builder.release_status
    payload['figure_data_links'] = builder.figure_data_coverage()
    payload['report_selection'] = source.get('report_selection', {'applied': False})
    (output/'report-manifest.json').write_text(json.dumps(payload, indent=2, default=str)+'\n', encoding='utf-8')
    sources = []
    for group in [source, *source.get('supplemental', [])]:
        sources.append({'project_id': group['project'].get('project_id'), 'metadata': group['project'].get('metadata', {}),
                        'parameters': [{'analysis_id': r.get('analysis_id'), 'artifact_id': r.get('artifact_id'), 'parameters': r.get('parameters', {}), 'method_metadata': r.get('metadata', {})} for r in group['artifacts']]})
    if hasattr(builder, 'preserved_source_projects'): sources = builder.preserved_source_projects
    provenance = {'schema': 'katana.report_provenance/v1', 'source_files': builder.audit, 'source_projects': sources,
                  'report_parameters': {'event_set': builder.event_set}, 'derived_views': builder.derived_audits,
                  'no_upstream_rerun': True, 'no_source_statistics_refit': True}
    provenance['figure_data_links'] = builder.figure_data_coverage()
    provenance['report_selection'] = source.get('report_selection', {'applied': False})
    (output/'provenance.json').write_text(json.dumps(provenance, indent=2, default=str)+'\n', encoding='utf-8')
    if builder.create_archive:
        with zipfile.ZipFile(output/'REPORT.zip', 'w', zipfile.ZIP_DEFLATED, compresslevel=4, strict_timestamps=False) as archive:
            for path in sorted(output.rglob('*')):
                if path.is_file() and path.name != 'REPORT.zip': archive.write(path, path.relative_to(output).as_posix())
    return payload


def _outputs(builder):
    output = builder.output
    groups = []
    for number, (title, _) in CHAPTERS.items():
        tables = [t for t in builder.tables if t.get('chapter') == number or (not t.get('chapter') and chapter_for(t.get('analysis_id', '')) == number)]
        rows = ''.join('<tr><td>'+esc(t['name'])+'</td><td>'+esc(t.get('rows', ''))+'</td><td>'+_link('../'+t.get('view', t['source']), 'Full table / browser')+' '+_link('../'+t['source'], 'Complete source')+'</td></tr>' for t in tables)
        figures = ''.join('<details><summary>'+esc(p['title'])+'</summary><p>'+esc(p.get('caption', ''))+'</p>'+_link('../'+p['image'], 'Figure')+_data_links(p, '../')+'<img loading="lazy" style="max-width:100%" src="../'+esc(p['image'])+'" alt="'+esc(p['title'])+'"></details>' for p in builder.panels if p['chapter'] == number)
        groups.append('<section id="chapter-'+str(number)+'"><h2>'+str(number)+'. '+esc(title)+'</h2><table class="metadata"><tr><th>Complete result</th><th>Rows</th><th>Access</th></tr>'+rows+'</table><h3>Complete figure gallery</h3>'+figures+'</section>')
    (output/'outputs').mkdir(exist_ok=True)
    (output/'outputs/index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Complete outputs · Katana</title><link rel="stylesheet" href="../assets/report.css"><main>'+_link('../REPORT.html', 'Report')+'<h1>Complete outputs</h1><p>Source tables retain every row, native P/FDR and effects. Top 10 and overview panels are summaries. Complete long tables, matrix browsers and figure galleries live here outside the report body.</p>'+''.join(groups)+'</main></html>', encoding='utf-8')
