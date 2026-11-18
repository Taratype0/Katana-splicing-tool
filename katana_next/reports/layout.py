"""Approved executive consulting report layout; presentation only, zero AI speculation."""
from pathlib import Path
from datetime import datetime, timezone
import json
import re
import zipfile
from .engine import CHAPTERS, CHAPTER_NUMBERS, esc, safe, chapter_for

_NOISY_KEY_RE = re.compile(
    r"(?:sha256|fingerprint|_path$|source_file|portable_root|local-reference|cluster|klemming|work_dir|^row_id$|^source_row_id$)",
    re.I,
)
_PATH_VAL_RE = re.compile(r"^(?:/cfs/|/mnt/|[A-Za-z]:[\\/]|local-reference:)")


def _scope(record):
    return (
        ' data-contrast="' + esc(record.get('contrast_id', record.get('comparison_id', 'all')))
        + '" data-celltype="' + esc(record.get('cell_type', 'all'))
        + '" data-policy="' + esc(record.get('policy') or 'common')
        + '" data-comparisons="' + esc(json.dumps(record.get('comparison_ids') or [])) + '"'
    )


def _policy_control(identity):
    return (
        '<label>DAS selection <select class="das-select" id="' + identity + '">'
        '<option value="stable_significant">PS-DAS · Precision-Supported (default)</option>'
        '<option value="significant">Ordinary DAS · Unfiltered</option>'
        '</select></label>'
    )


def _link(path, label):
    return '<a href="' + esc(str(path).replace('\\', '/')) + '">' + esc(label) + '</a>'


def _data_links(record, prefix=''):
    links = record.get('data_links', [])
    if not links:
        return '<span class="selection-note">Exact plot data not linked for this imported figure.</span>'
    return ''.join(
        '<a data-data-role="' + esc(link['role']) + '" href="' + esc((prefix + str(link['path'])).replace('\\', '/')) + '" title="'
        + esc('Unit: ' + str(link['unit']) + '; comparison: ' + str(link['contrast_id']) + '; policy: ' + str(link['policy'] or 'fixed/common'))
        + '">' + esc(link['label']) + '</a>'
        for link in links
    )


def _measurement_note(record):
    scope = record.get('metadata', {}).get('das_scope', {})
    if scope.get('mode') == 'rs_required':
        return '<span class="selection-note">PS-DAS required; fixed recorded precision scope.</span>'
    if scope.get('mode') == 'fixed_completed_result':
        name = {'stable_significant': 'PS-DAS', 'significant': 'Ordinary DAS'}.get(scope.get('recorded_event_set'), 'original collection')
        return '<span class="selection-note">Fixed completed ' + esc(name) + ' scope; unchanged by the DAS selector.</span>'
    return ''


def _cell(value):
    if value is None:
        return 'NA'
    if isinstance(value, float):
        return f'{value:.5g}'
    return esc(value)



def _image_dimensions(output, image):
    path = Path(output) / str(image).replace('\\', '/')
    if not path.is_file() or path.suffix.lower() not in {'.png', '.jpg', '.jpeg', '.webp'}:
        return ''
    try:
        from PIL import Image
        with Image.open(path) as asset:
            width, height = asset.size
        return f' width="{width}" height="{height}"'
    except Exception:
        return ''


def _analysis(record):
    title = record.get('title', '')
    for name in ('DEXSeq', 'DTU', 'DEG'):
        if re.search(r'\b' + name + r'\b', title, re.I):
            return name
    if re.search(r'\b(?:PS-DAS|RS-DAS|DAS|AS volcano)\b', title, re.I):
        return 'DAS'
    return 'integration'


def _ordered_panels(builder, chapter):
    records = [r for r in builder.panels if r['chapter'] == chapter and not r.get('output_only')]
    if chapter == 1:
        return sorted(records, key=lambda r: 0 if r['id'] == 'biological-design' else 1)
    if chapter == 2:
        def rank(record):
            title = record['title'].lower()
            if any(x in title for x in ('raw', 'read support', 'read-support', 'event composition', 'event types', 'utr', 'regional')):
                return 0
            if 'pca' in title or 'umap' in title or 'manifold' in title:
                return 1
            if 'heatmap' in title:
                return 2
            return 3
        return sorted(records, key=rank)
    if chapter == 3:
        contrasts = {c['contrast_id']: i for i, c in enumerate(builder.source['project'].get('contrasts', []))}
        levels = {'DEG': 0, 'DAS': 1, 'DTU': 2, 'DEXSeq': 3, 'integration': 4}
        return sorted(records, key=lambda r: (contrasts.get(r.get('contrast_id'), 10**6), levels[_analysis(r)]))
    if chapter == 7:
        def evidence_layer(record):
            if record.get('metadata', {}).get('motif_sequence_renderer'):
                return -4 + record['metadata'].get('sequence_priority', 3)
            title = record['title'].lower()
            if 'motiftrace' in title or 'splicetrace' in title:
                return 3 if any(x in title for x in ('window', 'cluster', 'statistics')) else 2
            return 1 if any(x in title for x in ('motif', 'clip', 'mapp')) else 0
        return sorted(records, key=evidence_layer)
    return records


def _clean_detail_fields(columns, row):
    """Filter out cluster path / hash noise from inline Top-10 detail drawers; full tables retain every byte."""
    items = []
    for k in columns:
        if _NOISY_KEY_RE.search(str(k)):
            continue
        v = row.get(k)
        if v is None or str(v).strip() == '':
            continue
        if isinstance(v, str) and (_PATH_VAL_RE.match(v) or len(v) > 180):
            continue
        items.append((k, v))
    return items[:24]


def _top10(record, output):
    columns = record.get('columns', [])
    raw_rows = record.get('rows', [])
    rows = [r if isinstance(r, dict) else dict(zip(columns, r)) for r in raw_rows]
    if not columns and rows:
        columns = list(rows[0])
    preferred = (
        'gene_symbol', 'gene_name', 'gene_id', 'event_id', 'transcript_id',
        'feature_id', 'featureID', 'event_type', 'log2fc', 'dpsi',
        'fdr', 'deg_fdr', 'dtu_fdr', 'effect', 'native_effect', 'padj',
        'module', 'consequence_class', 'domain_impact',
    )
    visible = [c for c in preferred if c in columns][:8] or [c for c in columns if not _NOISY_KEY_RE.search(str(c))][:7]
    num_cols = {'log2fc', 'dpsi', 'fdr', 'deg_fdr', 'dtu_fdr', 'effect', 'native_effect', 'padj'}
    body = ''
    for rank, row in enumerate(rows[:10], 1):
        clean_pairs = _clean_detail_fields(columns, row)
        details = '<dl>' + ''.join('<dt>' + esc(k) + '</dt><dd>' + _cell(v) + '</dd>' for k, v in clean_pairs) + '</dl>'
        cells = ''
        for k in visible:
            cls = ' class="num"' if k in num_cols else ''
            cells += f'<td{cls} title="{esc(row.get(k))}"><span class="top10-cell-value">{_cell(row.get(k))}</span></td>'
        body += f'<tr><td class="num">{rank}</td>{cells}<td><details><summary>Metrics</summary>{details}</details></td></tr>'
    if not rows:
        body = '<tr><td colspan="20">No qualifying result for this selection. No nonsignificant result has been promoted to fill ten rows.</td></tr>'
    full = record.get('full_path', record.get('full_table', ''))
    if full:
        full = Path(full)
        full = full.relative_to(output).as_posix() if full.is_absolute() and full.is_relative_to(output) else full.as_posix()
    headers = ''.join(f'<th{" class=\"num\"" if k in num_cols else ""}>{esc(k)}</th>' for k in visible)
    return (
        '<article class="exhibit top10 wide" id="' + esc(record['id']) + '"' + _scope(record) + '>'
        '<div class="exhibit-head"><span class="exhibit-no">EXECUTIVE TOP 10</span>'
        '<h3>' + esc(record.get('kind', record.get('analysis', 'Result'))) + ' — ' + esc(record.get('contrast_id', record.get('comparison_id', ''))) + '</h3></div>'
        '<p class="fig-notes">' + esc(record.get('caption', record.get('ranking', 'Native FDR ascending; absolute native effect descending when available; stable object identity.'))) + ' Unit: ' + esc(record.get('unit', 'native result')) + '.</p>'
        + (_policy_control(safe(record['id']) + '-selection') if record.get('policy') else '')
        + '<div class="table-scroll"><table class="metadata"><thead><tr><th class="num">Rank</th>' + headers + '<th>Record</th></tr></thead><tbody>' + body + '</tbody></table></div>'
        + ('<div class="fig-links">' + _link(full, 'Complete lossless result table') + '</div>' if full else '')
        + '</article>'
    )


def _rbp_summary(record):
    rows = ''
    for row in record['rows']:
        genes = ', '.join(row['hit_genes'])
        rows += '<tr><td>' + esc(row['term']) + '</td><td class="num">' + str(row['hits']) + '</td><td class="num">' + _cell(row['fold_enrichment']) + '</td><td class="num">' + _cell(row['pvalue']) + '</td><td class="num">' + _cell(row['fdr']) + '</td><td><details><summary>All ' + str(row['hits']) + ' hit genes</summary><p>' + esc(genes) + '</p></details></td></tr>'
    if not rows:
        rows = '<tr><td colspan="6">No recorded RBP set passes FDR ≤ 0.05. All tested sets remain in the complete results.</td></tr>'
    families = '; '.join(dict.fromkeys(row['fdr_family'] for row in record['rows']))
    return '<article class="exhibit wide rbp-enrichment" id="' + esc(record['id']) + '"' + _scope(record) + '><div class="exhibit-head"><h3>' + esc(record.get('title', 'RBP identity / declared-set enrichment')) + ' | ' + esc(record['contrast_id']) + '</h3></div><p class="fig-notes">' + esc(str(len(record['rows'])) + ' significant / ' + str(record['tested_sets']) + ' tested sets. Foreground: ' + record['foreground_source'] + '. Background: ' + record['background_source'] + '. ' + record['resource_semantics']) + '</p><div class="table-scroll"><table class="metadata"><thead><tr><th>RBP set</th><th class="num">Hit genes</th><th class="num">' + esc(record.get('effect_label', 'Fold enrichment')) + '</th><th class="num">P</th><th class="num">FDR</th><th>Recorded contributors</th></tr></thead><tbody>' + rows + '</tbody></table></div><p class="fig-notes">' + esc(families) + '</p><div class="fig-links">' + _link(record['full_path'], 'Complete native enrichment statistics') + _link(record.get('all_hits_path', record['full_path']), 'All tested sets and full hit genes') + '</div></article>'


def _motif_summary(record):
    if record.get('sequence_report'):
        from .motif_sequences import summary_html
        return '<article class="exhibit wide motif-significance" id="motif-status-' + esc(record.get('contrast_id', 'all')) + '"' + _scope(record) + '>' + summary_html(record['sequence_report']) + '</article>'
    status = record.get('status')
    heading = 'Significant motif sequences' if status == 'significant' else 'No significant motif' if status == 'no_significant_motif' else 'Motif tests not estimable'
    minimum = record.get('min_native_fdr')
    note = (' Minimum recorded test FDR: ' + _cell(minimum) + '.') if minimum is not None else ''
    full = _link(record['ranking_path'], 'All models and original P / FDR') if record.get('ranking_path') else ''
    notice_txt = record.get('notice', record.get('summary', 'Completed MotifTrace sequence summary; see complete output tables.'))
    fam_txt = record.get('native_test_family', '')
    return '<article class="exhibit wide motif-significance" id="motif-status-' + esc(record.get('contrast_id', 'all')) + '"' + _scope(record) + '><h3>SpliceTrace / MotifTrace: ' + heading + '</h3><p>' + esc(notice_txt) + note + '</p>' + ('<p class="fig-notes">' + esc(fam_txt) + '</p>' if fam_txt else '') + full + '</article>'


def _network_browser(record):
    if record.get('whole_study'):
        return '<article class="exhibit wide network-browser" id="' + esc(record['id']) + '"' + _scope(record) + '><div class="exhibit-head"><h3>CoSplice / SpliceNet | all study samples and conditions</h3></div><p class="fig-notes">One shared sample-level network across declared biological samples. Changing the comparison selector does not refit or subset this network. Each module has its own genes, events, sample response and enrichment.</p><div class="fig-links">' + _link(record['source'], 'Open whole-study network report') + '</div><iframe loading="lazy" data-src="' + esc(str(record['source']).replace('\\', '/')) + '" title="CoSplice whole-study modules" class="network-frame"></iframe></article>'
    return '<article class="exhibit wide network-browser" id="' + esc(record['id']) + '"' + _scope(record) + '><div class="exhibit-head"><h3>Module network and complete gene lists | ' + esc(record['contrast_id']) + '</h3></div><p class="fig-notes">Select a module below. The graph uses saved module membership and native adjacency; display-node and edge limits are stated in the browser. Gene symbols, IDs and all event memberships are searchable and downloadable.</p><div class="fig-links">' + _link(record['source'], 'Open network and module genes in a full window') + '</div><iframe loading="lazy" data-src="' + esc(str(record['source']).replace('\\', '/')) + '" title="Module network and gene lists" class="network-frame"></iframe></article>'


def render(builder):
    source, output = builder.source, builder.output
    project = source['project']
    from katana_next.flow.selection import active_chapters
    active_set = set(active_chapters(builder.analyses))
    # Only render Chapter 8 if integrative panels or tables were actually supplied
    has_ch8 = any(p.get('chapter') == 8 for p in builder.panels) or any(t.get('chapter') == 8 for t in builder.tables) or any(r.get('chapter') == 8 for r in builder.top10)
    if not has_ch8:
        active_set.discard(8)
    chapters = {n: v for n, v in CHAPTERS.items() if n in active_set}

    sections = []
    for number, (title, note) in chapters.items():
        blocks = []
        if number == 6:
            blocks.extend(_network_browser(r) for r in builder.network_browsers)
        if number == 4 and len(project.get('contrasts', [])) == 1 and not any(p.get('chapter') == 4 for p in builder.panels):
            blocks.append('<div class="status-panel"><strong>Between-comparison overlap: single-comparison study</strong>'
                          'One comparison was declared. Multi-comparison intersection requires at least two comparisons.</div>')
        if number == 7:
            blocks.extend(_rbp_summary(r) for r in builder.enrichment_summaries)
            blocks.extend(_motif_summary(r) for r in builder.motif_summaries)
        chapter_items = _ordered_panels(builder, number)
        inserted = set()
        display_number = CHAPTER_NUMBERS.get(number, number)
        for index, item in enumerate(chapter_items, 1):
            if item.get('metadata', {}).get('native_multiqc'):
                blocks.append(
                    '<article class="exhibit wide native-multiqc" id="' + esc(item['id']) + '">'
                    '<div class="exhibit-head"><h3>MultiQC sequencing quality</h3></div><p class="fig-notes">' + esc(item.get('caption', '')) + '</p>'
                    '<div class="fig-links">' + _link(item['image'], 'Open full MultiQC report') + _data_links(item) + '</div>'
                    '<iframe loading="lazy" data-src="' + esc(str(item['image']).replace('\\', '/')) + '" title="MultiQC sequencing quality" class="network-frame"></iframe></article>'
                )
                continue
            links = _link(item['image'], 'Figure (PNG)')
            if item.get('svg'):
                links += _link(item['svg'], 'Vector (SVG)')
            links += _data_links(item)
            links += _measurement_note(item)
            sequence_details = ''
            if item.get('metadata', {}).get('sequence_evidence'):
                from .motif_sequences import card_html
                sequence_details = card_html(item['metadata']['sequence_evidence'])
            if item.get('selection_scope'):
                links += '<span class="selection-note">Recorded result scope; see provenance</span>'
            subtitle_html = '<p class="exhibit-sub">' + esc(item['subtitle']) + '</p>' if item.get('subtitle') else ''
            img_pos = str(item['image']).replace('\\', '/')
            blocks.append(
                '<figure class="exhibit wide" id="' + esc(item['id']) + '"' + _scope(item) + '>'
                '<div class="exhibit-head"><span class="exhibit-no">EXHIBIT ' + str(display_number) + '.' + str(index) + '</span>'
                '<h3>' + esc(item['title']) + '</h3>' + subtitle_html + '</div>'
                '<button class="plot-button"><img loading="lazy"' + _image_dimensions(output, img_pos) + ' src="' + esc(img_pos) + '" alt="' + esc(item['title']) + '"></button>'
                '<figcaption><p class="fig-notes">' + esc(item.get('caption', '')) + '</p>' + sequence_details + '<div class="fig-links">' + links + '</div></figcaption></figure>'
            )
            if number == 3:
                group = (item.get('contrast_id'), _analysis(item))
                following = chapter_items[index] if index < len(chapter_items) else None
                if following is None or (following.get('contrast_id'), _analysis(following)) != group:
                    for record in builder.top10:
                        if record.get('chapter', 3) == 3 and (record.get('contrast_id', record.get('comparison_id')), record.get('analysis', record.get('kind'))) == group:
                            blocks.append(_top10(record, output))
                            inserted.add(record['id'])
        if number in (3, 8):
            blocks.extend(_top10(r, output) for r in builder.top10 if r.get('chapter', 3) == number and r['id'] not in inserted)

        missing = [r for r in builder.exclusions if r.get('chapter') == number]
        available_analyses = {p.get('metadata', {}).get('analysis_id') for p in builder.panels}
        if number == 6 and any(str(x).startswith('programs.') for x in available_analyses):
            missing = [r for r in missing if 'network' not in str(r.get('id', '')).lower()]
        if not blocks:
            enrichment_results = number == 5 and any(
                (r.get('chapter') or chapter_for(r.get('analysis_id', ''))) == 5
                and r.get('analysis_id') in {'prioritize.enrichment', 'interpretation.joint_pathways'}
                and {'term', 'pvalue', 'fdr'} <= set(r.get('columns', [])) for r in builder.tables
            )
            if enrichment_results:
                blocks.append('<div class="status-panel"><strong>No GO/KEGG panel for this selection</strong>'
                              'No GO/KEGG panel is available for this selection. See complete output tables '
                              'for recorded enrichment results, input genes and backgrounds.</div>')
            else:
                blocks.append('<div class="status-panel"><strong>No completed result available</strong>This component is not computed. No synthetic data or other study has been substituted.</div>')
        seen_missing = set()
        for entry in missing:
            text = entry.get('reason') or entry.get('caption') or 'No completed result supplied.'
            key = (entry.get('title'), text)
            if key in seen_missing:
                continue
            seen_missing.add(key)
            blocks.append('<div class="status-panel"><strong>' + esc(entry.get('title', 'Not computed')) + '</strong>' + esc(text) + '</div>')
        controls = '<div class="local-controls">' + _policy_control('chapter3-das-selection') + '</div>' if number == 3 and source['mode'] != 'single-cell' else ''
        callout = getattr(builder, 'chapter_callouts', {}).get(number, '')
        callout_html = '<div class="quant-callout">' + callout + '</div>' if callout else ''
        sections.append(
            '<section class="chapter" id="chapter-' + str(number) + '">'
            '<div class="chapter-head"><span class="chapter-number">' + str(display_number).zfill(2) + '</span><h2>' + esc(title) + '</h2></div>'
            '<p class="chapter-note">' + esc(note) + '</p>'
            + callout_html
            + controls
            + '<p class="fig-notes">' + _link('outputs/index.html#chapter-' + str(number), 'Complete data, figures and offline browsers for Chapter ' + str(display_number).zfill(2) + ' →') + '</p>'
            + '<div class="figure-grid">' + ''.join(blocks) + '</div>'
            + '<p class="empty-scope">No exhibit matches the active comparison filter in this section.</p></section>'
        )

    contrasts = project.get('contrasts', [])
    is_sc = str(source.get('mode', '')).lower().startswith('single')
    unit_label = 'Cells / Biological Donors' if is_sc else 'Biological Samples'
    sample_rows = ''.join(
        '<tr><td>' + esc(s['sample_id']) + '</td><td>' + esc(s.get('biological_sample_id', s['sample_id'])) + '</td><td>' + esc(s.get('condition', '')) + '</td></tr>'
        for s in project.get('samples', [])
    )
    contrast_rows = ''.join(
        '<tr><td><strong>' + esc(c.get('label', c['contrast_id'])) + '</strong></td><td>' + esc(c['numerator'].get('label', '')) + '</td><td>' + esc(c['denominator'].get('label', '')) + '</td><td class="num">' + str(len(c['numerator'].get('sample_ids', []))) + ' vs ' + str(len(c['denominator'].get('sample_ids', []))) + '</td></tr>'
        for c in contrasts
    )
    annotation = project.get('annotation', {})
    replicate_note = (
        'Single-cell / Smart-seq3 study: biological inference respects donor / biological replicate structure and Kish effective sample size (cells are not treated as independent biological organisms).'
        if is_sc else
        'Effects follow the explicitly declared numerator / denominator. Independent biological units are recorded per sample; no pairing is inferred from sample names.'
    )
    design = (
        '<h3>Declared Study Architecture &amp; Comparisons</h3>'
        '<p class="fig-notes">' + esc(replicate_note) + '</p>'
        '<details><summary>All recorded ' + unit_label.lower() + ' (' + str(len(project.get('samples', []))) + ')</summary>'
        '<div class="table-scroll"><table class="metadata"><thead><tr><th>ID</th><th>Biological unit / Donor</th><th>Condition / State</th></tr></thead><tbody>' + sample_rows + '</tbody></table></div></details>'
        '<details open><summary>Declared statistical comparisons (' + str(len(contrasts)) + ')</summary>'
        '<div class="table-scroll"><table class="metadata"><thead><tr><th>Comparison</th><th>Numerator</th><th>Denominator</th><th class="num">Group sizes</th></tr></thead><tbody>' + contrast_rows + '</tbody></table></div></details>'
        '<p class="fig-notes">Reference assembly: <strong>' + esc(annotation.get('genome_build', 'Not recorded')) + '</strong> · ' + esc(project.get('metadata', {}).get('annotation_release', 'Annotation release recorded in provenance')) + '.</p>'
    )
    if sections:
        sections[0] = sections[0].replace('<div class="figure-grid">', design + '<div class="figure-grid">', 1)

    # Clean consulting TOC: chapter number + title only (no KPI pill badges)
    toc = ''.join('<a href="#chapter-' + str(n) + '"><b>' + str(CHAPTER_NUMBERS.get(n, n)).zfill(2) + '</b>' + esc(value[0]) + '</a>' for n, value in chapters.items())
    options = ''.join('<option value="' + esc(c['contrast_id']) + '">' + esc(c.get('label', c['contrast_id'])) + '</option>' for c in contrasts)
    pills = ''.join(
        '<button type="button" class="comparison-pill" data-contrast-value="' + esc(c['contrast_id']) + '">' + esc(c.get('label', c['contrast_id'])) + '</button>'
        for c in contrasts
    ) if len(contrasts) > 1 else ''

    footnote = project.get('metadata', {}).get('gene_filter_policy', {}).get('footnote_en', 'Actual source filtering and fitting scope are retained in provenance. Missing filtering records are not inferred.')
    appendix = (
        '<section class="chapter" id="appendix"><div class="chapter-head"><span class="chapter-number">A</span><h2>Methods, parameters and complete results</h2></div>'
        '<p>' + esc(footnote) + '</p>'
        '<p>PS-DAS (Precision-Supported DAS): native caller FDR/effect criteria plus ≥75% replicate-quantile support on effective-length normalized counts (N_eff ≥ 10, or AS-NMD high-anchor rescue) and Jeffreys Binomial SE(ΔPSI) ≤ 0.16. Ordinary DAS: the same native caller criteria without the read-support gate. Original caller P/FDR and effects are preserved without refitting.</p>'
        '<p>Executive Top 10: native FDR ascending, absolute native effect descending when available, stable object identity. Complete machine-readable Parquet/TSV outputs retain every tested gene, event, module, motif test and background gene set with SHA256 verification.</p>'
        '<div class="fig-links">' + _link('outputs/index.html', 'Complete outputs index') + _link('provenance.json', 'Parameters & SHA256 provenance') + _link('REPORT.zip', 'Download portable REPORT.zip') + '</div>'
        '<p class="fig-notes">Katana Alternative Splicing Suite · ' + _link('https://github.com/Taratype0/Katana-splicing-tool', 'Repository & documentation') + '.</p></section>'
    )

    facts_html = (
        '<div class="facts">'
        f'<div class="fact"><strong>{len(project.get("samples", []))}</strong><span>{esc(unit_label)}</span></div>'
        f'<div class="fact"><strong>{len(contrasts)}</strong><span>Declared Comparisons</span></div>'
        f'<div class="fact"><strong>{esc(annotation.get("genome_build", "Recorded"))}</strong><span>Reference Genome</span></div>'
        f'<div class="fact"><strong>{len(builder.tables)}</strong><span>Verified Complete Tables</span></div>'
        '</div>'
    )

    das_ctrl = _policy_control('report-das-selection') if not is_sc else ''
    expand_btn = '<button type="button" id="expandAllScope">Show all comparisons</button>' if len(contrasts) > 1 else ''

    eyebrow_prefix = 'SYNTHETIC DISPLAY REVIEW · ' if source.get('dummy') else 'ANALYTICAL STUDY REPORT · '
    page = (
        '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
        '<title>' + esc(project.get('name', 'Katana report')) + '</title><link rel="stylesheet" href="assets/report.css"></head><body>'
        '<header class="topbar"><a class="wordmark" href="REPORT.html">KATANA</a><small>' + ('SINGLE-CELL SMART-SEQ3 SPLICING SUITE' if is_sc else 'ALTERNATIVE SPLICING &amp; ISOFORM SUITE') + '</small>'
        '<div class="top-actions">' + _link('outputs/index.html', 'Complete Outputs') + _link('REPORT.zip', 'Download REPORT.zip') + '<button id="printReport">Print / PDF</button></div></header>'
        '<section class="cover"><p class="eyebrow">' + eyebrow_prefix + esc(source['mode'].upper()) + '</p>'
        '<h1>' + esc(project.get('name', 'Katana report')) + '</h1>'
        '<p class="lead">Quantitative multi-layer transcriptomic, co-splicing network, cis-regulatory RBP motif, and mature-mRNA protein-domain / NMD report. Every visual exhibit links directly to its exact plotted coordinates and complete un-truncated statistical tables.</p>'
        + facts_html + '</section>'
        '<div class="shell"><nav class="toc"><p class="toc-title">Report Contents</p>' + toc + '<a href="#appendix"><b>A</b><span>Methods &amp; outputs</span></a></nav>'
        '<main><div class="scopebar"><label>Comparison <select id="contrastSelect">' + options + '</select></label>'
        + ('<div class="comparison-pills">' + pills + '</div>' if pills else '')
        + das_ctrl + expand_btn + '<button type="button" id="resetScope">Reset</button><span id="scopeSummary" class="scope-summary"></span></div>'
        + ''.join(sections) + appendix + '</main></div>'
        '<dialog id="figureModal"><div class="modal-top"><span id="modalTitle"></span><div class="modal-actions"><button id="fitFigure">Fit</button><button id="actualFigure">100%</button><a id="downloadFigure" download>Download PNG</a><button id="closeFigure">Close</button></div></div><div class="modal-scroll" id="modalScroll"><img id="modalImage" alt=""></div></dialog>'
        '<script id="report-settings" type="application/json">' + json.dumps({'event_set': builder.event_set, 'initial_comparison': contrasts[0]['contrast_id'] if contrasts else 'all'}) + '</script>'
        '<script src="assets/report.js"></script></body></html>'
    )
    (output / 'REPORT.html').write_text(page, encoding='utf-8')
    _outputs(builder, chapters)

    payload = {
        'schema': 'katana.approved_report/v1',
        'specification_version': '1.2',
        'presentation_revision': 'executive-20261118',
        'project_id': project.get('project_id'),
        'project_name': project.get('name'),
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'event_set': builder.event_set,
        'samples': len(project.get('samples', [])),
        'contrasts': [c['contrast_id'] for c in contrasts],
        'panels': builder.panels,
        'top10': builder.top10,
        'policy_data': builder.policy_data,
        'complete_results': builder.tables,
        'rbp_enrichment_summaries': builder.enrichment_summaries,
        'source_files': builder.audit,
        'unavailable': builder.exclusions,
        'no_upstream_rerun': True,
        'no_source_statistics_refit': True,
        'dummy': bool(source.get('dummy', False)),
    }
    from .presentation import portable_display_project
    payload['display_project'] = portable_display_project(project)
    payload['mode'] = source['mode']
    payload['analyses'] = {key: ('selected' if value else 'not_selected') for key, value in builder.analyses.items()}
    payload['motif_summaries'] = builder.motif_summaries
    payload['network_browsers'] = builder.network_browsers
    payload['presentation'] = {
        'default_comparison': contrasts[0]['contrast_id'] if contrasts else 'all',
        'view': 'one_comparison',
        'theme': 'executive-editorial-v2',
        'network_dendrogram_in_main': False,
        'chapter_order': list(chapters),
        'chapter_numbers': CHAPTER_NUMBERS,
        'motif_sequence_rule': 'sequence-first; regional FWER and window BH separated; at most 10 distinct qualifying models; affected RNA first; all models retained',
    }
    if hasattr(builder, 'release_status'):
        payload['release_status'] = builder.release_status
    payload['figure_data_links'] = builder.figure_data_coverage()
    payload['report_selection'] = source.get('report_selection', {'applied': False})
    (output / 'report-manifest.json').write_text(json.dumps(payload, indent=2, default=str) + '\n', encoding='utf-8')

    sources = []
    for group in [source, *source.get('supplemental', [])]:
        sources.append({
            'project_id': group['project'].get('project_id'),
            'metadata': group['project'].get('metadata', {}),
            'parameters': [{'analysis_id': r.get('analysis_id'), 'artifact_id': r.get('artifact_id'), 'parameters': r.get('parameters', {}), 'method_metadata': r.get('metadata', {})} for r in group.get('artifacts', [])],
        })
    if hasattr(builder, 'preserved_source_projects'):
        sources = builder.preserved_source_projects
    provenance = {
        'schema': 'katana.report_provenance/v1',
        'source_files': builder.audit,
        'source_projects': sources,
        'report_parameters': {'event_set': builder.event_set, 'analyses': payload['analyses']},
        'derived_views': builder.derived_audits,
        'no_upstream_rerun': True,
        'no_source_statistics_refit': True,
        'figure_data_links': builder.figure_data_coverage(),
        'report_selection': source.get('report_selection', {'applied': False}),
    }
    (output / 'provenance.json').write_text(json.dumps(provenance, indent=2, default=str) + '\n', encoding='utf-8')
    if builder.create_archive:
        with zipfile.ZipFile(output / 'REPORT.zip', 'w', zipfile.ZIP_DEFLATED, compresslevel=4, strict_timestamps=False) as archive:
            for path in sorted(output.rglob('*')):
                if path.is_file() and path.name != 'REPORT.zip':
                    archive.write(path, path.relative_to(output).as_posix())
    return payload


def _outputs(builder, chapters):
    output = builder.output
    groups = []
    for number, (title, _) in chapters.items():
        disp = CHAPTER_NUMBERS.get(number, number)
        tables = [t for t in builder.tables if t.get('chapter') == number or (not t.get('chapter') and chapter_for(t.get('analysis_id', '')) == number)]
        rows = ''.join('<tr><td>' + esc(t['name']) + '</td><td class="num">' + esc(t.get('rows', '')) + '</td><td>' + _link('../' + t.get('view', t['source']), 'Full table / browser') + ' · ' + _link('../' + t['source'], 'Raw file') + '</td></tr>' for t in tables)
        figures = ''.join('<details><summary>' + esc(p['title']) + '</summary><p class="fig-notes">' + esc(p.get('caption', '')) + '</p><div class="fig-links">' + _link('../' + p['image'], 'Figure') + _data_links(p, '../') + '</div><img loading="lazy" style="max-width:100%" src="../' + esc(str(p['image']).replace('\\', '/')) + '" alt="' + esc(p['title']) + '"></details>' for p in builder.panels if p['chapter'] == number)
        groups.append('<section class="chapter" id="chapter-' + str(number) + '"><div class="chapter-head"><span class="chapter-number">' + str(disp).zfill(2) + '</span><h2>' + str(disp) + '. ' + esc(title) + '</h2></div><table class="metadata"><thead><tr><th>Complete result</th><th class="num">Rows</th><th>Access</th></tr></thead><tbody>' + rows + '</tbody></table><h3>Complete figure gallery</h3>' + figures + '</section>')
    (output / 'outputs').mkdir(exist_ok=True)
    (output / 'outputs/index.html').write_text('<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Complete outputs · Katana</title><link rel="stylesheet" href="../assets/report.css"></head><body><header class="topbar"><a class="wordmark" href="../REPORT.html">← RETURN TO MAIN REPORT</a><small>COMPLETE DATA LEDGER</small></header><main style="max-width:1400px;margin:auto;padding:32px">' + '<h1>Complete Verified Outputs</h1><p class="lead">Source tables retain every row, native P/FDR, background gene universe and SHA256 digest.</p>' + ''.join(groups) + '</main></body></html>', encoding='utf-8')
