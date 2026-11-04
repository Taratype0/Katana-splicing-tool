"""Sequence-first views of recorded MotifTrace statistics; no scientific fits.

This module is shared by portable reports and native event-consequence reports.
Regional FWER and window BH values are read from their own original columns.
"""
from __future__ import annotations

import csv
import hashlib
import html
import io
import json
import math
import os
from pathlib import Path

VERSION = 'motif-sequence-report-v2-rbp-genes'
BASES = 'ACGU'
REGIONAL_P = 'cluster_P_FWER_new_family'
WINDOW_Q = 'q_new_window_family'
WINDOW_P = 'raw_P_new_paired_swaps'
IUPAC = {'A':'A','C':'C','G':'G','U':'U','AG':'R','CU':'Y','CG':'S','AU':'W',
         'GU':'K','AC':'M','CGU':'B','AGU':'D','ACU':'H','ACG':'V','ACGU':'N'}


def digest(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def atomic(path, content):
    """Never write through an immutable payload hard link on report refresh."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.sequence-report.tmp')
    temporary.write_bytes(content.encode('utf-8') if isinstance(content, str) else content)
    os.replace(temporary, path)


def number(value):
    if value is None or str(value).strip().lower() in {'', 'na', 'nan', 'none'}:
        return None
    value = float(value)
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError('Recorded motif significance must be in [0, 1] or missing')
    return value


def representation(model):
    matrix = model.get('matrix')
    literal = str(model.get('sequence', model.get('motif', model.get('pattern', '')))).upper().replace('T', 'U')
    if matrix is None:
        if literal and set(literal) <= set(BASES):
            matrix = [[float(base == char) for base in BASES] for char in literal]
        else:
            return {'sequence': literal or 'not available', 'iupac': literal or 'not available',
                    'matrix': None, 'kind': 'recorded pattern; base probabilities unavailable'}
    normalized = []
    for row in matrix:
        if len(row) != 4 or any(not math.isfinite(float(p)) or float(p) < 0 for p in row) or sum(row) <= 0:
            raise ValueError('Invalid recorded position-by-four motif matrix')
        normalized.append([float(p) / sum(row) for p in row])
    if not normalized:
        raise ValueError('An empty PWM cannot represent a motif')
    sequence, degenerate = [], []
    for row in normalized:
        order = sorted(range(4), key=lambda i: (-row[i], i))
        sequence.append(BASES[order[0]])
        total, selected = 0., []
        for i in order:
            total += row[i]; selected.append(BASES[i])
            if total >= .9: break
        degenerate.append(IUPAC[''.join(sorted(selected))])
    return {'sequence': ''.join(sequence), 'iupac': ''.join(degenerate), 'matrix': normalized,
            'kind': 'maximum-probability PWM consensus' if model.get('matrix') is not None else 'literal sequence'}


def _rows(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def _tsv(path, rows, columns):
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t', extrasaction='ignore')
    writer.writeheader(); writer.writerows(rows)
    atomic(path, stream.getvalue())


def _best_windows(path):
    """Select real rows, retaining paired raw P and q; never compute a model q."""
    import pandas as pd
    columns = ['model_id', WINDOW_Q, WINDOW_P, 'event_type', 'direction', 'region', 'method', 'window_nt', 'window_start0']
    best = {}
    for chunk in pd.read_csv(path, sep='\t', usecols=lambda c: c in columns,
                             dtype={'model_id': str}, chunksize=100000, float_precision='round_trip'):
        if not {WINDOW_Q, WINDOW_P, 'model_id'} <= set(chunk):
            if chunk.empty: continue
            raise ValueError('Native window P/q columns are required')
        valid = chunk[chunk[WINDOW_Q].notna() & chunk[WINDOW_P].notna()]
        if valid.empty: continue
        for col in (WINDOW_Q, WINDOW_P):
            if not valid[col].between(0, 1).all(): raise ValueError('Invalid native window P/q')
        sort = [c for c in columns if c != 'model_id']
        for row in valid.sort_values(sort, kind='stable').drop_duplicates('model_id').to_dict('records'):
            key = (float(row[WINDOW_Q]), float(row[WINDOW_P]), str(row))
            if row['model_id'] not in best or key < best[row['model_id']][0]: best[row['model_id']] = key, row
    return [v[1] for v in best.values()]


def evidence(models, regions, windows, *, regional_available=True, limit=10):
    """At most ten distinct cards; keep all original rows in complete files."""
    catalogue = {str(k): representation(v) for k,v in models.items()}
    sig_regions, sig_windows, estimable_windows = [], [], set()
    for row in regions:
        p = number(row.get(REGIONAL_P))
        if p is not None and p <= .05:
            if str(row.get('model_id')) not in catalogue: raise ValueError('Undeclared regional motif')
            sig_regions.append(dict(row))
    for row in windows:
        q, p = number(row.get(WINDOW_Q)), number(row.get(WINDOW_P))
        if q is not None and p is not None: estimable_windows.add(str(row.get('model_id')))
        if q is not None and p is not None and q <= .05:
            if str(row.get('model_id')) not in catalogue: raise ValueError('Undeclared window motif')
            sig_windows.append(dict(row))
    layer_rank = {'affected_sequence': 0, 'regulatory_pre_mRNA': 1}
    sig_regions.sort(key=lambda r: (layer_rank.get(r.get('layer'), 2), number(r[REGIONAL_P]), str(r['model_id']), str(r.get('region')), str(r.get('covered_region_start0'))))
    sig_windows.sort(key=lambda r: (number(r[WINDOW_Q]), number(r[WINDOW_P]), str(r['model_id'])))
    order = list(dict.fromkeys(str(r['model_id']) for r in sig_regions + sig_windows))
    cards = []
    for token in order[:limit]:
        regional = [r for r in sig_regions if str(r['model_id']) == token]
        window = next((r for r in sig_windows if str(r['model_id']) == token), None)
        cards.append({'model_id': token, **catalogue[token], 'regions': regional, 'window': window,
            'library_annotations': str(models[token].get('RBPs', 'not recorded')),
            'assignments': models[token].get('assignments', []),
            'priority': min([layer_rank.get(r.get('layer'), 2) for r in regional] or [3])})
    counts = {layer: sum(r.get('layer') == layer for r in sig_regions) for layer in layer_rank}
    counts['unclassified'] = len(sig_regions) - sum(counts.values())
    return {'version': VERSION, 'regional_available': regional_available, 'regional_counts': counts,
            'regional_significant_regions': len(sig_regions),
            'window_significant_models': len({r['model_id'] for r in sig_windows}),
            'window_estimable_models': len(estimable_windows),
            'eligible_models': len(order), 'displayed_models': len(cards), 'limit': limit,
            'cards': cards, 'catalogue': catalogue,
            'selection_rule': 'At most 10 distinct models: affected-RNA regional FWER, flank regional FWER, then window BH; each endpoint <=0.05; complete models and tests retained',
            'new_scans': 0, 'new_statistical_tests': 0, 'native_statistics_unchanged': True}


def position(row, regional=True):
    if regional:
        start, end = row.get('covered_region_start0'), row.get('covered_region_end0')
        span = f'[{start}, {end}) nt' if start not in (None, '') and end not in (None, '') else 'span not recorded'
    else:
        width = row.get('window_nt')
        span = 'full-sequence test' if str(width) in ('0', '0.0') else f"window start {row.get('window_start0')}; width {width} nt"
    coordinate = row.get('coordinate_system') or 'native region-relative coordinates'
    return ' / '.join(str(row.get(k, 'not recorded')) for k in ('event_type','direction','region')) + ' | ' + span + ' | ' + coordinate


def _logo(path, rep, title):
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.textpath import TextPath
    from matplotlib.transforms import Affine2D
    from matplotlib.patches import PathPatch
    fig = Figure(figsize=(9.5, 2.7), layout='constrained'); FigureCanvasAgg(fig)
    ax = fig.subplots(); matrix = rep['matrix']
    if matrix is None:
        ax.text(.5, .5, rep['sequence']+'\nRecorded pattern; no probability logo available', ha='center', va='center', transform=ax.transAxes)
        ax.axis('off')
    else:
        colors = ['#238B45','#2863A5','#D49A25','#C63F4A']
        for i, row in enumerate(matrix):
            bottom = 0.
            for j in sorted(range(4), key=lambda j: row[j]):
                height = row[j]
                if height <= 1e-9: continue
                glyph = TextPath((0,0), BASES[j], size=1, prop={'weight':'bold'})
                box = glyph.get_extents()
                transform = Affine2D().translate(-box.x0,-box.y0).scale(.82/box.width,height/box.height).translate(i+.09,bottom)
                ax.add_patch(PathPatch(glyph, transform=transform+ax.transData, color=colors[j], linewidth=0))
                bottom += height
        ax.set(xlim=(0,len(matrix)),ylim=(0,1.02),ylabel='Base probability',xlabel='Motif position (RNA, 5′ to 3′)')
        ax.set_xticks([i+.5 for i in range(len(matrix))], range(1,len(matrix)+1))
        ax.spines[['top','right']].set_visible(False)
    ax.set_title(title, color='#2863A5', fontweight='bold')
    for suffix in ('.png','.svg'):
        stream = io.BytesIO(); fig.savefig(stream, format=suffix[1:], dpi=150)
        atomic(path.with_suffix(suffix), stream.getvalue())
    fig.clear()


def write_assets(models, summary, directory, comparison_id):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    catalogue = []
    for token, rep in summary['catalogue'].items():
        catalogue.append({'comparison_id': comparison_id, 'model_id': token, 'motif_RNA': rep['sequence'],
            'IUPAC_RNA_90pct': rep['iupac'], 'representation':rep['kind'],
            'library_annotations': models[token].get('RBPs',''),
            'source_assignments': json.dumps(models[token].get('assignments',[]), ensure_ascii=False)})
    all_path = directory/'MOTIF_REPORT_ALL_SEQUENCES.tsv'
    _tsv(all_path, catalogue, ['comparison_id','model_id','motif_RNA','IUPAC_RNA_90pct','representation','library_annotations','source_assignments'])
    original_path = directory/'MOTIF_REPORT_ORIGINAL_MODELS.json'
    atomic(original_path, json.dumps(models, ensure_ascii=False, indent=2, allow_nan=False))
    files = [all_path, original_path]
    for card in summary['cards']:
        token = card['model_id']; name = 'MOTIF_REPORT_' + hashlib.sha256(token.encode()).hexdigest()[:16]
        base = directory/name
        matrix_path = base.with_suffix('.tsv')
        rows = [{'comparison_id':comparison_id,'model_id':token,'position_1based':i+1,**dict(zip(BASES,row))}
                for i,row in enumerate(card['matrix'] or [])]
        # Pattern-only models have an explicitly empty numeric matrix; never invent probabilities.
        _tsv(matrix_path, rows, ['comparison_id','model_id','position_1based',*BASES])
        region_rows = [{'comparison_id':comparison_id, 'endpoint':'regional FWER P', **r} for r in card['regions']]
        if card['window']: region_rows.append({'comparison_id':comparison_id, 'endpoint':'window/full-sequence BH FDR', **card['window']})
        evidence_path = base.with_name(name+'-evidence.tsv')
        columns = list(dict.fromkeys(k for row in region_rows for k in row))
        _tsv(evidence_path, region_rows, columns)
        _logo(base, card, '5′-'+card['sequence']+'-3′')
        card.update(image=str(base.with_suffix('.png')), svg=str(base.with_suffix('.svg')),
                    matrix_path=str(matrix_path), evidence_path=str(evidence_path))
        files.extend([matrix_path,evidence_path,base.with_suffix('.png'),base.with_suffix('.svg')])
    return files


def summary_html(summary):
    esc = lambda v: html.escape(str(v), quote=True)
    counts = summary['regional_counts']
    regional = (f"{counts['affected_sequence']} significant affected-RNA regions; {counts['regulatory_pre_mRNA']} significant regulatory-flank regions (regional FWER P ≤ 0.05)."
                if summary['regional_available'] else 'Regional FWER: not available from an explicit native correction column.')
    if counts.get('unclassified'): regional += f" {counts['unclassified']} significant regions have no recorded biological layer."
    title = ('Motif consequences: significant affected-RNA regions' if counts['affected_sequence'] else
             'Motif consequences: completed, no significant affected-RNA regions' if summary['regional_available'] else 'Motif consequences: regional results unavailable')
    if summary.get('window_test_count') == 0 and not summary['regional_available']:
        title = 'Motif consequences: not estimable for this comparison'
        regional = 'No eligible native window/full-sequence tests were recorded. Regional significance is unavailable; this is not a nonsignificant result.'
    window = f"{summary['window_significant_models']} motif models have a recorded window/full-sequence test with BH FDR ≤ 0.05."
    if not summary.get('window_estimable_models') and not summary.get('window_test_count'):
        window = 'No estimable native window/full-sequence P/q values are available; this is not a measured zero-discovery result.'
    return ('<h3>'+esc(title)+'</h3><p>'+esc(regional)+'</p><p><strong>Window / full-sequence BH FDR:</strong> '+esc(window)+'</p>'
        '<p class="fig-notes">Regional FWER and window BH FDR are separate endpoints. Regional significance does not establish significant individual windows or binding sites. '
        'Sequence and logo describe the supplied motif model; RBP names are downstream library annotations, not binding or causality evidence.</p>'
        '<p class="fig-notes">'+esc(f"Showing {summary['displayed_models']} of {summary['eligible_models']} qualifying distinct models. "+summary['selection_rule'])+'</p>')


def card_html(card):
    from .motif_gene_candidates import candidate_html
    esc = lambda v: html.escape(str(v), quote=True)
    rows = []
    for row in card['regions'][:20]:
        rows.append('<tr><td>'+esc(row.get('layer','not recorded'))+'</td><td>'+esc(position(row))+'</td><td>Regional FWER P = '+esc(f'{number(row[REGIONAL_P]):.6g}')+'</td></tr>')
    if card['window']:
        row=card['window']
        rows.append('<tr><td>Window / full-sequence test</td><td>'+esc(position(row,False))+'</td><td>BH FDR = '+esc(f'{number(row[WINDOW_Q]):.6g}')+'; paired raw P = '+esc(f'{number(row[WINDOW_P]):.6g}')+'</td></tr>')
    return ('<div class="table-scroll"><table class="metadata"><tr><th>Endpoint / biological layer</th><th>Position and event direction</th><th>Native significance</th></tr>'+''.join(rows)+'</table></div>'
        '<p class="fig-notes">Coordinates are zero-based and half-open when endpoints are supplied; aggregate region-relative coordinates are not individual genomic loci. '
        'At most 20 regional rows are shown here; every row is retained in the evidence download. '
        'No window q is inferred from regional P.</p>'+candidate_html(card)+
        '<details><summary>Downstream RBP annotation and model provenance</summary><p>'+esc(card['library_annotations'])+'</p><p>'+esc(json.dumps(card['assignments'],ensure_ascii=False))+'</p>'
        '<p>Model: '+esc(card['model_id'])+'; representation: '+esc(card['kind'])+'; 90% cumulative-probability IUPAC consensus: '+esc(card['iupac'])+'.</p></details>')


def attach_to_builder(builder):
    """One automatic hook for fresh reports, API builds and portable refreshes."""
    root = builder.output.resolve()
    tables = [t for t in builder.tables if not t.get('metadata',{}).get('motif_sequence_renderer')]
    def member(cid, basename):
        choices=[]
        for record in tables:
            if str(record.get('contrast_id')) != cid: continue
            path=(root/record.get('source','')).resolve()
            if path.name == basename or path.name.endswith('-'+basename):
                if not path.is_relative_to(root) or not path.is_file(): raise ValueError('Motif evidence escapes portable report')
                if path not in choices: choices.append(path)
        if len(choices)>1:
            if len({digest(p) for p in choices}) != 1: raise ValueError('Conflicting motif evidence for one comparison: '+basename)
        return choices[0] if choices else None
    all_cids = sorted({str(t.get('contrast_id')) for t in tables if t.get('contrast_id')})
    for cid in all_cids:
        plan_path=member(cid,'ANALYSIS_PLAN_BEFORE_SCORING.json')
        result_path=member(cid,'RESULTS.json')
        if not plan_path and not result_path: continue
        plan=json.loads(plan_path.read_text('utf-8-sig')) if plan_path else {}
        result=json.loads(result_path.read_text('utf-8-sig')) if result_path else {}
        # Existing installed pipeline outputs retain their validated legacy
        # panels and native-test acceptance. Opt into the new sequence renderer
        # only through the recorded EventGraph engine, never a filename guess.
        engine=plan.get('analysis_engine') or result.get('policy',{}).get('analysis_engine')
        if engine != 'eventgraph_20261030':
            continue
        models=plan.get('models') or result.get('policy',{}).get('models')
        if not models: continue
        for data in (plan,result):
            if data.get('comparison_id',cid)!=cid: raise ValueError('Motif models belong to another comparison')
        region_path=member(cid,'NEW_CLUSTER_FWER.tsv')
        regions=_rows(region_path) if region_path else []
        available=bool(region_path and (not regions or REGIONAL_P in regions[0]))
        if regions and not available: regions=[]
        ranking_path=member(cid,'all-motif-ranked-recorded-tests.tsv')
        window_path=member(cid,'NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz')
        windows=_rows(ranking_path) if ranking_path else (_best_windows(window_path) if window_path else [])
        if ranking_path and any(r.get('comparison_id',cid)!=cid for r in windows): raise ValueError('Window ranking has different comparison')
        summary=evidence(models,regions,windows,regional_available=available)
        summary['window_test_count']=result.get('eligible_new_window_tests')
        if result.get('eligible_new_window_tests') == 0:
            summary['regional_available']=False
        events_path=member(cid,'FROZEN_SELECTED_EVENTS.tsv')
        matching_path=member(cid,'MATCHING_ALL_SELECTED_WITH_LOSSES.tsv')
        changes_path=member(cid,'PAIRED_FORM_INSTANCE_CHANGES.tsv.gz')
        evidence_files=[p for p in (plan_path,result_path,region_path,ranking_path,window_path,events_path,matching_path,changes_path) if p]
        before={str(p):digest(p) for p in evidence_files}
        directory=root/'outputs/motif-sequences'/hashlib.sha256(cid.encode()).hexdigest()[:16]
        files=write_assets(models,summary,directory,cid)
        from katana_motiftrace.annotations import enrich, project_species
        target_species,species_basis=project_species(builder.source.get('project',{}))
        files.extend(enrich(summary,models,directory,cid,target_species=target_species,species_basis=species_basis,
            events_path=events_path,matching_path=matching_path,changes_path=changes_path))
        old_ids={p['id'] for p in builder.panels if p.get('metadata',{}).get('motif_sequence_renderer') and p.get('contrast_id')==cid}
        builder.panels=[p for p in builder.panels if p.get('id') not in old_ids]
        builder.tables=[t for t in builder.tables if not(t.get('metadata',{}).get('motif_sequence_renderer') and str(t.get('contrast_id'))==cid)]
        builder.data_link_audit=[r for r in builder.data_link_audit if r.get('panel_id') not in old_ids]
        # Prior sequence cards remain downloadable; the new card is sequence-first.
        for panel in builder.panels:
            if panel.get('contrast_id') != cid:
                continue
            metadata = panel.get('metadata', {})
            native_name = Path(metadata.get('relative_path', '')).name
            native_logo = (metadata.get('analysis_engine') == 'eventgraph_20261030'
                           and native_name.startswith('MOTIF_REPORT_')
                           and Path(native_name).suffix in {'.png', '.svg'})
            if metadata.get('motif_view_kind') == 'top_motif' or native_logo:
                # The native report supplies standalone logo files too. The
                # comparison-specific cards below already show these models,
                # with exact PWM links and original P/FDR. Preserve the native
                # files as attachments instead of displaying duplicate cards.
                panel['output_only'] = True
        views={'panels':[], 'tables':[]}
        recorded_scope=next((p.get('metadata',{}).get('das_scope') for p in builder.panels
            if p.get('contrast_id')==cid and p.get('metadata',{}).get('das_scope')), {})
        recorded_policy=recorded_scope.get('recorded_event_set') if recorded_scope.get('mode') in {'fixed_completed_result','rs_required'} else None
        for card in summary['cards']:
            links=[{'role':role,'path':str(path),'label':label,'unit':unit,'contrast_id':cid,'policy':recorded_policy}
                for role,path,label,unit in [('plot_data',card['matrix_path'],'Exact motif PWM','motif_position_probability'),
                    ('full_results',card['evidence_path'],'All native evidence for this motif','motif_region_or_window'),
                    ('methods',directory/'MOTIF_REPORT_ORIGINAL_MODELS.json','Original models and annotations','method_record')]]
            for key,label,unit in [('candidate_genes_path','Complete motif-to-RBP gene mapping','candidate_RBP_gene'),
                ('target_events_path','All changed target genes and events for this motif','descriptive_target_event'),
                ('gene_mapping_reference_path','Frozen gene identity and orthology evidence','annotation_reference')]:
                links.append({'role':'full_results' if key!='gene_mapping_reference_path' else 'methods','path':card[key],
                    'label':label,'unit':unit,'contrast_id':cid,'policy':recorded_policy})
            for basename,label,unit in [('ACTUAL_LOCUS_HITS.tsv.gz','Existing RNA-form candidate sequences and genomic blocks','descriptive_motif_site'),
                    ('PAIRED_FORM_INSTANCE_CHANGES.tsv.gz','Existing gained / lost / retained sequence instances','motif_instance')]:
                path=member(cid,basename)
                if path: links.append({'role':'full_results','path':str(path),'label':label,'unit':unit,'contrast_id':cid,'policy':recorded_policy})
            public={k:v for k,v in card.items() if k not in ('matrix','image','svg','matrix_path','evidence_path','candidate_genes_path','target_events_path','gene_mapping_reference_path')}
            views['panels'].append({'id':'motif-seq-'+hashlib.sha256((cid+'\0'+card['model_id']).encode()).hexdigest()[:16],
                'chapter':7,'title':'Motif sequence: 5′-'+card['sequence']+'-3′ | Comparison: '+cid,
                'path':card['image'],'contrast_id':cid,'policy':None,'output_only':False,
                'caption':'Original supplied motif model, not de novo motif discovery. Logo shows base probabilities; T is displayed as U for RNA. Exact native regional FWER and window BH statistics are listed separately below.',
                'metadata':{'component':'motiftrace_sequence','motif_view_kind':'sequence_evidence',
                    'das_scope':recorded_scope,
                    'motif_sequence_renderer':VERSION,'sequence_evidence':public,'sequence_priority':card['priority'],
                    'completion_verified':True,'completion_scope':'motiftrace_sequence','analysis_id':'prioritize.motiftrace_consequences'},
                'data_links':links})
        for path in files:
            views['tables'].append({'name':path.name+' | '+cid,'path':str(path),'chapter':7,'contrast_id':cid,'kind':'file',
                                    'metadata':{'motif_sequence_renderer':VERSION}})
        builder.add_views(views)
        record=next((r for r in builder.motif_summaries if str(r.get('contrast_id'))==cid),None)
        if record is None:
            record={'contrast_id':cid,'chapter':7}; builder.motif_summaries.append(record)
        public_summary={k:v for k,v in summary.items() if k not in ('cards','catalogue')}
        record.update(sequence_report=public_summary,regional_FWER_significant_regions=summary['regional_significant_regions'],regional_FWER_is_window_FDR=False)
        builder.derived_audits.append({'schema':VERSION,'contrast_id':cid,'source_sha256':before,
            'native_statistics_unchanged':all(digest(p)==h for p,h in before.items()),'new_scans':0,'new_statistical_tests':0})
        if not builder.derived_audits[-1]['native_statistics_unchanged']: raise RuntimeError('Native motif evidence changed during presentation')


def standalone_fragment(summary, clusters, destination):
    """Return native report HTML and files covered by its normal completion manifest."""
    root=Path(destination).parent
    models=summary.get('policy',{}).get('models',{})
    if not models: return ''
    window_path=root/'NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz'
    windows=_best_windows(window_path) if summary.get('new_window_q_lt_0_05',0) and window_path.is_file() else []
    region_rows=clusters.to_dict('records')
    report=evidence(models,region_rows,windows,regional_available=REGIONAL_P in clusters.columns)
    report['window_test_count']=summary.get('eligible_new_window_tests')
    if summary.get('eligible_new_window_tests') == 0: report['regional_available']=False
    write_assets(models,report,root,str(summary['comparison_id']))
    from katana_motiftrace.annotations import enrich, project_species
    provenance=summary.get('policy',{}).get('provenance',{})
    target_species,species_basis=project_species({'annotation':{'genome_build':provenance.get('genome_build')},'metadata':{'species':provenance.get('species')}})
    def existing(name):
        path=root/name
        return path if path.is_file() else None
    enrich(report,models,root,str(summary['comparison_id']),target_species=target_species,species_basis=species_basis,
        events_path=existing('FROZEN_SELECTED_EVENTS.tsv'),matching_path=existing('MATCHING_ALL_SELECTED_WITH_LOSSES.tsv'),
        changes_path=existing('PAIRED_FORM_INSTANCE_CHANGES.tsv.gz'))
    output=['<section id="motif-sequences">',summary_html(report)]
    for card in report['cards']:
        output.append('<article><h3>Motif sequence: 5′-'+html.escape(card['sequence'])+'-3′</h3><img style="max-width:100%" data-file="'+html.escape(Path(card['image']).name)+'" src="'+html.escape(Path(card['image']).name)+'" alt="Original motif PWM logo">'+card_html(card))
        for key,label in [('matrix_path','Exact PWM'),('evidence_path','All regional/window evidence'),
            ('candidate_genes_path','Complete motif-to-RBP gene mapping'),('target_events_path','All changed target genes and events'),
            ('gene_mapping_reference_path','Gene identity and orthology evidence')]:
            output.append('<p><a data-file="'+html.escape(Path(card[key]).name)+'" href="'+html.escape(Path(card[key]).name)+'">'+label+'</a></p>')
        output.append('</article>')
    output.append('<p><a data-file="MOTIF_REPORT_ALL_SEQUENCES.tsv" href="MOTIF_REPORT_ALL_SEQUENCES.tsv">All motif sequences</a> · <a data-file="MOTIF_REPORT_ORIGINAL_MODELS.json" href="MOTIF_REPORT_ORIGINAL_MODELS.json">Original models and annotations</a></p></section>')
    return ''.join(output)
