"""Check EventGraph scientific outputs and the approved sequence-card report."""
from pathlib import Path
import pandas as pd
from .model import read_json, file_sha256


def verify_eventgraph_report(receipt, transported, cid, report, report_dir):
    from .presentation import _producer_path
    from katana_next.reports.motif_sequences import evidence, _rows, _best_windows
    if receipt.get('status') not in {'completed','reused'} or receipt.get('job',{}).get('status') not in {'complete','cache_hit'}:
        raise ValueError('EventGraph has no completed scientific job')
    markers = [r for r in receipt['outputs'] if _producer_path(r['path']).name == 'COMPLETE.json']
    if len(markers) != 1:
        raise ValueError('EventGraph requires one unsplit complete family per comparison')
    marker = markers[0]; complete = read_json(transported(receipt,marker['path']))
    if complete.get('analysis_engine') != 'eventgraph_20261030' or complete.get('complete') is not True or complete.get('statistical_family_split') is not False:
        raise ValueError('Wrong or incomplete EventGraph scientific family')
    required = {'ANALYSIS_PLAN_BEFORE_SCORING.json','RESULTS.json','NEW_CLUSTER_FWER.tsv',
                'NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz','FROZEN_SELECTED_EVENTS.tsv',
                'MATCHING_ALL_SELECTED_WITH_LOSSES.tsv','ACTUAL_RNA_FORMS.tsv.gz',
                'PAIRED_FORM_INSTANCE_CHANGES.tsv.gz','MOTIF_REPORT_ALL_SEQUENCES.tsv'}
    if not required <= set(complete['members']):
        raise ValueError('EventGraph omits complete sequence, matching or inference outputs')
    base = _producer_path(marker['path']).parent
    paths = {}
    for name,record in complete['members'].items():
        if Path(name).name != name:
            raise ValueError('EventGraph member escapes the scientific family')
        path = transported(receipt,str(base/name))
        if file_sha256(path) != record['sha256'] or path.stat().st_size != record['bytes']:
            raise ValueError('EventGraph scientific member changed: '+name)
        paths[name] = path
    plan,summary = read_json(paths['ANALYSIS_PLAN_BEFORE_SCORING.json']),read_json(paths['RESULTS.json'])
    if any(r.get('comparison_id') != cid or r.get('analysis_engine') != 'eventgraph_20261030' for r in (plan,summary)):
        raise ValueError('EventGraph source belongs to a different comparison or engine')
    original = evidence(plan['models'],_rows(paths['NEW_CLUSTER_FWER.tsv']),
                        _best_windows(paths['NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz']))
    cards = [p for p in report.get('panels',[]) if p.get('contrast_id') == cid and p.get('metadata',{}).get('motif_view_kind') == 'sequence_evidence']
    if [p['metadata']['sequence_evidence']['model_id'] for p in cards] != [c['model_id'] for c in original['cards']]:
        raise ValueError('Report Top 10 differs from complete native EventGraph significance and ranking')
    for panel,card in zip(cards,original['cards']):
        shown = panel['metadata']['sequence_evidence']
        if shown.get('regions') != card['regions'] or shown.get('window') != card['window']:
            raise ValueError('Regional FWER or window BH statistics changed in the report')
    summaries = [s for s in report.get('motif_summaries',[]) if s.get('contrast_id') == cid]
    if len(summaries) != 1 or summaries[0].get('sequence_report',{}).get('displayed_models') != original['displayed_models']:
        raise ValueError('Report is missing its comparison-specific motif result or no-significance summary')
    report_dir = Path(report_dir).resolve()
    # Every source file must be downloadable, including non-significant models.
    available = set()
    for row in report.get('complete_results',[]):
        if row.get('contrast_id') != cid:
            continue
        path = (report_dir/row['source']).resolve()
        if not path.is_relative_to(report_dir) or not path.is_file():
            raise ValueError('Motif report download is absent or outside its report')
        available.add(file_sha256(path))
    # HTML may be rebased for local links; numerical and method files may not.
    for name,record in complete['members'].items():
        if Path(name).suffix not in {'.html','.png','.svg'} and record['sha256'] not in available:
            raise ValueError('Report omits complete EventGraph output: '+name)
    return {'comparison_id':cid,'models':len(plan['models']),'displayed_models':original['displayed_models']}
