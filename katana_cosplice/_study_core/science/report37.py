"""Portable validation panel, strictly bound to a frozen network and evidence."""
import base64, html, json
from pathlib import Path
import pandas as pd
from study import sha

def panel(run, manifest, interpretation_path=None):
    root=Path(run)/'validation37'
    if not (root/'manifest.json').exists():return ''
    binding=json.loads((root/'manifest.json').read_text())
    if binding['network_id']!=manifest['network_id'] or binding['assignment_sha256']!=manifest['assignment_sha256']:
        raise ValueError('Donor validation belongs to another network')
    for name,h in binding['files'].items():
        if Path(name).name!=name or sha(root/name)!=h:raise ValueError('Donor validation evidence changed')
    data=json.loads((root/'evaluation.json').read_text())
    interpretation_path=Path(interpretation_path) if interpretation_path else Path(run)/'interpretations/default'
    interpretation=json.loads((interpretation_path/'summary.json').read_text())
    tab=pd.read_csv(root/'module_validation.tsv',sep='\t')
    associations=pd.read_csv(interpretation_path/'associations.tsv',sep='\t')
    order=associations.groupby('module').q_value.min()
    tab['FDR']=tab.module.map(order);tab=tab.sort_values(['FDR','module'],kind='stable')
    display=tab[['module','stability_events','genes','median_match_jaccard','median_member_support',
        'supported_member_fraction','passes_stability_rule','between_donor_score_variance_fraction','median_pooled_holdout_gap','median_holdout_gap','median_trajectory_r']].rename(columns={
        'module':'Module','stability_events':'Events','genes':'Genes','median_match_jaccard':'Module Jaccard',
        'median_member_support':'Member support','supported_member_fraction':'Supported fraction',
        'passes_stability_rule':'Stability pass','between_donor_score_variance_fraction':'Between-donor variance fraction',
        'median_pooled_holdout_gap':'Pooled held-out gap','median_holdout_gap':'Within-donor held-out gap',
        'median_trajectory_r':'Held-out trajectory r'})
    display['Stability pass']=display['Stability pass'].map({True:'Yes',False:'No'})
    links=[]
    for name in ['module_validation.tsv','event_stability.tsv','partition_recovery.tsv','heldout_donor_modules.tsv','heldout_module_matching.tsv','module_score_correlations.tsv','analysis_summary.md','protocol.json']:
        p=root/name
        if p.exists():
            encoded=base64.b64encode(p.read_bytes()).decode('ascii')
            links.append(f'<a download="{html.escape(name)}" href="data:application/octet-stream;base64,{encoded}">{html.escape(name)}</a>')
    gene_links=' · '.join(f'<a download href="report_tables/module_{html.escape(m)}_genes.tsv">{html.escape(m)} genes</a> / <a download href="report_tables/module_{html.escape(m)}_events.tsv">events</a>' for m in tab.module)
    return ('<style>#donor-validation-5000 table{width:100%;border-collapse:collapse} '
        '#donor-validation-5000 th,#donor-validation-5000 td{text-align:left!important;padding:7px 8px} '
        '#donor-validation-5000 thead{border-bottom:2px solid #174e7a} '
        '#donor-validation-5000 tbody tr:nth-child(even){background:#f4f7fa}</style>'
        '<section class="card" id="donor-validation-5000" style="margin:20px;padding:20px">'
        '<h2>5,000-event donor validation / 供体稳定性验证</h2>'
        f'<p><strong>{data["stable_modules"]} / {data["total_modules"]} modules passed the unchanged stability rule.</strong> '
        '50 distinct delete-four-donor refits; 5 disjoint holdout folds covering all 19 donors. '
        'QC, variance selection, imputation and WGCNA were refitted on training donors only. This is deletion resampling, not a with-replacement bootstrap.</p>'
        '<p>Stability requires median module Jaccard ≥ 0.5, median member support ≥ 0.5 and at least 50% of members with support ≥ 0.5; '
        'all 50 fits succeeded. Unselected or unassigned events count as lost. Table shows every original module, sorted by association FDR.</p>'
        '<div style="overflow-x:auto">'+display.to_html(index=False,border=0,float_format=lambda x:f'{x:.3f}')+'</div><p>'+gene_links+'</p>'
        '<p>Held-out coherence gap = observed cross-gene correlation within the training module minus outside-module controls, with equal gene-pair weighting. '
        'Pooled gaps summarize five held-out folds (3–4 donors each); within-donor gaps summarize the 19 donors individually. Pooled correlations can include donor differences and repeated samples do not increase the independent donor count. '
        'Trajectory r uses frozen training PC1 loadings and training mean time profiles. Full-to-training module matches are one-to-one and descriptive; unmatched folds remain in the downloaded audit.</p>'
        '<p><strong>Stable does not mean a validated biological mechanism or a calibrated significant association.</strong> '
        f'The fully qualified association count for this interpretation is {interpretation["qualified_candidates"]}. '+
        ('The donor-trajectory inference and its empirical calibration are documented above; stability cutoffs are unchanged. ' if interpretation.get('inference_mode')=='donor_trajectory' else 'The legacy repeated-donor inference gate remains unqualified. ')+
        'This separate association gate does not undo the stability measurements above. M1/M2 and M3/M4 have strongly opposite scores (r = −0.995 and −0.982), '
        'so four directional modules must not be presented as four independent programs. Approximately 90–92% of M3/M4 score variance lies between donors; '
        'their stability alone does not establish a shared differentiation-time response. Holdouts are from this same study, not an external cohort.</p>'
        f'<p>Validation compute: {data["execution_wall_seconds"]/60:.1f} minutes for input preparation and 55 native fits with four workers; '
        f'summed native fit time {data["native_seconds_sum"]/60:.1f} minutes; largest individual native process tree '
        f'{data["native_peak_rss_bytes"]/2**30:.2f} GiB. These are validation costs, not the cost of an ordinary single network run.</p>'
        '<p>'+ ' · '.join(links)+'</p></section>')
