"""Portable, study-first report over a frozen network. No fitting or selection."""
import html, json, shutil
from pathlib import Path
import pandas as pd
from study import verify_frozen


def read_json(path, default=None):
    return json.loads(path.read_text(encoding='utf-8')) if path.exists() else (default or {})


def table(path):
    try:
        result=pd.read_csv(path, sep='\t', dtype={'sample_id':str, 'biological_sample_id':str, 'condition':str,'gene_id':str,'event_key':str,'source':str,'target':str})
        if 'status' in result:result['status']=result.status.fillna('NA')
        return result
    except (FileNotFoundError, pd.errors.EmptyDataError): return pd.DataFrame()


def records(tab):
    # pandas.to_json defaults to fixed decimal precision and can turn small
    # nonzero P/FDR values into 0. Preserve Python floats through JSON encoding.
    return tab.astype(object).where(pd.notna(tab),None).to_dict(orient='records')


def build(run, destination, interpretation=None, summary=None):
    run=Path(run); destination=Path(destination); destination.mkdir(parents=True, exist_ok=True)
    interpretation=Path(interpretation) if interpretation else run/'interpretations/default'
    manifest=verify_frozen(run)
    result=summary if summary is not None else read_json(run/'summary.json')
    analysis=read_json(interpretation/'summary.json')
    samples=table(interpretation/'analysis_samples.tsv')
    if samples.empty: samples=table(run/'samples.tsv')
    scores=table(run/'study/scores.tsv')
    if scores.empty: scores=pd.DataFrame({'sample_id':samples.sample_id})
    modules=table(run/'study/modules.tsv')
    factor_members=table(run/'model/membership_evidence.tsv')
    if not factor_members.empty:
        definition=factor_members[factor_members.module.gt(0)].groupby('module').strongest_factor.first()
        modules['latent_factor']=modules.module.map(lambda m:'F'+str(int(definition.get(int(m[1:]),0))))
        modules['loading_direction']=modules.module.map(lambda m:'positive' if int(m[1:])%2 else 'negative')
    support=table(run/'modules.tsv')
    if 'module' in support:
        support=support.copy(); support['module']=support.module.map(lambda x:'M'+str(int(float(x))))
        modules=modules.merge(support.drop(columns=['events'],errors='ignore'),on='module',how='left')
    diagnostics=table(interpretation/'co_response_diagnostic.tsv')
    if 'module' in diagnostics: modules=modules.merge(diagnostics,on='module',how='left')
    from module_qualification import qualification
    associations,confidence=qualification(run,table(interpretation/'associations.tsv'),diagnostics)
    events=table(run/'events.tsv')
    if not factor_members.empty:
        events=events.merge(factor_members[['event_key','strongest_factor','absolute_loading','loading_margin']],on='event_key',how='left')
    if 'module' in events:
        events=events[pd.to_numeric(events.module,errors='coerce').gt(0)].copy()
        events['module']=events.module.map(lambda x:'M'+str(int(float(x))))
    paths=table(run/'pathway_enrichment.tsv')
    validation=read_json(run/'scientific_validation.json')
    if validation and (validation.get('network_id')!=manifest['network_id'] or validation.get('assignment_sha256')!=manifest['assignment_sha256']):
        raise ValueError('Scientific validation is bound to a different network')
    resolution=read_json(run/'resolution_validation.json')
    if resolution and (resolution.get('network_id')!=manifest['network_id'] or resolution.get('assignment_sha256')!=manifest['assignment_sha256']):
        raise ValueError('Resolution validation is bound to a different network')
    if 'module' in paths:
        paths['module']=paths.module.map(lambda x:'M'+str(int(float(x))))
        paths=paths.sort_values('q_value',kind='stable')
    pathway_total=len(paths)
    # A report should not embed hundreds of thousands of repeated term records.
    # This is a disclosed display subset; the complete correction family and
    # all nonsignificant results remain in the downloadable table / application.
    if 'q_value' in paths:
        paths=paths.loc[paths.q_value.le(.05)]
    from study_report_data import member_evidence
    usage=table(run/'usage.tsv.gz')
    if not usage.empty:usage=usage.set_index(usage.columns[0])
    events,genes,event_usage,graphs=member_evidence(events,usage,scores,samples,table(run/'edges.tsv'))
    sources={
        'sample_manifest.tsv':run/'samples.tsv', 'module_scores.tsv':run/'study/scores.tsv',
        'module_coverage.tsv':run/'study/coverage.tsv', 'module_basis.tsv':run/'study/module_basis.tsv',
        'module_members.tsv':run/'events.tsv', 'event_qc.tsv':run/'event_audit.tsv',
        'associations.tsv':interpretation/'associations.tsv', 'condition_scores.tsv':interpretation/'condition_scores.tsv',
        'unit_influence.tsv':interpretation/'co_response_diagnostic.tsv', 'pathways.tsv':run/'pathway_enrichment.tsv',
        'study_manifest.json':run/'study/manifest.json', 'input_assembly.json':run/'inputs/assembly.json',
        'analysis_summary.json':interpretation/'summary.json', 'analysis_config.json':interpretation/'config.json',
        'network_edges.tsv':run/'edges.tsv',
        'module_stability.tsv':run/'stability/module_evidence.tsv',
        'stability_provenance.json':run/'stability/qualification_manifest.json',
        'statistical_audit.json':run/'study/statistical_audit.json',
        'statistical_audit.tsv':run/'study/statistical_audit.tsv',
        'factor_loadings.tsv':run/'model/loadings.tsv',
        'factor_membership_evidence.tsv':run/'model/membership_evidence.tsv',
        'scientific_validation.json':run/'scientific_validation.json',
        'matched_method_comparison.tsv':run/'validation/comparison.tsv',
        'matched_event_sensitivity.tsv':run/'validation/matched_coverage.tsv',
        'alternative_caller_corroboration.tsv':run/'validation/jum_corroboration.tsv',
        'sample_qc_correlations.tsv':run/'validation/sample_qc_associations.tsv',
        'gene_balanced_recovery.tsv':run/'validation/gene_balanced_recovery.tsv',
        'refit_gene_recovery.tsv':run/'validation/refit_gene_recovery.tsv'}
    sources['resolution_validation.json']=run/'resolution_validation.json'
    for name in ['enrichment_manifest.json','enrichment_sources.json','enrichment_backgrounds.tsv','enrichment_network_genes.tsv','enrichment_gene_mapping.tsv']:
        sources[name]=run/name
    if analysis.get('inference_mode')=='donor_trajectory':
        for name in ['trajectory_components.tsv','donor_coefficients.tsv','legacy_associations.tsv','legacy_summary.json']:
            sources[name]=interpretation/name
    for name in ['programs.tsv','program_scores.tsv','program_loadings.tsv','module_program_links.tsv','local_modules.tsv']:
        sources[name]=run/'model'/name
    for name in ['resolution_profile.tsv','resolution_panel.tsv','resolution_curves.tsv','resolution_candidates.tsv',
                 'native_granularity.tsv','megena_coverage_diagnosis.tsv','resolution_audit.png','resolution_evidence_manifest.json']:
        sources[name]=run/'validation29'/name
    downloads=[]; exports=destination/'report_tables'; exports.mkdir(exist_ok=True)
    for name, source in sources.items():
        if source.exists():
            shutil.copyfile(source,exports/name); downloads.append(name)
    for name,frame in [('gene_event_membership.tsv',events),('module_gene_summary.tsv',genes)]:
        frame.to_csv(exports/name,sep='\t',index=False);downloads.append(name)
    associations.to_csv(exports/'module_qualification.tsv',sep='\t',index=False);downloads.append('module_qualification.tsv')
    associations.loc[associations.qualified_candidate].to_csv(exports/'qualified_modules.tsv',sep='\t',index=False) if 'qualified_candidate' in associations else associations.to_csv(exports/'qualified_modules.tsv',sep='\t',index=False)
    downloads.append('qualified_modules.tsv')
    is31=read_json(run/'report_revision.json').get('schema')=='study-report-3.1'
    is29=read_json(run/'report_revision.json').get('schema')=='study-report-2.9' or is31
    module_files={}
    if is29:
        from report29 import module_exports
        module_files,extra=module_exports(exports,genes,events,associations);downloads.extend(extra)
    if len(events):
        usage.reindex(events.event_key).to_csv(exports/'module_event_psi.tsv.gz',sep='\t',index_label='event_key');downloads.append('module_event_psi.tsv.gz')
    payload={'schema':'study-report-2.7', 'manifest':manifest, 'result':result, 'analysis':analysis,
             'scientific_validation':validation,'resolution_validation':resolution,'module_files':module_files,
             'assembly':read_json(run/'inputs/assembly.json'), 'samples':records(samples), 'scores':records(scores),
             'modules':records(modules), 'associations':records(associations),'qualification':confidence,
             'statistical_audit':read_json(run/'study/statistical_audit.json'),
             'members':records(events),'genes':records(genes),'event_usage':event_usage,'graphs':graphs,
             'pathways':records(paths), 'pathway_rows_total':pathway_total,
             'module_display_rule':'Default: FDR <= 0.05 AND matched full-QC stability criteria AND no single-unit dominance. No fallback to unstable modules. Sort FDR, then P, then stable ID. Explicit exploratory view retains FDR-only associations. Audit ledger includes every test and exclusion reason; correction precedes stability filtering.',
             'pathway_display_rule':'Only FDR <= 0.05; all nonsignificant tests retained in complete downloads', 'downloads':downloads}
    if int(result.get('bootstrap_requested',0))==0:
        payload['module_display_rule']='Association FDR <= 0.05, FDR then P then module ID. Stability not assessed (0 refits); kME is membership, not reproducibility. Qualification flags remain unchanged in the audit.'
    # JSON embedded as inert data. Escape HTML delimiters, including adversarial IDs.
    is_core=manifest['discovery_identity']['method']=='factor_core'
    if is_core:payload['schema']='study-report-2.8'
    if is29:payload['schema']='study-report-2.9'
    if is31:
        payload['schema']='study-report-3.1'
        payload['global_local']={name:records(table(run/'model'/file)) for name,file in
            [('programs','programs.tsv'),('scores','program_scores.tsv'),('links','module_program_links.tsv'),('modules','local_modules.tsv')]}
    payload=json.loads(json.dumps(payload,default=lambda x:x.item()),parse_constant=lambda _:None)
    data=json.dumps(payload,ensure_ascii=False,allow_nan=False).replace('&','\\u0026').replace('<','\\u003c').replace('>','\\u003e')
    template=Path(__file__).with_name('study_report31.html' if is31 else 'study_report29.html' if is29 else 'study_report28.html' if is_core else 'study_report27.html').read_text(encoding='utf-8')
    if read_json(run/'report_revision.json').get('presentation_revision')==39:
        template=Path(__file__).with_name('study_report39.html').read_text(encoding='utf-8')
    css=Path(__file__).with_name('katana_report.css').read_text(encoding='utf-8')
    rendered=template.replace('__KATANA_CSS__',css).replace('__STUDY_TITLE__',html.escape(manifest['study_id'])).replace('__REPORT_DATA__',data)
    rendered=rendered.replace("['No global event cap','Full eligible matrix']",
        "['Network input',D.result.network_events??D.result.eligible_events],['User event limit',D.result.feature_cap==null?'Unlimited':D.result.feature_cap]")
    selection=result.get('feature_selection')
    if selection:
        cap=selection['requested_event_limit']
        selection_html='<aside class="callout" id="event-selection-summary"><strong>Event selection / 建网事件范围</strong><p>User limit: '+('Unlimited' if cap is None else f'{cap:,}')+' · QC eligible: '+f'{selection["qc_eligible_events"]:,}'+' · Entered network: '+f'{selection["network_events"]:,}'+' · Excluded by limit: '+f'{selection["excluded_by_event_limit"]:,}'+'</p><p>Selection follows QC and observed PSI variance; it does not use differential-splicing significance. Events outside the limit remain in the QC audit.</p></aside>'
        rendered=rendered.replace('</header>','</header>'+selection_html,1)
    if (run/'validation37/manifest.json').exists():
        from report37 import panel
        # The association gate also requires calibrated inference; a failure of
        # that gate must not be described as a failed stability measurement.
        rendered=rendered.replace('Stability-qualified','Fully qualified').replace('stability-qualified','fully qualified')
        rendered=rendered.replace('</header>','</header>'+panel(run,manifest,interpretation),1)
    if analysis.get('inference_mode')=='donor_trajectory':
        from report38 import panel, adapt_template
        rendered=adapt_template(rendered)
        rendered=rendered.replace('</header>','</header>'+panel(run,interpretation,associations),1)
    rendered=rendered.replace('<div id="downloads" class="file-list"></div>', '<div id="downloads" class="file-list">'+''.join('<a download href="report_tables/'+html.escape(name)+'">'+html.escape(name)+'</a>' for name in downloads)+'</div>')
    from no_refit_report import adapt as adapt_without_refits
    rendered=adapt_without_refits(rendered,result)
    (destination/'report.html').write_text(rendered,encoding='utf-8')
    (destination/'report_manifest.json').write_text(json.dumps({'schema':payload['schema'],'network_id':manifest['network_id'],
        'samples':len(samples),'conditions':len(manifest['conditions']),'modules':len(modules),'contrasts':analysis.get('contrasts',[]),
        'network_refitted':False,'source_run':str(run.resolve()),'interpretation':str(interpretation.resolve())},indent=2),encoding='utf-8')
    verify_frozen(run)
    if read_json(run/'report_revision.json').get('presentation_revision')==39:
        from chapter6_42 import build as build_chapter6
        build_chapter6(run,destination)
    return payload
