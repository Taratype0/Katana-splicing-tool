"""Evidence-bound donor trajectory inference, separate from legacy omnibus."""
from pathlib import Path
import json,time,zipfile
import numpy as np
import pandas as pd
from donor38 import infer,METHOD

def compatible(cfg):
    return (not cfg.get('traits_path') and not cfg.get('covariates') and cfg.get('contrasts',[])==[]
        and not cfg.get('permutation_strata') and cfg.get('traits',{'condition':'categorical'}) in
        [{'condition':'categorical'},{'time_day':'continuous'}])

def validate_certificate(run):
    from study import sha,verify_frozen
    run=Path(run);folder=run/'inference38';binding=json.loads((folder/'manifest.json').read_text())
    frozen=verify_frozen(run)
    if binding['network_id']!=frozen['network_id'] or binding['frozen_files']!=frozen['frozen_files']:raise ValueError('Inference calibration is bound to another network')
    if sha(Path(__file__).with_name('donor38.py'))!=binding['test_source_sha256']:raise ValueError('Calibrated test source changed')
    for name,h in binding['files'].items():
        if Path(name).name!=name or sha(folder/name)!=h:raise ValueError('Calibration artifact changed')
    for name,h in binding['legacy_files'].items():
        if Path(name).name!=name or sha(folder/'legacy_reference'/name)!=h:raise ValueError('Legacy reference changed')
    calibration=json.loads((folder/'calibration.json').read_text())
    scores=json.loads((folder/'score_calibration.json').read_text())
    if not calibration['completed'] or not calibration['passes'] or not scores['completed'] or not scores['passes']:raise ValueError('Calibration not passed')
    calibrated_sources={k.replace('\\','/'):v for k,v in calibration['sources'].items()}
    if calibration['network_id']!=frozen['network_id'] or calibrated_sources.get('science/donor38.py')!=binding['test_source_sha256']:raise ValueError('Calibration method or network mismatch')
    if calibration['native_null_experiments']!=128 or calibration['native_upper95']>.05 or calibration['native_tested_families']<116:raise ValueError('Insufficient native null calibration')
    if scores['experiments']!=16000 or any(s['upper95']>.05 for s in scores['scenarios']):raise ValueError('Score null calibration failed')
    native=pd.read_csv(folder/'native_calibration.tsv',sep='\t',keep_default_na=False)
    if len(native)!=144 or not native.number.is_unique or not native.status.eq('completed').all() or not native.selected_events.eq(5000).all():raise ValueError('Incomplete full-pipeline calibration')
    null=native[native.kind.eq('null')]
    if len(null)!=128 or int(null.positive_modules.gt(0).sum())!=calibration['native_false_positive_families']:raise ValueError('Null calibration ledger mismatch')
    meta=pd.read_csv(run/'samples.tsv',sep='\t',dtype={'sample_id':str,'biological_sample_id':str}).set_index('sample_id')
    if meta.biological_sample_id.nunique()!=19 or len(meta)!=297 or set(meta.time_day)!=set(range(16)):raise ValueError('Design outside this calibration certificate')
    if frozen['discovery_identity']['method']!='wgcna' or frozen['discovery_identity']['qc'].get('event_cap')!=5000:raise ValueError('Network outside calibrated scope')
    return binding

def calculate(run):
    run=Path(run);scores=pd.read_csv(run/'study/scores.tsv',sep='\t',index_col=0)
    meta=pd.read_csv(run/'samples.tsv',sep='\t',dtype={'sample_id':str,'biological_sample_id':str}).set_index('sample_id').loc[scores.index]
    return infer(scores,meta)

def verify_associations(run,tab):
    try:
        validate_certificate(run);expected,_,_=calculate(run)
        query=tab[tab.test.eq(METHOD)]
        if len(query)!=len(expected) or not query.module.is_unique or set(query.module)!=set(expected.module):raise ValueError('Incomplete calibrated hypothesis family')
        q=query.set_index('module').reindex(expected.module);e=expected.set_index('module')
        if not q.factor.eq('time_day').all() or not q.kind.eq('omnibus').all():raise ValueError('Wrong hypothesis scope')
        for col in ['p_value','q_value','p_fwer','effect']:
            if not np.allclose(q[col],e[col],rtol=1e-10,atol=1e-14,equal_nan=True):raise ValueError('Calibrated statistics mismatch')
        return True,'matching_donor_trajectory_calibration'
    except (ValueError,KeyError,OSError,TypeError) as exc:return False,str(exc)

def analyze(run,cfg,out):
    from study import verify_frozen,sha
    from module_qualification import qualification
    from cosplice_research import write_json
    from study_report import build
    if not compatible(cfg):raise ValueError('Donor trajectory calibration does not support extra traits, covariates, contrasts or replacement metadata')
    run=Path(run);out=Path(out);binding=validate_certificate(run);start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    results,parts,coef=calculate(run)
    results['inference_gate_eligible']=results.status.eq('estimated')
    results['inference_limit']='Within-study empirical calibration; fixed-score exact signs; not universal selective inference or external replication'
    results['calibration_id']=binding['calibration_id'];results['effect_definition']='Median donor linear Legendre coefficient; not total start-to-end change'
    results.to_csv(out/'associations.tsv',sep='\t',index=False);parts.to_csv(out/'trajectory_components.tsv',sep='\t',index=False)
    coef.to_csv(out/'donor_coefficients.tsv',sep='\t',index=False)
    legacy=run/'inference38/legacy_reference'
    for name in ['associations.tsv','summary.json']:
        (out/('legacy_'+name)).write_bytes((legacy/name).read_bytes())
    for name in ['condition_scores.tsv','co_response_diagnostic.tsv','sample_scores.tsv','analysis_samples.tsv']:
        (out/name).write_bytes((legacy/name).read_bytes())
    summary=json.loads((legacy/'summary.json').read_text())
    diagnostics=pd.read_csv(out/'co_response_diagnostic.tsv',sep='\t')
    candidates,gate=qualification(run,results,diagnostics)
    candidates.to_csv(out/'candidate_audit.tsv',sep='\t',index=False)
    summary.update(traits={'time_day':'continuous'},contrasts=[],covariates={},tests_requested=len(results),tests_estimable=int(results.p_value.notna().sum()),
        q05=int(results.q_value.le(.05).sum()),qualified_candidates=gate['qualified'],qualification=gate,
        inference_mode='donor_trajectory',calibration_id=binding['calibration_id'],
        inference_revision='3.8: per-donor cubic coefficients; exact two-sided signs; Bonferroni within three components; BH across all modules plus Holm FWER <=0.05 for qualification',
        multiplicity='Three planned components per module, then ALL module hypotheses. Holm gate additional to unchanged BH FDR <=0.05; no selected-only correction.',
        limits=['Typical-donor polynomial time-response null, not the legacy 15-df equality of all condition means',
            'Fixed-score exact binomial tests; data-adaptive selection calibrated empirically under the documented native simulation generators',
            'Independent donors and noninformative observation/missingness; no support for dependent donors, confounded time effects or causal interpretation',
            'Same-study internal validation, not external confirmation; four directional modules are not four independent programs'],
        interpretation_seconds=time.monotonic()-start,interpretation_source_sha256=sha(__file__),
        effect_definition='Median donor linear Legendre coefficient; component table contains quadratic/cubic effects and exact median intervals')
    write_json(out/'summary.json',summary);write_json(out/'config.json',dict(cfg,inference_mode='donor_trajectory'))
    verify_frozen(run);build(run,out,out)
    with zipfile.ZipFile(out/'results.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in out.rglob('*'):
            if p.is_file() and p.name!='results.zip':z.write(p,p.relative_to(out))
    return summary
