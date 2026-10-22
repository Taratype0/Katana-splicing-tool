from dataclasses import replace
from pathlib import Path
import json
import numpy as np
import pandas as pd
import pytest

from katana_next.science.hub_events import HubSettings, prepare_inputs, run_hub_analysis, summarize_hubs
from katana_next.science.module_dashboard import build_module_dashboard

RSCRIPT = Path(r'C:\Program Files\R\R-4.5.0\bin\Rscript.exe')


def planted(n=40,p=15):
    rng=np.random.default_rng(8)
    latent=rng.normal(0,.12,n)
    values=np.vstack([.5+latent]+[.5+latent+rng.normal(0,.05+i*.004,n) for i in range(p-1)])
    # An uncorrelated member contributes weak edges that must not be thresholded out of kWithin.
    values[-1]=rng.uniform(.2,.8,n)
    matrix=pd.DataFrame(np.clip(values,.01,.99),index=[f'e{i}' for i in range(p)],columns=[f's{i}' for i in range(n)])
    assignments=pd.DataFrame({'event_key':matrix.index,'module':'M1'})
    samples=pd.DataFrame({'sample_id':matrix.columns,'biological_sample_id':[f'd{i}' for i in range(n)],'condition':['A']*(n//2)+['B']*(n-n//2)})
    events=pd.DataFrame({'event_key':matrix.index,'gene_id':[f'G{i}' for i in range(p)]})
    return matrix,assignments,samples,events


@pytest.fixture(scope='module')
def fitted(tmp_path_factory):
    inputs=planted()
    settings=HubSettings(method='pearson',bootstrap_iterations=20,top_fraction=.2,minimum_selection_frequency=.5)
    directory=tmp_path_factory.mktemp('real_R_hub_analysis')
    metrics,boot,plan=run_hub_analysis(*inputs,settings,directory,RSCRIPT)
    return inputs,settings,directory,metrics,boot,plan


def test_complete_connectivity_equals_unthresholded_adjacency_and_kme_equals_PC1(fitted):
    (matrix,*_),settings,_,metrics,_,_=fitted
    actual=metrics.loc[metrics.condition.eq('all')].set_index('event_key').loc[matrix.index]
    cor=np.corrcoef(matrix.to_numpy()); adj=((1+cor)/2)**settings.power;np.fill_diagonal(adj,0)
    np.testing.assert_allclose(actual.kWithin,adj.sum(axis=1),rtol=1e-10,atol=1e-11)
    thresholded=np.where(adj>=settings.edge_threshold,adj,0).sum(axis=1)
    assert np.max(actual.kWithin.to_numpy()-thresholded)>1e-5
    z=(matrix.to_numpy()-matrix.to_numpy().mean(axis=1,keepdims=True))/matrix.to_numpy().std(axis=1,ddof=1,keepdims=True)
    _,_,vt=np.linalg.svd(z,full_matrices=False); score=vt[0]
    if np.corrcoef(score,z.mean(axis=0))[0,1]<0:score=-score
    expected=[np.corrcoef(row,score)[0,1] for row in z]
    np.testing.assert_allclose(actual.kME,expected,rtol=1e-9,atol=1e-10)


def test_numerical_benchmark_against_unmodified_official_WGCNA_functions(fitted):
    import subprocess
    (matrix,*_),settings,directory,metrics,_,_=fitted
    script=directory/'official_WGCNA_reference.R'
    script.write_text('''suppressPackageStartupMessages(library(WGCNA))
args<-commandArgs(TRUE)
x<-t(as.matrix(read.delim(args[1],row.names=1,check.names=FALSE)))
a<-adjacency(x,type="signed",power=6,corFnc="cor",corOptions=list(use="p"))
k<-intramodularConnectivity(a,rep("blue",ncol(x)))$kWithin
me<-moduleEigengenes(x,rep("blue",ncol(x)),impute=FALSE,subHubs=FALSE)$eigengenes
mm<-signedKME(x,me,corFnc="cor",corOptions="use = 'p'")
write.table(data.frame(event_key=colnames(x),kWithin=k,kME=mm[,1]),args[2],sep="\\t",quote=TRUE,row.names=FALSE)
''',encoding='utf-8')
    reference=directory/'official_WGCNA_reference.tsv'
    result=subprocess.run([str(RSCRIPT),'--vanilla',str(script),str(directory/'hub_input_psi.tsv'),str(reference)],capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    expected=pd.read_csv(reference,sep='\t').set_index('event_key')
    actual=metrics.loc[metrics.condition.eq('all')].set_index('event_key').loc[expected.index]
    np.testing.assert_allclose(actual[['kWithin','kME']],expected[['kWithin','kME']],rtol=1e-9,atol=1e-10)


def test_planted_hub_is_recovered_and_all_bootstraps_and_conditions_are_retained(fitted):
    _,settings,_,metrics,boot,_=fitted
    center=metrics.loc[metrics.event_key.eq('e0')]
    assert center.hub_candidate.all()
    assert set(center.condition)=={'all','A','B'}
    assert set(boot.bootstrap)==set(range(1,settings.bootstrap_iterations+1))
    assert metrics.bootstrap_frequency_all_edges.between(0,1).all()
    assert 'selection_reason' in metrics and metrics.event_key.nunique()==15


def test_official_WGCNA_bicor_path_is_finite_and_retains_missing_measurements(tmp_path):
    matrix,a,s,e=planted()
    matrix.iloc[1,0]=np.nan
    settings=HubSettings(bootstrap_iterations=20,top_fraction=.2,minimum_selection_frequency=.5)
    metrics,_,_=run_hub_analysis(matrix,a,s,e,settings,tmp_path,RSCRIPT,scope='all')
    assert np.isfinite(metrics.kWithin).all() and metrics.kME.notna().all()
    assert pd.read_csv(tmp_path/'hub_input_psi.tsv',sep='\t').isna().any().any()
    assert (tmp_path/'R_session.txt').read_text().find('WGCNA')>=0


def test_same_gene_cluster_loses_cross_gene_hub_status(tmp_path):
    m,a,s,e=planted(p=15)
    # All high-correlation members belong to one gene; unrelated members cannot rescue it.
    e.loc[:13,'gene_id']='same_host'
    settings=HubSettings(method='pearson',bootstrap_iterations=20,top_fraction=.2,minimum_selection_frequency=.5)
    metrics,_,_=run_hub_analysis(m,a,s,e,settings,tmp_path,RSCRIPT,scope='all')
    center=metrics.set_index('event_key').loc['e0']
    assert center.kWithin_no_same_gene < center.kWithin*.1
    assert center.rank_no_same_gene != center.rank_all_edges
    assert center.candidate_all_edges and not center.hub_candidate


def test_block_bootstrap_preserves_donor_condition_profiles_and_is_deterministic():
    m,a,s,e=planted()
    s['biological_sample_id']=[f'd{i%20}' for i in range(40)]
    settings=HubSettings(bootstrap_iterations=20)
    _,_,_,plan,audit=prepare_inputs(m,a,s,e,settings)
    _,_,_,again,_=prepare_inputs(m,a,s,e,settings)
    assert plan==again
    pooled=plan['scopes'][0]
    for boot in pooled['bootstraps']:
        assert len(boot['indices'])==40
        assert all(boot['units'].count(u)==2 for u in set(boot['units']))
    assert len(audit)==40


def test_technical_repeats_do_not_inflate_independent_units(tmp_path):
    m,a,s,e=planted(n=24)
    s['biological_sample_id']=[f'd{i//2}' for i in range(24)]
    metrics,boot,_=run_hub_analysis(m,a,s,e,HubSettings(bootstrap_iterations=20),tmp_path,RSCRIPT,scope='conditions')
    assert metrics.biological_units.eq(6).all() and not metrics.hub_candidate.any()
    assert boot.empty
    assert metrics.selection_reason.str.contains('insufficient_biological_units').all()


def test_no_edges_never_produces_zero_score_hubs(tmp_path):
    m,a,s,e=planted();m.loc[:,:]=.5
    metrics,_,_=run_hub_analysis(m,a,s,e,HubSettings(bootstrap_iterations=20),tmp_path,RSCRIPT,scope='all')
    assert not metrics.hub_candidate.any() and metrics.supported_degree.eq(0).all()
    d=build_module_dashboard(m,a,s,'M1',metadata=e)
    assert d['hub_events']==[]


def test_invalid_bootstraps_count_against_frequency(fitted):
    _,settings,directory,_,_,_=fitted
    raw=pd.read_csv(directory/'hub_raw_metrics.tsv',sep='\t')
    raw.loc[raw.bootstrap.gt(5),['kWithin','kWithin_no_same_gene','kME']]=np.nan
    metrics,_=summarize_hubs(raw,settings)
    assert metrics.bootstrap_frequency_all_edges.max()<=5/settings.bootstrap_iterations
    assert not metrics.hub_candidate.any()


def test_missing_gene_annotation_is_unassessable_not_cross_gene_evidence(tmp_path):
    m,a,s,e=planted();e=e[['event_key']]
    metrics,_,_=run_hub_analysis(m,a,s,e,HubSettings(method='pearson',bootstrap_iterations=20),tmp_path,RSCRIPT,scope='all')
    assert not metrics.hub_candidate.any()
    assert metrics.selection_reason.str.contains('cross_gene_comparison_unavailable').all()


def test_incomplete_modules_and_missing_biological_identity_are_rejected():
    m,a,s,e=planted();settings=HubSettings()
    with pytest.raises(ValueError,match='subset'):prepare_inputs(m.iloc[:-1],a,s,e,settings)
    with pytest.raises(ValueError,match='biological'):prepare_inputs(m,a,s.drop(columns='biological_sample_id'),e,settings)


def test_drawing_same_gene_identity_uses_host_id_not_display_symbol():
    m,a,s,e=planted();e['gene_symbol']=[f'label{i}' for i in range(len(e))]
    e.loc[1,'gene_id']=e.loc[0,'gene_id']
    edges=pd.DataFrame({'event_a':['e0'],'event_b':['e1'],'adjacency':[.99],'correlation':[.99]})
    d=build_module_dashboard(m,a,s,'M1',metadata=e,edges=edges)
    assert d['same_gene_edges']==1 and d['edges'][0]['same_gene']
