import numpy as np
import pandas as pd

from scripts import finish_motiftrace5712_comparison_20261016 as report


def setup_module():
    report.load_science()


def test_physical_coordinates_preserve_masked_reverse_window():
    assert report.physical_window('x', 'masked_acceptor_3prime_minus20', 0, 50)[:2] == (-70, -20)
    assert report.physical_window('x', 'masked_acceptor_3prime_minus20', 20, 50)[:2] == (-90, -40)
    assert report.physical_window('x', 'masked_donor_5prime_plus6', 0, 50)[:2] == (6, 56)
    assert report.physical_window('x', 'affected_RNA_5prime', -1, 0)[:2] == (None, None)


def test_information_flag_does_not_mislabel_every_other_window_as_powered():
    frame = pd.DataFrame({'nonzero_gene_pairs': [0, 20, 200], 'minimum_attainable_P': [1., 2.**-20, 2.**-200]})
    actual = report.annotate_information(frame, 4e-10)
    assert actual.information_limit.tolist() == ['no_nonzero_pair_differences',
        'cannot_reach_single_isolated_discovery_boundary', 'power_not_estimated_minimum_P_only']


def test_adjacent_windows_merge_only_for_display_same_width_method():
    base = dict(model_id='m',RBP_labels='RBP1;RBP2',event_type='SE',direction='inclusion_down',
        region='inclusion_junction_1_acceptor_intron',method='finite_nonoverlap_counts',window_nt=50,
        coordinate_system='masked_acceptor_3prime_minus20',P_exact=1e-12,q_BH=.01,n_matched_gene_pairs=250,
        nonzero_gene_pairs=200,mean_paired_score_difference=.2)
    frame = pd.DataFrame([{**base,'window_start0':p} for p in [0,1,3]])
    extra = {**base,'window_start0':0,'window_nt':20}
    frame = pd.concat([frame,pd.DataFrame([extra])],ignore_index=True)
    rows = report.merge_display_windows(frame,'q_BH','BH_FDR_05')
    assert len(rows)==3
    long=[r for r in rows if r['window_nt']==50]
    assert long[0]['first_window_start0']==0 and long[0]['last_window_start0']==1
    assert long[0]['displayed_RNA_start_nt']==-71 and long[0]['displayed_RNA_end_nt']==-20
    assert all(r['this_is_not_cluster_significance'] for r in rows)


def test_grouped_common114_adjustment_counts_alias_duplicate_rows():
    import duckdb
    from katana_next.science.motiftrace_exact import full_family_adjust
    values=np.r_[np.repeat(1e-8,4),np.linspace(.01,1,60)]
    frame=pd.DataFrame({'P_exact':values})
    con=duckdb.connect()
    con.execute('SET threads=1')
    con.register('raw',frame)
    report.grouped_adjustment(con,'raw','actual',len(values))
    actual=con.execute('SELECT * FROM actual').df().set_index('p')
    for column,expected in zip(['q114_BH','q114_BY','P114_Holm'],full_family_adjust(values)):
        np.testing.assert_allclose(actual.loc[values,column],expected,rtol=1e-13,atol=1e-15)
    con.close()


def test_full_result_reconstructs_decisions_from_raw_not_unsealed_cache(tmp_path,monkeypatch):
    up,out=tmp_path/'up',tmp_path/'out'
    case='test_case'
    (up/case/'models').mkdir(parents=True)
    (up/case/'layouts').mkdir()
    out.mkdir()
    monkeypatch.setattr(report,'UPSTREAM',up)
    monkeypatch.setattr(report,'OUT',out)
    numbers=pd.DataFrame({'within_method_row_id':[0,1],'method':['PWM_mean_odds']*2,
        'P_exact':[1e-10,.5],'paired_rank_Z':[6.,0.], 'nonzero_gene_pairs':[200,20],
        'minimum_attainable_P':[2.**-200,2.**-20],'mean_paired_score_difference':[.2,0.]})
    numbers.to_parquet(up/case/'models/m.numbers.parquet',index=False)
    layout=pd.DataFrame({'within_method_row_id':[0,1], 'event_type':['SE']*2,'direction':['inclusion_down']*2,
        'region':['inclusion_junction_1_acceptor_intron']*2,'window_nt':[50]*2,
        'coordinate_system':['masked_acceptor_3prime_minus20']*2,'window_start0':[0,1], 'window_end0':[50,51],
        'n_matched_gene_pairs':[250,250],'status':['estimable']*2})
    layout.to_parquet(up/case/'layouts/width_4.parquet',index=False)
    adjust=pd.DataFrame({'p':[1e-10,.5],'q_BH':[.001,1.],'q_BY':[.01,1.],'P_Holm':[.001,1.]})
    path=up/case/'ALL_FULL_FAMILY_ADJUSTMENTS.parquet'
    adjust.to_parquet(path,index=False)
    report.save(up/case/'COMPLETE.json',{'full_adjustments_SHA256':report.digest(path)})
    # No adjusted_discoveries cache is provided: report must rebuild selections.
    models=pd.DataFrame({'model_id':['m'],'kind':['PWM'],'width':[4],'RBPs':['R1;R2'],'sources':['native'], 'evidence':['reference']})
    result=report.full_results(case,models,{'all_estimable_hypotheses':1000000})
    assert len(result)==2
    selected=pd.read_csv(out/case/'SIGNIFICANT_WINDOWS_FULL5712.tsv.gz',sep='\t')
    assert len(selected)==1 and selected.window_start0.iloc[0]==0
    assert selected.q_BH.iloc[0]==.001
    assert selected.physical_RNA_start_nt.iloc[0]==-70
