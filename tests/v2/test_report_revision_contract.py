from types import SimpleNamespace as NS
import pandas as pd
import numpy as np
from katana_next.science.report_revision import module_contrasts, transcript_usage_reliability
from katana_next.science.bulk import STABLE_PRESETS, apply_stable_event_qc

def contrast():
    return NS(contrast_id="A_vs_B", numerator=NS(label="A",sample_ids=["a1","a2","a3"]),
              denominator=NS(label="B",sample_ids=["b1","b2","b3"]))

def test_module_contrast_is_matched_not_rest():
    scores=pd.DataFrame({"C01":[3,4,5,0,1,2,999]},index=["a1","a2","a3","b1","b2","b3","other"])
    result=module_contrasts(scores,[contrast()]).iloc[0]
    assert result.effect == 3
    assert result.ci_low < result.effect < result.ci_high
    assert result.n_numerator == 3

def test_usage_not_abundance_and_no_gene_fallback():
    records=[]
    for sample in ["a1","a2","a3","b1","b2","b3"]:
        for tx,tpm in [("t1",80),("t2",20)]:
            records.append(dict(sample_id=sample,gene_id="g",gene_symbol="G",transcript_id=tx,tpm=tpm*(2 if sample.startswith("a") else 1)))
    dtu=pd.DataFrame([dict(comparison_id="A_vs_B",source_kind="dtu_result",transcript_id="t1",dtu_fdr=float("nan"),gene_screen_fdr=.001)])
    _,r=transcript_usage_reliability(pd.DataFrame(records),[contrast()],dtu)
    assert r.delta_fraction.abs().max() < 1e-10
    assert not r.DTU_supported.any()

def test_balanced_reads_are_ten():
    assert STABLE_PRESETS["balanced"]["min_reads"] == 10

def test_approximate_se_does_not_exclude_by_default():
    row=dict(IncLevel1='0.3,0.3',IncLevel2='0.6,0.6',IJC_SAMPLE_1='5,5',SJC_SAMPLE_1='5,5',
        IJC_SAMPLE_2='5,5',SJC_SAMPLE_2='5,5',numerator_replicates=2,denominator_replicates=2,
        replicate_pairwise_direction_consistency=1.,leave_one_out_direction_consistency=1.,max_measurement_se=.4)
    frame=pd.DataFrame([row])
    assert apply_stable_event_qc(frame,{}).stable_event.iloc[0]
    assert not apply_stable_event_qc(frame,{'apply_measurement_se_gate':True}).stable_event.iloc[0]

def test_stager_gene_screen_is_not_transcript_confirmation(tmp_path):
    from katana_next.science.bulk import load_dtu_results
    path=tmp_path/'stage.tsv'
    pd.DataFrame([dict(geneID='G',txID='T',gene=.001,transcript=float('nan'))]).to_csv(path,sep='\t',index=False)
    context=NS(manifest=NS(metadata={},assets=[NS(kind='dtu_result',path=path,metadata={'contrast_id':'A_vs_B'})]))
    r=load_dtu_results(context,['A_vs_B']).iloc[0]
    assert pd.isna(r.dtu_fdr)
    assert r.gene_screen_fdr==.001

def test_motif_refuses_pooled_conditions(tmp_path,monkeypatch):
    import pytest
    import katana_next.workflows.motif as module
    path=tmp_path/'events.parquet'
    pd.DataFrame({'comparison_id':['A_vs_B','C_vs_D']}).to_parquet(path)
    monkeypatch.setattr(module,'selected_source',lambda *a:NS(path=path))
    monkeypatch.setattr(module,'selection_rows',lambda context,parameters,source,frame:frame)
    with pytest.raises(ValueError,match='exactly one comparison'):
        module.EventMotifWorkflow().run(None,None,'test',{})

def test_dtu_background_includes_screened_nonsignificant_genes():
    from katana_next.science.report_revision import dtu_background_genes
    frame=pd.DataFrame({'gene_id':['hit','screened','missing'], 'source_kind':['dtu_result']*3,
        'dtu_fdr':[.01,np.nan,np.nan],'gene_screen_fdr':[.01,1.,np.nan]})
    assert dtu_background_genes(frame)=={'hit','screened'}
