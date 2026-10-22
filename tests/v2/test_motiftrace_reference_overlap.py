import numpy as np
import pandas as pd
import pytest

from katana_next.science.motiftrace_reference_overlap import audit_windows
from katana_next.science.rmaps_sliding_comparison import window_design


def inputs(tmp_path):
    def events(n, start):
        return pd.DataFrame([dict(event_id=f"event_{start}_{i}",gene_id=f"g{i}",gene_symbol=f"G{i}",chr="chr1",strand="+",
            exonStart_0base=start,exonEnd=start+60,upstreamES=start-100,upstreamEE=start-50,
            downstreamES=start+300,downstreamEE=start+350) for i in range(n)])
    groups={"dn":events(20,300),"bg":events(21,1000)}
    sequences={"dn":{"DownstreamIntron":["TTT"+"A"*247]*20},"bg":{"DownstreamIntron":["A"*250]*21}}
    np.savez(tmp_path/"m000.npz",dn=np.zeros((20,1000)),bg=np.zeros((21,1000)))
    windows=pd.DataFrame([dict(contrast="TEST_ONLY",window_id="5:0:inclusion_down:0",RBP="Test",pattern="TTT",
        region_index=5,native_position=0,direction="inclusion_down",region="DownstreamIntron",
        window_start_RNA_relative=0,window_end_RNA_relative=50,n_foreground=20,n_background=21,motif_index=0,
        original_author_P=.001,original_author_q=.02,original_MotifTrace_P=.01,original_MotifTrace_q=.2)])
    return windows,groups,sequences


def test_no_reference_is_unknown_not_zero_binding(tmp_path):
    w,g,s=inputs(tmp_path)
    summary,events,hits=audit_windows(w,g,s,tmp_path,{})
    assert summary.foreground_CLIP.isna().all()
    assert summary.external_CLIP_raw_P.isna().all()
    assert events.CLIP_overlap.isna().all()
    assert hits.experimental_peak_overlaps_motif_in_window.isna().all()
    assert len(hits)==20 and set(hits.actual_motif_RNA)=={"UUU"}
    assert summary.iloc[0].original_MotifTrace_q==.2


def test_actual_strand_and_motif_overlap_are_separate(tmp_path):
    w,g,s=inputs(tmp_path)
    ref={"Test":{"label":"synthetic_not_real", "quality":"test_only", "rows":pd.DataFrame([dict(chr="chr1",start=360,end=361,strand="+")])}}
    summary,events,hits=audit_windows(w,g,s,tmp_path,ref)
    assert summary.iloc[0].foreground_CLIP==20 and summary.iloc[0].background_CLIP==0
    assert summary.iloc[0].foreground_motif_CLIP==20
    assert hits.experimental_peak_overlaps_motif_in_window.all()
    ref["Test"]["rows"].loc[0,"strand"]="-"
    summary,_,_=audit_windows(w,g,s,tmp_path,ref)
    assert summary.iloc[0].foreground_CLIP==0
    assert summary.iloc[0].external_CLIP_raw_P==1


@pytest.mark.parametrize("mutation",["denominator","geometry","eligibility"])
def test_frozen_opportunities_cannot_be_changed(tmp_path,mutation):
    w,g,s=inputs(tmp_path)
    if mutation=="denominator":w.loc[0,"n_foreground"]-=1
    elif mutation=="geometry":w.loc[0,"window_end_RNA_relative"]+=1
    else:
        design=window_design()
        column=int(design[design.region_index.eq(5)&design.native_position.eq(0)].index[0])
        scores=np.zeros((20,1000));scores[0,column]=np.nan
        np.savez(tmp_path/"m000.npz",dn=scores,bg=np.zeros((21,1000)))
    with pytest.raises(ValueError):audit_windows(w,g,s,tmp_path,{})
