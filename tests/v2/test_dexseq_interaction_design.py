from pathlib import Path
import pandas as pd
import pytest
from katana_next.science.dexseq_interaction import read_design
from katana_next.workflows.dexseq_interaction import DexseqInteractionWorkflow


def design(tmp_path):
    rows=[]
    for a in ['WT','MUT']:
        for b in ['BS','KA']:
            for i in range(2):
                sid=f'{a}_{b}_{i}'
                path=tmp_path/f'{sid}.txt'
                path.write_text('gene:001\t10\ngene:002\t20\n')
                rows.append(dict(sample_id=sid,factor_a=a,factor_b=b,count_file=path.name))
    path=tmp_path/'design.tsv'
    pd.DataFrame(rows).to_csv(path,sep='\t',index=False)
    return path


def test_explicit_design_and_absolute_counts(tmp_path):
    frame=read_design(design(tmp_path),'WT','BS')
    assert len(frame)==8
    assert all(Path(p).is_absolute() for p in frame.count_file)


@pytest.mark.parametrize('problem',['reference','missing_cell','duplicate','outside','reused_counts','paired'])
def test_design_rejects_ambiguous_inputs(tmp_path,problem):
    path=design(tmp_path); frame=pd.read_csv(path,sep='\t')
    ref='WT'; samples=None
    if problem=='reference': ref='unknown'
    if problem=='missing_cell': frame=frame.iloc[:-2]
    if problem=='duplicate': frame.loc[1,'sample_id']=frame.loc[0,'sample_id']
    if problem=='outside': samples=['not_a_real_sample']
    if problem=='reused_counts': frame.loc[1,'count_file']=frame.loc[0,'count_file']
    if problem=='paired': frame['biological_sample_id']='one_donor'
    frame.to_csv(path,sep='\t',index=False)
    with pytest.raises(ValueError): read_design(path,ref,'BS',samples)


def test_official_workflow_has_no_automatic_pair_or_event_selection():
    definition=DexseqInteractionWorkflow.definition
    assert not definition.supports_contrast_selection
    keys={p.key for p in definition.parameters}
    assert 'selection_id' not in keys
    assert {'design_path','reference_a','reference_b'} <= keys
    assert definition.analysis_unit.startswith('DEXSeq exon')
