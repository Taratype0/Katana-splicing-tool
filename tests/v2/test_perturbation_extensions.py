from types import SimpleNamespace
import numpy as np
import pandas as pd
import pytest
from katana_next.science.perturbation import alignment_evidence, gene_exons
from katana_next.workflows.joint_pathways import joint_pathways
from katana_next.workflows.registry import WorkflowRegistry
from katana_next.workflows.perturbation import PerturbationWorkflow


def read(flag=0,nh=1):
    return SimpleNamespace(flag=flag,mapq=255,pos=10,cigar=[(0,10),(3,20),(0,10)],get_tags=lambda:[('NH',nh)])


def test_cigar_coverage_and_duplicate_sensitivity():
    depth,j,q=alignment_evidence([read(),read(1024),read(256),read(nh=2)],0,60)
    assert depth[10:20].tolist()==[2]*10
    assert depth[20:40].sum()==0
    assert depth[40:50].tolist()==[2]*10
    assert j[(20,40,'primary')]==2
    assert j[(20,40,'nonduplicate_anchored')]==1
    assert q['retained_alignments']==2


def test_zero_coverage_not_deletion_call():
    depth,j,q=alignment_evidence([],0,10)
    assert not depth.any() and not j
    assert 'knockout' not in q


def test_annotation_exact_gene_and_unquoted_exon(tmp_path):
    p=tmp_path/'a.gtf'
    p.write_text('chr1\tx\texon\t11\t20\t.\t+\t.\tgene_id "G.1"; gene_name "ABC"; transcript_id "T"; exon_number 5;\n')
    result=gene_exons(p,'ABC')
    assert result.start.iloc[0]==10 and result.exon_number.iloc[0]==5
    with pytest.raises(ValueError):gene_exons(p,'AB')


def test_joint_backgrounds_and_missing_are_distinct():
    f=pd.DataFrame(dict(comparison_id=['x']*4,gene_id=list('abcd'),gene_symbol=list('ABCD'),deg_available=[True]*4,as_available=[True,True,False,False],deg_supported=[True,False,False,False],as_supported=[True,False,False,False]))
    r,a,_,_=joint_pathways(f,{'lib / one':{'A','B'},'lib / two':{'C','D'}})
    assert a.set_index('layer').loc['splicing','tested_genes']==2
    assert r[r.layer.eq('expression')].background_genes.eq(4).all()
    assert r[r.layer.eq('splicing')].background_genes.eq(2).all()
    assert r[r.layer.eq('expression')].shape[0]==2
    assert r[r.term.eq('lib / two')].pvalue.eq(1).all()
    f['as_available']=False
    _,a,_,_=joint_pathways(f,{'lib / one':{'A','B'}})
    assert a[a.layer.eq('splicing')].status.str.startswith('unavailable').all()


def test_production_registry_new_modules():
    registry=WorkflowRegistry()
    for name in ['validate.perturbation','interpretation.joint_pathways','validate.replicate_robustness']:
        assert registry.get(name).definition.enabled


def test_single_project_contrast_scope_and_missing_bam(tmp_path):
    contrast=SimpleNamespace(contrast_id='a',numerator=SimpleNamespace(sample_ids=['s1']),denominator=SimpleNamespace(sample_ids=['s2']))
    context=SimpleNamespace(manifest=SimpleNamespace(contrasts=[contrast],samples=[SimpleNamespace(sample_id=s,biological_sample_id=s,assets={}) for s in ['s1','s2']]))
    # Workflow service canonicalizes all-selected to an omitted contrast_ids key.
    with pytest.raises(ValueError,match='Missing BAM/index'):
        PerturbationWorkflow()._inputs(context,{})


def test_all_multimappers_do_not_yield_false_junction_zero_evidence():
    depth,j,q=alignment_evidence([read(nh=2)],0,60)
    assert not j and not depth.any() and q.get('retained_alignments',0)==0
