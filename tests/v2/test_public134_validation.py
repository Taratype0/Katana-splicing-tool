import json
import shutil

import numpy as np
import pandas as pd

from scripts import public_motiftrace134_worker_20261016 as worker
from scripts import public_motiftrace5712_20261016 as common
from katana_next.science.motiftrace_full_execution import null_cluster_max


def functions(tmp_path):
    shutil.copytree(common.LOCAL.with_name('public_kd_clip_5712_v4')/'bundle/snapshots',tmp_path/'snapshots')
    shutil.copyfile(common.LOCAL.with_name('public_kd_clip_5712_v4')/'bundle/public_worker.py',tmp_path/'frozen_public.py')
    for name in ['motiftrace_events.py','motiftrace_clusters.py','se_rbp.py']:
        shutil.copyfile(common.REPO/'katana_next/science'/name,tmp_path/'snapshots'/name)
    return worker.core(tmp_path)[1]


def test_compressed_cluster_null_matches_dense_with_NA_and_shared_signs(tmp_path):
    api=functions(tmp_path)
    rng=np.random.default_rng(77)
    diff=rng.integers(-3,4,size=(24,37)).astype(float)
    diff[:,6:9]=np.nan
    diff[:,20:23]=0
    diff[:,13:17]=diff[:,12:13]
    signs=rng.choice([-1,1],size=(129,24))
    z,n,actual=worker.cluster_values(diff,signs,api)
    ranks,expected_z,expected_n=api['signed_ranks'](diff)
    dense=signs@ranks
    dense[:,~np.isfinite(expected_z)]=np.nan
    expected=null_cluster_max(dense,api['CLUSTER_Z'])
    np.testing.assert_allclose(actual,expected,rtol=1e-13,atol=1e-13)
    np.testing.assert_allclose(z,expected_z,equal_nan=True)
    np.testing.assert_array_equal(n,expected_n)


def test_actual_core_minus_strand_junction_gain_is_projected_losslessly(tmp_path):
    api=functions(tmp_path)
    rna=pd.DataFrame([{'event_id':'real_contract','event_type':'SE','gene_id':'g','comparison_id':'case',
        'region':form,'layer':'mature_local_form','chr':'chr1','strand':'-',
        'blocks_json':json.dumps(blocks),'junctions_json':json.dumps(joins),
        'sequence':word,'sequence_RNA':word.replace('T','U'),'length':len(word)}
        for form,blocks,joins,word in [('inclusion',[[200,206],[100,106]],[6],'ACTTCCACTTCC'),
                                      ('exclusion',[[100,106]],[],'ACTTCC')]])
    events=pd.DataFrame([{'event_id':'real_contract','event_type':'SE','gene_id':'g','dpsi':.2,
                          'IncLevel1':'.6,.7','IncLevel2':'.2,.3','reverse_direction':False}])
    model={'kind':'consensus','pattern':'CAC','width':3,'RBPs':'fixture'}
    hits,changes,weighted=api['trace_instances'](rna,events,{'m':model})
    assert len(hits)==worker.assert_hit_projection(hits,rna)==1
    assert hits.crosses_splice_join.all()
    assert changes.instance_state.tolist()==['inclusion_only']
    assert changes.delta_sequence_opportunity.tolist()==[.2]


def test_below20_independent_pairs_is_not_a_negative_result(tmp_path):
    api=functions(tmp_path)
    z,n,maximum=worker.cluster_values(np.ones((19,60)),np.ones((20,19)),api)
    assert np.isnan(z).all()
    assert not maximum.any()
    np.testing.assert_array_equal(n,np.repeat(19,60))
