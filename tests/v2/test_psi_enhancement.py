import numpy as np
import pandas as pd
import pytest
from katana_next.science.frozen_psi import _digest, project_psi_modules, load_frozen_model, save_frozen_model
from katana_next.science.psi_prediction import fit_multimodule_predictor, predict_multimodule


def test_observed_loading_score_recovers_rank_one_under_missingness():
    t = np.linspace(-1, 1, 31)
    p = pd.DataFrame(.5 + .2*np.tile(t, (20, 1)), index=[f'e{i}' for i in range(20)], columns=[f'n{i}' for i in range(31)])
    model = dict(schema='katana-frozen-psi-1.0',feature_space_id='x',training_samples=['train'],training_biological_ids=['train'],trait_models=[],
        modules=[dict(module='M1',events=list(p.index),imputation_means=[.5]*20,centers=[.5]*20,scales=[.2]*20,
            loading=[1/np.sqrt(20)]*20,singular_value=np.sqrt(20),training_mean_pairwise_r=1)])
    model['sha256']=_digest(model)
    p.iloc[:3]=np.nan
    meta=pd.DataFrame(dict(sample_id=p.columns,biological_sample_id=p.columns))
    old=project_psi_modules(model,p,meta,feature_space_id='x')
    new=project_psi_modules(model,p,meta,feature_space_id='x',projection_method='observed_loading_ls')
    assert np.allclose(new['scores'].M1,t,atol=1e-12)
    assert np.allclose(old['scores'].M1,.85*t,atol=1e-12)
    # A new patient's score does not depend on other patients in the target batch.
    one=project_psi_modules(model,p.iloc[:,[4]],meta.iloc[[4]],feature_space_id='x',projection_method='observed_loading_ls')
    assert one['scores'].iloc[0,0] == pytest.approx(new['scores'].iloc[4,0])


def test_loading_gate_rejects_loss_of_a_dominant_event():
    p=pd.DataFrame(.5,index=[f'e{i}' for i in range(20)],columns=['new'])
    p.iloc[0]=np.nan
    load=np.sqrt([.9]+[.1/19]*19)
    model=dict(schema='katana-frozen-psi-1.0',feature_space_id='x',training_samples=[],training_biological_ids=[],trait_models=[],
        modules=[dict(module='M1',events=list(p.index),imputation_means=[.5]*20,centers=[.5]*20,scales=[.2]*20,
            loading=load.tolist(),singular_value=1,training_mean_pairwise_r=None)])
    model['sha256']=_digest(model)
    meta=pd.DataFrame(dict(sample_id=['new'],biological_sample_id=['new']))
    out=project_psi_modules(model,p,meta,feature_space_id='x',projection_method='observed_loading_ls')
    assert np.isnan(out['scores'].iloc[0,0])
    assert out['sample_coverage'].score_status.iloc[0]=='insufficient_loading_information'


def test_joint_prediction_recovers_two_independent_components_and_roundtrips(tmp_path):
    rng=np.random.default_rng(12)
    x=pd.DataFrame(rng.normal(size=(250,3)),columns=['M1','M2','M3'],index=[f't{i}' for i in range(250)])
    y=pd.DataFrame({'joint':2*x.M1-3*x.M2+rng.normal(0,.1,len(x))},index=x.index)
    entries=fit_multimodule_predictor(x,y)
    z=pd.DataFrame(rng.normal(size=(100,3)),columns=x.columns)
    pred=predict_multimodule(entries[0],z)
    assert np.mean((pred-(2*z.M1-3*z.M2))**2)<.05
    assert entries[0]['alpha'] is not None
    model=dict(schema='katana-frozen-psi-1.1',multimodule_trait_models=entries)
    model['sha256']=_digest(model)
    save_frozen_model(model,tmp_path/'model.json')
    assert load_frozen_model(tmp_path/'model.json')==model
    z.loc[0,'M1']=np.nan
    assert np.isnan(predict_multimodule(entries[0],z)[0])


def test_predictor_does_not_use_validation_traits_or_modify_inputs():
    rng=np.random.default_rng(1)
    x=pd.DataFrame(rng.normal(size=(60,2)),index=[f't{i}' for i in range(60)])
    y=pd.DataFrame({'y':rng.normal(size=60)},index=x.index)
    before=x.copy(deep=True)
    first=fit_multimodule_predictor(x,y)
    assert first==fit_multimodule_predictor(x,y)
    pd.testing.assert_frame_equal(x,before)
    with pytest.raises(ValueError,match='identical unique'):
        fit_multimodule_predictor(x,y.iloc[:-1])
