from types import SimpleNamespace
import json
import numpy as np
import pandas as pd
import pytest
from katana_next.science.module_expression import correlate_modules
from katana_next.workflows.registry import WorkflowRegistry
from katana_next.workflows.expression_response import ExpressionResponseWorkflow

def test_registered_primary_and_exploratory_are_distinct():
    registry=WorkflowRegistry()
    assert 'DEXSeq' in registry.get('differential.replicate_interaction').definition.method
    assert 'Beta-binomial' in registry.get('differential.exploratory_event_interaction').definition.method
    assert registry.get('differential.expression_response')
    assert registry.get('programs.explain_expression')

def test_condition_adjustment_and_sample_identity():
    rng=np.random.default_rng(12); names=[f's{i}' for i in range(24)]
    groups=pd.Series(np.repeat(['a','b','c','d'],6),index=names)
    shift=np.repeat([0,10,20,30],6)
    scores=pd.DataFrame({'C01':shift+rng.normal(size=24)},index=names)
    features=pd.DataFrame({'expression':shift+rng.normal(size=24)},index=names)
    result=correlate_modules(scores,features,groups).set_index('association')
    assert result.loc['overall','rho']>.8
    assert abs(result.loc['condition_adjusted','rho'])<.8
    assert result.loc['condition_adjusted','residual_df']==19
    with pytest.raises(ValueError):correlate_modules(scores,features.iloc[:-1],groups)

def test_expression_rejects_repeated_biological_samples(tmp_path):
    d=tmp_path/'design.tsv'
    pd.DataFrame({'sample_id':['s1','s2'],'biological_sample_id':['same','same'],
        'factor_a':['reference','target'],'factor_b':['reference','target']}).to_csv(d,sep='\t',index=False)
    config={k:str(d) for k in ('design','tximport','gene_names','gmt')}
    p=tmp_path/'config.json';p.write_text(json.dumps(config))
    with pytest.raises(ValueError,match='biological'):
        ExpressionResponseWorkflow()._config(SimpleNamespace(),{'config_path':str(p)})


def test_expression_requires_readable_design_labels(tmp_path):
    names=[f's{i}' for i in range(8)]
    d=tmp_path/'design.tsv'
    pd.DataFrame({'sample_id':names,'biological_sample_id':names,
        'factor_a':['reference']*4+['target']*4,
        'factor_b':['reference','reference','target','target']*2}).to_csv(d,sep='\t',index=False)
    config={k:str(d) for k in ('design','tximport','gene_names','gmt')}
    p=tmp_path/'config.json';p.write_text(json.dumps(config))
    context=SimpleNamespace(manifest=SimpleNamespace(samples=[SimpleNamespace(sample_id=s) for s in names]))
    with pytest.raises(ValueError,match='labels required'):
        ExpressionResponseWorkflow()._config(context,{'config_path':str(p)})
