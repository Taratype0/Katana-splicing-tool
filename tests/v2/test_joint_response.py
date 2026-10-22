import pandas as pd
import pytest
from katana_next.workflows.joint_response import classify_response

def test_missing_is_not_negative():
    data = pd.DataFrame(dict(comparison_id=['A']*5, gene_id=list('abcde'),
        as_available=[True]*4+[False], deg_available=[True]*5,
        as_supported=[True,True,False,False,False], deg_supported=[True,False,True,False,True]))
    out, summary = classify_response(data)
    assert out.joint_class.nunique() == 5
    assert out.iloc[-1].joint_class == 'not jointly assessable'
    assert summary.jointly_assessable_genes.eq(4).all()
    assert summary.genes.sum() == 5
    with pytest.raises(ValueError):
        classify_response(pd.concat([data,data]))
