import hashlib
import pandas as pd
import pytest
from katana_next.reports.sources import read_table

@pytest.mark.parametrize('filters',[None,[('comparison_id','=','B')],[[('comparison_id','=','A'),('value','>=',0)],[('comparison_id','=','B')]]])
def test_bounded_parquet_scan_preserves_native_values_types_order_and_projection(tmp_path,monkeypatch,filters):
    import pyarrow.dataset as ds
    frame=pd.DataFrame({'comparison_id':pd.Categorical(['A','B','A','B','C']),
                        'value':[float('inf'),5e-310,-1.2,float('nan'),float('-inf')],
                        'gene':pd.Series(['g1',None,'g3','g4','g5'],dtype='string')})
    path=tmp_path/'native.parquet';frame.to_parquet(path,index=False,row_group_size=2)
    before=path.read_bytes()
    expected=pd.read_parquet(path,columns=['gene','value'],filters=filters)
    actual=read_table(path,columns=['gene','value'],filters=filters)
    pd.testing.assert_frame_equal(actual,expected)
    assert path.read_bytes()==before

def test_report_scan_disables_parallel_readahead(tmp_path,monkeypatch):
    import pyarrow.dataset as ds
    p=tmp_path/'native.parquet';pd.DataFrame({'comparison_id':['A','B'],'fdr':[0.1,5e-310]}).to_parquet(p,index=False)
    original=ds.dataset; observed={}
    class Bounded:
        def __init__(self,source):self.source=source
        def scanner(self,**kwargs):
            observed.update(kwargs);return self.source.scanner(**kwargs)
    monkeypatch.setattr(ds,'dataset',lambda *a,**kw:Bounded(original(*a,**kw)))
    assert len(read_table(p,filters=[('comparison_id','=','B')]))==1
    assert observed['use_threads'] is False
    assert observed['batch_readahead']==observed['fragment_readahead']==1
    assert observed['batch_size']==16384
