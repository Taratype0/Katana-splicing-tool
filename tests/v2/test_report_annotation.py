import pandas as pd
from katana_next.science.report_annotation import annotate_exons


def test_composite_exons_keep_identity_and_statistics():
    data=pd.DataFrame({'groupID':['ENSMUSG00001.2+ENSMUSG00002.1','ENSMUSG00003.1'], 'interaction_FDR':[.02,.8]})
    names=pd.DataFrame({'gene_id':['ENSMUSG00001.3','ENSMUSG00002.1'],'gene_name':['A','B']})
    out=annotate_exons(data,names)
    assert out.gene_name.tolist()==['A; B','ENSMUSG00003.1']
    pd.testing.assert_frame_equal(out[data.columns],data)
