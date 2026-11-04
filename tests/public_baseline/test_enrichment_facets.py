import pandas as pd
import pytest
from katana_next.reports.enrichment_facets import render_enrichment_facets


def test_facets_use_saved_counts_ratios_and_adjusted_p(tmp_path):
    rows=pd.DataFrame([dict(category=c,term=f'{c} term {i}',Count=3+i,GeneRatio=(3+i)/40,
        **{'p.adjust':.001*(i+1) if i<6 else .6}) for c in ['BP','CC','MF','KEGG'] for i in range(7)])
    original=rows.copy(deep=True)
    selected=render_enrichment_facets(rows,tmp_path/'facets.png')
    assert len(selected)==20
    assert selected.category.drop_duplicates().tolist()==['BP','CC','MF','KEGG']
    assert (selected['p.adjust']<=.05).all()
    pd.testing.assert_frame_equal(rows,original)
    for row in selected.to_dict('records'):
        source=original.loc[(original.category==row['category'])&(original.term==row['term'])].iloc[0]
        assert all(row[key]==source[key] for key in ['Count','GeneRatio','p.adjust'])
    svg=(tmp_path/'facets.svg').read_text('utf8')
    assert all(x in svg for x in ['Gene Ratio','p.adjust','Counts','BP','CC','MF','KEGG'])


def test_unknown_ontology_is_not_guessed(tmp_path):
    frame=pd.DataFrame([{'category':'Unclassified','term':'RNA processing','Count':4,'GeneRatio':.1,'p.adjust':.01}])
    with pytest.raises(ValueError,match='explicitly recorded'):
        render_enrichment_facets(frame,tmp_path/'bad.png')


def test_empty_facet_is_reported_without_inventing_terms(tmp_path):
    frame=pd.DataFrame([{'category':'BP','term':'nonsignificant','Count':4,'GeneRatio':.1,'p.adjust':.8}])
    result=render_enrichment_facets(frame,tmp_path/'empty.png')
    assert result.empty
    assert 'No significant terms' in (tmp_path/'empty.svg').read_text('utf8')
