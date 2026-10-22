import pandas as pd

from katana_next.science.enrichment import canonical_gene_mapping, fill_gene_symbols, map_gene_sets
from katana_next.science.enrichment import enrichment_background


def test_gtf_fills_id_only_rows_without_changing_identity_or_background():
    source = pd.DataFrame({'row_id':['tx1','tx2','tx3'],
        'gene_id':['ENSMUSG000001.2','ENSMUSG000002','ENSMUSG000003'],
        'gene_symbol':['',None,'Existing']}, index=[2,4,7])
    annotation = pd.DataFrame({'gene_id':['ENSMUSG000001','ENSMUSG000001','ENSMUSG000002','ENSMUSG999999'],
        'gene_symbol':['Ptbp1','Ptbp1','Srsf7','NotInBackground']})
    filled = fill_gene_symbols(source, annotation)
    assert filled.row_id.equals(source.row_id) and filled.gene_id.equals(source.gene_id)
    assert filled.gene_symbol.tolist() == ['Ptbp1','Srsf7','Existing']
    assert filled.gene_symbol_source.tolist() == ['project_gtf','project_gtf','source']
    canonical, aliases, _ = canonical_gene_mapping(filled)
    mapped, _ = map_gene_sets({'splicing':{'Ptbp1','Srsf7','NotInBackground'}}, aliases)
    assert mapped['splicing'] == {'ENSMUSG000001','ENSMUSG000002'}
    assert len(set(canonical)) == 3
    assert source.gene_symbol.iloc[0] == ''


def test_ambiguous_project_annotation_is_not_guessed():
    source = pd.DataFrame({'gene_id':['G1','G2'], 'gene_symbol':['G1','Authoritative']})
    annotation = pd.DataFrame({'gene_id':['G1','G1','G2'], 'gene_symbol':['One','Two','Other']})
    filled = fill_gene_symbols(source, annotation)
    assert filled.gene_symbol.tolist() == ['G1','Authoritative']
    assert filled.gene_symbol_source.tolist() == ['unresolved','source']


def test_background_excludes_failed_models_not_nonsignificant_events():
    frame = pd.DataFrame({'gene_id':['A','B','C','D'], 'model_status':['fit','fit','optimizer_failed','fit'],
        'pvalue':[0.001, 0.9, float('nan'), float('inf')]})
    kept, audit = enrichment_background(frame)
    assert kept.gene_id.tolist() == ['A','B']
    assert audit.background_eligible.tolist() == [True,True,False,False]
    assert len(audit) == len(frame)


def test_dexseq_only_evidence_remains_in_background():
    frame = pd.DataFrame({'as_available':[False,False], 'dexseq_available':[True,False]})
    kept, audit = enrichment_background(frame)
    assert kept.index.tolist() == [0]
    assert audit.background_reason.iloc[1] == 'no_available_evidence'
