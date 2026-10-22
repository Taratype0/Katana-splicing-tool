"""Opt-in fit/project interface for the experimental count copula backend."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np,pandas as pd
from katana_next.science.splice_copula import fit,predict,edge_scores,SCHEMA

def digest(m):return hashlib.sha256(json.dumps(m,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='action',required=True)
    for name in ['fit','project']:
        p=sub.add_parser(name);p.add_argument('--counts',required=True);p.add_argument('--totals',required=True);p.add_argument('--output',required=True);p.add_argument('--feature-space',required=True)
        p.add_argument('--context',help='Sample-by-covariate TSV, identically ordered named columns; training scaling must already be frozen')
        if name=='fit':
            p.add_argument('--annotation',required=True,help='TSV indexed by event_key with donor_group');p.add_argument('--robust',action='store_true');p.add_argument('--nested',action='store_true');p.add_argument('--midpoint',action='store_true')
        else:p.add_argument('--model',required=True)
    a=parser.parse_args();out=Path(a.output);out.mkdir(parents=True,exist_ok=False)
    Y=pd.read_csv(a.counts,sep='\t',index_col=0);N=pd.read_csv(a.totals,sep='\t',index_col=0)
    if not Y.index.equals(N.index) or not Y.columns.equals(N.columns):raise ValueError('Counts/totals require identical feature and sample order')
    if not Y.index.is_unique or not Y.columns.is_unique:raise ValueError('Unique feature and biological sample IDs required; aggregate technical runs first')
    X=pd.read_csv(a.context,sep='\t',index_col=0).loc[Y.columns] if a.context else None
    if a.action=='fit':
        A=pd.read_csv(a.annotation,sep='\t',index_col=0).loc[Y.index]
        model,audit=fit(Y,N,Y.index,A.donor_group,context=None if X is None else X.to_numpy(),exact=not a.midpoint,robust=a.robust,nested=a.nested)
        model.pop('sha256');model.update(feature_space=a.feature_space,training_sample_ids=list(Y.columns),context_names=[] if X is None else list(X.columns));model['sha256']=digest(model)
        (out/'frozen_model.json').write_text(json.dumps(model,indent=2));(out/'fit_audit.json').write_text(json.dumps(audit,indent=2))
        pd.DataFrame(audit['raw_correlation'],index=Y.index,columns=Y.index).to_csv(out/'conditional_co_usage.tsv.gz',sep='\t')
        pd.DataFrame(edge_scores(model),index=Y.index,columns=Y.index).to_csv(out/'partial_associations.tsv.gz',sep='\t')
        if model.get('core'):pd.DataFrame(edge_scores(model,residual=True),index=Y.index,columns=Y.index).to_csv(out/'residual_partial_associations.tsv.gz',sep='\t')
    else:
        model=json.loads(Path(a.model).read_text());sha=model.pop('sha256')
        if digest(model)!=sha:raise ValueError('Frozen model checksum mismatch')
        if model['schema']!=SCHEMA or model['feature_space']!=a.feature_space:raise ValueError('Frozen feature definition mismatch')
        if set(model['training_sample_ids'])&set(Y.columns):raise ValueError('Projection sample IDs overlap training; supply true biological identities')
        if model['feature_ids']!=list(Y.index):raise ValueError('Frozen feature identity/order mismatch')
        if model['context_names']!=([] if X is None else list(X.columns)):raise ValueError('Frozen context columns mismatch')
        prediction,audit=predict(model,Y.to_numpy(),N.to_numpy(),context=None if X is None else X.to_numpy())
        pd.DataFrame(prediction,index=Y.index,columns=Y.columns).to_csv(out/'usage_predictions.tsv.gz',sep='\t');(out/'projection_audit.json').write_text(json.dumps(audit,indent=2))
    (out/'COMPLETE.json').write_text(json.dumps(dict(status='complete',action=a.action,schema=SCHEMA,experimental=True)))
if __name__=='__main__':main()
