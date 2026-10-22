"""Experimental count-based splice-competition programs; preserved PSI workflow stays available."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np,pandas as pd
from katana_next.science.splice_structure import fit_competition_model,project_competition_model,co_usage_edges,SCHEMA
from katana_next.science.psi_robust import fit_reliable_predictors,predict_reliable

def seal(model):
    payload={k:v for k,v in model.items() if k!='sha256'}
    payload['sha256']=hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
    return payload

def load(path):
    model=json.loads(Path(path).read_text(encoding='utf-8'))
    if model.get('schema')!=SCHEMA or seal(model)['sha256']!=model['sha256']:raise ValueError('Frozen competition model schema/checksum mismatch')
    return model

def main():
    parser=argparse.ArgumentParser(description=__doc__);subs=parser.add_subparsers(dest='command',required=True)
    for command in ['fit','project']:
        p=subs.add_parser(command)
        for name in ['counts','donor-totals','samples','feature-space-id','out']:p.add_argument('--'+name,required=True)
        p.add_argument('--traits',help='Optional numeric sample-by-trait TSV; training-only selection or held-out evaluation')
        if command=='fit':p.add_argument('--annotation',required=True);p.add_argument('--rank',type=int,default=8)
        else:p.add_argument('--model',required=True)
    args=parser.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    Y=pd.read_csv(args.counts,sep='\t',index_col=0);N=pd.read_csv(args.donor_totals,sep='\t',index_col=0)
    if not Y.index.is_unique or not Y.columns.is_unique or not N.index.equals(Y.index) or not N.columns.equals(Y.columns):raise ValueError('Counts and actual donor totals must match exactly; no invented totals')
    samples=pd.read_csv(args.samples,sep='\t').set_index('sample_id')
    if not samples.index.is_unique or set(samples.index)!=set(Y.columns) or not samples.biological_sample_id.is_unique:raise ValueError('Require one row and unique independent biological unit per sample')
    traits=pd.read_csv(args.traits,sep='\t',index_col=0) if args.traits else None
    if traits is not None and (not traits.index.is_unique or set(traits.index)!=set(Y.columns)):raise ValueError('Trait rows must exactly match samples')
    if args.command=='fit':
        A=pd.read_csv(args.annotation,sep='\t').set_index('event_key').loc[Y.index]
        if 'donor_group' not in A:raise ValueError('Require an explicit structural donor_group annotation')
        model,audit=fit_competition_model(Y,N,Y.index,A.donor_group,rank=args.rank)
        model.update(feature_space_id=args.feature_space_id,training_biological_ids=samples.biological_sample_id.tolist(),
            experimental=True,causal_regulation='not_tested',method_origin='conditional factorial-moment adaptation; measurement/latent separation inspired by CS-CORE and CONCISE')
        pred,z,coverage=project_competition_model(model,Y,N);scores=pd.DataFrame(z,index=Y.columns,columns=['CP%03d'%(i+1) for i in range(z.shape[1])])
        model['trait_predictors']=fit_reliable_predictors(scores,traits.loc[Y.columns]) if traits is not None else []
        model=seal(model);audit.to_csv(out/'moment_audit.tsv',sep='\t',index_label='event_key')
        loadings=pd.DataFrame(model['loadings'],index=Y.index,columns=scores.columns);loadings.to_csv(out/'overlapping_program_loadings.tsv',sep='\t',index_label='event_key')
        edges=co_usage_edges(model);rows=[]
        for i,j in zip(*np.triu_indices(len(Y),1)):
            if 'gene_id' in A and A.gene_id.iloc[i]==A.gene_id.iloc[j]:continue
            rows.append(dict(source=Y.index[i],target=Y.index[j],latent_correlation=edges[i,j],statistical_significance='not_computed'))
        pd.DataFrame(rows).to_csv(out/'cross_gene_associations.tsv.gz',sep='\t',index=False)
    else:
        model=load(args.model)
        if args.feature_space_id!=model.get('feature_space_id'):raise ValueError('Feature definition mismatch')
        if list(Y.index)!=model['feature_ids']:raise ValueError('Counts must match the exact frozen feature order')
        if set(samples.biological_sample_id)&set(model.get('training_biological_ids',[])):raise ValueError('Projection samples overlap training biological units')
        pred,z,coverage=project_competition_model(model,Y,N);scores=pd.DataFrame(z,index=Y.columns,columns=['CP%03d'%(i+1) for i in range(z.shape[1])])
    pd.DataFrame(pred,index=Y.index,columns=Y.columns).to_csv(out/'predicted_competing_usage.tsv.gz',sep='\t')
    scores.to_csv(out/'program_scores.tsv',sep='\t',index_label='sample_id')
    coverage.insert(0,'sample_id',Y.columns);coverage.to_csv(out/'sample_coverage.tsv',sep='\t',index=False)
    records=[]
    for predictor in model.get('trait_predictors',[]):
        prediction,audit=predict_reliable(predictor,scores)
        for i,sample in enumerate(Y.columns):
            records.append(dict(sample_id=sample,trait=predictor['trait'],prediction=float(prediction[i]),
                truth=float(traits.loc[sample,predictor['trait']]) if traits is not None else None,**audit.iloc[i].to_dict()))
    pd.DataFrame(records).to_csv(out/'phenotype_predictions.tsv',sep='\t',index=False)
    (out/'frozen_competition_model.json').write_text(json.dumps(model,indent=2,allow_nan=False),encoding='utf-8')
    summary=dict(status='complete',experimental=True,command=args.command,features=len(Y),samples=len(Y.columns),programs=model['rank'],
        model_sha256=model['sha256'],module_refit=args.command=='fit',regulatory_accuracy='not_tested')
    (out/'COMPLETE.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary))

if __name__=='__main__':main()
