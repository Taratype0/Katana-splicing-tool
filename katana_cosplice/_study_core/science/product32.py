"""Production count-aware adaptive factor engine for whole-study trajectories."""
from pathlib import Path
import json,time,traceback,shutil,zipfile
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
import pandas as pd
from common32 import ROOT,prepare,write,sha
import masks32
from run_native33 import fit
from graph33 import run as partition
from statistics32 import balanced_scores,trajectory_tests,gene_recovery
from module_qualification import summarize_modules

def input_source(cfg,out):
 from product import read_matrix
 c=read_matrix(cfg['counts']);s=read_matrix(cfg['totals'])
 if not c.index.equals(s.index) or not c.columns.equals(s.columns):raise ValueError('Counts/total identities and order must match')
 for t in [c,s]:
  a=t.to_numpy()
  if not np.isfinite(a).all() or (a<0).any() or (a!=np.floor(a)).any():raise ValueError('Counts/totals must be finite nonnegative integers; use zero total for unobserved cells')
 if (c>s).any().any():raise ValueError('Junction numerator exceeds competing-cluster total')
 m=pd.read_csv(cfg['samples'],sep='\t',dtype={'sample_id':str,'biological_sample_id':str})
 required={'sample_id','biological_sample_id','condition','time_day'}
 if not required<=set(m) or m[list(required)].isna().any().any() or not m.sample_id.is_unique:raise ValueError('Unique samples require sample_id, biological_sample_id, condition and numeric time_day')
 if set(m.sample_id)!=set(c.columns):raise ValueError('Sample manifest and count columns must match exactly')
 m=m.set_index('sample_id').loc[c.columns].reset_index().rename(columns={'index':'sample_id'});m.time_day=pd.to_numeric(m.time_day,errors='raise')
 if not np.isfinite(m.time_day).all():raise ValueError('Time values must be finite')
 if m.groupby('biological_sample_id').time_day.nunique().lt(4).any():raise ValueError('This trajectory engine requires at least four distinct times per donor')
 if m.duplicated(['biological_sample_id','time_day']).any():raise ValueError('Merge technical replicates before fitting; donor/time pairs must be unique')
 e=pd.read_csv(cfg['annotation'],sep='\t',dtype=str)
 if not {'event_key','gene_id'}<=set(e) or not e.event_key.is_unique:raise ValueError('One annotation row per event with event_key and gene_id is required')
 if not set(c.index)<=set(e.event_key):raise ValueError('Every event needs a unique gene annotation')
 e=e.set_index('event_key').loc[c.index].reset_index().rename(columns={'index':'event_key'})
 if e.gene_id.isna().any() or e.gene_id.str.contains('[;,|]',regex=True).any():raise ValueError('Ambiguous/missing gene mappings must be resolved before analysis')
 if 'ensembl_gene_id' not in e:e['ensembl_gene_id']=e.gene_id
 if e.ensembl_gene_id.isna().any():raise ValueError('Missing stable gene identifiers')
 u=np.divide(c.to_numpy(),s.to_numpy(),out=np.full(c.shape,np.nan),where=s.to_numpy()>=10)
 out.mkdir(exist_ok=True);e.to_csv(out/'events.tsv',sep='\t',index=False);m.to_csv(out/'samples.tsv',sep='\t',index=False)
 for name,a in [('usage',u),('counts',c.to_numpy()),('support',s.to_numpy())]:np.save(out/(name+'.npy'),a.astype(np.float32 if name=='usage' else np.float64))
 write(out/'manifest.json',dict(input_files={k:dict(path=str(cfg[k]),sha256=sha(cfg[k])) for k in ['counts','totals','samples','annotation']},events=len(e),samples=len(m)))
 return e,m,u,c.to_numpy(),s.to_numpy()

def qualify(out,input,model,refits,annotation):
 ev=pd.read_csv(input/'events.tsv',sep='\t');meta=pd.read_csv(input/'samples.tsv',sep='\t',dtype={'biological_sample_id':str})
 lab=pd.read_csv(model/'graph/hierarchical/target200/assignments.tsv',sep='\t',index_col=0).module
 mods,members=summarize_modules(lab,refits);mods=mods.merge(gene_recovery(lab,refits,annotation),on='module')
 x=np.fromfile(input/'input.bin',dtype='<f8').reshape(len(ev),len(meta))/np.sqrt(np.load(input/'observation_weights.npy'))[None,:]
 scores=balanced_scores(x,ev.ensembl_gene_id.to_numpy(),lab.reindex(ev.event_key).to_numpy(int));scores.index=meta.sample_id
 tests,slopes=trajectory_tests(scores,meta);tab=tests.merge(mods,on='module')
 tab['stability_pass']=tab.median_match_jaccard.ge(.5)&tab.median_gene_jaccard.ge(.5)&tab.median_member_support.ge(.5)&tab.supported_member_fraction.ge(.5)&tab.successful_refits.ge(18)&tab.requested_refits.ge(20)
 tab['qualified_candidate']=tab.stability_pass&tab.q_value.le(.05)&tab.maximum_unit_variance_fraction.lt(.5)
 from recovery_details32 import recovery_details,verify_aggregates
 detail=recovery_details(lab,refits,annotation);verify_aggregates(detail,tab)
 detail.to_csv(out/'module_refit_recovery.tsv.gz',sep='\t',index=False)
 tab=tab.sort_values(['q_value','p_value','module'],kind='stable')
 tab.to_csv(out/'modules.tsv',sep='\t',index=False)
 scores.to_csv(out/'module_scores.tsv',sep='\t');slopes.to_csv(out/'donor_time_slopes.tsv',sep='\t')
 events=ev.merge(lab.rename('module'),left_on='event_key',right_index=True).merge(members[['event_key','within_module_coassignment']],on='event_key',how='left')
 events.to_csv(out/'events.tsv',sep='\t',index=False);events[['event_key','module','gene_id']].to_csv(out/'module_genes.tsv',sep='\t',index=False)
 shutil.copyfile(input/'event_qc.tsv.gz',out/'event_audit.tsv.gz');shutil.copyfile(model/'loadings.tsv',out/'factor_loadings.tsv')
 shutil.copyfile(model/'factor_scores.tsv',model/'native_factor_scores.tsv');scores.to_csv(model/'scores.tsv',sep='\t')
 return tab,events,members

def finalize(cfg,out,input,model,refits,annotation,input_count,benchmark_pathways=None,native_refit_models=None):
 tab,events,members=qualify(out,input,model,refits,annotation)
 tab['gene_count']=tab.module.map({f'M{m}':g.ensembl_gene_id.nunique() for m,g in events[events.module.gt(0)].groupby('module')})
 if cfg.get('gene_sets'):
  from enrichment32 import enrich
  if benchmark_pathways is None:bio=enrich(events,cfg['gene_sets'],out/'pathway_enrichment.tsv')
  else:
   bio=benchmark_pathways.copy();bio['module']=bio.module.str.removeprefix('M').astype(int);bio['genes']=bio.contributing_genes;bio['pathway_background_genes']=bio.term_genes
   bio.to_csv(out/'pathway_enrichment.tsv',sep='\t',index=False)
  sig=bio[bio.q_value.le(.05)&bio.overlap.ge(3)]
  tab['functional_support']=tab.module.isin('M'+sig.module.astype(str));tab['qualified_with_functional_support']=tab.qualified_candidate&tab.functional_support
 else:tab['functional_support']=False;tab['qualified_with_functional_support']=False
 tab.to_csv(out/'modules.tsv',sep='\t',index=False)
 edges=pd.read_csv(model/'graph/edges.tsv.gz',sep='\t');keys=events.event_key.to_numpy();edges['source']=keys[edges.a];edges['target']=keys[edges.b]
 edges[['source','target','weight']].to_csv(out/'edges.tsv',sep='\t',index=False)
 fullmeta=pd.read_csv(input/'samples.tsv',sep='\t',dtype={'biological_sample_id':str});shape=json.loads((model/'graph/hierarchical/target200/selection.json').read_text())
 summary=dict(version='3.2.0-development',method='factor_adaptive',analysis_scope='whole_study',value_kind='usage',
  input_events=input_count,eligible_events=len(events),assigned_events=int(events.module.gt(0).sum()),samples=len(fullmeta),biological_units=fullmeta.biological_sample_id.nunique(),
  selected_factors=len(pd.read_csv(model/'factors.tsv',sep='\t')),local_modules=len(tab),qualified_candidates=int(tab.qualified_candidate.sum()),
  qualified_with_functional_support=int(tab.qualified_with_functional_support.sum()),shape=shape,feature_cap=None,
  stability=dict(requested=len(refits),successful=sum(r is not None for r in refits),median_support=float(members.within_module_coassignment.median()),events_support_ge05=int(members.within_module_coassignment.ge(.5).sum())),
  interpretation='Candidate modules; gene-based pathway annotation is not independent or causal validation. Whole-donor cubic trajectory inference requires symmetric donor errors.')
 summary['pathway_correction_scope']='All reported benchmark methods and resolutions' if benchmark_pathways is not None else 'Every module and eligible pathway in this analysis'
 scores=pd.read_csv(out/'module_scores.tsv',sep='\t',index_col=0).to_numpy();F=pd.read_csv(model/'factor_scores.tsv',sep='\t',index_col=0).to_numpy()/np.sqrt(np.load(input/'observation_weights.npy'))[:,None]
 F=np.column_stack([np.ones(len(F)),F]);pred=F@np.linalg.lstsq(F,scores,rcond=None)[0]
 cor=np.atleast_2d(np.corrcoef(scores,rowvar=False));eig=np.linalg.eigvalsh(cor)
 summary['effective_local_score_rank']=float(eig.sum()**2/(eig*eig).sum())
 summary['median_local_score_variance_explained_by_global_factors']=float(np.median(1-((scores-pred)**2).sum(axis=0)/(scores*scores).sum(axis=0)))
 summary['local_pairs_abs_r_ge08_fraction']=float(np.mean(np.abs(cor[np.triu_indices(len(cor),1)])>=.8)) if len(cor)>1 else None
 from interpret_modules32 import interpret
 scoreframe=pd.read_csv(out/'module_scores.tsv',sep='\t',index_col=0);factorframe=pd.read_csv(model/'factor_scores.tsv',sep='\t',index_col=0)
 factorframe=factorframe.div(np.sqrt(np.load(input/'observation_weights.npy')),axis=0)
 factorframe.to_csv(out/'program_scores.tsv',sep='\t',index_label='sample_id')
 tab,links=interpret(tab,bio if cfg.get('gene_sets') else None,scoreframe,factorframe)
 tab=tab.sort_values(['q_value','p_value','module'],kind='stable')
 tab.to_csv(out/'modules.tsv',sep='\t',index=False);links.to_csv(out/'module_program_links.tsv',sep='\t',index=False)
 from program_recovery32 import program_recovery
 program_recovery(input,model,native_refit_models or [],out)
 from program_gene_scope32 import program_gene_scope
 program_gene_scope(input,model,out)
 write(out/'summary.json',summary)
 from report32 import product_report
 product_report(cfg,out,summary)
 return summary

def run(cfg,out):
 out=Path(out);start=time.time();status=dict(status='running',stage='validate_counts_and_donor_design',method='factor_adaptive',started=start);write(out/'status.json',status)
 try:
  source=input_source(cfg,out/'source');e,m,u,c,s=source;donors=sorted(m.biological_sample_id.unique());requested=int(cfg.get('bootstrap',20))
  retain=int(np.floor(.8*len(donors)))
  if requested and retain<15:raise ValueError('At least 19 donors required for 80 percent whole-donor stability refits; fit-only mode needs 15')
  masks32.DATA=out/'source';input=out/'input';prepare(donors,input,source=source);masks32.prepare_masks(input)
  status['stage']='adaptive_factor_selection_and_fit';write(out/'status.json',status)
  model=out/'model';result=fit('ebmf_counts',input,model)
  if result['status']!='completed':raise RuntimeError('Native adaptive factor model failed; see model/engine.log')
  partition('factor_adaptive',input,model,model/'graph');refs=[];rng=np.random.default_rng(32041);draws=[]
  while len(draws)<requested:
   selected=sorted(rng.choice(donors,retain,replace=False).tolist())
   if selected not in draws:draws.append(selected)
  write(out/'stability_plan.json',dict(draws=draws,seed=32041,kind='whole donor deletion with full event QC and rank reselection'))
  refs=[None]*len(draws)
  def refit_one(i,selected):
   dest=out/'stability'/f'refit_{i:03}';inp=dest/'input';prepare(selected,inp,source=source);masks32.prepare_masks(inp)
   r=fit('ebmf_counts',inp,dest/'model')
   if r['status']=='completed':
    partition('factor_adaptive',inp,dest/'model',dest/'model/graph');return pd.read_csv(dest/'model/graph/hierarchical/target200/assignments.tsv',sep='\t',index_col=0).module
   return None
  status.update(stage='whole_donor_refits',refits_completed=0,refits_requested=requested);write(out/'status.json',status)
  with ThreadPoolExecutor(max_workers=2) as pool:
   futures={pool.submit(refit_one,i,d):i for i,d in enumerate(draws)}
   for completed,future in enumerate(as_completed(futures),1):
    refs[futures[future]]=future.result()
    status['refits_completed']=completed;write(out/'status.json',status)
  # Keep actual un-imputed observed usage for event inspection.
  selected=pd.read_csv(input/'events.tsv',sep='\t').source_event_row.to_numpy()
  pd.DataFrame(u[selected],index=e.event_key.iloc[selected],columns=m.sample_id).to_csv(out/'usage.tsv.gz',sep='\t',index_label='event_key');m.to_csv(out/'samples.tsv',sep='\t',index=False)
  status['stage']='inference_and_report';write(out/'status.json',status);finalize(cfg,out,input,model,refs,e,len(e),native_refit_models=[out/'stability'/f'refit_{i:03}/model' for i in range(len(draws))])
  status.update(status='completed',stage='complete',seconds=time.time()-start);write(out/'status.json',status)
 except Exception as exc:
  traceback.print_exc();status.update(status='failed',error=f'{type(exc).__name__}: {exc}',seconds=time.time()-start);write(out/'status.json',status);raise
