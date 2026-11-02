"""Katana-format English report with separate programs and module exports."""
from pathlib import Path
import html,json,shutil,zipfile,re
import numpy as np
import pandas as pd
from common32 import ROOT,sha,write

def table(t,escape=True):
 labels={'biological_sample_id':'Donor','gene_id':'Gene','ensembl_gene_id':'Ensembl gene ID','factor':'Global factor','pve':'Variance fraction',
  'q_value':'FDR','p_value':'P value','scope':'Resolution','median_match_jaccard':'Median event Jaccard','median_module_jaccard':'Median event Jaccard',
  'median_gene_jaccard':'Median gene Jaccard','median_member_support':'Median coassignment','maximum_unit_variance_fraction':'Largest donor variance share',
  'maximum_gene_loading_energy_fraction':'Largest gene loading share','effective_loading_genes':'Effective loading genes',
  'posterior_sign_confident_events':'Sign-confident events','posterior_sign_confident_genes':'Sign-confident genes',
  'qualified_candidate':'Passes candidate gate','qualified_candidates':'Candidates','qualified_with_functional_support':'Candidates with pathway support'}
 display=t.rename(columns={c:labels.get(c,str(c).replace('_',' ').capitalize()) for c in t.columns})
 return '<div class="table-scroll">'+display.to_html(index=False,border=0,classes='data-table',float_format=lambda v:f'{v:.4g}',na_rep='not computed',escape=escape)+'</div>'
def chapter(i,title,body):return f'<section class="chapter" id="chapter-{i}"><div class="chapter-head"><span class="chapter-number">{i:02}</span><h2>{html.escape(title)}</h2></div>{body}</section>'
def document(title,chapters):
 css=(ROOT/'science/katana_report.css').read_text(encoding='utf8')
 toc=''.join(f'<a href="#chapter-{i}"><b>{i:02}</b>{re.search(r"<h2>(.*?)</h2>",c).group(1)}</a>' for i,c in enumerate(chapters,1))
 return '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+html.escape(title)+'</title><style>'+css+'''
 .module-card{padding:18px;border:1px solid #d8e1eb;margin:12px 0;background:white}.module-card summary{cursor:pointer;font-weight:650;color:#164e81}.gene-list{line-height:1.9;overflow-wrap:anywhere}.figure{max-width:100%;height:auto}details{margin:14px 0}.table-scroll{overflow:auto;max-height:600px}.data-table{border-collapse:collapse;width:100%;font-size:12px}.data-table th,.data-table td{text-align:left;padding:9px;border-bottom:1px solid #e3e7eb;max-width:450px;overflow-wrap:anywhere}.data-table th{color:#5b6470;font-size:11px}.notice{border-left:3px solid #c63f4a;padding:15px;background:#fff4f4}.signal{font-size:18px;line-height:1.6;color:#2863a5}
 </style><body><div class="topbar"><span class="wordmark">KATANA</span><small>CoSplice Next · 3.2</small></div><header class="cover"><p class="eyebrow">Whole-study splicing analysis</p><h1>'''+html.escape(title)+'</h1><p class="lead">Count-aware global factor programs, bounded local modules, and measured biological evidence.</p></header><div class="shell"><nav class="toc"><p class="toc-title">REPORT CONTENTS</p>'+toc+'</nav><main>'+''.join(chapters)+'</main></div></body></html>'
def link(name,label,prefix='report_tables/'):
 return '<a href="'+prefix+html.escape(name,quote=True)+'" download>'+html.escape(label)+'</a>'

def membership_exports(events,qualification,dest,prefix='report_tables/'):
 dest=Path(dest);dest.mkdir(parents=True,exist_ok=True);cards=[];manifest=[]
 q=qualification.sort_values(['q_value','p_value','module'],kind='stable')
 for r in q.itertuples():
  mid=str(r.module);m=int(mid[1:]);part=events[events.module.eq(m)].copy()
  gs=part.groupby(['gene_id','ensembl_gene_id'],dropna=False).agg(events=('event_key','size'),event_keys=('event_key',lambda s:';'.join(s))).reset_index()
  ef=f'module_{mid}_events.tsv';gf=f'module_{mid}_genes.tsv';part.to_csv(dest/ef,sep='\t',index=False);gs.to_csv(dest/gf,sep='\t',index=False)
  manifest.append(dict(module=mid,genes=len(gs),events=len(part),genes_file=gf,events_file=ef,q_value=r.q_value,qualified_candidate=bool(r.qualified_candidate)))
  fmt=lambda v:f'{v:.3f}' if np.isfinite(v) else 'not computed'
  metric=f'FDR {r.q_value:.3g} · event Jaccard {fmt(r.median_match_jaccard)} · gene Jaccard {fmt(r.median_gene_jaccard)}'
  reason='Passes operational stability + time-association gate' if r.qualified_candidate else 'Exploratory: does not pass the complete stability + time-association gate'
  biology='<p>Pathway annotation: '+html.escape(str(getattr(r,'leading_pathways','See complete pathway table')))+'.</p>'
  if hasattr(r,'largest_correlated_factor'):biology+=f'<p>Largest score association: global factor {r.largest_correlated_factor}, r={r.factor_score_correlation:+.3f}. All global factors explain {r.variance_explained_by_all_factors:.1%} of this module score. This is an association, not a regulatory assignment.</p>'
  card=f'<details class="module-card"><summary>{mid} · {len(gs)} genes · {len(part)} events · {metric}</summary><p>{reason}</p>'+biology+'<p>'+link(gf,'Complete gene list',prefix)+' · '+link(ef,'Complete event list',prefix)+'</p><div class="gene-list">'+html.escape(', '.join(sorted(gs.gene_id.astype(str).unique())))+'</div></details>'
  cards.append((bool(r.qualified_candidate),bool(r.q_value<=.05),card))
 pd.DataFrame(manifest).to_csv(dest/'module_membership_manifest.tsv',sep='\t',index=False)
 for name,qualified in [('all_module_lists.zip',False),('qualified_module_lists.zip',True)]:
  with zipfile.ZipFile(dest/name,'w',zipfile.ZIP_DEFLATED) as z:
   for row in manifest:
    if not qualified or row['qualified_candidate']:
     for key in ['genes_file','events_file']:z.write(dest/row[key],row[key])
 return cards,manifest

def rank_figure(model,dest):
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 if (model/'parallel_rank.tsv').exists():
  p=pd.read_csv(model/'parallel_rank.tsv',sep='\t');k=int(p.retained.sum());fig,ax=plt.subplots(figsize=(8,3.3))
  ax.plot(p['rank'],p.observed_eigenvalue,color='#175b94',label='Observed spectrum')
  ax.plot(p['rank'],p.null_q95,color='#b43d50',label='Block-sign null, 95th percentile')
  ax.axvline(k,color='#888',ls='--',label=f'Selected budget: {k}')
  ax.set(xlabel='Ordered component',ylabel='Eigenvalue',yscale='log',title='Shared-program rank from structured FlipPA; 999 randomizations')
  ax.legend(fontsize=8);fig.tight_layout();fig.savefig(dest,dpi=150);plt.close(fig);return
 p=pd.read_csv(model/'rank_selection.tsv',sep='\t');fig,ax=plt.subplots(figsize=(8,3.3))
 ax.plot(p['rank'],p['mean'],'o-',color='#175b94',ms=4);r=p[p.selected].iloc[0];ax.axvline(r['rank'],color='#b43c50',ls='--',label=f'Selected budget: {int(r["rank"])}')
 ax.set(xlabel='Number of latent factors',ylabel='Heldout donor-by-gene loss',title='Internal rank selection — no fixed four-factor setting');ax.legend();fig.tight_layout();fig.savefig(dest,dpi=150);plt.close(fig)

def program_figure(model,meta,dest,observation_weights=None):
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 f=pd.read_csv(model/'factor_scores.tsv',sep='\t',index_col=0).reindex(meta.sample_id);k=f.shape[1]
 if observation_weights is not None:
  weights=np.asarray(observation_weights)
  if weights.shape!=(len(meta),) or (weights<=0).any():raise ValueError('Invalid observation weights for program display')
  f=f.div(np.sqrt(weights),axis=0)
 fig,axes=plt.subplots(int(np.ceil(k/3)),min(3,k),figsize=(11,2.7*int(np.ceil(k/3))),squeeze=False)
 for i,ax in enumerate(axes.flat):
  if i>=k:ax.axis('off');continue
  values=f.iloc[:,i].to_numpy();sd=np.std(values);values/=max(sd,1e-12)
  for d,ix in meta.groupby('biological_sample_id').indices.items():
   order=ix[np.argsort(meta.time_day.iloc[ix].to_numpy())];ax.plot(meta.time_day.iloc[order],values[order],alpha=.22,color='#19598b',lw=.8)
  avg=pd.Series(values).groupby(meta.time_day.to_numpy()).mean();ax.plot(avg.index,avg.values,color='#b53d52',lw=2)
  ax.set(title=f'Global factor {i+1}',xlabel='Time',ylabel='Standardized score')
 fig.tight_layout();fig.savefig(dest,dpi=130);plt.close(fig)

def product_report(cfg,out,summary):
 out=Path(out);files=out/'report_tables';files.mkdir(exist_ok=True);model=out/'model'
 events=pd.read_csv(out/'events.tsv',sep='\t');q=pd.read_csv(out/'modules.tsv',sep='\t');meta=pd.read_csv(out/'samples.tsv',sep='\t')
 cards,manifest=membership_exports(events,q,files)
 if (out/'module_program_links.tsv').exists():shutil.copyfile(out/'module_program_links.tsv',files/'module_program_links.tsv')
 if (out/'program_scores.tsv').exists():shutil.copyfile(out/'program_scores.tsv',files/'program_scores.tsv')
 for name in ['global_program_summary.tsv','program_refit_matches.tsv','program_subspace_recovery.tsv','module_refit_recovery.tsv.gz','global_program_gene_scope.tsv','global_program_gene_loadings.tsv.gz']:
  if (out/name).exists():shutil.copyfile(out/name,files/name)
 for name in ['parallel_rank.tsv','parallel_rank.json','parallel_null_eigenvalues.npy']:
  if (model/name).exists():shutil.copyfile(model/name,files/name)
 for src,name in [(out/'modules.tsv','module_qualification.tsv'),(out/'module_scores.tsv','module_scores.tsv'),(out/'donor_time_slopes.tsv','donor_time_slopes.tsv'),
   (out/'events.tsv','gene_event_membership.tsv'),(out/'summary.json','analysis_summary.json'),(model/'rank_selection.tsv','rank_selection.tsv'),(model/'factors.tsv','global_factors.tsv'),
   (model/'factor_scores.tsv','factor_scores.tsv'),(model/'loadings.tsv','factor_loadings.tsv'),(out/'samples.tsv','sample_manifest.tsv')]:shutil.copyfile(src,files/name)
 rank_figure(model,files/'rank_selection.png');program_figure(model,meta,files/'global_trajectories.png',np.load(out/'input/observation_weights.npy'))
 if (out/'pathway_enrichment.tsv').exists():shutil.copyfile(out/'pathway_enrichment.tsv',files/'pathway_enrichment.tsv')
 stable=int(q.stability_pass.sum());sig=int(q.q_value.le(.05).sum());qualified=int(q.qualified_candidate.sum())
 rank_description=('Native flashier EBMF with sparse point-normal loadings. Rank comes from 999 gene-by-donor block-sign randomizations inspired by FlipPA. Entire donor trajectories and all events of a gene share a sign; the original independent-entry theoretical guarantees are not assumed for this adaptation. Native null-check can remove unsupported factors. Local covariance centers fitted scores within each donor, matching input baseline removal and whole-donor sign invariance.' if (model/'parallel_rank.tsv').exists() else 'Historical selector: native flashier EBMF, with rank selected by internal gene-by-donor heldout loss and a paired one-standard-error rule.')
 chapters=[chapter(1,'Study question and measured conclusion',f'<p class="signal">{summary["selected_factors"]} learned global factors; {len(q)} bounded local modules; {qualified} operational candidates pass both stability and time FDR.</p><p>{summary["eligible_events"]:,} QC events across {summary["biological_units"]} independent donor trajectories ({summary["samples"]} observations). All comparable conditions are analysed together.</p><p class="notice">Candidate status is not external biological confirmation. Local loading-pattern communities may share the same global program; module count is not the number of independent biological processes.</p>'),
 chapter(2,'Input, experimental units and count uncertainty',f'<p>{html.escape(cfg.get("study_id",""))}: {html.escape(cfg.get("cohort_scope",""))}</p><p>Competing-junction counts / cluster totals; minimum support 10; coverage and within-donor variation QC; no significance filter or event cap. Within-donor centering removes stable between-donor baselines. Donors receive equal aggregate weight.</p><p>Binomial delta sampling uncertainty enters native EBMF, with an additional learned common noise term. The diagonal likelihood does not fully model within-donor error correlation; validation deletes whole donor trajectories.</p>'+table(meta.groupby('biological_sample_id').agg(observations=('sample_id','size'),first_time=('time_day','min'),last_time=('time_day','max')).reset_index())),
 chapter(3,'Adaptive global programs',f'<p>{rank_description} The safety cap is 64, not a requested output count.</p><img class="figure" src="report_tables/rank_selection.png" alt="Rank selection evidence"><img class="figure" src="report_tables/global_trajectories.png" alt="Every global factor across all donors">'+table(pd.read_csv(out/'global_program_summary.tsv',sep='\t'))+'<p>Program recovery uses one-to-one loading matches across complete donor-deletion fits. It is distinct from local-module membership recovery; no program is certified from this diagnostic alone. Subspace recovery remains available in the complete downloads. Loading-energy concentration reports whether one gene dominates a factor. Posterior sign confidence (local false-sign rate ≤0.05) is not a gene-level association FDR.</p>'+f'<p>Effective rank of local module scores: {summary["effective_local_score_rank"]:.2f}; median score variance explained by global factors: {summary["median_local_score_variance_explained_by_global_factors"]:.1%}. Module count must not be interpreted as this many independent programs.</p>'),
 chapter(4,'Local modules, ordered by association significance',f'<p>{stable} stable modules; {sig} time-FDR-significant modules; {qualified} pass both. Default candidates are filtered by evidence, not by a Top-N display limit. Each module has separate complete gene and event lists.</p>'+(''.join(c for good,s,c in cards if good) or '<p class="notice">No module currently passes the complete candidate gate. Exploratory memberships remain downloadable.</p>')+'<details><summary>All FDR-significant exploratory modules that fail the candidate gate</summary>'+''.join(c for good,s,c in cards if s and not good)+'</details><details><summary>Remaining modules and all independent lists</summary>'+''.join(c for good,s,c in cards if not s)+'</details>'),
 chapter(5,'Whole-donor stability and statistical design',f'<p>{summary["stability"]["successful"]}/{summary["stability"]["requested"]} full-QC refits completed. Missing QC members count as losses. Candidate thresholds: event and gene Jaccard ≥0.5; median member coassignment ≥0.5; at least half of members supported; ≥18/20 successful refits; donor variance share &lt;0.5; time FDR ≤0.05.</p><p>The cubic time-response omnibus signs entire donor trajectories together (49,999 draws). Its null assumes independent, centrally symmetric donor trajectories with no mean response in the three-degree time basis. Days are not independently permuted. BH includes all modules before stability filtering.</p>'+table(q[['module','stability_events','q_value','median_match_jaccard','median_gene_jaccard','median_member_support','maximum_unit_variance_fraction','qualified_candidate']])),
 chapter(6,'Biological interpretation',(f'<p>Correction family: {summary.get("pathway_correction_scope","This analysis")}. Unique-gene Reactome overlap uses all QC genes as background and conditions on event multiplicity (1, 2, 3–4, 5–8, 9+ events/gene), with BH across every module × eligible pathway. Annotation support is not an independent validation cohort or causal edge accuracy. Bulk co-splicing can also reflect changing cell composition; these results do not establish regulation within a single cell type.</p>'+table(pd.read_csv(out/'pathway_enrichment.tsv',sep='\t').query('q_value <= .05 and overlap >= 3').sort_values('q_value'))) if (out/'pathway_enrichment.tsv').exists() else '<p>Pathway analysis: not computed; no gene-set file supplied.</p>'),
 chapter(7,'Reproducibility and complete downloads','<p>'+link('all_module_lists.zip','Every module: separate genes and events')+' · '+link('qualified_module_lists.zip','Qualified candidates only')+'</p>'+''.join('<p>'+link(p.name,p.name)+'</p>' for p in sorted(files.iterdir()) if p.is_file() and not p.name.startswith('module_M'))+'<p>Native algorithms: <a href="https://willwerscheid.github.io/flashier/">flashier / EBMF</a>. The sealed CoSplice release and earlier reports remain separate.</p>')]
 (out/'report.html').write_text(document(cfg.get('name','Adaptive factor CoSplice'),chapters),encoding='utf8')
 with zipfile.ZipFile(out/'results.zip','w',zipfile.ZIP_DEFLATED) as z:
  z.write(out/'report.html','report.html')
  for p in files.iterdir():
   if p.is_file():z.write(p,'report_tables/'+p.name)
