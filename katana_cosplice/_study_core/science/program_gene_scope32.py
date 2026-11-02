"""Describe whether a global factor is broad or concentrated in a few genes."""
from pathlib import Path
import numpy as np,pandas as pd

def program_gene_scope(input,model,out):
 input,model,out=map(Path,[input,model,out]);out.mkdir(parents=True,exist_ok=True)
 ev=pd.read_csv(input/'events.tsv',sep='\t').set_index('event_key')
 l=pd.read_csv(model/'loadings.tsv',sep='\t',index_col=0)
 sign=pd.read_csv(model/'loading_lfsr.tsv',sep='\t',index_col=0)
 assert l.index.equals(sign.index) and set(l.index)<=set(ev.index)
 genes=ev.loc[l.index,['gene_id','ensembl_gene_id']];records=[];gene_rows=[]
 for j in range(l.shape[1]):
  frame=genes.copy();frame['loading_energy']=l.iloc[:,j].to_numpy()**2
  frame['events']=1;frame['posterior_sign_confident_events']=sign.iloc[:,j].le(.05).astype(int).to_numpy()
  g=frame.groupby(['gene_id','ensembl_gene_id']).sum().reset_index();total=g.loading_energy.sum()
  g['energy_fraction']=g.loading_energy/total if total>0 else 0.;g['factor']=j+1
  g=g.sort_values(['energy_fraction','ensembl_gene_id'],ascending=[False,True]);lead=g.iloc[0]
  records.append(dict(factor=j+1,leading_loading_gene=lead.gene_id,maximum_gene_loading_energy_fraction=lead.energy_fraction,
    effective_loading_genes=float(1/(g.energy_fraction**2).sum()) if total>0 else 0.,
    posterior_sign_confident_events=int(g.posterior_sign_confident_events.sum()),
    posterior_sign_confident_genes=int(g.posterior_sign_confident_events.gt(0).sum())))
  gene_rows.append(g)
 result=pd.DataFrame(records);result.to_csv(out/'global_program_gene_scope.tsv',sep='\t',index=False)
 pd.concat(gene_rows,ignore_index=True).to_csv(out/'global_program_gene_loadings.tsv.gz',sep='\t',index=False)
 summary=out/'global_program_summary.tsv'
 if summary.exists():
  s=pd.read_csv(summary,sep='\t');s=s.drop(columns=[n for n in result if n!='factor' and n in s])
  s.merge(result,on='factor',validate='one_to_one').to_csv(summary,sep='\t',index=False)
 return result
