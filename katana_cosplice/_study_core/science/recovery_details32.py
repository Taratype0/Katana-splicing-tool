"""Export auditable per-refit membership recovery; never select or modify modules."""
import numpy as np
import pandas as pd

def recovery_details(base,refits,annotation):
 genes=annotation.set_index('event_key').ensembl_gene_id
 groups={int(m):set(v.index) for m,v in base[base.gt(0)].groupby(base[base.gt(0)])}
 gene_groups={m:set(genes.reindex(list(v)).dropna()) for m,v in groups.items()}
 rows=[]
 for i,ref in enumerate(refits):
  if ref is None:
   rows.extend(dict(refit=i,module=f'M{m}',fit_success=False,event_jaccard=np.nan,gene_jaccard=np.nan) for m in groups)
   continue
  r={int(m):set(v.index) for m,v in ref[ref.gt(0)].groupby(ref[ref.gt(0)])}
  rg={m:set(genes.reindex(list(v)).dropna()) for m,v in r.items()};available=set(ref.index)
  for m,events in groups.items():
   def best(g,others):
    values=[(len(g&v)/len(g|v),-n,len(g&v),len(v)) for n,v in others.items()]
    return max(values) if values else (0.,0,0,0)
   je,em,eo,es=best(events,r);jg,gm,go,gs=best(gene_groups[m],rg)
   rows.append(dict(refit=i,module=f'M{m}',fit_success=True,events=len(events),genes=len(gene_groups[m]),
    events_missing_refit_qc=len(events-available),event_jaccard=je,event_best_match=f'M{-em}' if em else None,
    event_overlap=eo,event_matched_module_size=es,gene_jaccard=jg,gene_best_match=f'M{-gm}' if gm else None,
    gene_overlap=go,gene_matched_module_size=gs))
 columns=['refit','module','fit_success','events','genes','events_missing_refit_qc','event_jaccard','event_best_match',
          'event_overlap','event_matched_module_size','gene_jaccard','gene_best_match','gene_overlap','gene_matched_module_size']
 return pd.DataFrame(rows,columns=columns)

def verify_aggregates(detail,qualification):
 if detail.empty:return
 med=detail[detail.fit_success].groupby('module')[['event_jaccard','gene_jaccard']].median()
 q=qualification.set_index('module').reindex(med.index)
 np.testing.assert_allclose(med.event_jaccard,q.median_match_jaccard,atol=1e-12,rtol=1e-10)
 np.testing.assert_allclose(med.gene_jaccard,q.median_gene_jaccard,atol=1e-12,rtol=1e-10)
