"""User-selected WGCNA event budget, after QC and before imputation.

QC eligibility and network inclusion remain separate audit fields. Ranking is
recomputed using training observations in every refit; no condition P values.
"""
from pathlib import Path
import numpy as np
import pandas as pd

METHODS = {'wgcna', 'wgcna_robust', 'wgcna_enhanced'}

def select_events(measured, audit, feature_cap=None):
    if feature_cap is not None and (isinstance(feature_cap, bool) or
            not isinstance(feature_cap, (int, np.integer)) or feature_cap < 40):
        raise ValueError('Event limit must be an integer >= 40 or null for unlimited')
    audit = audit.copy()
    eligible = audit.index[audit.eligible]
    ranks = pd.DataFrame({'variance': measured.reindex(eligible).var(axis=1),
                          'event_key': eligible.astype(str)}, index=eligible)
    ranks = ranks.sort_values(['variance', 'event_key'], ascending=[False, True], kind='stable')
    selected = ranks.index if feature_cap is None else ranks.index[:feature_cap]
    audit['qc_eligible'] = audit.eligible
    audit['variance_rank'] = pd.Series(np.arange(1, len(ranks)+1), index=ranks.index).reindex(audit.index)
    audit['network_selected'] = audit.index.isin(selected)
    audit['selection_reason'] = np.where(~audit.eligible, 'failed_qc',
        np.where(audit.network_selected, 'selected', 'outside_user_event_limit'))
    # Preserve the input order; changing the budget must not silently reorder
    # already eligible events or change the unlimited engine's input ordering.
    return measured.loc[measured.index.isin(selected)], audit

def selection_summary(audit, cap):
    return dict(requested_event_limit=cap, qc_eligible_events=int(audit.eligible.sum()),
        network_events=int(audit.network_selected.sum()),
        excluded_by_event_limit=int((audit.eligible & ~audit.network_selected).sum()),
        ranking='Observed sample PSI variance (ddof=1) after existing QC; event ID breaks ties; no significance selection')

def refit_stability(cfg, raw, meta, support, assignment, keys, out, iterations):
    from product import load_inputs
    from validation24 import deletion_draws, summarize
    from validation22 import labels_at
    from cosplice_research import fit, write_json
    out=Path(out);out.mkdir(exist_ok=False)
    draws=deletion_draws(meta,iterations,24031)
    write_json(out/'planned_draws.json',dict(training_sample_draws=draws,seed=24031,
        feature_cap=cfg.get('feature_cap'),selection_refitted=True))
    records=[]
    for i,cols in enumerate(draws):
        d=out/f'refit_{i:03}';d.mkdir()
        raw.loc[:,cols].to_csv(d/'raw.tsv.gz',sep='\t',index_label='event_key')
        meta.set_index('sample_id').loc[cols].reset_index().to_csv(d/'samples.tsv',sep='\t',index=False)
        sub=dict(cfg,psi=str(d/'raw.tsv.gz'),samples=str(d/'samples.tsv'),counts=None,totals=None,psi_sources=[],read_support=None)
        if support is not None:
            support.loc[:,cols].to_csv(d/'support.tsv.gz',sep='\t',index_label='event_key');sub['read_support']=str(d/'support.tsv.gz')
        try:
            observed,complete,audit,means,m,genes,expression=load_inputs(sub)
            audit.to_csv(d/'event_audit.tsv.gz',sep='\t',index_label='event_key')
            result=fit(cfg['method'],complete,m,d/'model')
            result.update(selection_summary(audit,cfg.get('feature_cap')))
        except Exception as exc:
            result=dict(status='failed',error=f'{type(exc).__name__}: {exc}')
        write_json(d/'status.json',result)
        labels=labels_at(d/'model',keys) if result['status']=='completed' else None
        records.append((i,labels,result))
    return summarize(assignment.reindex(keys).module.to_numpy(),keys,records,out)
