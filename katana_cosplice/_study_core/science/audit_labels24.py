"""Clarify non-estimable QC labels without changing masks or eligibility.

Primary benchmark uses the frozen numerical unit_audit implementation. Its
undefined-share rows must not blame an arbitrary first subject in product UI.
"""
import numpy as np

def explain_audit(audit):
    out=audit.copy();known=np.isfinite(out.maximum_unit_variance_fraction)
    out.loc[~known,'highest_contribution_unit']=None
    reasons=[]
    for r in out.itertuples():
        why=[]
        if r.observed_fraction<.8:why.append('missingness_or_read_support')
        if not np.isfinite(r.variance):why.append('variance_not_estimable')
        elif r.variance<=1e-12:why.append('constant')
        if not np.isfinite(r.maximum_unit_variance_fraction):why.append('unit_influence_not_estimable')
        elif r.maximum_unit_variance_fraction>=.5:why.append('unit_dominated_variation')
        reasons.append(';'.join(why) or 'passed')
    out['reason']=reasons
    return out
