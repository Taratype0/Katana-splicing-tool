"""Donor-level cubic trajectory sign tests; no observation-wise replication.

H0: every prespecified within-donor polynomial coefficient has median zero
(balanced positive/negative signs conditional on being nonzero). This is not
the old 15-df equality-of-condition-means null. Exactness is for fixed scores;
the adaptive network procedure additionally needs empirical null calibration.
"""
import numpy as np
import pandas as pd
from scipy.stats import binom, beta
from statsmodels.stats.multitest import multipletests

METHOD='Donor-level exact sign tests of three prespecified cubic trajectory coefficients'
MIN_DONORS=15
DEGREES=3

def coefficients(scores,meta):
    if not scores.index.equals(meta.index):raise ValueError('Scores and metadata must have identical sample order')
    days=pd.to_numeric(meta.time_day,errors='raise').to_numpy(float)
    if not np.isfinite(days).all() or np.ptp(days)<=0:raise ValueError('Finite varying time_day required')
    z=2*(days-days.min())/np.ptp(days)-1
    design=np.polynomial.legendre.legvander(z,DEGREES)
    rows=[]
    for donor,ix in meta.groupby('biological_sample_id',sort=True).indices.items():
        for module in scores:
            y=scores[module].iloc[ix].to_numpy(float);ok=np.isfinite(y)
            if ok.sum()<6 or np.linalg.matrix_rank(design[ix][ok])<4:continue
            coef=np.linalg.lstsq(design[ix][ok],y[ok],rcond=None)[0]
            for degree,value in enumerate(coef[1:],1):
                rows.append(dict(donor=str(donor),module=module,degree=degree,coefficient=float(value),observations=int(ok.sum())))
    return pd.DataFrame(rows,columns=['donor','module','degree','coefficient','observations'])

def test_coefficients(coef,modules):
    rows=[];components=[]
    for module in modules:
        ps=[];valid=True
        for degree in range(1,4):
            values=coef.loc[coef.module.eq(module)&coef.degree.eq(degree),'coefficient'].to_numpy(float)
            values=values[np.isfinite(values)]
            # A fixed numerical zero band does not depend on significance.
            positives=int((values>1e-12).sum());negatives=int((values< -1e-12).sum());n=positives+negatives
            p=float(min(1.,2*binom.cdf(min(positives,negatives),n,.5))) if n else 1.
            enough=n>=MIN_DONORS;valid &= enough
            ps.append(p if enough else 1.)
            median=float(np.median(values)) if len(values) else np.nan
            # Distribution-free central order-statistic interval for a median.
            ordered=np.sort(values);nn=len(ordered);cut=int(binom.ppf(.025,nn,.5)) if nn else 0
            low=float(ordered[cut-1]) if cut>=1 else -np.inf
            high=float(ordered[nn-cut]) if cut>=1 else np.inf
            components.append(dict(module=module,degree=degree,donors=len(values),nonzero_donors=n,
                positive_donors=positives,negative_donors=negatives,zero_donors=len(values)-n,p_value=p,
                median_coefficient=median,median_ci_low=low,median_ci_high=high,estimable=enough))
        p_module=min(1.,3*min(ps)) if valid else np.nan
        linear=components[-3]
        rows.append(dict(module=module,kind='omnibus',factor='time_day',a=None,b=None,
            effect=linear['median_coefficient'],ci_low=linear['median_ci_low'],ci_high=linear['median_ci_high'],
            p_value=p_module,status='estimated' if valid else 'NA',biological_units=linear['donors'],
            test=METHOD,df_num=3,within_module_correction='Bonferroni over all three planned polynomial coefficients',
            inference_gate_eligible=False,inference_limit='Requires matching completed empirical calibration receipt',
            null_hypothesis='All three donor-specific time coefficients have balanced signs about zero; typical-donor trajectory, not equality of all condition means',
            reason='' if valid else 'Fewer than 15 nonzero donor coefficients in a planned component'))
    tab=pd.DataFrame(rows);tab['q_value']=np.nan;tab['p_fwer']=np.nan
    ok=tab.p_value.notna()
    if ok.any():
        # Unestimable hypotheses stay in the family as p=1, not removed.
        p=tab.p_value.fillna(1).to_numpy()
        tab.loc[ok,'q_value']=multipletests(p,method='fdr_bh')[1][ok]
        tab.loc[ok,'p_fwer']=multipletests(p,method='holm')[1][ok]
    return tab.sort_values(['q_value','p_value','module'],kind='stable',na_position='last'),pd.DataFrame(components)

def infer(scores,meta):
    coef=coefficients(scores,meta);tab,parts=test_coefficients(coef,list(scores))
    tab['samples']=tab.module.map(scores.notna().sum())
    return tab,parts,coef

def upper_binomial(k,n):return float(beta.ppf(.95,k+1,n-k)) if k<n else 1.
