"""Experimental depth-conditioned splice copula; opt-in, never the default.

Count margins use the measured donor total, not ranks of raw junction counts.
Pair association is Gaussian rectangle likelihood (heterocop's DD mechanism)
with deterministic quadrature. Robust precision and a nested residual stage
are explicit ablations; neither is labelled a reproduction of MVTLASSO or
Nested-WGCNA. Associations do not identify causal regulators.
"""
from __future__ import annotations
import hashlib, json
import warnings
from functools import lru_cache
from collections import OrderedDict
import numpy as np
from scipy.special import expit, logit, ndtr, ndtri
from scipy.stats import betabinom, binom, beta
from scipy.optimize import minimize, minimize_scalar
from sklearn.covariance import graphical_lasso
from sklearn.exceptions import ConvergenceWarning

SCHEMA = 'katana-splice-conditional-copula-0.1'
_INTERVAL_CACHE=OrderedDict()

@lru_cache(maxsize=4)
def _quadrature(nodes):
    return np.polynomial.legendre.leggauss(nodes)

def validate(y, n):
    y, n = np.asarray(y, float), np.asarray(n, float)
    if y.ndim != 2 or y.shape != n.shape:
        raise ValueError('Count and actual donor-total matrices must have equal feature-by-sample shapes')
    if not np.array_equal(np.isfinite(y), np.isfinite(n)):
        raise ValueError('Missing counts and totals must be paired')
    finite=np.isfinite(y)&np.isfinite(n)
    seen = finite & (n > 0)
    if np.any(finite & ((y < 0) | (y > n) | (y != np.round(y)) | (n != np.round(n)))) or np.any(np.isfinite(n) & (n < 0)):
        raise ValueError('Require integer counts 0 <= y <= actual donor total n')
    return y, n, seen

def fit_margins(y, n, context=None):
    y,n,seen = validate(y,n)
    x = np.ones((y.shape[1],1)) if context is None else np.column_stack([np.ones(y.shape[1]),np.asarray(context,float)])
    if not np.isfinite(x).all(): raise ValueError('Context must be observed, with training-frozen scaling')
    coefficients=[]; dispersion=[]
    for j in range(len(y)):
        use=seen[j]; a,b=y[j,use],n[j,use]; xx=x[use]
        if use.sum()<15: raise ValueError('Each feature needs at least 15 training observations')
        initial=np.zeros(x.shape[1]);initial[0]=logit(np.clip((a.sum()+.5)/(b.sum()+1),1e-5,1-1e-5))
        def objective(c):
            eta=xx@c; mu=expit(eta)
            return (np.sum(b*np.logaddexp(0,eta)-a*eta)+.5*np.sum(c[1:]**2),xx.T@(b*mu-a)+np.r_[0,c[1:]])
        fit=minimize(objective,initial,jac=True,method='L-BFGS-B')
        if not fit.success and np.linalg.norm(fit.jac)>1e-3: raise RuntimeError(str(fit.message))
        mu=np.clip(expit(xx@fit.x),1e-5,1-1e-5); p=a/b
        numerator=np.mean((p-mu)**2-mu*(1-mu)/b)
        denominator=np.mean(mu*(1-mu)*(1-1/b))
        d=float(np.clip(numerator/max(denominator,1e-10),0,.95))
        coefficients.append(fit.x.tolist());dispersion.append(d)
    return dict(coefficients=coefficients,dispersion=dispersion,context_columns=x.shape[1]-1)

def intervals(margins,y,n,context=None):
    y,n,seen=validate(y,n)
    x=np.ones((y.shape[1],1)) if context is None else np.column_stack([np.ones(y.shape[1]),np.asarray(context,float)])
    key=hashlib.sha256(json.dumps(margins,sort_keys=True).encode()+str(y.shape).encode()+y.tobytes()+n.tobytes()+x.tobytes()).hexdigest()
    if key in _INTERVAL_CACHE:
        _INTERVAL_CACHE.move_to_end(key)
        return _INTERVAL_CACHE[key]
    c=np.asarray(margins['coefficients']);mu=np.clip(expit(c@x.T),1e-5,1-1e-5)
    lower=np.zeros_like(y);upper=np.ones_like(y)
    for j,d in enumerate(margins['dispersion']):
        use=seen[j]; yy=y[j,use];nn=n[j,use];mm=mu[j,use]
        if d<1e-8:
            lower[j,use]=binom.cdf(yy-1,nn,mm);upper[j,use]=binom.cdf(yy,nn,mm)
        else:
            concentration=1/d-1
            lower[j,use]=betabinom.cdf(yy-1,nn,mm*concentration,(1-mm)*concentration)
            upper[j,use]=betabinom.cdf(yy,nn,mm*concentration,(1-mm)*concentration)
    result=(lower,upper,mu,seen)
    for value in result:value.setflags(write=False)
    _INTERVAL_CACHE[key]=result
    if len(_INTERVAL_CACHE)>4:_INTERVAL_CACHE.popitem(last=False)
    return result

def rectangle_loss(rho,lo1,hi1,lo2,hi2,*,nodes=24,floor=1e-9):
    """Integrate conditional Gaussian CDF over first marginal's U interval."""
    q,w=_quadrature(nodes)
    u=lo1[:,None]+(hi1-lo1)[:,None]*(q+1)/2
    z=ndtri(np.clip(u,1e-14,1-1e-14)); sd=np.sqrt(1-rho*rho)
    a=ndtri(hi2)[:,None];b=ndtri(lo2)[:,None]
    probability=(hi1-lo1)*((ndtr((a-rho*z)/sd)-ndtr((b-rho*z)/sd))@w)/2
    return -float(np.log(np.maximum(probability,floor)).sum())

def pair_correlation(lower,upper,seen,*,exact=True):
    z=ndtri(np.clip((lower+upper)/2,1e-6,1-1e-6));z=np.where(seen,z,0)
    R=np.eye(len(z));audit=[]
    for i in range(len(z)-1):
        for j in range(i+1,len(z)):
            use=seen[i]&seen[j]
            if use.sum()<15:rho=0.
            elif exact:
                fit=minimize_scalar(rectangle_loss,bounds=(-.95,.95),method='bounded',args=(lower[i,use],upper[i,use],lower[j,use],upper[j,use]),options={'xatol':2e-3})
                if not fit.success:raise RuntimeError('Pair likelihood failed')
                rho=float(fit.x)
            else:
                a,b=z[i,use],z[j,use]
                rho=float(np.corrcoef(a,b)[0,1]) if a.std()>1e-10 and b.std()>1e-10 else 0.
            R[i,j]=R[j,i]=rho
            audit.append([i,j,int(use.sum()),rho])
    return R,z,audit

def positive_correlation(R):
    e,v=np.linalg.eigh((R+R.T)/2);R=(v*np.maximum(e,1e-4))@v.T
    s=np.sqrt(np.diag(R));return R/s[:,None]/s[None,:]

def estimate_precision(R,z,*,robust=False,alpha=.15):
    cov=positive_correlation(R)
    if robust:
        # Student-t reweighting is an ablation, not MVTLASSO's ICA/multi-study model.
        for _ in range(12):
            inv=np.linalg.pinv(cov);d=np.einsum('is,ij,js->s',z,inv,z)
            w=(3+len(z))/(3+d);w=np.clip(w,.05,5)
            c=(z*w)@z.T/w.sum();sd=np.sqrt(np.maximum(np.diag(c),1e-8))
            cov=positive_correlation(.5*cov+.5*c/sd[:,None]/sd[None,:])
    # A fixed numerical ridge, shared by all comparison graphs. Do not treat
    # a returned nonconverged precision matrix as a successful scientific fit.
    cov=.95*cov+.05*np.eye(len(cov))
    with warnings.catch_warnings():
        warnings.simplefilter('error',ConvergenceWarning)
        covariance,precision=graphical_lasso(cov,alpha=alpha,max_iter=2000,tol=1e-4,enet_tol=1e-6)
    return covariance,precision

def fit(y,n,feature_ids,donor_groups,*,context=None,exact=True,robust=False,nested=False,alpha=.15):
    ids=list(map(str,feature_ids));groups=list(map(str,donor_groups));y,n,seen=validate(y,n)
    if len(ids)!=len(y) or len(set(ids))!=len(ids) or len(groups)!=len(y) or len(set(groups))!=len(groups):
        raise ValueError('Unique ordered focal junction and donor IDs are required')
    margins=fit_margins(y,n,context);lo,hi,mu,seen=intervals(margins,y,n,context)
    R,z,pairs=pair_correlation(lo,hi,seen,exact=exact)
    cov,precision=estimate_precision(R,z,robust=robust,alpha=alpha)
    core=None
    if nested:
        # A frozen broad factor; residual edges are reported separately, never
        # reinterpret PSI/core division as a probability. No phenotype is used.
        e,v=np.linalg.eigh(cov);loading=v[:,-1]*np.sqrt(max(e[-1]-1,0))
        score=loading@z/max(loading@loading,1e-8)
        residual=z-loading[:,None]*score
        rr=np.corrcoef(residual);rr=np.nan_to_num(rr);np.fill_diagonal(rr,1)
        res_cov,res_precision=estimate_precision(rr,residual,robust=robust,alpha=alpha)
        core=dict(loading=loading.tolist(),residual_precision=res_precision.tolist(),scope='experimental Gaussian residual hierarchy, not author GenFocus')
        combined=positive_correlation(res_cov+np.outer(loading,loading))
        cov=combined;precision=np.linalg.inv(combined)
    model=dict(schema=SCHEMA,feature_ids=ids,donor_groups=groups,margins=margins,
        covariance=cov.tolist(),precision=precision.tolist(),core=core,
        config=dict(exact_rectangle_likelihood=exact,robust=robust,nested=nested,alpha=alpha),
        training_samples=y.shape[1],limits=['Beta-binomial margins are approximations','Training context must be independently observed in projection','Conditional associations are not causal regulator edges','No network FDR claim'])
    model['sha256']=hashlib.sha256(json.dumps(model,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return model,dict(pair_audit=pairs,raw_correlation=R.tolist())

def edge_scores(model,*,residual=False):
    p=np.asarray(model['core']['residual_precision'] if residual else model['precision'])
    d=np.sqrt(np.diag(p));e=-p/d[:,None]/d[None,:];np.fill_diagonal(e,0);return e

def predict(model,y,n,*,context=None):
    """Frozen leave-one-site-out network prediction; own reads are excluded.

    All comparator graphs can use this same head; marginal/graph effects must
    be ablated separately. Unobserved neighbours revert to training margins.
    """
    lo,hi,mu,seen=intervals(model['margins'],y,n,context)
    z=ndtri(np.clip((lo+hi)/2,1e-6,1-1e-6));z=np.where(seen,z,0)
    omega=np.asarray(model['precision']);coeff=-omega/np.diag(omega)[:,None];np.fill_diagonal(coeff,0)
    latent=coeff@z;result=mu.copy()
    covariance=np.asarray(model['covariance'])
    for sample in range(y.shape[1]):
        if seen[:,sample].all():continue
        # Integrate absent neighbours out of the Gaussian graph, rather than
        # treating a missing coordinate as an observed zero residual.
        available=seen[:,sample]
        if not available.any():latent[:,sample]=0;continue
        inv=np.linalg.inv(covariance[np.ix_(available,available)])
        values=z[available,sample];weighted=inv@values
        latent[available,sample]=values-weighted/np.diag(inv)
        latent[~available,sample]=covariance[np.ix_(~available,available)]@weighted
    for j,d in enumerate(model['margins']['dispersion']):
        if d>1e-8:
            concentration=1/d-1
            result[j]=beta.ppf(np.clip(ndtr(latent[j]),1e-5,1-1e-5),mu[j]*concentration,(1-mu[j])*concentration)
    return np.clip(result,1e-5,1-1e-5),dict(observed_neighbours=seen.sum(0).tolist(),missing_route='training_margins_zero_latent')
