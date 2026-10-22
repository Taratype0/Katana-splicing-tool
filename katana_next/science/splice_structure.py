"""Experimental splice-competition moment model.

This is a conditional-binomial adaptation of measurement/latent-expression
separation, not a claim to invent GLM-PCA or to establish causal regulation.
One focal junction per donor competition is required. Counts and actual donor
totals are mandatory; PSI alone cannot identify sampling noise.
"""
from __future__ import annotations
import hashlib,json
import numpy as np
import pandas as pd
from scipy.linalg import eigh

SCHEMA='katana-splice-competition-moments-0.1'
_NULL_CACHE={}

def _counts(y,n):
    y=np.asarray(y,dtype=float);n=np.asarray(n,dtype=float)
    if y.ndim!=2 or y.shape!=n.shape:raise ValueError('Junction counts and donor totals must have identical feature-by-sample shapes')
    finite=np.isfinite(y)&np.isfinite(n)
    if np.any(finite&((y<0)|(n<0)|(y>n)|(abs(y-np.round(y))>1e-8)|(abs(n-np.round(n))>1e-8))):
        raise ValueError('Require observed integer junction counts 0 <= y <= actual donor total n')
    return y,n,finite&(n>=2)

def moments(y,n):
    """Unbiased per-observation binomial variance subtraction, before PSD/rank regularization.

    E[y(n-y)/(n^2(n-1))|p,n] = p(1-p)/n.
    Cross-donor read sampling is assumed conditionally independent. Missing
    pairs are centered on their common observations. Nonpositive latent
    variances remain in the audit and become unestimable for correlations.
    """
    y,n,seen=_counts(y,n)
    p=np.divide(y,n,out=np.zeros_like(y),where=seen);mask=seen.astype(float)
    observed=mask.sum(1);mean=np.divide(p.sum(1),observed,out=np.zeros(len(y)),where=observed>0)
    pairs=mask@mask.T;sums=p@mask.T
    cross=p@p.T-np.divide(sums*sums.T,pairs,out=np.zeros_like(pairs),where=pairs>0)
    covariance=np.divide(cross,pairs-1,out=np.zeros_like(cross),where=pairs>1)
    sampling=np.divide(p*(1-p),n-1,out=np.zeros_like(y),where=seen)
    noise=np.divide(sampling.sum(1),observed,out=np.zeros(len(y)),where=observed>0)
    latent=covariance.copy();latent.flat[::len(y)+1]-=noise
    raw_var=np.diag(covariance).copy();true_var=np.diag(latent).copy()
    reliability=np.divide(np.maximum(true_var,0),raw_var,out=np.zeros_like(raw_var),where=raw_var>1e-12)
    return dict(mean=mean,covariance=covariance,latent_covariance=latent,mean_sampling_variance=noise,
                observed_variance=raw_var,latent_variance=true_var,reliability=reliability,
                observed_fraction=observed/y.shape[1],pair_observations=pairs)

def covariance_to_correlation(covariance):
    diag=np.diag(covariance);d=np.sqrt(np.maximum(diag,1e-12))
    result=np.clip(covariance/d[:,None]/d[None,:],-1,1)
    result[diag<=1e-12,:]=0;result[:,diag<=1e-12]=0
    np.fill_diagonal(result,1)
    return result

def _parallel_analysis(y,n,correct_noise,indices=None):
    """Shared training-only reference, not an analytic network/FDR test.

    Independently permuting each site's count/total pair preserves its
    marginal read-depth distribution and variation. A candidate eigenvalue
    must exceed all 32 null maximum eigenvalues. This generic safeguard is
    supplied identically to corrected/uncorrected and hard-module controls.
    Fixed, data-selected hard modules do not have a post-selection FDR claim.
    """
    key=hashlib.sha256(y.tobytes()+n.tobytes()).hexdigest()
    if key not in _NULL_CACHE:
        if len(_NULL_CACHE)>=2:_NULL_CACHE.pop(next(iter(_NULL_CACHE)))
        rng=np.random.default_rng(20261019);raw=[];latent=[]
        for repeat in range(32):
            order=np.argsort(rng.random(y.shape),axis=1)
            mm=moments(np.take_along_axis(y,order,axis=1),np.take_along_axis(n,order,axis=1))
            raw.append(mm['covariance']);latent.append(mm['latent_covariance'])
        _NULL_CACHE[key]=dict(raw=raw,latent=latent,thresholds={})
    cached=_NULL_CACHE[key];ix=tuple(indices) if indices is not None else None;tag=(bool(correct_noise),ix)
    if tag not in cached['thresholds']:
        values=[]
        for cov in cached['latent' if correct_noise else 'raw']:
            a=cov if indices is None else cov[np.ix_(indices,indices)]
            e=eigh(a,subset_by_index=(len(a)-1,len(a)-1),eigvals_only=True,check_finite=False)
            values.append(float(e[-1]))
        cached['thresholds'][tag]=max(0.,max(values))
    return cached['thresholds'][tag]

def fit_competition_model(y,n,feature_ids,donor_groups,*,rank=8,correct_noise=True,labels=None):
    """Freeze a global overlapping program model or a comparator's hard-module model.

    All methods use the same count-aware projection below. labels constrains
    each feature to its original WGCNA/DJExpress module (grey stays reference).
    The sole new-estimator contrast is moment-corrected covariance and global
    overlapping loadings, measured against both uncorrected/global and
    corrected/hard-module controls.
    """
    y,n,seen=_counts(y,n);ids=list(map(str,feature_ids));groups=list(map(str,donor_groups))
    if len(ids)!=len(y) or len(set(ids))!=len(ids):raise ValueError('Unique feature IDs required')
    if len(groups)!=len(ids) or len(set(groups))!=len(groups):raise ValueError('Exactly one focal junction per donor competition is required')
    if y.shape[1]<15:raise ValueError('At least 15 independent training samples required')
    if len(y)<2 or (seen.sum(1)<2).any():raise ValueError('Each of at least two focal sites requires measured support in at least two training samples')
    mm=moments(y,n);cov=mm['latent_covariance'] if correct_noise else mm['covariance']
    thresholds=[]
    if labels is None:
        k=min(int(rank),len(y)-1,y.shape[1]-2)
        if k<1:loadings=np.zeros((len(y),0));eigen=np.array([])
        else:
            eigen,vectors=eigh(cov,subset_by_index=(len(y)-k,len(y)-1),check_finite=False)
            threshold=_parallel_analysis(y,n,correct_noise);thresholds.append(threshold)
            keep=eigen>max(1e-10,threshold);eigen=eigen[keep][::-1];vectors=vectors[:,keep][:,::-1]
            loadings=vectors*np.sqrt(eigen)[None,:]
    else:
        labels=np.asarray(labels).astype(str);parts=[];ev=[]
        for label in sorted(set(labels)-{'grey','gray','0','0.0','nan','None'}):
            ix=np.flatnonzero(labels==label)
            if len(ix)<2:continue
            e,v=eigh(cov[np.ix_(ix,ix)],subset_by_index=(len(ix)-1,len(ix)-1),check_finite=False)
            threshold=_parallel_analysis(y,n,correct_noise,ix);thresholds.append(threshold)
            if e[-1]<=max(1e-10,threshold):continue
            a=np.zeros(len(y));a[ix]=v[:,-1]*np.sqrt(e[-1]);parts.append(a);ev.append(e[-1])
        loadings=np.column_stack(parts) if parts else np.zeros((len(y),0));eigen=np.array(ev)
    # Stable orientation; no test values or labels determine signs.
    for k in range(loadings.shape[1]):
        anchor=np.argmax(abs(loadings[:,k]))
        if loadings[anchor,k]<0:loadings[:,k]*=-1
    residual=np.maximum(np.maximum(mm['latent_variance'],0)-np.square(loadings).sum(1),1e-6)
    model=dict(schema=SCHEMA,feature_ids=ids,donor_groups=groups,training_samples=y.shape[1],
        correction=bool(correct_noise),rank=int(loadings.shape[1]),requested_rank=int(rank),
        partition='overlapping_global' if labels is None else 'frozen_comparator_hard_modules',
        means=mm['mean'].tolist(),loadings=loadings.tolist(),residual_variance=residual.tolist(),
        training_mean_sampling_variance=mm['mean_sampling_variance'].tolist(),
        minimum_projection_depth=2,feature_definition='focal annotated junction / all annotated same-donor competitors',
        shared_parallel_analysis=dict(training_permutations=32,seed=20261019,criterion='eigenvalue exceeds maximum of 32 independently permuted-site null maxima',thresholds=thresholds,network_FDR='not_established'),
        limits=['Conditional independence of read sampling across donor sites','Binomial sampling variance excludes mapping uncertainty and PCR dependence','Latent programs are not identified causal splicing factors'])
    model['sha256']=hashlib.sha256(json.dumps(model,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    audit=pd.DataFrame({k:v for k,v in mm.items() if np.ndim(v)==1},index=ids)
    return model,audit

def project_competition_model(model,y,n,*,exclude_features=None):
    """Common Gaussian conditional score using observed binomial precision.

    The empirical moment model is fitted on proportions; the observation
    variance is estimated using Jeffreys probabilities only for stable new
    sample precision. Missing data are integrated out. Intercepts are explicit.
    """
    y,n,seen=_counts(y,n);mu=np.asarray(model['means']);L=np.asarray(model['loadings'])
    if len(y)!=len(mu):raise ValueError('Feature order/number must match the frozen model')
    if exclude_features is not None:seen[np.asarray(exclude_features,dtype=int),:]=False
    ratio=np.divide(y,n,out=np.zeros_like(y),where=seen)
    posterior=(y+.5)/(n+1)
    noise=np.divide(posterior*(1-posterior),n+1,out=np.ones_like(y),where=np.isfinite(n)&(n>=0))
    residual=np.asarray(model['residual_variance']);pred=np.repeat(mu[:,None],y.shape[1],axis=1)
    scores=np.zeros((y.shape[1],L.shape[1]));coverage=[]
    for j in range(y.shape[1]):
        valid=seen[:,j];a=L[valid];d=noise[valid,j]+residual[valid]
        if L.shape[1] and valid.any():
            precision=np.eye(L.shape[1])+(a.T/d)@a
            z=np.linalg.solve(precision,(a.T/d)@(ratio[valid,j]-mu[valid]));scores[j]=z;pred[:,j]+=L@z
        coverage.append(dict(observed_features=int(valid.sum()),total_features=len(y),uses_splice_information=bool(L.shape[1] and valid.any()),
                             route='posterior_program' if L.shape[1] and valid.any() else 'training_reference'))
    return np.clip(pred,1e-5,1-1e-5),scores,pd.DataFrame(coverage)

def co_usage_edges(model):
    L=np.asarray(model['loadings']);cov=L@L.T
    return covariance_to_correlation(cov+np.diag(model['residual_variance']))
