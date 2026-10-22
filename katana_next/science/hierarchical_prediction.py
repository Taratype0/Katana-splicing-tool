"""Training-only BB moment estimate and quadrature over the event mean.

This empirical-Bayes diagnostic assumes independent replicate probabilities given
an event mean and a shared rho. It is not a universal calibration guarantee.
"""
from functools import lru_cache
from itertools import combinations

import numpy as np
from scipy.special import betaln, logsumexp, roots_jacobi


def validated_counts(inclusion, exclusion):
    inclusion, exclusion = np.asarray(inclusion, float), np.asarray(exclusion, float)
    if inclusion.ndim != 2 or inclusion.shape != exclusion.shape:
        raise ValueError("Counts require matching event-by-replicate matrices")
    for counts in (inclusion, exclusion):
        if np.any(np.isinf(counts) | (counts < 0) | (np.isfinite(counts) & (counts != np.floor(counts)))):
            raise ValueError("Counts must be nonnegative integers or missing")
    paired = np.isfinite(inclusion+exclusion) & ((inclusion+exclusion) > 0)
    return np.where(paired, inclusion, np.nan), np.where(paired, exclusion, np.nan)


def paired_moment_rho(inclusion, exclusion):
    inclusion, exclusion = validated_counts(inclusion, exclusion)
    n = inclusion+exclusion
    q = np.divide(inclusion,n,out=np.full_like(n,np.nan,dtype=float),where=n>0)
    numerator, denominator, pairs = 0., 0., 0
    for a,b in combinations(range(n.shape[1]),2):
        valid = np.isfinite(q[:,a]+q[:,b]) & (n[:,a]>1) & (n[:,b]>1)
        x,y = q[valid,a],q[valid,b]
        sampling = 1/n[valid,a]+1/n[valid,b]
        # E[q_a(1-q_b)] = mu(1-mu), without fitting an event mean to two repeats.
        u = (x*(1-y)+y*(1-x))/2
        numerator += np.sum((x-y)**2-u*sampling)
        denominator += np.sum(u*(2-sampling))
        pairs += int(valid.sum())
    if pairs < 10 or denominator <= 0:
        return None
    raw = float(numerator/denominator)
    return {"rho":float(np.clip(raw,0,1-1e-6)),"unclipped_rho":raw,"training_pairs":pairs,
            "boundary":bool(raw<0 or raw>=1)}


@lru_cache(maxsize=4)
def mean_quadrature(points=256):
    if not isinstance(points, int) or points < 16:
        raise ValueError("At least 16 quadrature points are required")
    nodes,weights=roots_jacobi(points,-.5,-.5)
    return (nodes+1)/2,np.log(weights/weights.sum())


def hierarchical_draws(inclusion, exclusion, held_n, *, rho, rng, draws, points=256):
    if not 0 <= rho < 1:
        raise ValueError("rho must be in [0,1)")
    inclusion, exclusion = validated_counts(inclusion, exclusion)
    held_n = np.asarray(held_n, float)
    if held_n.shape != (len(inclusion),) or np.any(~np.isfinite(held_n) | (held_n <= 0) | (held_n != np.floor(held_n))):
        raise ValueError("Held depths require one positive integer per event")
    if not isinstance(draws, int) or draws < 1:
        raise ValueError("Predictive draws must be a positive integer")
    n=inclusion+exclusion
    if np.any(np.isfinite(n).sum(axis=1) < 2):
        raise ValueError("At least two observed training replicates per event are required")
    if not len(n):
        return np.empty((0, draws), dtype=np.int64)
    if rho <= 1e-12:
        mean=rng.beta(np.nansum(inclusion,axis=1)[:,None]+.5,
                      np.nansum(exclusion,axis=1)[:,None]+.5,(len(n),draws))
    else:
        mu,log_weights=mean_quadrature(points)
        k=1/rho-1
        a,b=mu*k,(1-mu)*k
        logposterior=np.broadcast_to(log_weights,(len(n),len(mu))).copy()
        for j in range(n.shape[1]):
            valid=np.isfinite(n[:,j]) & (n[:,j]>0)
            logposterior[valid] += betaln(inclusion[valid,j,None]+a,exclusion[valid,j,None]+b)-betaln(a,b)
        weights=np.exp(logposterior-logsumexp(logposterior,axis=1)[:,None])
        cumulative=weights.cumsum(axis=1); cumulative[:,-1]=1.
        mean=np.stack([mu[np.searchsorted(cdf,rng.random(draws))] for cdf in cumulative])
        mean=rng.beta(mean*k,(1-mean)*k)
    return rng.binomial(np.asarray(held_n,dtype=np.int64)[:,None],mean)
