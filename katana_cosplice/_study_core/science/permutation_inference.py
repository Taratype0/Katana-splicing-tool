"""Finite Monte Carlo label test for independent, unadjusted one-way designs.

Null: scores are exchangeable across condition labels (identical distributions).
This is not a heteroscedastic equal-means-only null. Discovery must be label-free;
never use this route for repeated units, batch-adjusted or paired designs.
"""
from functools import lru_cache
import numpy as np

PERMUTATIONS=9999
SEED=27041


@lru_cache(maxsize=8)
def permutation_indices(n, permutations=PERMUTATIONS, seed=SEED):
    rng=np.random.default_rng(seed)
    result=np.stack([rng.permutation(n) for _ in range(permutations)])
    result.flags.writeable=False
    return result


def permutation_anova(values, labels, permutations=PERMUTATIONS, seed=SEED):
    y=np.asarray(values,float);labels=np.asarray(labels)
    if y.ndim!=1 or len(y)!=len(labels) or not np.isfinite(y).all():raise ValueError('Finite individual scores required')
    levels,inverse=np.unique(labels,return_inverse=True);k=len(levels);n=len(y)
    if k<2 or n<=k:raise ValueError('At least two nonempty groups and residual degrees of freedom required')
    y=y-y.mean();total=float(y@y)
    if total<1e-12:raise ValueError('Constant score; association is not estimable')
    design=np.eye(k)[inverse];sizes=design.sum(axis=0)
    between=float(np.sum((y@design)**2/sizes))
    # Compare between-group SS: F is monotone in it with fixed n/groups/SST.
    shuffled=y[permutation_indices(n,permutations,seed)]
    null=np.sum((shuffled@design)**2/sizes,axis=1)
    exceed=int(np.sum(null>=between-1e-12*max(total,1)))
    statistic=(between/(k-1))/max((total-between)/(n-k),np.finfo(float).tiny)
    return dict(test='Monte Carlo condition-label permutation ANOVA (exchangeable independent units)',
                statistic=float(statistic),p_value=float((exceed+1)/(permutations+1)),
                df=float(n-k),df_num=float(k-1),permutations=permutations,permutation_exceedances=exceed,
                permutation_seed=seed,p_resolution=1/(permutations+1),
                null_hypothesis='Exchangeable condition-label distributions, not equal means under arbitrary unequal variances')
