"""Equivalent cached implementation; no scientific parameters are changed."""
from __future__ import annotations

import re
import numpy as np
from scipy.stats import mannwhitneyu
from katana_next.science.rbp_short_extension import PositionIndex

_INDEX = {}
_OPPORTUNITY = {}


def consensus_scores(sequences, pattern, junctions=None):
    identity=(tuple(sequences),None if junctions is None else tuple(tuple(x) for x in junctions))
    if identity not in _INDEX:
        _INDEX[identity]=PositionIndex(sequences,25,junctions)
    index=_INDEX[identity]
    width=len(re.findall(r"\[[ACGT]+\]|[ACGT]",pattern))
    key=(id(index),width)
    if key not in _OPPORTUNITY:
        _,owners,local=index.positions(width)
        contribution=np.full(len(owners),width,dtype=float)
        if len(owners)>1:
            same=owners[1:]==owners[:-1]
            contribution[1:]=np.where(same,np.minimum(local[1:]-local[:-1],width),width)
        denominator=np.bincount(owners,weights=contribution,minlength=len(sequences))
        _OPPORTUNITY[key]=denominator
    denominator=_OPPORTUNITY[key]
    regex=re.compile("(?=("+pattern+"))")
    values=np.full(len(sequences),np.nan)
    hits,masks=[],[]
    for i,sequence in enumerate(sequences):
        covered=np.zeros(len(sequence),dtype=bool)
        present=[]
        if denominator[i]>0:
            for match in regex.finditer(sequence):
                a,b=match.start(1),match.end(1)
                if junctions is None or any(a<j<b for j in junctions[i]):
                    present.append((a,b))
                    covered[a:b]=True
            values[i]=covered.sum()/denominator[i]
        hits.append(present)
        masks.append(covered)
    return values,hits,masks


def rank_test_matrix(foreground,background,minimum=20):
    a,b=np.asarray(foreground,float),np.asarray(background,float)
    if a.ndim==1:
        a,b=a[:,None],b[:,None]
    n1,n2=np.isfinite(a).sum(axis=0),np.isfinite(b).sum(axis=0)
    eligible=(n1>=minimum)&(n2>=minimum)
    minimum_value=np.minimum(np.min(np.where(np.isfinite(a),a,np.inf),axis=0,initial=np.inf),
                             np.min(np.where(np.isfinite(b),b,np.inf),axis=0,initial=np.inf))
    maximum_value=np.maximum(np.max(np.where(np.isfinite(a),a,-np.inf),axis=0,initial=-np.inf),
                             np.max(np.where(np.isfinite(b),b,-np.inf),axis=0,initial=-np.inf))
    constant=minimum_value==maximum_value
    p=np.full(a.shape[1],np.nan)
    p[eligible&constant]=1.
    active=eligible&~constant
    if active.any():
        p[active]=mannwhitneyu(a[:,active],b[:,active],axis=0,alternative="greater",method="asymptotic",
                              use_continuity=True,nan_policy="omit").pvalue
    mean_a=np.divide(np.nansum(a,axis=0),n1,out=np.full(len(n1),np.nan),where=n1>0)
    mean_b=np.divide(np.nansum(b,axis=0),n2,out=np.full(len(n2),np.nan),where=n2>0)
    return p,mean_a,mean_b,n1,n2
