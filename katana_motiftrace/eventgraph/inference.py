"""Event signed ranks with whole-gene sign randomization.

Validity is conditional on joint central symmetry within independent gene blocks.
This is a sequence association test, not a model of biological sample counts.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import rankdata, norm
from scipy.ndimage import label

THRESHOLD=float(norm.isf(.01))

def bh(p):
    p=np.asarray(p,float);out=np.full(len(p),np.nan);valid=np.flatnonzero(np.isfinite(p))
    if len(valid):
        order=valid[np.argsort(p[valid],kind='stable')];q=p[order]*len(order)/np.arange(1,len(order)+1)
        out[order]=np.minimum(1,np.minimum.accumulate(q[::-1])[::-1])
    return out

def signed_block_ranks(differences,blocks,minimum_blocks=20):
    d=np.asarray(differences,float)
    if d.ndim!=2 or len(blocks)!=len(d):raise ValueError('Event matrix/block identities disagree')
    keys,codes=np.unique(np.asarray(blocks,str),return_inverse=True)
    valid=np.isfinite(d);n_events=valid.sum(axis=0)
    values=np.where(valid&(d!=0),abs(d),np.nan)
    ranks=np.nan_to_num(rankdata(values,axis=0,method='average',nan_policy='omit'))*np.sign(np.nan_to_num(d))
    sums=np.zeros((len(keys),d.shape[1]));presence=np.zeros((len(keys),d.shape[1]),dtype=np.int32)
    np.add.at(sums,codes,ranks);np.add.at(presence,codes,valid.astype(np.int32))
    n_blocks=(presence>0).sum(axis=0)
    denom=np.sqrt(np.square(sums).sum(axis=0));denom[denom==0]=1
    normalized=sums/denom
    observed=normalized.sum(axis=0);observed[n_blocks<minimum_blocks]=np.nan
    return keys,normalized,observed,n_events,n_blocks

def runs(z,threshold=THRESHOLD):
    idx=np.flatnonzero(np.isfinite(z)&(z>threshold))
    if not len(idx):return []
    return [(int(p[0]),int(p[-1]+1),float((z[p]-threshold).sum())) for p in np.split(idx,np.flatnonzero(np.diff(idx)!=1)+1)]

def max_cluster_rows(z,threshold=THRESHOLD):
    # Label only horizontally: each row is a separate permutation.
    mask=np.isfinite(z)&(z>threshold)
    labs,n=label(mask,structure=np.array([[0,0,0],[1,1,1],[0,0,0]]))
    result=np.zeros(len(z))
    if n:
        rr,cc=np.nonzero(mask);masses=np.bincount(labs[rr,cc],weights=z[rr,cc]-threshold,minlength=n+1)
        np.maximum.at(result,rr,masses[labs[rr,cc]])
    return result

class BlockFamily:
    def __init__(self,blocks,permutations=999,seed=20261015,minimum_blocks=20):
        if permutations<99:raise ValueError('At least 99 randomizations required')
        self.keys=np.array(sorted(set(map(str,blocks))));self.lookup={x:i for i,x in enumerate(self.keys)}
        self.signs=np.random.default_rng(seed).choice(np.array([-1.,1.]),size=(permutations,len(self.keys)))
        self.null_max=np.zeros(permutations);self.B=permutations;self.minimum=minimum_blocks

    def evaluate(self,differences,blocks,meta,compute_clusters=True):
        keys,ranks,z,n,ng=signed_block_ranks(differences,blocks,self.minimum)
        indices=np.array([self.lookup[x] for x in keys]);signs=self.signs[:,indices]
        p=np.full(len(z),np.nan);good=np.isfinite(z)
        # Chunk permutations, not positions: preserve entire contiguous spatial runs.
        exceed=np.zeros(len(z),dtype=np.int32)
        for start in range(0,self.B,64):
            null=signs[start:start+64]@ranks
            null[:,~good]=np.nan
            exceed+=np.sum(null>=z-1e-12,axis=0)
            if compute_clusters and meta.get('window_nt',0)>0:
                self.null_max[start:start+len(null)]=np.maximum(self.null_max[start:start+len(null)],max_cluster_rows(null))
        p[good]=(1+exceed[good])/(self.B+1)
        sums=np.nansum(differences,axis=0);means=np.divide(sums,n,out=np.full(len(n),np.nan),where=n>0)
        rows=pd.DataFrame({'window_start0':np.arange(len(z)) if meta.get('window_nt',0)>0 else -np.ones(len(z),int),
                           'n_matched_events':n,'n_independent_gene_blocks':ng,'paired_rank_Z':z,
                           'mean_paired_score_difference':means,'raw_P':p})
        for k,v in meta.items():rows[k]=v
        clusters=[]
        if compute_clusters and meta.get('window_nt',0)>0:
            for a,b,mass in runs(z):
                clusters.append({**meta,'first_window_start0':a,'last_window_start0':b-1,'cluster_mass':mass,
                                 'minimum_events':int(n[a:b].min()),'minimum_gene_blocks':int(ng[a:b].min())})
        return rows,clusters

    def adjust_clusters(self,clusters):
        c=pd.DataFrame(clusters)
        if len(c):c['cluster_P_FWER']=[(1+np.sum(self.null_max>=x-1e-12))/(self.B+1) for x in c.cluster_mass]
        return c

def fast_pwm_positions(sequence,model):
    if model['kind']!='PWM':raise ValueError('Explicit PWM required; no silent consensus conversion')
    seq=sequence.upper().replace('U','T');w=len(model['matrix']);n=max(0,len(seq)-w+1)
    if not n:return np.empty(0)
    table=np.full(256,4,np.int16);table[np.frombuffer(b'ACGT',dtype=np.uint8)]=np.arange(4)
    encoded=table[np.frombuffer(seq.encode('ascii'),dtype=np.uint8)]
    invalid=np.r_[0,np.cumsum(encoded==4)];valid=(invalid[w:]-invalid[:-w])==0
    starts=np.flatnonzero(valid);odds=(20*np.asarray(model['matrix'])+.025)/20.1/.25
    values=np.ones(len(starts))
    for k in range(w):values*=odds[k,encoded[starts+k]]
    out=np.full(n,np.nan);out[starts]=values
    return out

def pwm_profiles(sequence,model,widths=(0,10,20,50)):
    pos=fast_pwm_positions(sequence,model);w=len(model['matrix']);valid=np.isfinite(pos)
    out={}
    for width in widths:
        if width==0:
            out[0]=np.array([np.log2(1+np.nanmean(pos)) if valid.any() else np.nan]);continue
        span=max(0,len(sequence)-width+1);values=np.full(span,np.nan);k=width-w+1
        if k>0 and span:
            # Match the original per-window reduction rather than subtracting
            # long cumulative sums, which can perturb exact ties in rank tests.
            view=np.lib.stride_tricks.sliding_window_view(pos,k)
            mask=np.isfinite(view).any(axis=1)
            values[mask]=np.log2(1+np.nanmean(view[mask],axis=1))
        out[width]=values
    return out
