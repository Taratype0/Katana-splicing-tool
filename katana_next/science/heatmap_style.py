"""Shared presentation defaults; scaling never changes inference or clustering."""
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


def row_scaled(values):
    values=np.asarray(values,dtype=float)
    finite=np.isfinite(values)
    count=finite.sum(axis=1,keepdims=True)
    mean=np.divide(np.where(finite,values,0).sum(axis=1,keepdims=True),count,out=np.zeros_like(count,dtype=float),where=count>0)
    centered=np.where(finite,values-mean,0)
    variance=np.divide((centered**2).sum(axis=1,keepdims=True),count-1,out=np.zeros_like(mean),where=count>1)
    sd=np.sqrt(variance)
    return np.where(finite,centered/np.where(sd>0,sd,1),np.nan)


def pastel_diverging():
    cmap=LinearSegmentedColormap.from_list('katana_pastel',['#8bb9d7','#ffffff','#e8a0a5'])
    cmap.set_bad('#d5d9dd')
    return cmap
