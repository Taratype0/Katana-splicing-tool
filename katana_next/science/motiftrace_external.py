"""External assay validation helpers; no quantification or motif-model fitting."""
from collections import Counter
import hashlib
import re

import numpy as np


def interval_from_rna(block, strand, start, stop):
    if not 0 <= start < stop <= block[1]-block[0]:
        raise ValueError("RNA interval does not fit the actual source sequence")
    return (block[0]+start, block[0]+stop) if strand == "+" else (block[1]-stop, block[1]-start)


class PeakIndex:
    def __init__(self, rows):
        self.groups = {}
        for (chrom, strand), part in rows.groupby(["chr", "strand"], observed=True):
            part = part.sort_values(["start", "end"])
            starts, ends = part.start.to_numpy(), part.end.to_numpy()
            self.groups[(str(chrom), str(strand))] = (starts, ends, np.maximum.accumulate(ends))

    def overlap(self, chrom, strand, start, stop):
        values = self.groups.get((str(chrom), str(strand)))
        if values is None:
            return False
        starts, _, prefix = values
        i = int(np.searchsorted(starts, stop, side="left"))-1
        return i >= 0 and bool(prefix[i] > start)


def dinucleotide_shuffle(sequence, rng):
    """Random Euler trail preserves every dinucleotide and both terminal bases."""
    adjacent = {base: [] for base in "ACGT"}
    for first, second in zip(sequence, sequence[1:]):
        adjacent[first].append(second)
    for values in adjacent.values():
        rng.shuffle(values)
    stack, path = [sequence[0]], []
    while stack:
        if adjacent[stack[-1]]:
            stack.append(adjacent[stack[-1]].pop())
        else:
            path.append(stack.pop())
    shuffled = "".join(path[::-1])
    if len(shuffled) != len(sequence) or Counter(zip(sequence, sequence[1:])) != Counter(zip(shuffled, shuffled[1:])):
        raise ValueError("Dinucleotide-preserving null generation failed")
    return shuffled


def tpr_at_fpr(positive, negative, alpha=.05):
    """Tie-randomized ROC sensitivity at an exactly specified empirical null FPR."""
    positive, negative = np.asarray(positive), np.asarray(negative)
    if not len(positive) or not len(negative):
        return np.nan
    boundary = np.sort(negative)[::-1][min(int(alpha*len(negative)), len(negative)-1)]
    higher = int((negative > boundary).sum())
    tied = int((negative == boundary).sum())
    fraction = np.clip((alpha*len(negative)-higher)/tied, 0, 1)
    return float(((positive > boundary).sum()+fraction*(positive == boundary).sum())/len(positive))


def overlapping_hit_positions(sequence, pattern):
    return [(match.start(1), match.end(1)) for match in re.finditer("(?=("+pattern+"))", sequence)]


def stable_hash(value):
    return hashlib.sha256(str(value).encode()).hexdigest()
