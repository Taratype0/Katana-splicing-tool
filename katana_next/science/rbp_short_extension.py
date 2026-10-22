"""Full RNA-position scoring, including short exons and genuine splice joins."""
from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np
from scipy.stats import mannwhitneyu


@dataclass
class PositionIndex:
    sequences: list[str]
    max_width: int
    junctions: list[tuple[int, ...]] | None = None

    def __post_init__(self):
        lookup = np.full(256, 4, dtype=np.uint8)
        lookup[np.frombuffer(b"ACGT", dtype=np.uint8)] = np.arange(4)
        self.lengths = np.array([len(s) for s in self.sequences])
        self.offsets = np.r_[0, np.cumsum(self.lengths + self.max_width)[:-1]]
        self.encoded = np.concatenate([np.r_[lookup[np.frombuffer(s.encode("ascii"), dtype=np.uint8)],
                                             np.full(self.max_width, 4, dtype=np.uint8)] for s in self.sequences])
        self.cache = {}

    def positions(self, width):
        if width not in self.cache:
            starts, owners, local = [], [], []
            for i, length in enumerate(self.lengths):
                positions = np.arange(max(0, length - width + 1))
                if self.junctions is not None:
                    crossing = np.zeros(len(positions), dtype=bool)
                    for junction in self.junctions[i]:
                        crossing |= (positions < junction) & (positions + width > junction)
                    positions = positions[crossing]
                starts.append(positions + self.offsets[i])
                local.append(positions)
                owners.append(np.full(len(positions), i, dtype=np.int32))
            starts, owners, local = np.concatenate(starts), np.concatenate(owners), np.concatenate(local)
            valid = np.ones(len(starts), dtype=bool)
            for j in range(width):
                valid &= self.encoded[starts+j] < 4
            self.cache[width] = (starts[valid], owners[valid], local[valid])
        return self.cache[width]

    def pwm(self, matrix, background):
        matrix = np.asarray(matrix, dtype=float)
        matrix = matrix / matrix.sum(axis=1, keepdims=True)
        matrix = (20 * matrix + .1 * np.asarray(background)) / 20.1
        starts, owners, local = self.positions(len(matrix))
        logodds = np.log2(matrix / background)
        relative = np.log2(matrix / matrix.max(axis=1, keepdims=True))
        score, rel = np.zeros(len(starts)), np.zeros(len(starts))
        for j in range(len(matrix)):
            score += logodds[j, self.encoded[starts+j]]
            rel += relative[j, self.encoded[starts+j]]
        values = np.exp2(score)
        total = np.bincount(owners, weights=values, minlength=len(self.sequences))
        number = np.bincount(owners, minlength=len(self.sequences))
        mean = np.divide(total, number, out=np.full(len(number), np.nan), where=number > 0)
        return np.log2(1+mean), total, (owners, local, values, np.exp2(rel))


def window_validity(sequences, window, span):
    valid = np.zeros((len(sequences), span), dtype=bool)
    for i, sequence in enumerate(sequences):
        invalid = np.r_[0, np.cumsum([base not in "ACGT" for base in sequence])]
        number = min(span, max(0, len(sequence)-window+1))
        valid[i, :number] = invalid[window:window+number] == invalid[:number]
    return valid


def pwm_windows(site_data, n_sequences, motif_width, window, span, validity):
    if motif_width > window:
        return np.full((n_sequences, span), np.nan)
    owners, positions, values, _ = site_data
    number = window - motif_width + 1
    dense = np.zeros((n_sequences, span+number-1))
    chosen = positions < dense.shape[1]
    dense[owners[chosen], positions[chosen]] = values[chosen]
    cumulative = np.c_[np.zeros(n_sequences), np.cumsum(dense, axis=1)]
    score = np.log2(1+(cumulative[:, number:] - cumulative[:, :-number]) / number)
    score[~validity] = np.nan
    return score


def consensus_scores(sequences, pattern, junctions=None):
    regex = re.compile("(?=(" + pattern + "))")
    scores, hits, masks = [], [], []
    for i, sequence in enumerate(sequences):
        covered = np.zeros(len(sequence), dtype=bool)
        present = []
        width = None
        for match in regex.finditer(sequence):
            a, b = match.start(1), match.end(1)
            width = b-a
            if junctions is not None and not any(a < j < b for j in junctions[i]):
                continue
            covered[a:b] = True
            present.append((a,b))
        # Determine physical eligibility even for a sequence with zero hits.
        tokens = re.findall(r"\[[ACGT]+\]|[ACGT]", pattern)
        width = len(tokens)
        possible = [a for a in range(max(0,len(sequence)-width+1))
                    if set(sequence[a:a+width]) <= set("ACGT")
                    and (junctions is None or any(a < j < a+width for j in junctions[i]))]
        opportunity = np.zeros(len(sequence), dtype=bool)
        for start in possible:
            opportunity[start:start+width] = True
        scores.append(covered.sum()/opportunity.sum() if possible else np.nan)
        masks.append(covered)
        hits.append(present)
    return np.asarray(scores), hits, masks


def consensus_windows(masks, window, span, validity, motif_width):
    scores = np.full((len(masks), span), np.nan)
    if motif_width > window:
        return scores
    for i, mask in enumerate(masks):
        count = min(span, max(0,len(mask)-window+1))
        prefix = np.r_[0,np.cumsum(mask)]
        scores[i,:count] = (prefix[window:window+count] - prefix[:count])/window
    scores[~validity] = np.nan
    return scores


def rank_test_matrix(foreground, background, minimum=20):
    a = np.asarray(foreground, float)
    b = np.asarray(background, float)
    if a.ndim == 1:
        a = a[:,None]
        b = b[:,None]
    n1, n2 = np.isfinite(a).sum(axis=0), np.isfinite(b).sum(axis=0)
    keep = (n1 >= minimum) & (n2 >= minimum)
    p = np.full(a.shape[1], np.nan)
    if keep.any():
        p[keep] = mannwhitneyu(a[:,keep], b[:,keep], axis=0, alternative="greater",
                              method="asymptotic", use_continuity=True, nan_policy="omit").pvalue
    mean_a = np.divide(np.nansum(a,axis=0),n1,out=np.full(len(n1),np.nan),where=n1>0)
    mean_b = np.divide(np.nansum(b,axis=0),n2,out=np.full(len(n2),np.nan),where=n2>0)
    return p, mean_a, mean_b, n1, n2
