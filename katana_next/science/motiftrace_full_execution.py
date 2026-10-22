"""Execution-only vectorization of the frozen MotifTrace paired analysis."""
from __future__ import annotations

import re

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from katana_next.science.motiftrace_custom import fixed_pattern
from katana_next.science.motiftrace_events import psi_values


def fast_match(events, selected_ids, rna):
    features = {}
    for token, group in rna.groupby("event_id", sort=False):
        word = "".join(group.sequence)
        n = len(word)
        valid = sum(word.count(b) for b in "ACGT")
        features[token] = (n, (word.count("G")+word.count("C"))/valid if valid else np.nan, 1-valid/n if n else 1)
    significant = events.stable_significant.fillna(False) | events.fdr.le(.05) | events.event_id.isin(selected_ids)
    excluded = set(events.loc[significant, "gene_id"])
    pool = events[events.measurement_eligible.fillna(False) & events.fdr.gt(.1) & events.dpsi.abs().lt(.05) & ~events.gene_id.isin(excluded)]
    def info(row):
        n, gc, amb = features.get(row["event_id"], (0, np.nan, 1))
        p = [p["PSI"] for p in psi_values(row) if p["group"] == "denominator" and np.isfinite(p["PSI"])]
        return n, gc, amb, np.mean(p) if p else np.nan, min(float(row["numerator_mean_reads"]), float(row["denominator_mean_reads"]))
    controls = pool.to_dict("records")
    values = np.asarray([info(r) for r in controls], float).reshape(-1, 5)
    kinds = np.array([r["event_type"] for r in controls])
    genes = np.array([str(r["gene_id"]) for r in controls])
    ids = np.array([r["event_id"] for r in controls])
    available = np.array([pd.notna(r["gene_id"]) and bool(str(r["gene_id"]).strip()) for r in controls])
    used, rows = set(), []
    for row in events[events.event_id.isin(selected_ids)].sort_values("event_id").to_dict("records"):
        gene = row["gene_id"]
        result = {"foreground_event": row["event_id"], "foreground_gene": gene, "event_type": row["event_type"],
            "direction": "inclusion_up" if row["dpsi"] > 0 else "inclusion_down", "background_event": "", "background_gene": "", "status": "unmatched"}
        if pd.isna(gene) or not str(gene).strip():
            result["status"] = "missing_gene_id"
        elif gene in used:
            result["status"] = "same_gene_extra_event_trace_only"
        else:
            used.add(gene)
            n, gc, amb, psi, depth = info(row)
            if n > 0 and amb <= .05 and np.isfinite([gc, psi, depth]).all() and depth > 0:
                with np.errstate(divide="ignore", invalid="ignore"):
                    distance = np.array([np.abs(values[:, 1]-gc), np.abs(values[:, 0]-n)/n,
                        np.abs(values[:, 3]-psi), np.abs(np.log2(depth/values[:, 4]))])
                eligible = available & (kinds == row["event_type"]) & (values[:, 2] <= .05) & (values[:, 4] > 0)
                eligible &= np.isfinite(distance).all(axis=0) & (distance <= np.array([.05, .2, .15, 1])[:, None]).all(axis=0)
                chosen = np.flatnonzero(eligible)
                if len(chosen):
                    order = np.lexsort((ids[chosen], distance[3, chosen], distance[2, chosen], distance[1, chosen], distance[0, chosen]))
                    index = chosen[order[0]]
                    result.update(background_event=ids[index], background_gene=controls[index]["gene_id"], status="matched")
                    available[genes == genes[index]] = False
        rows.append(result)
    return pd.DataFrame(rows)


def placements(word, pattern):
    regex, width = fixed_pattern(pattern)
    n = max(0, len(word)-width+1)
    bases = np.frombuffer(word.encode("ascii"), np.uint8)
    bad = ~np.isin(bases, [65, 67, 71, 84])
    cumulative = np.r_[0, np.cumsum(bad)]
    valid = np.flatnonzero(cumulative[width:width+n] == cumulative[:n])
    hits = np.fromiter((m.start(1) for m in re.finditer("(?=("+regex+"))", word)), dtype=np.int64)
    return valid, hits, width


def _union_lengths(positions, starts, ends, width):
    if not len(positions):
        return np.zeros(len(starts))
    sums = np.r_[0, np.cumsum(np.minimum(np.diff(positions), width))]
    left = np.searchsorted(positions, starts)
    stop = np.searchsorted(positions, ends-width, side="right")
    values = np.zeros(len(starts), float)
    chosen = stop > left
    values[chosen] = width+sums[stop[chosen]-1]-sums[left[chosen]]
    return values


def profiles(word, pattern, widths=(0, 10, 20, 50)):
    valid, hits, motif_width = placements(word, pattern)
    output = {}
    next_hit = np.full(len(word)+motif_width+1, len(word)+motif_width, dtype=np.int64)
    if len(hits):
        next_hit[hits] = hits
        next_hit = np.minimum.accumulate(next_hit[::-1])[::-1]
    for width in widths:
        start = np.arange(max(0, len(word)-width+1)) if width else np.array([0])
        end = start+width if width else np.array([len(word)])
        possible = np.searchsorted(valid, end-motif_width, side="right")-np.searchsorted(valid, start)
        count = np.zeros(len(start), float)
        current = next_hit[start]
        for _ in range((width if width else len(word))//motif_width):
            active = current <= end-motif_width
            if not active.any():
                break
            count[active] += 1
            current[active] = next_hit[current[active]+motif_width]
        coverage = _union_lengths(hits, start, end, motif_width)
        opportunity = _union_lengths(valid, start, end, motif_width)
        with np.errstate(invalid="ignore", divide="ignore"):
            coverage /= opportunity
        count[possible <= 0] = np.nan
        coverage[possible <= 0] = np.nan
        output[width] = (count, coverage)
    return output


def signed_ranks(differences, minimum=20):
    n = np.isfinite(differences).sum(axis=0)
    values = np.abs(differences)
    values[(values == 0) | ~np.isfinite(values)] = np.nan
    ranks = rankdata(values, axis=0, method="average", nan_policy="omit")
    ranks = np.nan_to_num(ranks)*np.sign(np.nan_to_num(differences))
    denominator = np.sqrt(np.sum(ranks*ranks, axis=0))
    denominator[denominator == 0] = 1
    observed = ranks.sum(axis=0)/denominator
    ranks /= denominator
    observed[n < minimum] = np.nan
    return ranks, observed, n


def null_cluster_max(values, threshold):
    positive = np.isfinite(values) & (values > threshold)
    sums = np.cumsum(np.where(positive, values-threshold, 0), axis=1)
    reset = np.maximum.accumulate(np.where(~positive, sums, 0), axis=1)
    return np.max(sums-reset, axis=1, initial=0)


def model_datasets(rna, pairs, pattern):
    regions = rna[rna.layer.isin(["affected_sequence", "regulatory_pre_mRNA"])]
    lookup = {(r.event_id, r.region): r.sequence for r in regions.itertuples()}
    cache = {}
    for kind, direction in sorted(set(zip(pairs.event_type, pairs.direction))):
        group = pairs[pairs.event_type.eq(kind) & pairs.direction.eq(direction)]
        for region in sorted(set(regions.loc[regions.event_type.eq(kind), "region"])):
            words = [(lookup.get((r.foreground_event, region), ""), lookup.get((r.background_event, region), "")) for r in group.itertuples()]
            for word in {w for pair in words for w in pair if w}:
                if word not in cache:
                    cache[word] = profiles(word, pattern)
            for width in (0, 10, 20, 50):
                span = max(1, max([max(0, min(len(a), len(b))-width+1) for a, b in words], default=0)) if width else 1
                for m, method in enumerate(("finite_nonoverlap_counts", "overlapping_union_coverage")):
                    differences = np.full((len(group), span), np.nan)
                    for i, (a, b) in enumerate(words):
                        if not a or not b:
                            continue
                        x, y = cache[a][width][m], cache[b][width][m]
                        if width and region.endswith("acceptor_intron"):
                            x, y = x[::-1], y[::-1]
                        n = min(span, len(x), len(y))
                        differences[i, :n] = x[:n]-y[:n]
                    coordinate = "masked_acceptor_3prime_minus20" if region.endswith("acceptor_intron") else "masked_donor_5prime_plus6" if region.endswith("donor_intron") else "affected_RNA_5prime"
                    yield {"meta": {"event_type": kind, "direction": direction, "region": region, "method": method,
                        "window_nt": width, "coordinate_system": coordinate}, "differences": differences,
                        "pair_ids": group.pair_id.tolist(), "foreground_events": group.foreground_event.tolist(),
                        "foreground_genes": group.foreground_gene.tolist()}
