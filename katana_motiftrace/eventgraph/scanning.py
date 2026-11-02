"""Packaging adapter for the public-validation scan; frozen numerical kernels.

The traversal, positional truncation, gene blocks and one complete correction
family follow scripts/scan.py from the completed 2026-10-31 validation. Only I/O,
cancellation and descriptive output columns are added here.
"""
from pathlib import Path
import gzip
import json
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from .inference import BlockFamily, pwm_profiles, bh

WINDOW_COLUMNS = ["event_type", "direction", "region", "model_id", "method", "window_nt", "layer",
                  "window_start0", "n_matched_events", "n_independent_gene_blocks", "paired_rank_Z",
                  "mean_paired_score_difference", "raw_P", "RBP", "coordinate_system",
                  "n_matched_gene_pairs", "raw_P_new_paired_swaps", "q_new_window_family", "status"]
CLUSTER_COLUMNS = ["event_type", "direction", "region", "model_id", "method", "window_nt", "layer",
                   "first_window_start0", "last_window_start0", "cluster_mass", "minimum_events",
                   "minimum_gene_blocks", "cluster_P_FWER", "cluster_P_FWER_new_family", "RBP",
                   "covered_region_start0", "covered_region_end0", "coordinate_system",
                   "number_correlated_windows", "supporting_independent_gene_pairs"]


def coordinate_system(region):
    if region.endswith('acceptor_intron'):
        return 'distance upstream from masked acceptor edge (-20nt), 3prime-aligned'
    if region.endswith('donor_intron'):
        return 'distance from masked donor edge (+6nt)'
    return 'affected RNA 5prime-relative'


def scan_family(rna, matching, models, destination, *, permutations=999, seed=20261015,
                check_cancelled=lambda: None):
    """Scan all declared PWMs, then apply one global BH and one max-cluster FWER.

``n_matched_gene_pairs`` is deliberately missing for this engine. Independent
foreground gene blocks, including every matched event, are recorded separately.
The legacy raw-P column is a value-preserving export alias for the joint swaps
of matched differences within each gene; it never contains a regional P value.
"""
    if not models or any(m.get('kind') != 'PWM' for m in models.values()):
        raise ValueError('EventGraph requires PWM models; select legacy analysis_engine for consensus models. No models were discarded.')
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    pairs = matching[matching.status.eq('matched')].copy()
    pairs['block'] = pairs.foreground_gene
    family = BlockFamily(pairs.block, permutations, seed, 20)
    regions = rna[rna.layer.isin(['affected_sequence', 'regulatory_pre_mRNA'])]
    lookup = {(r.event_id, r.region): r.sequence for r in regions.itertuples()}
    writer = None
    clusters, counts = [], []
    raw_path = destination / 'WINDOWS_RAW.parquet'
    if raw_path.exists():
        raise FileExistsError(raw_path)
    try:
        for kind in sorted(regions.event_type.unique()):
            for direction in ['inclusion_down', 'inclusion_up']:
                group = pairs[pairs.event_type.eq(kind) & pairs.direction.eq(direction)].copy()
                if group.block.nunique() < 20:
                    counts.append(dict(event_type=kind, direction=direction, status='fewer_than_20_independent_gene_blocks',
                                       estimable_positions=0, max_events=len(group), max_genes=group.block.nunique()))
                    continue
                events = set([*group.foreground_event, *group.background_event])
                for region in sorted(regions.loc[regions.event_type.eq(kind), 'region'].unique()):
                    seqs = {ev: lookup.get((ev, region), '') for ev in events}
                    for token, model in models.items():
                        check_cancelled()
                        profiles = {ev: pwm_profiles(seq, model) for ev, seq in seqs.items()}
                        for width in [0, 10, 20, 50]:
                            check_cancelled()
                            lens = np.array([min(len(profiles[r.foreground_event][width]), len(profiles[r.background_event][width])) for r in group.itertuples()])
                            lengths = pd.DataFrame({'gene': group.block.to_numpy(), 'length': lens}).groupby('gene').length.max().sort_values(ascending=False)
                            span = int(lengths.iloc[19])
                            meta = dict(event_type=kind, direction=direction, region=region, model_id=token,
                                        method='PWM_mean_odds', window_nt=width,
                                        layer='affected_sequence' if region.startswith('affected_') else 'regulatory_pre_mRNA')
                            if not span:
                                counts.append(dict(**meta, estimable_positions=0, max_events=len(group),
                                                   max_genes=group.block.nunique(), status='insufficient_shared_sequence_length'))
                                continue
                            diff = np.full((len(group), span), np.nan)
                            for j, row in enumerate(group.itertuples()):
                                a = profiles[row.foreground_event][width]
                                b = profiles[row.background_event][width]
                                if region.endswith('acceptor_intron') and width:
                                    a, b = a[::-1], b[::-1]
                                n = min(len(a), len(b), span)
                                if n:
                                    diff[j, :n] = a[:n] - b[:n]
                            rows, found = family.evaluate(diff, group.block.tolist(), meta)
                            rows = rows[rows.raw_P.notna()]
                            if len(rows):
                                tab = pa.Table.from_pandas(rows, preserve_index=False)
                                if writer is None:
                                    writer = pq.ParquetWriter(raw_path, tab.schema, compression='zstd')
                                writer.write_table(tab)
                            clusters.extend(found)
                            counts.append(dict(**meta, estimable_positions=len(rows),
                                               max_events=int(rows.n_matched_events.max()) if len(rows) else 0,
                                               max_genes=int(rows.n_independent_gene_blocks.max()) if len(rows) else 0,
                                               status='computed' if len(rows) else 'not_estimable'))
    finally:
        if writer is not None:
            writer.close()
    check_cancelled()
    raw = pq.read_table(raw_path, columns=['raw_P']).to_pandas().raw_P.to_numpy() if raw_path.exists() else np.array([])
    q = bh(raw)
    with gzip.open(destination / 'NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz', 'wt', encoding='utf-8', newline='') as stream:
        offset = 0
        if raw_path.exists():
            for batch in pq.ParquetFile(raw_path).iter_batches(batch_size=100000):
                check_cancelled()
                frame = batch.to_pandas()
                frame['RBP'] = frame.model_id.map(lambda k: models[k].get('RBPs', ''))
                frame['coordinate_system'] = frame.region.map(coordinate_system)
                frame['n_matched_gene_pairs'] = np.nan  # Not a valid unit for EventGraph.
                frame['raw_P_new_paired_swaps'] = frame.raw_P
                frame['q_new_window_family'] = q[offset:offset+len(frame)]
                frame['status'] = 'computed_conditional_on_gene_block_symmetry'
                frame[WINDOW_COLUMNS].to_csv(stream, sep='\t', index=False, header=offset == 0)
                offset += len(frame)
        else:
            pd.DataFrame(columns=WINDOW_COLUMNS).to_csv(stream, sep='\t', index=False)
    assert offset == len(q)
    corrected = family.adjust_clusters(clusters)
    if len(corrected):
        corrected['cluster_P_FWER_new_family'] = corrected.cluster_P_FWER
        corrected['RBP'] = corrected.model_id.map(lambda k: models[k].get('RBPs', ''))
        corrected['covered_region_start0'] = corrected.first_window_start0
        corrected['covered_region_end0'] = corrected.last_window_start0 + corrected.window_nt
        corrected['coordinate_system'] = corrected.region.map(coordinate_system)
        corrected['number_correlated_windows'] = corrected.last_window_start0 - corrected.first_window_start0 + 1
        corrected['supporting_independent_gene_pairs'] = np.nan
    corrected = corrected.reindex(columns=CLUSTER_COLUMNS)
    corrected.to_csv(destination / 'NEW_CLUSTER_FWER.tsv', sep='\t', index=False)
    np.save(destination / 'MAX_CLUSTER_NULL.npy', family.null_max, allow_pickle=False)
    pd.DataFrame(counts).to_csv(destination / 'TEST_FAMILY_ESTIMABILITY.tsv', sep='\t', index=False)
    result = dict(eligible_new_window_tests=len(raw), new_window_q_lt_0_05=int((q <= .05).sum()),
                  cluster_P_FWER_lt_0_05=int(corrected.cluster_P_FWER.le(.05).sum()),
                  cluster_rows=len(corrected), matched_events=len(pairs),
                  matched_independent_gene_blocks=int(pairs.block.nunique()),
                  minimum_q=float(q.min()) if len(q) else None)
    (destination / 'STATISTICAL_COLUMNS.json').write_text(json.dumps({
        'analysis_engine': 'eventgraph_20261030', 'window_FDR_column': 'q_new_window_family',
        'window_FDR_family': 'BH across all estimable scalar/window tests in this complete comparison',
        'regional_P_column': 'cluster_P_FWER_new_family',
        'regional_P_family': 'maximum cluster mass across all models/types/directions/regions/widths; FWER, not FDR',
        'raw_P_new_paired_swaps': 'Export alias of raw_P; matched differences swap jointly within foreground gene blocks',
        'n_matched_gene_pairs': 'Not applicable; missing. Use n_matched_events and n_independent_gene_blocks.',
        'region_significance_implies_window_significance': False,
        'nonestimable_groups': 'TEST_FAMILY_ESTIMABILITY.tsv; not included in BH denominator',
    }, indent=2), encoding='utf-8')
    return result, corrected
