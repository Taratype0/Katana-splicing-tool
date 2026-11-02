"""Input/output integration of the frozen event-level matching implementation."""
import numpy as np
import pandas as pd
from katana_next.science.motiftrace_events import psi_values
from katana_next.science.motiftrace_full_execution import fast_match
from .matching import global_match


def match_controls(events, selected_ids, rna):
    features = []
    for token, group in rna[rna.layer.eq('affected_sequence')].groupby('event_id', sort=False):
        seq = ''.join(group.sequence)
        valid = sum(seq.count(base) for base in 'ACGT')
        features.append(dict(event_id=token, length=len(seq),
                             GC=(seq.count('G')+seq.count('C'))/valid if valid else np.nan,
                             ambiguous_fraction=1-valid/len(seq) if seq else 1))
    data = events.merge(pd.DataFrame(features, columns=['event_id', 'length', 'GC', 'ambiguous_fraction']),
                        on='event_id', how='left', validate='one_to_one')
    data['length'] = data.length.fillna(0)
    data['ambiguous_fraction'] = data.ambiguous_fraction.fillna(1)
    data['baseline_PSI'] = [np.mean(values) if values else np.nan for values in
                            [[p['PSI'] for p in psi_values(r) if p['group']=='denominator' and np.isfinite(p['PSI'])]
                             for r in data.to_dict('records')]]
    data['depth'] = np.minimum(data.numerator_mean_reads, data.denominator_mean_reads)
    data['direction'] = np.where(data.dpsi > 0, 'inclusion_up', 'inclusion_down')
    forbidden = set(data.loc[data.stable_significant.fillna(False) | data.fdr.le(.05) | data.event_id.isin(selected_ids), 'gene_id'])
    foreground = data[data.event_id.isin(selected_ids)].copy()
    background = data[data.measurement_eligible.fillna(False) & data.fdr.gt(.1) & data.dpsi.abs().lt(.05)
                      & ~data.gene_id.isin(forbidden)].copy()
    legacy = fast_match(events, selected_ids, rna)
    protected = set(legacy.loc[legacy.status.eq('matched'), 'foreground_event'])
    foreground['retain_original_match'] = foreground.event_id.isin(protected)
    matching, certificate, _ = global_match(foreground, background)
    assert protected <= set(matching.loc[matching.status.eq('matched'), 'foreground_event'])
    assert len(matching) == len(selected_ids)
    return matching, certificate, legacy
