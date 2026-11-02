"""Assemble comparable experiments into one individual-sample discovery matrix.

The manifest, not a contrast, defines the cohort. Duplicate sample measurements
must agree; NA can be filled by another source but finite values are never averaged.
"""
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from study import sha, validate_contract
from cosplice_research import write_json


def assemble(cfg, output=None):
    from product import read_matrix
    sources = cfg.get('psi_sources') or []
    if not sources:
        sources = [{'psi': cfg['psi'], 'read_support': cfg.get('read_support')}]
    elif cfg.get('psi') or cfg.get('counts') or cfg.get('read_support'):
        raise ValueError('Choose a study input bundle OR one complete matrix; support belongs to each source')
    if cfg.get('value_kind', 'usage') != 'usage':
        raise ValueError('Study bundles require individual-sample PSI, not contrast response profiles')
    paths = [str(Path(s['psi']).resolve()) for s in sources]
    if len(paths) != len(set(paths)):
        raise ValueError('Duplicate input file; list each experiment matrix once')
    meta = pd.read_csv(cfg['samples'], sep=',' if str(cfg['samples']).endswith('.csv') else '\t', dtype=str)
    contract = validate_contract(cfg, meta)
    expected = meta.sample_id.tolist()
    supported = [bool(s.get('read_support')) for s in sources]
    if any(supported) and not all(supported):
        raise ValueError('Read support must be provided for every source or none')
    merged = pd.DataFrame(); support = pd.DataFrame(); records = []; occurrences = {}

    def join(old, new, label):
        ix = old.index.intersection(new.index); cols = old.columns.intersection(new.columns)
        if len(ix) and len(cols):
            a = old.loc[ix, cols].to_numpy(); b = new.loc[ix, cols].to_numpy()
            # Allow only numerical serialization noise, not measurement differences.
            bad = np.isfinite(a) & np.isfinite(b) & ~np.isclose(a, b, rtol=0, atol=1e-12)
            if bad.any():
                i, j = np.argwhere(bad)[0]
                raise ValueError(f'Conflicting {label} for shared sample {cols[j]}, event {ix[i]}; reconcile upstream')
        rows = old.index.union(new.index, sort=False); columns = old.columns.union(new.columns, sort=False)
        return old.reindex(index=rows, columns=columns).combine_first(new).reindex(index=rows, columns=columns)

    for source in sources:
        x = read_matrix(source['psi']); values = x.to_numpy()
        if np.isinf(values).any() or ((values < 0) | (values > 1)).any():
            raise ValueError('Study PSI values must be within [0,1] or NA')
        unknown = sorted(set(x.columns) - set(expected))
        if unknown: raise ValueError(f'Input has samples outside the whole-study manifest: {unknown[:5]}')
        for sample in x.columns: occurrences[sample] = occurrences.get(sample, 0) + 1
        merged = join(merged, x, 'PSI')
        record = {'psi': str(Path(source['psi']).resolve()), 'psi_sha256': sha(source['psi']),
                  'events': len(x), 'samples': list(x.columns)}
        if source.get('read_support'):
            r = read_matrix(source['read_support'])
            if set(r.index) != set(x.index) or set(r.columns) != set(x.columns):
                raise ValueError('Each support matrix must match its source PSI event/sample identities')
            v = r.to_numpy()
            if np.isinf(v).any() or (v < 0).any(): raise ValueError('Support must be nonnegative or NA')
            support = join(support, r, 'read support')
            record.update(read_support=str(Path(source['read_support']).resolve()), support_sha256=sha(source['read_support']))
        records.append(record)
    missing = sorted(set(expected) - set(merged.columns))
    if missing: raise ValueError(f'Whole-study samples absent from all sources: {missing[:5]}')
    merged = merged.loc[:, expected]
    audit = {'schema': 1, 'analysis_scope': 'whole_study', 'study': contract, 'sources': records,
             'source_count': len(sources), 'source_sample_occurrences': sum(occurrences.values()),
             'unique_samples': len(expected), 'duplicate_sample_occurrences_removed': sum(occurrences.values()) - len(expected),
             'shared_samples': {s: n for s, n in occurrences.items() if n > 1}, 'union_events': len(merged),
             'missing_fraction': float(merged.isna().to_numpy().mean()), 'read_support_available': all(supported),
             'metadata_sha256': sha(cfg['samples']), 'sample_order': expected,
             'reconciliation': 'Union of event IDs; individual sample IDs once; shared finite values agree within 1e-12; NA filled only from observed source values; no averaging',
             'scope_limit': 'All eligible samples/unscreened events are user-declared. Events absent from every source cannot be recovered; matching IDs must mean the same quantification.'}
    if output is not None:
        output = Path(output); output.mkdir(parents=True, exist_ok=False)
        merged.to_csv(output/'psi.tsv.gz', sep='\t', index_label='event_key', float_format='%.17g')
        if all(supported): support.reindex(index=merged.index, columns=expected).to_csv(output/'read_support.tsv.gz', sep='\t', index_label='event_key', float_format='%.17g')
        write_json(output/'assembly.json', audit)
        resolved = {**cfg, 'psi_sources': [], 'psi': str(output/'psi.tsv.gz'),
                    'read_support': str(output/'read_support.tsv.gz') if all(supported) else None}
        write_json(output/'resolved_config.json', resolved)
        return resolved, audit
    return audit


if __name__ == '__main__':
    try:
        print(json.dumps(assemble(json.loads(sys.stdin.read())), ensure_ascii=False))
    except Exception as exc:
        print(json.dumps({'error': str(exc)}, ensure_ascii=False)); sys.exit(1)
