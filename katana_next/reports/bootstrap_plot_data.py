"""Recover exact plotted bootstrap values from a completed native family."""
import hashlib
import json
from pathlib import Path


def _event_keys(path):
    """Read native NumPy string/object keys without unrestricted pickle globals."""
    import numpy as np
    import pickle
    import zipfile
    class StringArrayReader(pickle.Unpickler):
        def find_class(self, module, name):
            if module in ('numpy.core.multiarray', 'numpy._core.multiarray') and name == '_reconstruct':
                from importlib import import_module
                return import_module(module)._reconstruct
            if module == 'numpy' and name in ('ndarray', 'dtype'):
                return getattr(np, name)
            raise ValueError('Unapproved object in native bootstrap event keys')
    with zipfile.ZipFile(path) as archive, archive.open('event_keys.npy') as stream:
        version = np.lib.format.read_magic(stream)
        if version == (1, 0):
            shape, fortran, dtype = np.lib.format.read_array_header_1_0(stream)
        elif version == (2, 0):
            shape, fortran, dtype = np.lib.format.read_array_header_2_0(stream)
        else:
            raise ValueError('Unsupported bootstrap event-key encoding')
        if dtype.hasobject:
            values = StringArrayReader(stream).load()
        else:
            with np.load(path, allow_pickle=False) as saved:
                values = saved['event_keys']
    if not isinstance(values, np.ndarray) or values.ndim != 1 or values.shape != shape or any(not isinstance(v, str) for v in values):
        raise ValueError('Bootstrap event keys must be a one-dimensional string array')
    return values.astype(str)


def bootstrap_plot_links(artifact, panel, output, original_links):
    metadata = artifact.get('metadata', {})
    key = metadata.get('flow_artifact_key')
    if metadata.get('flow_stage') != 'network.bootstrap' or key not in ('heatmap', 'sizes'):
        return original_links
    import numpy as np
    import pandas as pd
    image = Path(artifact['path']).resolve()
    expected = {'heatmap': 'bootstrap_coassignment_heatmap.png', 'sizes': 'consensus_module_sizes.png'}[key]
    if image.name != expected or image.parent.name != 'plots' or image.parent.parent.name != '04c_consensus_modules':
        raise ValueError('Bootstrap plot does not match its declared native artifact key')
    root = image.parent.parent
    sources = [root / 'consensus_event_module_assignments.tsv', root / 'bootstrap_coassignment_matrix.npz'] if key == 'heatmap' else [root / 'consensus_module_summary.tsv']
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    token = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()[:24]
    target = Path(output) / 'outputs/bootstrap-plot-data' / token
    target.mkdir(parents=True, exist_ok=True)
    if key == 'heatmap':
        assignments = pd.read_csv(sources[0], sep='\t')
        visible = assignments.sort_values(['consensus_module', 'consensus_membership_confidence'], ascending=[True, False])
        visible = visible.loc[visible.consensus_module.ne('grey')].head(250)
        with np.load(sources[1], allow_pickle=False) as saved:
            events = _event_keys(sources[1])
            matrix = saved['coassignment']
            if len(set(events)) != len(events) or matrix.shape != (len(events), len(events)):
                raise ValueError('Bootstrap matrix identity or dimensions are invalid')
            positions = pd.Series(np.arange(len(events)), index=events).reindex(visible.event_key)
            if positions.isna().any():
                raise ValueError('Bootstrap plotted event is missing from its saved matrix')
            idx = positions.astype(int).to_numpy()
            data = {'event_keys': events[idx].tolist(), 'coassignment': matrix[np.ix_(idx, idx)].tolist(),
                    'selection': 'Native plot order: module ascending, confidence descending; non-grey first 250'}
        plotted = target / 'plotted-matrix.json'
        plotted.write_text(json.dumps(data, allow_nan=False), encoding='utf-8')
    else:
        data = pd.read_csv(sources[0], sep='\t').sort_values('n_events', ascending=False)
        plotted = target / 'plotted-module-sizes.tsv'
        data.to_csv(plotted, sep='\t', index=False)
    common = {'contrast_id': panel['contrast_id'], 'policy': panel.get('policy'), 'unit': 'native bootstrap event/module'}
    links = [{**common, 'role': 'plot_data', 'path': str(plotted), 'label': 'Exact plotted bootstrap values'}]
    links += [{**common, 'role': 'full_results', 'path': str(p), 'label': 'Complete native ' + p.name} for p in sources]
    record = {'native_sources_sha256': hashes, 'original_links': original_links, 'scientific_recomputation': False,
              'plot_values_sha256': hashlib.sha256(plotted.read_bytes()).hexdigest()}
    (target / 'binding.json').write_text(json.dumps(record, indent=2), encoding='utf-8')
    if hashes != {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}:
        raise ValueError('Bootstrap source changed during report preparation')
    return links
