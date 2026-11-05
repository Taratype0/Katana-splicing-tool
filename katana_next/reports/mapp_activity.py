"""Native sample activities with explicit experimental-condition labels."""
from pathlib import Path


def activity_axis(columns, project):
    samples = {str(s['sample_id']): s for s in project.get('samples', [])}
    axis = []
    for column in columns:
        sample_id = str(column).removeprefix('activities_')
        condition = str(samples.get(sample_id, {}).get('condition') or '')
        axis.append({'column': column, 'sample_id': sample_id, 'condition': condition,
                     'label': condition + '\n' + sample_id if condition else sample_id})
    return axis


def render_activity(frame, project, comparison, destination):
    import numpy as np
    import pandas as pd
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.colors import LinearSegmentedColormap
    columns = [c for c in frame if str(c).startswith('activities_')]
    axis = activity_axis(columns, project)
    values = frame[columns].apply(pd.to_numeric, errors='coerce').to_numpy()
    finite = np.abs(values[np.isfinite(values)])
    bound = float(finite.max()) if len(finite) else 1.
    fig = Figure(figsize=(max(10, .85*len(columns)), max(4, .28*len(frame)) + .8))
    FigureCanvasAgg(fig)
    ax = fig.subplots()
    cmap = LinearSegmentedColormap.from_list('katana-maei', ['#2863A5', '#FFFFFF', '#C63F4A'])
    image = ax.imshow(np.ma.masked_invalid(values), cmap=cmap, vmin=-bound or -1,
                      vmax=bound or 1, aspect='auto')
    ax.set_xticks(range(len(axis)), [x['label'] for x in axis], rotation=45, ha='right')
    for i in range(1, len(axis)):
        if axis[i]['condition'] and axis[i-1]['condition'] and axis[i]['condition'] != axis[i-1]['condition']:
            ax.axvline(i-.5, color='#1e2630', linewidth=1)
    ax.set_xlabel('Experimental condition / individual sample')
    ax.set_yticks(range(len(frame)), [str(s)+' | '+str(m) for s,m in zip(frame.site, frame.native_model_id)], fontsize=7)
    ax.set_title('Native MAEI motif activity | '+comparison)
    fig.colorbar(image, ax=ax, label='Original native sample activity')
    fig.tight_layout()
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=150)
    fig.clear()
    return axis


def refresh_activity_panel(panel, project, report_root):
    """Redraw a copied report from the exact saved display rows, without fits."""
    from hashlib import sha256
    import pandas as pd
    metadata = panel.get('metadata', {})
    if metadata.get('component') != 'mapp_native' or not project.get('samples'):
        return panel, None
    expected = metadata.get('native_plot_data_sha256')
    rows = [link for link in panel.get('data_links', []) if link.get('role') == 'plot_data']
    if not expected or len(rows) != 1:
        return panel, None
    path = Path(rows[0]['path'])
    if not path.is_absolute():
        path = Path(report_root) / path
    if sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError('Native MAPP plotted activity rows changed')
    cid = panel['contrast_id']
    destination = Path(report_root)/'figures/mapp'/(sha256(cid.encode()).hexdigest()[:16]+'-activities.png')
    axis = render_activity(pd.read_csv(path, sep='\t'), project, cid, destination)
    rendered_sha = sha256(destination.read_bytes()).hexdigest()
    result = {**panel, 'path': str(destination), 'metadata': {
        **metadata, 'activity_axis': axis,
        'activity_rendering': 'condition_sample_labels/v1',
        'rendered_activity_image_sha256': rendered_sha}}
    result['caption'] = panel.get('caption', '') + ' Columns identify experimental condition and individual sample; biological replicates are not averaged.'
    return result, {'component': 'mapp_condition_labels', 'comparison_id': cid,
                    'source_plot_data_sha256': expected, 'activity_axis': axis,
                    'rendered_activity_image_sha256': rendered_sha,
                    'native_values_unchanged': True, 'new_fits': 0}


def validate_activity_image(panel, project, report_root):
    """Authenticate the displayed image without replacing native provenance."""
    from hashlib import sha256
    import pandas as pd
    root = Path(report_root).resolve()
    metadata = panel['metadata']

    def member(relative):
        path = (root / relative).resolve()
        if Path(relative).is_absolute() or not path.is_relative_to(root):
            raise ValueError('Native MAEI image/data escapes its report')
        return path

    expected = metadata.get('native_activity_image_sha256')
    if metadata.get('activity_rendering') is not None:
        if metadata['activity_rendering'] != 'condition_sample_labels/v1':
            raise ValueError('Unknown native MAEI activity rendering')
        links = [row for row in panel.get('data_links', []) if row.get('role') == 'plot_data']
        if len(links) != 1:
            raise ValueError('Rendered native MAEI activity lacks its exact source rows')
        data = member(links[0]['path'])
        if sha256(data.read_bytes()).hexdigest() != metadata.get('native_plot_data_sha256'):
            raise ValueError('Rendered native MAEI activity source rows changed')
        columns = [c for c in pd.read_csv(data, sep='\t', nrows=0) if str(c).startswith('activities_')]
        if not columns or metadata.get('activity_axis') != activity_axis(columns, project):
            raise ValueError('Rendered native MAEI activity axis differs from the project samples')
        expected = metadata.get('rendered_activity_image_sha256')
    elif metadata.get('rendered_activity_image_sha256') is not None:
        raise ValueError('Rendered native MAEI image has no declared rendering')
    image = member(panel.get('image', ''))
    if not image.is_file() or sha256(image.read_bytes()).hexdigest() != expected:
        raise ValueError('Native MAEI activity panel is absent or changed')
