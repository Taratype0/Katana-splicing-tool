"""Display assay identity from recorded design, never from filenames."""

def assay_identity(metadata, assets=(), recorded_mode=None):
    raw = recorded_mode or metadata.get('assay_mode') or metadata.get('assay')
    text = str(raw or '').lower().replace('_', '-').replace(' ', '')
    if text.startswith(('single', 'scrna', 'snrna')) or text in {'sc', 'sn'}:
        mode = 'single-cell'
    elif text.startswith('bulk'):
        mode = 'bulk'
    elif any(str(a.get('kind', '') if isinstance(a, dict) else getattr(a, 'kind', '')).startswith('single_cell_') for a in assets):
        mode = 'single-cell'
    else:
        mode = 'unknown'
    return {'mode': mode, 'recorded': raw or None,
            'label': {'single-cell': 'Single-cell / single-nucleus', 'bulk': 'Bulk RNA-seq',
                      'unknown': 'Assay not recorded'}[mode],
            'replicate_unit': 'donor / biological sample; cells are not independent replicates' if mode == 'single-cell' else 'biological sample'}
