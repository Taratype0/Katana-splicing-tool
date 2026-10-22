"""Source-faithful single-cell result adapters, not replacements for source methods."""
from pathlib import Path
import base64
import hashlib
import html
import json
import numpy as np
import pandas as pd

def read_config(path):
    if not isinstance(path, (str, Path)) or (isinstance(path, str) and (not path.strip())):
        raise ValueError('Select an import configuration JSON file (config_path is required).')
    try:
        path = Path(path).resolve()
        if not path.is_file():
            raise ValueError(f'Import configuration must be an existing regular file: {path}')
        text = path.read_text(encoding='utf-8-sig')
    except (OSError, UnicodeError, RuntimeError) as exc:
        raise ValueError(f'Cannot read import configuration file: {path} ({exc})') from exc
    try:
        config = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError(f'Invalid import configuration JSON: {path} (line {exc.lineno}, column {exc.colno})') from exc
    if not isinstance(config, dict):
        raise ValueError('Import configuration JSON must contain an object, not a list or scalar.')
    for key in ('metadata_path', 'results_path', 'psi_path', 'counts_path', 'annotation_path'):
        value = config.get(key)
        if value is not None and value != '':
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f'Import configuration {key} must be a nonempty file path string.')
            try:
                config[key] = str((path.parent / value).resolve())
            except (OSError, ValueError, RuntimeError) as exc:
                raise ValueError(f'Invalid import configuration file path for {key}: {exc}') from exc
    return config

def input_files(config):
    return [Path(config[k]) for k in ('metadata_path', 'results_path', 'psi_path', 'counts_path', 'annotation_path') if config.get(k)]

def read_table(path):
    path = Path(path)
    return pd.read_parquet(path) if path.suffix == '.parquet' else pd.read_csv(path, sep='\t', dtype={'cell_id': str, 'gene_id': str, 'transcript_id': str, 'tran_id': str})

def metadata(config):
    frame = read_table(config['metadata_path'])
    for col in ('cell_id', 'donor_id', 'cell_type', 'condition'):
        if col not in frame or frame[col].isna().any() or frame[col].astype(str).str.strip().eq('').any():
            raise ValueError(f'Missing cell identity/design column: {col}')
    if frame.cell_id.duplicated().any():
        raise ValueError('Duplicate cell_id; prefix sample-specific barcodes before import')
    return frame

def bounded(series, lower, upper, name):
    values = pd.to_numeric(series, errors='raise')
    valid = values.dropna()
    if not np.isfinite(valid).all() or not valid.between(lower, upper).all():
        raise ValueError(f'{name} outside [{lower}, {upper}]')
    return values

def import_marvel(config):
    """Import a full CompareValues.PSI table and optional native wide PSI export."""
    meta = metadata(config)
    for key in ('method', 'tool_version', 'reference', 'target', 'cell_type', 'adjustment'):
        if not config.get(key):
            raise ValueError(f'MARVEL import requires {key}')
    if config['reference'] == config['target']:
        raise ValueError('Reference and target must differ')
    selected = meta.loc[meta.cell_type.eq(config['cell_type']) & meta.condition.isin([config['reference'], config['target']])]
    if set(selected.condition) != {config['reference'], config['target']}:
        raise ValueError('Both declared conditions must exist in the declared cell type')
    native = read_table(config['results_path'])
    required = {'tran_id', 'mean.g1', 'mean.g2', 'mean.diff', 'p.val', 'p.val.adj'}
    if not required.issubset(native):
        raise ValueError(f'Missing MARVEL native columns: {sorted(required - set(native))}')
    if native.tran_id.isna().any() or native.tran_id.duplicated().any():
        raise ValueError('MARVEL event IDs must be present and unique per imported contrast')
    result = native.copy()
    result['event_id'] = native.tran_id.astype(str)
    for col in ('mean.g1', 'mean.g2'):
        result[col] = bounded(native[col], 0, 100, col)
    diff = bounded(native['mean.diff'], -100, 100, 'mean.diff')
    if not np.allclose(diff, result['mean.g2'] - result['mean.g1'], equal_nan=True, atol=1e-05):
        raise ValueError('MARVEL mean.diff must equal group2 minus group1 in percentage points')
    result['delta_psi'] = diff / 100
    result['source_pvalue'] = bounded(native['p.val'], 0, 1, 'p.val')
    result['source_adjusted_value'] = bounded(native['p.val.adj'], 0, 1, 'p.val.adj')
    if config['adjustment'] not in ('BH', 'none', 'unknown'):
        raise ValueError('Declare adjustment as BH, none or unknown; no automatic relabeling')
    if config['method'] in ('dts', 'permutation') and config['adjustment'] == 'BH':
        raise ValueError('Native MARVEL dts/permutation p.val.adj is not assumed BH; import as none/unknown')
    result['source_fdr'] = result.source_adjusted_value if config['adjustment'] == 'BH' else np.nan
    result['passes_source_fdr_0_05'] = result.source_fdr.le(0.05)
    for key in ('reference', 'target', 'cell_type', 'method', 'tool_version'):
        result[key] = config[key]
    result['inferential_scope'] = 'source cell-group test; not donor-validated by Katana'
    tables = {'marvel_complete_results': result}
    if config.get('psi_path'):
        wide = read_table(config['psi_path'])
        if 'tran_id' not in wide or wide.tran_id.isna().any() or wide.tran_id.duplicated().any():
            raise ValueError('Native PSI matrix requires unique tran_id rows')
        cells = [c for c in wide if c != 'tran_id']
        if not set(cells).issubset(set(meta.cell_id)):
            raise ValueError('PSI cells absent from metadata')
        long = wide.melt(id_vars='tran_id', var_name='cell_id', value_name='psi').rename(columns={'tran_id': 'event_id'})
        long['psi'] = bounded(long.psi, 0, 1, 'native PSI')
        long = long.merge(meta, on='cell_id', validate='many_to_one')
        tables['marvel_cell_psi'] = long
        tables['marvel_donor_descriptive_psi'] = long.groupby(['event_id', 'donor_id', 'cell_type', 'condition'], as_index=False).agg(mean_cell_psi=('psi', 'mean'), measured_cells=('psi', 'count'))
    tables['donor_design'] = selected.groupby(['condition', 'donor_id', 'cell_type'], as_index=False).agg(cells=('cell_id', 'size'))
    return tables

def import_longread(config):
    """Join quantified molecules to SQANTI/Pigeon classification and existing donor DIU."""
    from katana_next.science.single_cell_isoforms import validate_isoform_inputs, analyze_isoform_usage
    meta = metadata(config)
    counts = read_table(config['counts_path'])
    counts, meta = validate_isoform_inputs(counts, meta)
    annotation = read_table(config['annotation_path'])
    id_col = config.get('annotation_id_column', 'isoform')
    gene_col = config.get('annotation_gene_column', 'associated_gene')
    if not {id_col, gene_col, 'structural_category'}.issubset(annotation):
        raise ValueError('Classification requires isoform ID, associated gene and structural_category')
    annotation = annotation.rename(columns={id_col: 'transcript_id', gene_col: 'annotation_gene_id'})
    if annotation.transcript_id.isna().any() or annotation.transcript_id.duplicated().any():
        raise ValueError('Classification IDs must be unique; harmonize across samples first')
    if not set(counts.transcript_id).issubset(set(annotation.transcript_id)):
        raise ValueError('Unmatched transcript IDs; do not join PB IDs across samples without harmonization')
    joined = counts.merge(annotation, on='transcript_id', validate='many_to_one')
    if not joined.gene_id.astype(str).eq(joined.annotation_gene_id.astype(str)).all():
        raise ValueError('Count/annotation gene mismatch; explicitly reconcile gene identifiers')
    if not config.get('tool_version') or not config.get('reference_build'):
        raise ValueError('Declare source tool_version and reference_build')
    measured = joined.loc[joined['count'].gt(0)].merge(meta, on='cell_id', validate='many_to_one')
    summary = measured.groupby(['condition', 'cell_type', 'structural_category'], as_index=False).agg(molecules=('count', 'sum'), transcripts=('transcript_id', 'nunique'), genes=('gene_id', 'nunique'))
    analysis = analyze_isoform_usage(counts, meta, numerator=config['target'], denominator=config['reference'])
    return {'isoform_counts': counts, 'isoform_annotations': annotation, 'structural_landscape': summary, 'donor_usage': analysis.usage, 'isoform_diversity': analysis.diversity, 'gene_DIU': analysis.differential, 'transcript_effects': analysis.effects, 'measurement_audit': analysis.diagnostics}

def run_import(config, output):
    """Portable English summary plus all results; no upstream tool execution."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    mode = config['mode']
    if mode not in ('marvel', 'longread_isoforms'):
        raise ValueError('Unsupported single-cell import mode')
    tables = import_marvel(config) if mode == 'marvel' else import_longread(config)
    for name, frame in tables.items():
        frame.to_csv(output / f'{name}.tsv', sep='\t', index=False)
    fig, ax = plt.subplots(figsize=(10, 5))
    if mode == 'marvel':
        df = tables['marvel_complete_results']
        ax.scatter(df.delta_psi, -np.log10(df.source_pvalue.clip(lower=1e-300)), s=12, c=np.where(df.passes_source_fdr_0_05, '#bb7886', '#8ba9bc'))
        ax.set(xlabel='PSI: target minus reference', ylabel='-log10(source P)', title=f"MARVEL | {config['cell_type']} | {config['target']} vs {config['reference']}")
        caution = ('Original statistics retained. Red requires declared source BH FDR<=0.05. Cell-group tests are not donor-level confirmation. Mean cell PSI is descriptive, not count-pooled pseudobulk.', '保留原版统计值；红色要求已声明的原方法BH FDR≤0.05。细胞群检验不等于动物层面验证；细胞PSI均值仅作描述，不是计数合并的pseudobulk。')
    else:
        df = tables['structural_landscape']
        labels = df.condition.astype(str) + ' | ' + df.cell_type.astype(str) + ' | ' + df.structural_category.astype(str)
        ax.barh(labels, df.molecules, color='#8ba9bc')
        ax.set(xlabel='Observed molecules', title='Long-read isoform structure by condition and cell type')
        caution = ('Processed long-read counts and annotations, not raw-read reconstruction. DIU uses Katana donor-level modeling, not the paper cell-level LRT. Gene FDR is not transcript-specific; NMD/coding annotations are predictions.', '输入为已处理的长读长计数和注释，不执行原始reads重建。DIU沿用Katana动物层面模型，不是论文细胞层面LRT。基因FDR不是独立转录本FDR；NMD和编码注释属于预测。')
    fig.tight_layout()
    image = output / 'overview.png'
    fig.savefig(image, dpi=160, bbox_inches='tight')
    fig.savefig(output / 'overview.svg', bbox_inches='tight')
    plt.close(fig)
    encoded = base64.b64encode(image.read_bytes()).decode()
    links = ''.join((f'<details><summary>{html.escape(name)}</summary><a href="{name}.tsv">Complete table</a>{frame.head(10).to_html(index=False, escape=True)}</details>' for name, frame in tables.items()))
    body = f"""<h1>{html.escape(mode)} | {html.escape(config['target'])} vs {html.escape(config['reference'])}</h1><p>{caution[0]}</p><img style="max-width:100%" src="data:image/png;base64,{encoded}">{links}"""
    (output / f"REPORT_{'EN'}.html").write_text('<!doctype html><meta charset="utf-8"><style>body{max-width:1100px;margin:30px auto;font:16px sans-serif}details{overflow:auto;margin:20px 0}td,th{padding:6px}</style>' + body, encoding='utf-8')
    provenance = {'config': config, 'inputs': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in input_files(config)}, 'source_execution': 'not run; imported results', 'validation': 'format and identity checks only; not biological validation'}
    (output / 'provenance.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
    return tables
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    run_import(read_config(args.config), args.output)
