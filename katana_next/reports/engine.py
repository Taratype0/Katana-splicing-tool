"""Portable seven-section reports from completed results and declared metadata.

No project service is opened, no source artifact is rewritten, and no scientific
workflow is implicitly executed. Report-only counts/rankings retain native P/FDR.
"""
from pathlib import Path
import hashlib
from datetime import datetime, timezone
import html, json, re, shutil, zipfile
from .sources import load_source, digest, read_table, table_schema

CHAPTERS = {
    1: ('Design, inputs and quality control', 'Declared biological samples, comparisons, references and available QC.'),
    2: ('AS landscape', 'Raw/read-support event composition, sample structure and PSI patterns.'),
    3: ('Differential molecular results', 'Gene expression, exact splice events, transcript and exon usage; original measurement units retained.'),
    4: ('Cross-comparison overlap and interaction', 'Exact memberships and effect directions; interaction tests keep their own model families.'),
    5: ('GO, KEGG and pathway enrichment', 'Comparison-specific input sets and detected/testable backgrounds; original P/FDR retained.'),
    6: ('Katana-Splice-Network', 'All available modules, composition, activity, enrichment, hubs and diagnostics.'),
    7: ('RBP and MotifTrace', 'RBP identity, matched motif/CLIP evidence and RNA sequence consequences are separate evidence layers.'),
}
ASSET_ROOT = Path(__file__).with_name('assets')


def esc(value):
    return html.escape(str(value), quote=True)


def safe(value):
    return re.sub(r'[^A-Za-z0-9._-]+', '-', str(value)).strip('-')[:150] or 'result'


def chapter_for(analysis):
    if analysis.startswith(('project.', 'quality.')):
        return 2 if analysis == 'quality.stable_as' else 1
    if analysis.startswith('landscape.') or analysis in ('single_cell.landscape', 'single_cell.cell_embedding'):
        return 2
    if analysis in ('differential.as', 'differential.deg', 'differential.dtu', 'prioritize.evidence', 'single_cell.differential', 'single_cell.isoform_usage'):
        return 3
    if analysis.startswith('differential.') or analysis == 'interpretation.joint_response':
        return 4
    if analysis in ('prioritize.enrichment', 'interpretation.joint_pathways'):
        return 5
    if analysis.startswith('programs.') or analysis == 'single_cell.programs':
        return 6
    if any(word in analysis for word in ('rbp', 'motif', 'binding')):
        return 7
    return None


class Builder:
    def __init__(self, source, output, event_set):
        from katana_next.science.das_selection import normalize_event_set
        self.source, self.output = source, Path(output).expanduser().resolve()
        self.event_set = normalize_event_set(event_set)
        self.panels, self.tables, self.audit, self._copies = [], [], [], {}
        self.top10, self.policy_data, self.derived_audits = [], {}, []
        self.exclusions = []

    def copy(self, source, *, prefix='results'):
        source = Path(source)
        if not source.is_file():
            raise FileNotFoundError(source)
        key = str(source.resolve())
        if key in self._copies:
            return self._copies[key]
        before = source.stat()
        sha = digest(source)
        destination = self.output / 'outputs' / prefix / (sha[:12] + '-' + source.name)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        if source.stat().st_size != before.st_size or source.stat().st_mtime_ns != before.st_mtime_ns or digest(destination) != sha:
            raise RuntimeError('Source changed during report collection: ' + str(source))
        relative = destination.relative_to(self.output).as_posix()
        self.audit.append({'source_name': source.name, 'source_sha256': sha, 'bytes': before.st_size, 'report_path': relative})
        self._copies[key] = relative
        return relative

    def table(self, artifact):
        path = Path(artifact['path'])
        record = {'name': artifact.get('name', path.stem), 'analysis_id': artifact.get('analysis_id', ''),
                  'artifact_id': artifact.get('artifact_id'), 'source': self.copy(path), **table_schema(path)}
        record.update({k: artifact[k] for k in ('chapter', 'contrast_id', 'policy', 'unit') if k in artifact})
        # Small complete tables receive a self-contained browser. Large result
        # tables retain complete Parquet/TSV access, never a false complete Top N.
        if record['rows'] <= 50000:
            frame = read_table(path)
            record['view'] = self.table_browser(frame, record)
        self.tables.append(record)
        return record

    def table_browser(self, frame, record):
        # All rows of this table, with original scalar values; no new statistics.
        values = frame.astype(object).where(frame.notna(), None).values.tolist()
        data = json.dumps({'columns': list(frame.columns), 'rows': values}, default=str, ensure_ascii=False, allow_nan=False).replace('<', '\\u003c')
        identity = str(record.get('artifact_id') or record['name'])
        context = str(record.get('contrast_id', '')) + '|' + str(record.get('policy', '')) + '|' + record['source']
        suffix = hashlib.sha256(context.encode('utf-8')).hexdigest()[:16]
        filename = safe(identity)[:110] + '-' + suffix + '.html'
        view = self.output / 'outputs/views' / filename
        view.parent.mkdir(parents=True, exist_ok=True)
        raw = Path('..') / Path(record['source']).relative_to('outputs')
        content = (ASSET_ROOT / 'table-template.html').read_text('utf-8').replace('__TITLE__', esc(record['name'])).replace('__RAW__', esc(raw.as_posix())).replace('__DATA__', data)
        view.write_text(content, encoding='utf-8')
        return view.relative_to(self.output).as_posix()

    def figure(self, artifact, chapter):
        path = Path(artifact['path'])
        exports = artifact.get('metadata', {}).get('exports', {})
        image = Path(exports.get('png', ''))
        svg = Path(exports.get('svg', ''))
        if path.suffix.lower() in ('.png', '.svg', '.jpg', '.jpeg'):
            image = path
        if not image.is_file() and path.name.endswith('.vl.json'):
            import vl_convert as vlc
            spec = json.loads(path.read_text('utf-8'))
            target = self.output / 'figures' / (safe(artifact.get('artifact_id', path.stem)) + '.svg')
            target.parent.mkdir(parents=True, exist_ok=True)
            # Rendering a saved specification does not fit a statistical model.
            target.write_text(vlc.vegalite_to_svg(spec), encoding='utf-8')
            image = target
        if not image.is_file():
            return None
        relative = image.relative_to(self.output).as_posix() if image.is_relative_to(self.output) else self.copy(image, prefix='figures')
        title = artifact.get('name', image.stem)
        output_only = any(word in title.lower() for word in ('pairwise', 'candidate rank', 'similarity', 'correlation matrix'))
        # Completed historical selection results remain accessible in outputs.
        # New report RS flags must never relabel their earlier input population.
        if not artifact.get('metadata', {}).get('report_matched_selection') and artifact.get('analysis_id', '') in ('landscape.bulk', 'project.audit', 'quality.stable_as', 'differential.as', 'differential.deg', 'differential.dtu', 'prioritize.evidence', 'differential.cross_comparison', 'prioritize.enrichment'):
            output_only = True
        source_run = artifact.get('run_fingerprint', '')
        if path.is_file() and path.suffix.lower() == '.json': self.copy(path, prefix='figure-specifications')
        record = {'id': safe(artifact.get('artifact_id', image.stem)), 'chapter': chapter, 'title': title,
                  'image': relative, 'source_artifact_id': artifact.get('artifact_id'),
                  'source_run': source_run, 'contrast_id': (artifact.get('contrast_ids') or ['all'])[0] if len(artifact.get('contrast_ids', [])) <= 1 else 'all',
                  'cell_type': artifact.get('metadata', {}).get('cell_type', 'all'),
                  'caption': artifact.get('metadata', {}).get('method', 'Saved completed result; source statistics retained.'),
                  'output_only': output_only or artifact.get('output_only', False) or artifact.get('metadata', {}).get('display_in_main') is False, 'metadata': artifact.get('metadata', {}),
                  'policy': artifact.get('policy') or artifact.get('metadata', {}).get('event_set')}
        if svg.is_file():
            record['svg'] = self.copy(svg, prefix='figures')
        self.panels.append(record)
        return record

    def collect(self):
        for artifact in [*self.source['artifacts'], *self.source.get('extra_assets', [])]:
            path = Path(artifact['path'])
            if not path.is_file():
                self.exclusions.append({'name': artifact.get('name'), 'reason': 'Source file missing; not computed/imported into this report.'})
                continue
            chapter = chapter_for(artifact.get('analysis_id', ''))
            if not chapter:
                continue
            if artifact.get('kind') == 'table' or path.suffix.lower() in ('.parquet', '.tsv', '.csv'):
                self.table(artifact)
            elif path.suffix.lower() == '.html' and chapter not in (6, 7):
                self.tables.append({'name': artifact.get('name', path.stem), 'source': self.copy(path, prefix='native-browsers'), 'analysis_id': artifact.get('analysis_id'), 'chapter': chapter})
            elif chapter in (6, 7) and path.suffix.lower() not in ('.html', '.png', '.svg', '.jpg', '.jpeg') and not path.name.endswith('.vl.json'):
                self.tables.append({'name': artifact.get('name', path.stem), 'source': self.copy(path, prefix='native-attachments'), 'analysis_id': artifact.get('analysis_id'), 'chapter': chapter})
            elif chapter not in (6, 7):
                self.figure(artifact, chapter)
        self.collect_rich({**self.source, 'artifacts': [*self.source['artifacts'], *self.source.get('extra_assets', [])]})
        for supplemental in self.source.get('supplemental', []):
            for artifact in supplemental['artifacts']:
                p = Path(artifact['path'])
                if p.is_file() and artifact.get('kind') == 'table':
                    self.table({**artifact, 'chapter': chapter_for(artifact.get('analysis_id', ''))})
            self.collect_rich(supplemental)
        from .landscape import build_landscape_views
        self.add_views(build_landscape_views(self.source['project'], self.source['artifacts'], self.output/'outputs/landscape'))
        from .comparisons import build_comparison_views
        self.add_views(build_comparison_views(self.source['project'], self.source['artifacts'], self.output, event_set=self.event_set))
        for views in self.source.get('extra_views', []):
            if isinstance(views, str):
                views = json.loads(Path(views).read_text('utf-8-sig'))
            self.add_views(views)

    def collect_rich(self, source):
        from .rich_results import collect_rich_results
        rich = collect_rich_results(source['artifacts'], search_roots=[Path(p) for p in source.get('search_roots', [])])
        for asset in rich:
            if asset.get('status') == 'not_computed':
                self.exclusions.append(asset)
                continue
            path = Path(asset['path'])
            if asset['kind'] in ('image', 'vega'):
                raw = {**asset, 'artifact_id': asset['id'], 'name': asset['title'], 'analysis_id': asset.get('metadata', {}).get('analysis_id', ''), 'kind': 'plot_spec', 'contrast_ids': [] if asset.get('contrast_id', 'all') == 'all' else [asset['contrast_id']]}
                panel = self.figure(raw, asset['chapter'])
                if panel:
                    panel['caption'] = asset['caption']
                    if not asset.get('metadata', {}).get('display_in_main', True): panel['output_only'] = True
                    scope = asset.get('metadata', {}).get('das_scope', {})
                    panel['selection_scope'] = scope
                    if scope.get('mode') == 'switchable_recorded_result':
                        selected = scope.get('recorded_event_set')
                        supported = selected in ('significant', 'stable_significant')
                        if selected == 'stable_significant':
                            from katana_next.science.bulk import READ_SUPPORT_POLICY_VERSION
                            supported = supported and scope.get('read_support_policy') == READ_SUPPORT_POLICY_VERSION
                        if supported:
                            panel['policy'] = selected
                        else:
                            panel['output_only'] = True
                            panel['caption'] += ' Historical or unrecorded foreground; retained in outputs and not substituted for the current DAS/RS-DAS selection.'
            elif asset['kind'] == 'html':
                from .rich_results import copy_rich_asset
                receipt = copy_rich_asset(asset, self.output / 'outputs/rich' / safe(asset['id']))
                target = Path(receipt['path'])
                self.tables.append({'name': asset['title'], 'source': target.relative_to(self.output).as_posix(), 'rows': None, 'analysis_id': 'rich_browser', 'chapter': asset['chapter']})
                self.derived_audits.append(receipt)

    def add_views(self, results):
        for raw in results.get('panels', []):
            panel = dict(raw)
            path = Path(panel.pop('path', panel.get('image', '')))
            if not path.is_absolute(): path = self.output/path
            panel['image'] = path.relative_to(self.output).as_posix() if path.is_relative_to(self.output) else self.copy(path, prefix='figures')
            panel.setdefault('id', safe(path.stem))
            if not panel.get('contrast_id'): panel['contrast_id'] = 'all'
            panel.setdefault('cell_type', 'all'); panel.setdefault('caption', 'Completed results; native effects and P/FDR retained.')
            svg = path.with_suffix('.svg')
            if svg.is_file(): panel['svg'] = svg.relative_to(self.output).as_posix() if svg.is_relative_to(self.output) else self.copy(svg, prefix='figures')
            self.panels.append(panel)
        for raw in results.get('tables', []):
            record = dict(raw)
            path = Path(record.pop('path', record.get('source', '')))
            if not path.is_absolute(): path = self.output/path
            record.setdefault('chapter', 3)
            if path.suffix.lower() == '.html' or record.get('kind') == 'file' or path.suffix.lower() not in ('.parquet', '.csv', '.tsv', '.txt', '.gz'):
                record['source'] = path.relative_to(self.output).as_posix() if path.is_relative_to(self.output) else self.copy(path)
            else:
                record.update(self.table({**record, 'path': str(path)}))
                # table() already appended the same record.
                continue
            self.tables.append(record)
        for raw in results.get('top10', []):
            record = dict(raw)
            for key in ('full_path', 'full_table'):
                if record.get(key):
                    complete = Path(record[key])
                    if complete.is_absolute():
                        record[key] = complete.relative_to(self.output).as_posix() if complete.is_relative_to(self.output) else self.copy(complete)
            self.top10.append(record)
        self.policy_data.update(results.get('policy_data', {}))
        if results.get('gene_sets'):
            (self.output/'gene-set-inputs.json').write_text(json.dumps(results['gene_sets'], indent=2), encoding='utf-8')
        self.exclusions.extend(results.get('unavailable', []))
        if results.get('audit'): self.derived_audits.append(results['audit'])

    def write(self):
        from .layout import render
        return render(self)


def build_report(project, output_dir, *, catalog=None, event_set='stable_significant'):
    output = Path(output_dir).expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError('Choose a new report output directory; completed reports are preserved.')
    source=load_source(project,catalog)
    # Never permit report outputs inside an original source result tree.
    roots=[Path(source['project']['output_root']).resolve()] if source['project'].get('output_root') else []
    if any(output==root or output.is_relative_to(root) for root in roots):
        raise ValueError('Report output must be separate from the source result directory.')
    output.mkdir(parents=True,exist_ok=True)
    shutil.copytree(ASSET_ROOT,output/'assets',dirs_exist_ok=True)
    builder=Builder(source,output,event_set)
    builder.collect()
    return builder.write()
