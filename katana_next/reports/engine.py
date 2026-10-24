"""Portable seven-section reports from completed results and declared metadata.

No project service is opened, no source artifact is rewritten, and no scientific
workflow is implicitly executed. Report-only counts/rankings retain native P/FDR.
"""
from pathlib import Path
import hashlib
from datetime import datetime, timezone
import copy
import math
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
        self.enrichment_summaries = []
        self.motif_summaries = []
        self.network_browsers = []
        self.create_archive = True
        self.result_store = None
        self.compact_tables = False
        self.exclusions = []
        self.data_link_audit = []
        self._extra_asset_ids = {str(a.get('artifact_id')) for a in source.get('extra_assets', []) if a.get('artifact_id')}

    def _comparison_scope(self, panel):
        selected = {str(c['contrast_id']) for c in self.source['project'].get('contrasts', [])}
        contrast = str(panel.get('contrast_id') or 'all')
        if selected and contrast != 'all' and contrast not in selected:
            panel['output_only'] = True
        return panel

    def bind_data_links(self, panel, links):
        """Validate explicit producer bindings, then collect immutable files.

        Associations are declared by the producer, never guessed from a name.
        Exact plotted tables additionally check recorded comparison/policy
        columns when present. Full masters may be explicitly labelled as a
        full_collection, but they cannot stand in for a plot_data binding.
        """
        result = []
        roles = {'plot_data', 'full_results', 'foreground', 'background', 'methods'}
        contrast = str(panel.get('contrast_id') or 'all')
        display_policy = panel.get('policy') or None
        recorded_scope = panel.get('metadata', {}).get('das_scope', {})
        policy = display_policy
        if recorded_scope.get('mode') in {'rs_required', 'fixed_completed_result'} and recorded_scope.get('recorded_event_set') in {'significant', 'stable_significant'}:
            policy = recorded_scope['recorded_event_set']
        def group_ids(values):
            if not isinstance(values, (list, tuple)) or any(not isinstance(v, str) or not v for v in values) or len(values) != len(set(values)):
                raise ValueError('Figure comparison group requires unique declared comparison IDs.')
            return list(values)
        declared_ids = group_ids(panel.get('comparison_ids') or [])
        project = self.source.get('native_project', self.source['project']) if panel.get('output_only') else self.source['project']
        known_ids = {str(c['contrast_id']) for c in project.get('contrasts', [])}
        if known_ids and not set(declared_ids) <= known_ids:
            raise ValueError('Figure comparison group contains an unknown or unselected comparison.')
        comparisons = set(declared_ids)
        if not comparisons and contrast == 'all':
            comparisons = known_ids
        for original in links or []:
            link = copy.deepcopy(original)
            missing = {'role', 'path', 'label', 'contrast_id', 'policy', 'unit'} - link.keys()
            if missing or link.get('role') not in roles or not str(link.get('label', '')).strip() or not str(link.get('unit', '')).strip():
                raise ValueError('Figure data link requires a declared role, path, label, comparison, policy and unit: ' + str(sorted(missing)))
            role = link['role']
            link_contrast = str(link['contrast_id'] or 'all')
            link_policy = None if link['policy'] in (None, '', 'common') else link['policy']
            full_collection = role == 'full_results' and link.get('scope') == 'full_collection'
            shared_method = role == 'methods' and link_contrast == 'all' and link_policy is None
            if not shared_method and not full_collection:
                valid_contrast = link_contrast == contrast or (contrast == 'all' and role in {'foreground', 'background', 'full_results'} and link_contrast in comparisons)
                if not valid_contrast:
                    raise ValueError('Figure data link belongs to a different comparison: ' + link_contrast)
                if link_policy != policy:
                    raise ValueError('Figure data link belongs to a different DAS policy.')
            if link.get('comparison_ids') is not None:
                linked_ids = group_ids(link['comparison_ids'])
                if known_ids and not set(linked_ids) <= known_ids:
                    raise ValueError('Figure data link contains an unknown or unselected comparison.')
                if declared_ids and (set(linked_ids) != set(declared_ids) or (role == 'plot_data' and linked_ids != declared_ids)):
                    raise ValueError('Figure data link belongs to a different comparison group or plot axis order.')
            if link.get('project_id') and link['project_id'] != self.source['project'].get('project_id'):
                raise ValueError('Figure data link belongs to a different project.')
            if '://' in str(link['path']):
                raise ValueError('Figure data links require locally available files; remote URIs are not fetched.')
            path = Path(link['path']).expanduser()
            if not path.is_absolute():
                path = self.output / path
                if not path.resolve().is_relative_to(self.output):
                    raise ValueError('Relative figure data path escapes the report; declare external source paths explicitly.')
            path = path.resolve()
            if not path.is_file():
                raise FileNotFoundError('Declared figure data is unavailable: ' + str(path))
            before = path.stat()
            if not before.st_size:
                raise ValueError('Declared figure data file is empty: ' + str(path))
            checks = {}
            if role == 'plot_data' and path.suffix.lower() in ('.parquet', '.tsv', '.csv', '.txt'):
                if path.suffix.lower() == '.parquet':
                    import pyarrow.parquet as pq
                    columns = pq.ParquetFile(path).schema_arrow.names
                else:
                    import pandas as pd
                    columns = list(pd.read_csv(path, sep='\t' if path.suffix.lower() in ('.tsv', '.txt') else ',', nrows=0))
                comparison_column = next((c for c in ('comparison_id', 'contrast_id') if c in columns), None)
                policy_column = next((c for c in ('report_policy', 'report_event_set') if c in columns), None)
                requested = [c for c in (comparison_column, policy_column) if c]
                scope = read_table(path, columns=requested) if requested else None
                if comparison_column:
                    observed = set(scope[comparison_column].dropna().astype(str))
                    if contrast != 'all' and ('all' in observed or (len(scope) and not observed)):
                        raise ValueError('Plot data does not establish the exact single comparison.')
                    found = observed - {'all'}
                    expected = {contrast} if contrast != 'all' else comparisons
                    if expected and not found <= expected:
                        raise ValueError('Plot data contains a foreign comparison despite its declared binding.')
                    checks['comparison_values'] = sorted(found)
                if policy_column:
                    found = {None if p in ('common', '') else p for p in scope[policy_column].dropna().astype(str)}
                    if (found and found != {policy}) or (len(scope) and policy is not None and not found):
                        raise ValueError('Plot data contains a different DAS policy despite its declared binding.')
                    checks['policy_values'] = sorted(found, key=str)
            sha = digest(path)
            if link.get('source_sha256') and link['source_sha256'] != sha:
                raise ValueError('Figure data source SHA does not match its producer binding.')
            if path.is_relative_to(self.output):
                relative = path.relative_to(self.output).as_posix()
            elif path.suffix.lower() == '.html':
                root = Path(link.get('portable_root', path.parent)).resolve()
                if not path.is_relative_to(root):
                    raise ValueError('Figure data browser must remain inside its declared portable root.')
                origin = str(root) + '\0' + path.relative_to(root).as_posix()
                origin_sha = hashlib.sha256(origin.encode('utf-8')).hexdigest()[:16]
                # Identical HTML can refer to different comparison-local files.
                # Its content SHA alone cannot identify a portable bundle.
                cache_key = 'data-browser:' + sha + ':' + origin
                relative = self._copies.get(cache_key)
                if relative is None:
                    from .rich_results import _dependencies, copy_rich_asset
                    deps = _dependencies(path, root, {})
                    if deps['portability_warnings'] or deps['external_urls']:
                        raise ValueError('Figure data browser has unavailable or external dependencies; provide a complete offline bundle.')
                    destination = self.output / 'outputs/figure-data-browsers' / (sha[:12] + '-' + origin_sha + '-' + safe(path.stem))
                    receipt = copy_rich_asset({'id': 'figure-data-' + sha[:12] + '-' + origin_sha, 'status': 'available', 'path': str(path),
                                               'metadata': {'portable_root': str(root), 'relative_path': path.relative_to(root).as_posix(), **deps}}, destination)
                    relative = Path(receipt['path']).relative_to(self.output).as_posix()
                    self._copies[cache_key] = relative
                    self.derived_audits.append(receipt)
            else:
                relative = self.copy(path, prefix='figure-data')
            if (path.stat().st_size, path.stat().st_mtime_ns) != (before.st_size, before.st_mtime_ns) or digest(path) != sha:
                raise RuntimeError('Figure data changed during collection: ' + str(path))
            portable_sha = digest(self.output / relative)
            link.update(path=relative, contrast_id=link_contrast, policy=link_policy, source_sha256=sha, portable_sha256=portable_sha, bytes=before.st_size)
            result.append(link)
            self.data_link_audit.append({'panel_id': panel.get('id'), 'role': role, 'path': relative,
                                        'source_sha256': sha, 'portable_sha256': portable_sha, 'bytes': before.st_size, 'contrast_id': link_contrast,
                                        'policy': link_policy, 'display_policy': display_policy, 'measurement_policy': policy,
                                        'unit': link['unit'], 'scope_checks': checks,
                                        'original_preserved': True})
        return result

    def figure_data_coverage(self):
        main = [p for p in self.panels if not p.get('output_only')]
        missing = [{'id': p['id'], 'title': p['title']} for p in main if not any(d['role'] == 'plot_data' for d in p.get('data_links', []))]
        return {'schema': 'katana.figure_data_links/v1', 'main_figures': len(main),
                'main_figures_with_plot_data': len(main) - len(missing), 'unbound_main_figures': missing,
                'strict_plot_data_required': bool(self.source.get('report_require_plot_data', False)),
                'bindings': self.data_link_audit}

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
        from .storage import copy_payload
        copy_payload(source, destination, store=self.result_store, expected_sha=sha)
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
        # JSON has no signed infinity. Preserve it explicitly in the browser;
        # missing values stay null and the copied raw result remains unchanged.
        values = [[('Infinity' if value > 0 else '-Infinity')
                   if isinstance(value, (int, float)) and math.isinf(value) else value
                   for value in row] for row in values]
        data = json.dumps({'columns': list(frame.columns), 'rows': values}, default=str, ensure_ascii=False, allow_nan=False).replace('<', '\\u003c')
        identity = str(record.get('artifact_id') or record['name'])
        context = str(record.get('contrast_id', '')) + '|' + str(record.get('policy', '')) + '|' + record['source']
        suffix = hashlib.sha256(context.encode('utf-8')).hexdigest()[:16]
        filename = safe(identity)[:110] + '-' + suffix + '.html'
        view = self.output / 'outputs/views' / filename
        view.parent.mkdir(parents=True, exist_ok=True)
        raw = Path('..') / Path(record['source']).relative_to('outputs')
        content = (ASSET_ROOT / 'table-template.html').read_text('utf-8').replace('__TITLE__', esc(record['name'])).replace('__RAW__', esc(raw.as_posix())).replace('__DATA__', data)
        if self.compact_tables:
            from .compact import _compress_table
            content, _ = _compress_table(content, 32768)
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
        source_run = artifact.get('run_fingerprint') or artifact.get('metadata', {}).get('run_fingerprint', '')
        if path.is_file() and path.suffix.lower() == '.json': self.copy(path, prefix='figure-specifications')
        record = {'id': safe(artifact.get('artifact_id', image.stem)), 'chapter': chapter, 'title': title,
                  'image': relative, 'source_artifact_id': artifact.get('artifact_id'),
                  'source_run': source_run, 'contrast_id': (artifact.get('contrast_ids') or ['all'])[0] if len(artifact.get('contrast_ids', [])) <= 1 else 'all',
                  'cell_type': artifact.get('metadata', {}).get('cell_type', 'all'),
                  'caption': artifact.get('metadata', {}).get('method', 'Saved completed result; source statistics retained.'),
                  'output_only': output_only or artifact.get('output_only', False) or artifact.get('metadata', {}).get('display_in_main') is False, 'metadata': artifact.get('metadata', {}),
                  'policy': artifact.get('policy') or artifact.get('metadata', {}).get('event_set')}
        record['report_origin'] = artifact.get('report_origin', 'native')
        if record['metadata'].get('das_scope', {}).get('mode') in {'rs_required', 'fixed_completed_result'}:
            record['policy'] = None
        if record['report_origin'] == 'native' and not self.source.get('report_native_figures_in_main', True):
            record['output_only'] = True
        self._comparison_scope(record)
        links = list(artifact.get('data_links', artifact.get('metadata', {}).get('data_links', [])) or [])
        from .bootstrap_plot_data import bootstrap_plot_links
        links = bootstrap_plot_links(artifact, record, self.output, links)
        if not any(link.get('role') == 'plot_data' for link in links):
            from .saved_plot_data import recorded_plot_data_links
            links.extend(recorded_plot_data_links(artifact, record, self.output))
        record['data_links'] = self.bind_data_links(record, links)
        if svg.is_file():
            record['svg'] = self.copy(svg, prefix='figures')
        self.panels.append(record)
        return record

    def collect(self):
        collections = [('native', self.source['artifacts']), ('extra', self.source.get('extra_assets', []))]
        for origin, artifact in [(origin, {**raw, 'report_origin': origin}) for origin, records in collections for raw in records]:
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
        from .rbp_enrichment import build_rbp_summaries, build_motif_hit_summaries
        all_artifacts = [*self.source['artifacts'], *self.source.get('extra_assets', [])]
        for supplemental in self.source.get('supplemental', []):
            all_artifacts.extend(supplemental['artifacts'])
        self.enrichment_summaries = build_rbp_summaries(all_artifacts) + build_motif_hit_summaries(all_artifacts, self.output)
        replaced = {(r['source_run'], r['contrast_id']) for r in self.enrichment_summaries}
        for record in self.enrichment_summaries:
            record['full_path'] = self.copy(Path(record['full_path']))
            if record.get('all_hits_path'):
                record['all_hits_path'] = Path(record['all_hits_path']).relative_to(self.output).as_posix()
        for panel in self.panels:
            if panel.get('metadata', {}).get('analysis_id') in {'prioritize.rbp_support', 'prioritize.event_motifs', 'prioritize.rbp_mechanism'} and (panel.get('source_run', ''), panel.get('contrast_id')) in replaced:
                panel['output_only'] = True
                panel['metadata']['main_summary'] = 'Significant RBP enrichment with complete recorded hit genes'

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
                raw['report_origin'] = 'extra' if str(asset.get('metadata', {}).get('source_artifact_id')) in self._extra_asset_ids else 'native'
                intended_scope = asset.get('metadata', {}).get('das_scope', {})
                if intended_scope.get('mode') == 'switchable_recorded_result' and intended_scope.get('recorded_event_set') in ('significant', 'stable_significant'):
                    raw['policy'] = intended_scope['recorded_event_set']
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
                self.tables.append({'name': asset['title'], 'source': target.relative_to(self.output).as_posix(),
                    'rows': None, 'analysis_id': 'rich_browser', 'chapter': asset['chapter'],
                    'contrast_id': asset.get('contrast_id', 'all'),
                    'source_artifact_id': asset.get('metadata', {}).get('source_artifact_id'),
                    'metadata': {**asset.get('metadata', {}), 'source_sha256': digest(path),
                                 'portable_sha256': digest(target)}})
                self.derived_audits.append(receipt)

    def add_views(self, results):
        for raw in results.get('panels', []):
            panel = dict(raw)
            from .mapp_catalogs import prepare_portable_mapp_panel
            panel, mapp_audit = prepare_portable_mapp_panel(panel, report_root=self.output, copy_file=self.copy)
            if mapp_audit: self.derived_audits.append(mapp_audit)
            path = Path(panel.pop('path', panel.get('image', '')))
            if not path.is_absolute(): path = self.output/path
            panel['image'] = path.relative_to(self.output).as_posix() if path.is_relative_to(self.output) else self.copy(path, prefix='figures')
            panel.setdefault('id', safe(path.stem))
            if not panel.get('contrast_id'): panel['contrast_id'] = 'all'
            panel.setdefault('cell_type', 'all'); panel.setdefault('caption', 'Completed results; native effects and P/FDR retained.')
            self._comparison_scope(panel)
            panel['data_links'] = self.bind_data_links(panel, panel.get('data_links', []))
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
        for raw in results.get('motif_summaries', []):
            summary = copy.deepcopy(raw)
            if summary.get('ranking_path'):
                ranking = Path(summary['ranking_path'])
                if not ranking.is_absolute(): ranking = self.output / ranking
                summary['ranking_path'] = (ranking.relative_to(self.output).as_posix()
                    if ranking.is_relative_to(self.output) else self.copy(ranking))
            self.motif_summaries.append(summary)
        if results.get('gene_sets'):
            (self.output/'gene-set-inputs.json').write_text(json.dumps(results['gene_sets'], indent=2), encoding='utf-8')
        self.exclusions.extend(results.get('unavailable', []))
        if results.get('audit'): self.derived_audits.append(results['audit'])

    def write(self):
        from .layout import render
        from .presentation import prepare_presentation
        from .motif_display import apply_motif_significance_display
        from .mapp_catalogs import validate_portable_mapp_panel
        for panel in self.panels:
            if panel.get('metadata', {}).get('component') == 'mapp_native':
                self.derived_audits.append(validate_portable_mapp_panel(panel, report_root=self.output))
        current, motif_audit = apply_motif_significance_display(
            {'panels':self.panels,'complete_results':self.tables},self.output)
        self.panels = current['panels']
        if current['motif_summaries']:
            self.motif_summaries = current['motif_summaries']
            self.derived_audits.append(motif_audit)
        prepare_presentation(self)
        coverage = self.figure_data_coverage()
        if coverage['strict_plot_data_required'] and coverage['unbound_main_figures']:
            raise ValueError('Every main figure requires exact plot_data; missing bindings: ' + ', '.join(p['id'] for p in coverage['unbound_main_figures']))
        return render(self)


def build_report(project, output_dir, *, catalog=None, event_set='stable_significant',
                 compact=False, create_archive=True, archive_path=None, result_store=None,
                 path_mappings=None):
    output = Path(output_dir).expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError('Choose a new report output directory; completed reports are preserved.')
    source=load_source(project,catalog,path_mappings=path_mappings)
    # Never permit report outputs inside an original source result tree.
    roots=[Path(source['project']['output_root']).resolve()] if source['project'].get('output_root') else []
    if any(output==root or output.is_relative_to(root) for root in roots):
        raise ValueError('Report output must be separate from the source result directory.')
    output.mkdir(parents=True,exist_ok=True)
    shutil.copytree(ASSET_ROOT,output/'assets',dirs_exist_ok=True)
    builder=Builder(source,output,event_set)
    builder.create_archive=False
    builder.compact_tables=compact
    if result_store:
        store=Path(result_store).expanduser().resolve()
        if store==output or store.is_relative_to(output) or output.is_relative_to(store) or any(store==root or store.is_relative_to(root) for root in roots):
            raise ValueError('Result store must be separate from source results and report output')
        builder.result_store=store
    if source.get('assembly_audit'):
        builder.derived_audits.append(source['assembly_audit'])
        for path in source['source_files']: builder.copy(path,prefix='provenance')
    builder.collect()
    result=builder.write()
    if compact:
        from .compact import compact_report
        compact_report(output)
        result=json.loads((output/'report-manifest.json').read_text('utf-8'))
    if create_archive:
        from .refresh import write_archive
        write_archive(output,archive_path=archive_path)
    return result
