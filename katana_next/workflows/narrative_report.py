from pathlib import Path
import json
import zipfile
from katana_next.workflows.base import WorkflowBase, WorkflowDefinition, ParameterDefinition, ArtifactWriter
from katana_next.science.narrative_report import build_report, package_report, is_delivery_file
from katana_next.science.report_language import copy_english_reader, is_retired_report

class NarrativeReportWorkflow(WorkflowBase):
    version = '1.3.0+english.20261021'
    definition = WorkflowDefinition(workflow_id='reports.interpreted_bundle', workspace_id='reports', label='English interpreted study report', question='What did this experiment show, how strong is the evidence, and what remains unresolved?', method='Question-led synthesis of completed result delivery folders; complete non-significant attachments and SHA256 inventory.', output_summary='English reports, source manifest and portable ZIP; no hidden upstream analysis.', analysis_unit='experiment and its declared contrasts', parameters=(ParameterDefinition('config_path', 'Delivery configuration JSON', 'text', '', 'Keys: study, result_root, manifest_path. Requires standardized numbered delivery directories.'),), limitations='Requires completed delivery tables. This synthesis is not causal validation or a claim of journal suitability.')

    def _config(self, parameters):
        path = Path(parameters.get('config_path', ''))
        if not path.is_file():
            raise ValueError('Provide report configuration')
        c = json.loads(path.read_text())
        if c.get('approved_reader_root'):
            source = Path(c['approved_reader_root'])
            for name in ('REPORT_EN.html', 'tables', 'figures'):
                if not (source / name).exists():
                    raise ValueError(f'Approved report input missing: {name}')
            return c
        root = Path(c['result_root'])
        if not (root / '20_Presentation_report/comparison_totals.tsv').is_file():
            raise ValueError('Standardized delivery summary missing')
        if not list((root / '23_Module_explanation').rglob('module_cards.tsv')):
            raise ValueError('Run module explanation first')
        return c

    def input_paths(self, context, parameters):
        c = self._config(parameters)
        if c.get('approved_reader_root'):
            return [Path(parameters['config_path']), Path(c['manifest_path']), Path(__file__),
                    Path(__file__).parents[1] / 'science' / 'report_language.py',
                    *[p for p in Path(c['approved_reader_root']).rglob('*') if p.is_file() and not is_retired_report(p)]]
        root = Path(c['result_root'])
        return [Path(parameters['config_path']), Path(c['manifest_path']), *[p for p in root.rglob('*') if is_delivery_file(p, root) and p.suffix in ('.tsv', '.json', '.png') and ('24_Narrative_report' not in p.parts)]]

    def run(self, context, job, run_fingerprint, parameters):
        c = self._config(parameters)
        if json.loads(Path(c['manifest_path']).read_text())['project_id'] != context.manifest.project_id:
            raise ValueError('Report project identity mismatch')
        if context.manifest.metadata.get('gene_filter'):
            source = Path(c.get('approved_reader_root') or c['result_root']).resolve()
            if not source.is_relative_to(context.manifest.output_root.resolve()):
                raise ValueError('Filtered reports require result tables and figures rebuilt in this filtered project scope; an unfiltered delivery cannot be relabeled as filtered.')
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version)
        if c.get('approved_reader_root'):
            from bs4 import BeautifulSoup
            source = Path(c['approved_reader_root']).resolve()
            out = writer.directory / 'reader_report'
            path = copy_english_reader(source, out)
            soup = BeautifulSoup(path.read_text(encoding='utf-8'), 'html.parser')
            if not soup.find('section') or not soup.find('img'):
                raise ValueError('Reader report lacks sections or figures')
            for a in soup.find_all('a', href=True):
                href = a['href']
                if href.startswith('#'):
                    if not soup.find(id=href[1:]):
                        raise ValueError(f'Broken report anchor: {href}')
                elif not (out / href).is_file():
                    raise ValueError(f'Broken report attachment: {href}')
            if context.manifest.metadata.get('gene_filter'):
                from katana_next.science.gene_filter import gene_filter_footnote
                footer = soup.new_tag('footer', id='gene-exclusion-footnote')
                footer.string = gene_filter_footnote(context.manifest)
                (soup.body or soup).append(footer)
                path.write_text(str(soup), encoding='utf-8')
            writer.existing_file(f"Complete reader report {'EN'}", path, entity_type='report')
            archive = writer.directory / (c['study'] + '_English_complete_results.zip')
            with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED, compresslevel=1) as z:
                for p in out.rglob('*'):
                    if p.is_file():
                        z.write(p, Path(c['study']) / p.relative_to(out))
            with zipfile.ZipFile(archive) as z:
                if z.testzip() is not None:
                    raise ValueError('Report ZIP failed CRC check')
            writer.existing_file('Portable English results', archive, entity_type='report_package')
            return writer.finish()
        root = Path(c['result_root'])
        out = root / '24_Narrative_report'
        build_report(c['study'], root, c['manifest_path'], out)
        for name in ('REPORT_EN.html', 'SOURCE_MANIFEST_SHA256.tsv', 'report_status.json'):
            writer.existing_file(name, out / name, entity_type='report')
        archive = writer.directory / (c['study'] + '_English_complete_results.zip')
        package_report(c['study'], root, out, archive)
        writer.existing_file('Portable English results', archive, entity_type='report_package')
        return writer.finish()
