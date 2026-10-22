from pathlib import Path
import hashlib
import json
import math
import zipfile
from katana_next.workflows.base import WorkflowBase, WorkflowDefinition, ParameterDefinition, ArtifactWriter
from katana_next.science.narrative_report import build_report, package_report, is_delivery_file
from katana_next.science.report_language import copy_english_reader, is_retired_report


def _recorded_artifact_file(context, artifact, source, *, allow_export_copy=False):
    """Match an unchanged, completed artifact or its equivalent TSV export."""
    from katana_next.infrastructure.fingerprints import fingerprint_file, fingerprint_value

    if (artifact.project_id != context.manifest.project_id
            or str(artifact.status) != 'complete'
            or artifact.metadata.get('artifact_contract_version') != '3.0-correctness'
            or not artifact.path.is_file()):
        return False
    try:
        source = Path(source).resolve()
        exports = artifact.metadata.get('exports', {})
        if not isinstance(exports, dict):
            return False
        exported = exports.get('tsv')
        export_source = Path(exported).resolve() if exported else None
        if source != artifact.path.resolve() and source != export_source:
            if not (allow_export_copy and export_source and export_source.is_file() and source.is_file()
                    and hashlib.sha256(source.read_bytes()).digest()
                    == hashlib.sha256(export_source.read_bytes()).digest()):
                return False
        marker = json.loads((artifact.path.parent / 'parameters.json').read_text(encoding='utf-8'))
        expected = fingerprint_value({'run': artifact.run_fingerprint, 'name': artifact.name,
            'path': artifact.path.name,
            'content': fingerprint_file(artifact.path, force_full=True).model_dump(mode='json')})
        if (not marker.get('completed') or expected != artifact.fingerprint
                or marker.get('workflow_id') != artifact.analysis_id
                or marker.get('run_fingerprint') != artifact.run_fingerprint
                or not any(a.get('artifact_id') == artifact.artifact_id
                           and a.get('fingerprint') == artifact.fingerprint
                           for a in marker.get('generated_artifacts', []))):
            return False
        if source != artifact.path.resolve():
            import pandas as pd
            expected_export = pd.read_parquet(artifact.path).to_csv(sep='\t', index=False).encode('utf-8')
            if source.read_bytes() != expected_export:
                return False
    except (OSError, ValueError, TypeError, KeyError, AssertionError):
        return False
    return True


def _network_module_parameter_record(context, config):
    """Follow the exact module-card inputs; never select a latest network run."""
    missing = {'status': 'not_recorded', 'reason':
        'Actual network module parameters are not recorded for the network results used by this report. '
        'Current defaults are not substituted for historical settings.'}
    # An approved reader is an immutable delivery. Its underlying network input
    # chain is not represented by this configuration, even if another run exists.
    if config.get('approved_reader_root'):
        return missing
    cards = list((Path(config['result_root']) / '23_Module_explanation').rglob('module_cards.tsv'))
    if not cards:
        return missing
    selected = max(cards, key=lambda path: path.stat().st_mtime)
    explanations = [a for a in context.catalog.list_artifacts(context.manifest.project_id,
                    analysis_id='programs.explain_expression')
                    if _recorded_artifact_file(context, a, selected, allow_export_copy=True)]
    if len(explanations) != 1:
        return missing
    inputs = explanations[0].metadata.get('input_files', {})
    if not isinstance(inputs, dict) or not all(inputs.get(key) for key in ('scores', 'confidence')):
        return missing
    networks = context.catalog.list_artifacts(context.manifest.project_id, analysis_id='programs.cosplice')
    linked = []
    for key in ('scores', 'confidence'):
        matches = [a for a in networks if _recorded_artifact_file(context, a, inputs[key])]
        if len(matches) != 1:
            return missing
        linked.append(matches[0])
    if len({a.run_fingerprint for a in linked}) != 1:
        return missing
    actuals = [a.metadata.get('network_module_parameters') for a in linked]
    if not isinstance(actuals[0], dict) or any(value != actuals[0] for value in actuals[1:]):
        return missing
    size = actuals[0].get('min_module_size')
    height = actuals[0].get('merge_cut_height')
    if (isinstance(size, bool) or not isinstance(size, int) or size < 1
            or isinstance(height, bool) or not isinstance(height, (int, float))
            or not math.isfinite(height) or not 0 <= height <= 1):
        return missing
    label_policies = [a.metadata.get('network_label_policy') for a in linked]
    label_policy = (dict(label_policies[0]) if isinstance(label_policies[0], dict)
                    and all(value == label_policies[0] for value in label_policies[1:]) else None)
    return {'status': 'recorded', 'network_module_parameters': dict(actuals[0]),
            'network_label_policy': label_policy,
            'network_run_fingerprint': linked[0].run_fingerprint,
            'network_source_artifact_ids': [a.artifact_id for a in linked],
            'module_cards_artifact_id': explanations[0].artifact_id,
            'module_cards_file': str(selected.resolve())}


def _append_network_module_methods(path, record, *, mark_delivery_references=False):
    """Append recorded settings to a newly generated report or delivery copy."""
    from bs4 import BeautifulSoup

    path = Path(path)
    soup = BeautifulSoup(path.read_text(encoding='utf-8'), 'html.parser')
    if mark_delivery_references:
        for delivery_section in soup.find_all('section'):
            heading = delivery_section.find(['h2', 'h3'], recursive=False)
            title = heading.get_text(' ', strip=True).lower() if heading else ''
            section_key = delivery_section.get('id', '')
            if (section_key in {'module-15_CoSplice', 'module-23_Module_explanation'}
                    or any(word in title for word in ('cosplice', 'katana-splice-network'))):
                reference = soup.new_tag('p', attrs={'class': 'network-delivery-reference'})
                reference.string = (
                    'Delivery reference: figures and table previews in this section are retained from their original delivery. '
                    'Their individual threshold-label policies have not been authenticated by the report generator. '
                    'They remain historical/internal references and must not be relabelled by the current policy switch. '
                    'The recorded policy appendix applies only to its explicitly linked source run; '
                    'use that run\'s complete current table views for label status and continuous metrics.')
                if heading:
                    heading.insert_after(reference)
                else:
                    delivery_section.insert(0, reference)
    appendix = path.parent / 'NETWORK_MODULE_PARAMETERS.json'
    section_id = 'network-module-parameters'
    # A previously approved reader may already contain historical provenance.
    # Preserve it in the delivery copy as well as in the original source.
    copy_number = 0
    while appendix.exists() or soup.find(id=section_id) is not None:
        copy_number += 1
        suffix = 'COPY' if copy_number == 1 else f'COPY_{copy_number}'
        appendix = path.parent / f'NETWORK_MODULE_PARAMETERS_{suffix}.json'
        section_id = 'network-module-parameters-' + suffix.lower().replace('_', '-')
    appendix.write_text(json.dumps(record, indent=2), encoding='utf-8')
    section = soup.new_tag('section', id=section_id)
    heading = soup.new_tag('h2')
    heading.string = 'Recorded network module parameters'
    section.append(heading)
    if record['status'] == 'recorded':
        actual = record['network_module_parameters']
        text = (f"Minimum module size: {actual['min_module_size']} splice events. "
                f"Module merge height: {actual['merge_cut_height']:g}. "
                'These are the actual recorded settings for the network run linked to this report, '
                'not significance thresholds. ')
    else:
        text = 'Not recorded. ' + record['reason']
    paragraph = soup.new_tag('p')
    paragraph.string = text
    section.append(paragraph)
    policy = record.get('network_label_policy')
    policy = policy if isinstance(policy, dict) else {}
    enabled = policy.get('apply_threshold_labels')
    filtering = policy.get('hard_filter')
    policy_text = (f"Network label policy: {policy.get('policy_id') or 'not recorded'}. "
                   f"Threshold labels: {'on' if enabled is True else 'off' if enabled is False else 'not recorded'}. "
                   f"Hard filtering: {'yes' if filtering is True else 'no' if filtering is False else 'not recorded'}. "
                   'Off means labels were not applied, not that members failed a threshold. '
                   'Historical policy values are never inferred from current defaults.')
    policy_paragraph = soup.new_tag('p', attrs={'class': 'network-label-policy'})
    policy_paragraph.string = policy_text
    section.append(policy_paragraph)
    if isinstance(policy.get('label_thresholds'), dict):
        rules = soup.new_tag('pre', attrs={'class': 'network-label-rules'})
        rules.string = json.dumps(policy['label_thresholds'], indent=2, sort_keys=True)
        section.append(rules)
    from katana_next.workflows.reporting import NETWORK_LABEL_METHODS
    basis = soup.new_tag('p', attrs={'class': 'network-label-methods'})
    basis.string = NETWORK_LABEL_METHODS
    section.append(basis)
    link = soup.new_tag('a', href=appendix.name)
    link.string = 'Network parameter provenance'
    section.append(link)
    (soup.find('main') or soup.body or soup).append(section)
    path.write_text(str(soup), encoding='utf-8')
    return appendix


class NarrativeReportWorkflow(WorkflowBase):
    version = '1.5.0+network_labels.20261022'
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
        return [Path(parameters['config_path']), Path(c['manifest_path']),
                *[p for p in root.rglob('*') if is_delivery_file(p, root) and p.suffix in ('.tsv', '.json', '.png') and ('24_Narrative_report' not in p.parts)],
                *[a.path for a in self.parent_artifacts(context, parameters)]]

    def parent_artifacts(self, context, parameters):
        record = _network_module_parameter_record(context, self._config(parameters))
        if record['status'] != 'recorded':
            return []
        ids = [record['module_cards_artifact_id'], *record['network_source_artifact_ids']]
        return [context.catalog.get_artifact(artifact_id) for artifact_id in dict.fromkeys(ids)]

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
            parameter_provenance = _append_network_module_methods(path, _network_module_parameter_record(context, c))
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
            writer.existing_file('Recorded network module parameters', parameter_provenance, entity_type='report_provenance')
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
        out = root / '24_Narrative_report' / run_fingerprint[:16]
        build_report(c['study'], root, c['manifest_path'], out)
        parameter_provenance = _append_network_module_methods(out / 'REPORT_EN.html', _network_module_parameter_record(context, c), mark_delivery_references=True)
        for name in ('REPORT_EN.html', 'SOURCE_MANIFEST_SHA256.tsv', 'report_status.json'):
            writer.existing_file(name, out / name, entity_type='report')
        writer.existing_file('Recorded network module parameters', parameter_provenance, entity_type='report_provenance')
        archive = writer.directory / (c['study'] + '_English_complete_results.zip')
        package_report(c['study'], root, out, archive)
        writer.existing_file('Portable English results', archive, entity_type='report_package')
        return writer.finish()
