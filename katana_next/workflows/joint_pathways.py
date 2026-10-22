"""Separate expression and splicing ORA with explicit assay backgrounds."""
import pandas as pd
from katana_next.science.enrichment import canonical_gene_mapping, map_gene_sets, declared_gene_sets, library_overrepresentation
from katana_next.workflows.base import WorkflowBase, WorkflowDefinition, ArtifactWriter
from katana_next.workflows.sources import resolve_table

def joint_pathways(frame, sets):
    if frame.duplicated(['comparison_id', 'gene_id']).any():
        raise ValueError('Expected one gene per comparison')
    canonical, aliases, mapping_audit = canonical_gene_mapping(frame)
    frame = frame.assign(canonical_gene_id=canonical)
    sets, set_audit = map_gene_sets(sets, aliases)
    rows = []
    availability = []
    for comparison, group in frame.groupby('comparison_id'):
        for layer, available, supported in [('expression', 'deg_available', 'deg_supported'), ('splicing', 'as_available', 'as_supported')]:
            tested = set(group.loc[group[available].fillna(False), 'canonical_gene_id']) - {''}
            selected = set(group.loc[group[available].fillna(False) & group[supported].fillna(False), 'canonical_gene_id']) - {''}
            mapped = set().union(*sets.values()) if sets else set()
            status = 'available' if tested & mapped else 'unavailable: no annotated assay background'
            availability.append(dict(comparison_id=comparison, layer=layer, tested_genes=len(tested), selected_genes=len(selected), status=status))
            if status != 'available':
                continue
            usable = {t: m for t, m in sets.items() if m & tested}
            part = library_overrepresentation(selected, tested, usable, minimum_overlap=0)
            part['comparison_id'] = comparison
            part['layer'] = layer
            part['FDR_family'] = 'comparison x assay x annotation library; all eligible terms'
            rows.append(part)
    result = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    return (result, pd.DataFrame(availability), mapping_audit, set_audit)

class JointPathwayWorkflow(WorkflowBase):
    version = '1.0.0+english.20261021'
    definition = WorkflowDefinition('interpretation.joint_pathways', 'prioritize', 'Expression-splicing pathway evidence', 'Do expression and splicing changes converge on the same annotated processes?', 'Separate unique-gene overrepresentation tests per assay, comparison and library with respective testable backgrounds; no combined P-value.', 'All terms, contribution genes, mapping and availability audits, paired evidence table and English report.', 'pathway within comparison and assay', supports_contrast_selection=True, required_asset_kinds=('gene_sets_gmt',), limitations='Assays are dependent. Co-enrichment does not establish causality or pathway activation. ORA may retain gene/event-detectability bias; not a matched-bias calibration. Non-significant expression is not unchanged expression.')

    def parent_artifacts(self, context, parameters):
        source = resolve_table(context, 'prioritize.evidence', 'Candidate evidence matrix', parameters)
        if source is None:
            raise ValueError('Run candidate evidence first; no implicit analysis')
        return [source]

    def input_paths(self, context, parameters):
        return [a.path for a in self.parent_artifacts(context, parameters)] + [a.path for a in context.manifest.assets if a.kind == 'gene_sets_gmt']

    def run(self, context, job, run_fingerprint, parameters):
        source = self.parent_artifacts(context, parameters)[0]
        frame = pd.read_parquet(source.path)
        if parameters.get('contrast_ids'):
            frame = frame[frame.comparison_id.isin(parameters['contrast_ids'])]
        sets, libraries = declared_gene_sets(context.manifest, 'gene_sets_gmt')
        result, availability, mapping, terms = joint_pathways(frame, sets)
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={'limitations': self.definition.limitations, 'method': self.definition.method})
        for name, table in [('pathway_evidence', result), ('assay_availability', availability), ('gene_mapping', mapping), ('term_mapping', terms), ('libraries', libraries)]:
            table = table.copy()
            table['row_id'] = table.index.astype(str)
            writer.table(name, name + '.parquet', table, entity_type='pathway_evidence', parents=(source.artifact_id,))
        paired = pd.DataFrame()
        if not result.empty:
            paired = result.pivot(index=['comparison_id', 'term'], columns='layer', values='fdr').reset_index()
            for layer in ('expression', 'splicing'):
                if layer not in paired:
                    paired[layer] = float('nan')
            paired['both_FDR_le_0_05'] = paired.expression.le(0.05) & paired.splicing.le(0.05)
            paired['row_id'] = paired.index.astype(str)
            writer.table('Paired pathway evidence', 'paired_pathways.parquet', paired, entity_type='pathway_evidence', parents=(source.artifact_id,))
        message = self.definition.limitations
        p = writer.directory / f"REPORT_{'EN'}.html"
        p.write_text('<meta charset="utf-8"><h1>Expression / splicing pathway evidence</h1><p>' + message + '</p>' + availability.to_html(index=False) + paired.head(30).to_html(index=False) + '<p>Preview: first 30 terms, not a significance-filtered list. <a href="pathway_evidence.tsv">All results and contributor genes</a></p>', encoding='utf-8')
        writer.existing_file(p.name, p, entity_type='report', parents=(source.artifact_id,))
        return writer.finish()
