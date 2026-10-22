"""Reusable locus-wide perturbation evidence from declared project inputs."""
import html
from pathlib import Path
import numpy as np
import pandas as pd
from katana_next.science.perturbation import gene_exons, audit_gene
from katana_next.workflows.base import WorkflowBase, WorkflowDefinition, ParameterDefinition, ArtifactWriter
from katana_next.workflows.sources import selected_source, selection_rows

class PerturbationWorkflow(WorkflowBase):
    version = '1.0.1+english.20261021'
    definition = WorkflowDefinition('validate.perturbation', 'validate', 'Target-gene perturbation evidence', 'Is target-gene RNA structure altered consistently across the declared experimental groups?', 'Whole annotated locus: primary BAM base coverage, junction alignments and nonduplicate anchored sensitivity counts; no assumed target exon.', 'All-sample exon and junction tables, annotation, QC, coverage figure and English evidence report.', 'gene locus and biological sample', parameters=(ParameterDefinition('gene', 'Target gene name or ID', 'text', '', 'One unambiguous project-GTF gene; alternatively select gene/event rows for one gene.'), ParameterDefinition('selection_id', 'Saved gene/event selection', 'selection', '', 'Selection must identify exactly one gene.'), ParameterDefinition('min_mapq', 'Minimum MAPQ', 'number', 20, 'NH>1 alignments excluded; missing NH audited.', minimum=0, maximum=255), ParameterDefinition('min_anchor', 'Sensitivity junction anchor bp', 'number', 8, 'Both immediately adjacent aligned blocks; not a significance threshold.', minimum=1, maximum=150)), supports_contrast_selection=True, limitations='RNA evidence cannot confirm genomic recombination, protein loss or knockout efficiency. Read counts are not independent biological replicates. Within-gene scaling can hide total expression loss. All samples in one declared contrast must have BAM plus index.')

    def parent_artifacts(self, context, parameters):
        source = selected_source(context, parameters)
        return [source] if source else []

    def _inputs(self, context, parameters):
        ids = parameters.get('contrast_ids') or [c.contrast_id for c in context.manifest.contrasts]
        if len(ids) != 1:
            raise ValueError('Choose exactly one matched-background comparison')
        contrast = next((c for c in context.manifest.contrasts if c.contrast_id == ids[0]), None)
        if contrast is None:
            raise ValueError('Unknown comparison')
        required = set(contrast.numerator.sample_ids) | set(contrast.denominator.sample_ids)
        samples = {s.sample_id: s for s in context.manifest.samples if s.sample_id in required}
        if set(samples) != required:
            raise ValueError('Comparison samples missing from project')
        if len({s.biological_sample_id for s in samples.values()}) != len(samples):
            raise ValueError('Repeated biological units need a paired/repeated-measures analysis, not independent tracks')
        paths = []
        for s in samples.values():
            bam = s.assets.get('bam')
            index = s.assets.get('bai') or s.assets.get('csi')
            if bam and (not index):
                index = next((p for p in (Path(str(bam) + '.bai'), Path(str(bam) + '.csi'), bam.with_suffix('.bai')) if p.is_file()), None)
            if not bam or not bam.is_file() or (not index) or (not index.is_file()):
                raise ValueError(f'Missing BAM/index for {s.sample_id}; missing evidence is not a negative result')
            if index not in (Path(str(bam) + '.bai'), Path(str(bam) + '.csi'), bam.with_suffix('.bai')):
                raise ValueError('Index must be adjacent to BAM for this reader')
            paths.extend([bam, index])
        gtf = context.manifest.annotation.gtf
        if not gtf or not gtf.is_file():
            raise ValueError('Matching project GTF required')
        return (contrast, samples, [gtf, *paths])

    def input_paths(self, context, parameters):
        return self._inputs(context, parameters)[2] + [a.path for a in self.parent_artifacts(context, parameters)]

    def run(self, context, job, run_fingerprint, parameters):
        contrast, samples, paths = self._inputs(context, parameters)
        gene = str(parameters.get('gene', '')).strip()
        parents = self.parent_artifacts(context, parameters)
        if parents:
            frame = selection_rows(context, parameters, parents[0], pd.read_parquet(parents[0].path))
            genes = frame.gene_id.dropna().astype(str).str.split('.').str[0].unique()
            if len(genes) != 1:
                raise ValueError('Select exactly one gene for the locus audit')
            if gene and gene not in set(genes) | set(frame.get('gene_symbol', [])):
                raise ValueError('Typed gene conflicts with saved selection')
            gene = genes[0]
        if not gene:
            raise ValueError('Specify a target gene or saved gene selection')
        exons = gene_exons(paths[0], gene)
        segments, junctions, qc, depth = audit_gene(exons, {k: s.assets['bam'] for k, s in samples.items()}, int(parameters['min_mapq']), int(parameters['min_anchor']), lambda i, n, s: job.report(0.1 + 0.65 * i / n, f'Inspecting {s} ({i + 1}/{n})'))
        groups = {s: contrast.numerator.label if s in contrast.numerator.sample_ids else contrast.denominator.label for s in samples}
        evidence_columns = ['chrom', 'intron_start_1based', 'intron_end_1based', 'filter']
        joint = pd.DataFrame()
        if not junctions.empty:
            joint = junctions.pivot_table(index=evidence_columns, columns='sample_id', values='alignment_count', fill_value=0).reindex(columns=list(samples), fill_value=0)
            for sample in samples:
                retained = qc.loc[qc.sample_id.eq(sample), 'retained_alignments'] if 'retained_alignments' in qc else pd.Series(dtype=float)
                if retained.empty or not retained.fillna(0).gt(0).any():
                    joint[sample] = float('nan')
            joint = joint.reset_index()
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={'limitations': self.definition.limitations, 'gene': gene, 'comparison': contrast.model_dump(mode='json'), 'coordinate_convention': 'tables 1-based inclusive; internal GTF segments 0-based half-open', 'genome_build': context.manifest.annotation.genome_build, 'missing_proof': ['genomic recombination', 'protein loss']})
        parent_ids = tuple((a.artifact_id for a in parents))
        for name, table in [('exon_segments', segments), ('junction_support', junctions), ('sample_QC', qc), ('transcript_annotation', exons)]:
            table = table.copy()
            if 'sample_id' in table:
                table['group'] = table.sample_id.map(groups)
            table['row_id'] = table.index.astype(str)
            writer.table(name, name + '.parquet', table, entity_type='gene_validation', parents=parent_ids)
        if not joint.empty:
            joint['row_id'] = joint.index.astype(str)
            writer.table('Junction comparison across all samples', 'junction_comparison.parquet', joint, entity_type='gene_validation', parents=parent_ids)
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(len(samples), 1, figsize=(12, max(3, 2 * len(samples))), sharex=True, squeeze=False)
        x = np.arange(int(exons.start.min()), int(exons.end.max())) + 1
        for i, s in enumerate(samples):
            axes[i, 0].fill_between(x, depth[i], color='#c77d86' if s in contrast.numerator.sample_ids else '#7babc9', lw=0)
            axes[i, 0].set_ylim(0, max(1, depth.max()) * 1.05)
            axes[i, 0].set_ylabel(s, fontsize=8)
        axes[-1, 0].set_xlabel(f'{exons.chrom.iloc[0]} | 1-based coordinate')
        fig.suptitle(f'{gene}: raw coverage, not knockout efficiency')
        fig.tight_layout()
        for ext in ('png', 'svg', 'pdf'):
            p = writer.directory / f'coverage.{ext}'
            fig.savefig(p, dpi=180)
            writer.existing_file(p.name, p, entity_type='plot', parents=parent_ids)
        plt.close(fig)
        message = self.definition.limitations
        p = writer.directory / f"REPORT_{'EN'}.html"
        p.write_text('<meta charset="utf-8"><h1>' + html.escape(gene) + '</h1><p>' + html.escape(message) + '</p><img style="max-width:100%" src="coverage.png">' + qc.to_html(index=False) + '<p><a href="junction_support.tsv">All junctions</a> | <a href="exon_segments.tsv">All exon segments</a></p>', encoding='utf-8')
        writer.existing_file(p.name, p, entity_type='report', parents=parent_ids)
        return writer.finish()
