from pathlib import Path
import json, html
import numpy as np
import pandas as pd
from katana_next.workflows.base import WorkflowBase, WorkflowDefinition, ParameterDefinition, ArtifactWriter
from katana_next.science.module_expression import correlate_modules

class ModuleExplanationWorkflow(WorkflowBase):
    version = '1.2.0+english.20261021'
    definition = WorkflowDefinition(workflow_id='programs.explain_expression', workspace_id='programs', label='Module identity, function and expression context', question='What is in each module, and which expression programs co-vary with it?', method='Fixed consensus modules; GMT mean standardized log-expression scores; overall and condition-adjusted rank correlations; BH separately per association family.', output_summary='Module cards, hub genes, enrichment contributors, all pathway/RBP associations, score heatmap and English explanation.', analysis_unit='fixed module, pathway/RBP, biological sample', parameters=(ParameterDefinition('config_path', 'Module explanation configuration JSON', 'text', '', 'Explicit paths to module scores, hubs, enrichment, confidence, expression, design, GMT and RBP membership.'),), limitations='Same-cohort exploratory association, not external validation. Condition adjustment may remove real treatment signal. Approximate small-sample P-values; RBP correlation is not causality. Module score sign is arbitrary.')

    def _config(self, parameters):
        p = Path(parameters.get('config_path', ''))
        if not p.is_file():
            raise ValueError('Provide module explanation configuration')
        c = json.loads(p.read_text())
        for k in ('scores', 'hubs', 'enrichment', 'confidence', 'expression', 'design', 'gmt', 'rbp'):
            if not Path(c[k]).is_file():
                raise ValueError('Missing ' + k)
        return c

    def input_paths(self, context, parameters):
        c = self._config(parameters)
        return [Path(parameters['config_path']), *[Path(c[k]) for k in ('scores', 'hubs', 'enrichment', 'confidence', 'expression', 'design', 'gmt', 'rbp')]]

    def run(self, context, job, run_fingerprint, parameters):
        c = self._config(parameters)
        read = lambda key: pd.read_csv(c[key], sep='\t', low_memory=False)
        scores = read('scores').set_index('sample_id')
        design = read('design').set_index('sample_id')
        if not set(scores.index) <= {s.sample_id for s in context.manifest.samples}:
            raise ValueError('Unknown score samples')
        expression = read('expression')
        hubs = read('hubs')
        enrich = read('enrichment')
        confidence = read('confidence')
        if set(hubs.module) - set(scores.columns):
            raise ValueError('Hub and score module IDs differ; do not infer M/C mapping')
        samples = list(scores.index)
        if not set(samples) <= set(expression):
            raise ValueError('Expression samples missing')
        values = expression.dropna(subset=['gene_symbol']).groupby('gene_symbol')[samples].mean()
        z = values.sub(values.mean(axis=1), axis=0).div(values.std(axis=1, ddof=1).replace(0, np.nan), axis=0).dropna()
        features = {}
        members = []
        for line in Path(c['gmt']).read_text().splitlines():
            term, _, *genes = line.split('\t')
            genes = sorted(set(genes) & set(z.index))
            if 15 <= len(genes) <= 500:
                key = 'pathway::' + term
                features[key] = z.loc[genes].mean(axis=0)
                members.append(dict(feature=key, genes=';'.join(genes), n_genes=len(genes)))
        rbp = read('rbp')
        rbp_ids = set(rbp.gene_id.dropna().astype(str).str.split('.').str[0])
        annotated = expression[expression.gene_id.str.split('.').str[0].isin(rbp_ids)]
        for symbol in sorted(set(annotated.gene_symbol.dropna()) & set(z.index)):
            features['RBP::' + symbol] = z.loc[symbol]
        if not features:
            raise ValueError('No mapped expression features')
        features = pd.DataFrame(features, index=samples)
        associations = correlate_modules(scores, features, design.condition)
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={'method': self.definition.method, 'limitations': self.definition.limitations, 'input_files': c})
        long = scores.reset_index().melt(id_vars='sample_id', var_name='module', value_name='score').merge(design[['condition']], left_on='sample_id', right_index=True)
        group = long.groupby(['module', 'condition']).score.agg(['count', 'mean', 'std']).reset_index()
        group['row_id'] = group.module + '::' + group.condition
        writer.table('Module condition score summary', 'module_condition_summary.parquet', group, entity_type='module')
        for name, table in [('module_expression_associations', associations), ('expression_feature_members', pd.DataFrame(members)), ('module_hub_genes', hubs), ('module_function_enrichment', enrich), ('module_consensus_confidence', confidence)]:
            table = table.copy()
            table['row_id'] = table.index.astype(str)
            writer.table(name, name + '.parquet', table, entity_type='module_evidence')
        cards = []
        for module in scores:
            eh = enrich[enrich.module.eq(module)].sort_values('FDR')
            significant = eh[eh.FDR.le(0.05)]
            hg = hubs[hubs.module.eq(module)].sort_values('hub_gene_rank')
            cf = confidence[confidence.consensus_module.eq(module)]
            cards.append(dict(module=module, function_label=significant.iloc[0].term if len(significant) else 'No FDR-supported functional label', strongest_tested_term=eh.iloc[0].term if len(eh) else '', strongest_term_FDR=float(eh.iloc[0].FDR) if len(eh) else None, contributing_genes=eh.iloc[0].overlap_genes if len(eh) else '', significant_terms=len(significant), top_hub_genes=';'.join(hg.gene_symbol.fillna(hg.gene_id).head(8)), member_gene_count=hg.gene_id.nunique(), n_events=int(cf.iloc[0].n_events) if len(cf) else None, mean_coassignment=float(cf.iloc[0].mean_within_module_coassignment) if len(cf) else None))
        cards = pd.DataFrame(cards)
        cards['row_id'] = cards.module
        writer.table('Module interpretation cards', 'module_cards.parquet', cards, entity_type='module')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        order = design.loc[samples].sort_values(['condition']).index
        fig, ax = plt.subplots(figsize=(max(9, len(samples) * 0.5), max(3, len(scores.columns) * 0.5)))
        im = ax.imshow(scores.loc[order].T, aspect='auto', cmap='RdBu_r', vmin=-scores.abs().max().max(), vmax=scores.abs().max().max())
        ax.set_xticks(range(len(order)), [s + '\n' + str(design.loc[s, 'condition']) for s in order], rotation=90, fontsize=8)
        ax.set_yticks(range(len(scores.columns)), scores.columns)
        ax.set_title('Fixed consensus modules across individual samples\nScore sign is arbitrary; this is not an interaction test')
        fig.colorbar(im, ax=ax, label='Module score')
        fig.tight_layout()
        for ext in ('png', 'svg', 'pdf'):
            p = writer.directory / ('module_sample_context.' + ext)
            fig.savefig(p, dpi=300, bbox_inches='tight')
            writer.existing_file(p.name, p, entity_type='plot')
        plt.close(fig)
        p = writer.directory / ('REPORT_' + 'EN' + '.html')
        sections = []
        for module in scores:
            terms = enrich[enrich.module.eq(module)].sort_values('FDR')
            members = hubs[hubs.module.eq(module)].sort_values('hub_gene_rank')
            sections.append('<section><h2>' + html.escape(str(module)) + '</h2><p>' + 'Ten lowest-FDR terms, not necessarily significant. Whole-module enrichment does not imply every member responds to condition.' + '</p>' + terms.head(10).to_html(index=False, escape=True) + '<h3>Hub genes</h3>' + members.head(10).to_html(index=False, escape=True) + '<h3>Condition scores (descriptive, not a test)</h3>' + group[group.module.eq(module)].to_html(index=False, escape=True) + '</section>')
        links = '<p><a href="module_function_enrichment.tsv">All functional tests and contributors</a> | <a href="module_hub_genes.tsv">All hub genes</a> | <a href="module_consensus_confidence.tsv">Consensus confidence</a> | <a href="module_expression_associations.tsv">All expression/RBP associations</a></p>'
        p.write_text('<meta charset="utf-8"><style>body{margin:30px;font:16px sans-serif}section{margin-top:30px;overflow:auto}table{font-size:12px;border-collapse:collapse}td,th{padding:6px;border:1px solid #ddd}</style><h1>' + 'Module identity and interpretation' + '</h1><p>' + 'Modules group co-varying splice events; they are not condition-exclusive gene lists. Functional labels require enrichment FDR support. Hubs are representative network members, not causal regulators. Overall association can be condition-driven, so adjusted associations are reported separately. Small-sample associations remain exploratory. Low co-assignment confidence indicates uncertain boundaries despite significant associations.' + '</p><img style="max-width:100%" src="module_sample_context.png">' + cards.to_html(index=False, escape=True) + links + ''.join(sections), encoding='utf-8')
        writer.existing_file('Module identity and interpretation', p, entity_type='report', metadata={'mime_type': 'text/html'})
        return writer.finish()
