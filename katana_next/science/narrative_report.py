"""Portable English, question-led reports with traceable result attachments."""
from pathlib import Path
import hashlib, html, json, re, zipfile, os
import pandas as pd
from katana_cosplice.branding import display_module_name
from katana_next.science.report_language import is_retired_report
STYLE = 'body{font:17px Georgia,"Microsoft YaHei",serif;color:#20313b;background:#f1f3f4;margin:0}main{max-width:1100px;margin:auto;padding:40px}section{background:white;padding:28px;margin:24px 0;border-top:3px solid #28617a}h1{font-size:36px}h2{font-size:25px}p,li{line-height:1.8}table{border-collapse:collapse;font:13px sans-serif;width:100%}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:left}.scroll{overflow:auto}img{max-width:100%}figcaption{font:14px sans-serif;line-height:1.6;color:#40535f}a{color:#146484}details{margin:15px 0}summary{cursor:pointer;font-weight:bold}nav{line-height:2}@media print{body{background:white}section{break-inside:avoid}details{display:block}}'

def is_delivery_file(path, base):
    """Exclude active application bookkeeping from scientific attachments."""
    path = Path(path)
    relative = path.relative_to(base)
    return ('.katana' not in relative.parts
            and path.suffix.lower() != '.lock'
            and path.is_file())

def build_report(study, base, manifest_path, output):
    base = Path(base)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    from katana_next.science.report_complete import polished_heatmaps, complete_sections
    polished_heatmaps(base)
    manifest = json.loads(Path(manifest_path).read_text())
    report_metadata = manifest.get('metadata', {})
    design_labels = report_metadata.get('report_design_labels', {})
    design_notes = report_metadata.get('report_design_notes', [])
    if isinstance(design_notes, str):
        design_notes = [design_notes]
    if not isinstance(design_labels, dict) or not isinstance(design_notes, list) or not all(isinstance(note, str) for note in design_notes):
        raise ValueError('Report design labels must be a mapping and design notes must be text or a list of text')
    e = lambda x: html.escape(str(x))

    def latest(pattern, root=None):
        paths = list((root or base).rglob(pattern))
        return max(paths, key=lambda p: p.stat().st_mtime) if paths else None
    totals = pd.read_csv(base / '20_Presentation_report/comparison_totals.tsv', sep='\t')
    intersections = pd.read_csv(base / '20_Presentation_report/DEG_DAS_gene_intersection.tsv', sep='\t')
    designs = pd.DataFrame([{'Sample': s['sample_id'], 'Condition': s['condition'], 'Biological replicate': s.get('biological_sample_id', 'unspecified')} for s in manifest['samples']])
    modules_path = latest('module_cards.tsv', base / '23_Module_explanation')
    cards = pd.read_csv(modules_path, sep='\t') if modules_path else pd.DataFrame()
    expression_roots = sorted({p.parent for p in (base / '22_Expression_interaction_pathways').rglob('interaction.tsv')})
    exon_roots = sorted({p.parent for p in (base / '16_Interaction').rglob('dexseq_interaction.tsv')})
    from katana_next.science.report_annotation import annotate_exons
    annotation_frames = [pd.read_csv(d / 'interaction.tsv', sep='\t') for d in expression_roots]
    annotation = pd.concat([f[['gene_id', 'gene_symbol']] for f in annotation_frames if {'gene_id', 'gene_symbol'} <= set(f)], ignore_index=True) if annotation_frames else pd.DataFrame(columns=['gene_id', 'gene_symbol'])
    annotated_exons = {}
    for number, d in enumerate(exon_roots):
        annotated_exons[d] = annotate_exons(pd.read_csv(d / 'dexseq_interaction.tsv', sep='\t'), annotation)
        annotated_exons[d].to_csv(output / f'dexseq_interaction_annotated_{number + 1}.tsv', sep='\t', index=False)
    if not expression_roots or not exon_roots:
        raise ValueError('Formal expression/exon results missing; absent analysis must not be reported as zero significant results')
    expression_counts = []
    exon_counts = []
    for d in expression_roots:
        x = pd.read_csv(d / 'interaction.tsv', sep='\t')
        g = pd.read_csv(d / 'ranked_pathways.tsv', sep='\t')
        expression_counts.append({'Design': d.parent.name, 'Tested genes': int(x.padj.notna().sum()), 'Interaction FDR <= 0.05': int(x.padj.le(0.05).sum()), 'Pathway tests FDR <= 0.05': int(g.study_family_FDR.le(0.05).sum())})
    for d in exon_roots:
        x = pd.read_csv(d / 'dexseq_interaction.tsv', sep='\t')
        design_key = d.relative_to(base).as_posix()
        design_label = design_labels.get(design_key, d.name if d.parent.name == '16_Interaction' else d.parent.name)
        exon_counts.append({'Design': design_label, 'Tested bins': int(x.interaction_FDR.notna().sum()), 'Interaction FDR <= 0.05': int(x.interaction_FDR.le(0.05).sum())})
    status = {'study': study, 'expression_designs': expression_counts, 'exon_designs': exon_counts, 'modules': len(cards), 'missing_inputs': ['BAM for Sashimi', 'independent cohort for preservation', 'matched CLIP for binding', 'perturbation/protein assays for causal mechanism']}
    (output / 'report_status.json').write_text(json.dumps(status, indent=2))
    all_files = [p for p in base.rglob('*') if is_delivery_file(p, base) and output not in p.parents and (p.suffix.lower() not in ('.zip', '.html'))]
    source_inventory = []
    for p in all_files:
        h = hashlib.sha256()
        with p.open('rb') as f:
            for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
                h.update(chunk)
        source_inventory.append({'path': p.relative_to(base).as_posix(), 'bytes': p.stat().st_size, 'sha256': h.hexdigest()})
    pd.DataFrame(source_inventory).to_csv(output / 'SOURCE_MANIFEST_SHA256.tsv', sep='\t', index=False)
    sections = []

    def table(frame):
        frame = frame.drop(columns=['row_id'], errors='ignore')
        return '<div class="scroll">' + frame.to_html(index=False, escape=True) + '</div>'

    def figure(path, caption):
        if not path or not Path(path).is_file():
            return ''
        rel = Path(os.path.relpath(path, output)).as_posix()
        return '<figure><img loading="lazy" src="' + e(rel) + '"><figcaption>' + e(caption) + '</figcaption></figure>'

    def add(title, text, extra=''):
        sections.append('<section id="s' + str(len(sections) + 1) + '"><h2>' + e(display_module_name(title)) + '</h2><p>' + e(display_module_name(text)) + '</p>' + extra + '</section>')
    n_expr = sum((r['Interaction FDR <= 0.05'] for r in expression_counts))
    n_exon = sum((r['Interaction FDR <= 0.05'] for r in exon_counts))
    add('1. Main findings and strength of evidence', f'The {len(totals)} reported contrasts identify {totals.significant_events.min()}–{totals.significant_events.max()} significant AS events and {totals.DEG.min()}–{totals.DEG.max()} genes passing the original DEG criteria per contrast. These are different units, not a direct comparison of biological effect strength. Added designs contain {n_expr} significant gene–design interaction records and {n_exon} significant exon-bin–design records; counts are not deduplicated across designs.', table(totals))
    add('2. Experimental question and design', 'Biological replicate and donor identities follow the project manifest; technical measurements do not create additional biological replicates. Contrast direction follows the declared numerator relative to the denominator. Simple effects ask whether a change occurs within one background; interactions ask whether responses depend on background. Shared controls can induce correlation and are not independent mechanistic evidence.' + (' ' + ' '.join(design_notes) if design_notes else ''), table(designs))
    overview = base / '20_Presentation_report/manuscript_assets/Figure_1_evidence_overview.png'
    add('3. Abundance and splicing: complementary responses', 'DEG measures abundance, DAS specific splice choices, DTU transcript proportions and ordinary DEXSeq relative exon usage. Joint classification requires assessability in both layers; unavailable evidence is separate. Non-significance does not establish no change. Attachments retain all non-significant rows.', figure(overview, 'Each contrast is counted separately. Event and gene units differ; bar heights are not statistical comparisons.') + table(intersections))
    heatmaps = sorted((output / 'figures').glob('*/significant/rowZ.png'))
    add('4. Patterns across biological replicates', 'The default is event-wise scaled PSI: (PSI minus the event sample mean)/sample standard deviation. Pale red indicates above-mean values, pale blue below-mean, white near-mean and gray missing. Colors are neither raw PSI nor significance. Only the display is clipped at +/-2.5; data and existing clustering remain unchanged. Separation after selection is not independent validation.', ''.join(('<h3>' + e(p.parent.parent.name) + '</h3>' + figure(p, 'Default: row-scaled relative splicing patterns.') + '<details><summary>' + 'Raw PSI and stable-event heatmaps' + '</summary>' + figure(p.with_name('PSI.png'), 'Supplementary view: raw PSI, range 0-1.') + figure(p.parent.parent / 'stable_significant/rowZ.png', 'Row-scaled stable-significant events; thresholds unchanged.') + '</details>' for p in heatmaps)))
    interaction_extra = table(pd.DataFrame(expression_counts)) + table(pd.DataFrame(exon_counts))
    for d in expression_roots:
        interaction_extra += figure(d / 'expression_effects.png', d.parent.name)
    for d in exon_roots:
        full_exons = annotated_exons[d]
        leading = full_exons.dropna(subset=['interaction_FDR']).nsmallest(12, 'interaction_FDR')
        cols = [c for c in ['gene_name', 'bin_id', 'pvalue', 'interaction_FDR', 'descriptive_fraction_interaction'] if c in leading]
        interaction_extra += '<p><a href="dexseq_interaction_annotated_' + str(exon_roots.index(d) + 1) + '.tsv">' + 'Download all DEXSeq interaction results with gene names' + '</a></p>'
        interaction_extra += '<h3>DEXSeq: ' + e(d.name if d.parent.name == '16_Interaction' else d.parent.name) + '</h3>' + figure(d / 'four_group_exon_usage.png', 'Four-group replicate exon-bin fractions; not rMATS PSI.') + table(leading[cols]) + '<p>' + 'The 12 smallest adjusted P-values are shown; they are not necessarily all significant. Full significant and non-significant table: ' + '<a href="' + e(Path(os.path.relpath(d / 'dexseq_interaction.tsv', output)).as_posix()) + '">dexseq_interaction.tsv</a></p>'
    add('5. Interaction: does the response depend on background?', 'Expression interactions use official DESeq2 with tximport effective-length handling; exon-usage interactions use official DEXSeq full/reduced models. Expression coefficients are differences of log2 effects; exon fractions are descriptive. Legacy beta-binomial results remain exploratory and their FDR is not merged. Non-significant interaction does not establish equivalence.', interaction_extra)
    pathway_extra = ''
    pathway_themes = []
    for d in expression_roots:
        gp = pd.read_csv(d / 'ranked_pathways.tsv', sep='\t')
        hits = gp[gp.contrast.eq('interaction') & gp.study_family_FDR.le(0.05)].sort_values('study_family_FDR')
        concise = hits[hits.pathway.str.startswith('KEGG')]
        concise = (concise if len(concise) else hits).head(5)
        if len(concise):
            show = concise[['pathway', 'NES', 'study_family_FDR', 'leadingEdge']].copy()
            show['leadingEdge'] = show.leadingEdge.fillna('').map(lambda s: ';'.join(s.split(';')[:12]))
            pathway_extra += '<h3>' + e(d.parent.name) + '</h3>' + table(show)
            pathway_themes.extend(concise.pathway.tolist())
        pathway_extra += figure(d / 'ranked_pathways.png', 'NES sign indicates enrichment at opposite ends of the ranking, not direct pathway activation/inhibition. Stars use the declared study-wide FDR family across contrasts and tested pathways.')
    add('6. Functional interpretation and contributing genes', 'Ranked analysis uses signed Wald statistics from all testable genes rather than relaxing DEG thresholds. Leading-edge genes contribute to enrichment; all pathways, including non-significant results, are retained in ranked_pathways.tsv. Expression and splicing enrichment use different backgrounds; their P-values are not combined.', pathway_extra)
    mp = latest('module_sample_context.png', base / '23_Module_explanation')
    module_functions = ''
    enrichment_path = latest('module_function_enrichment.tsv', base / '23_Module_explanation')
    if enrichment_path:
        enrich = pd.read_csv(enrichment_path, sep='\t')
        for module in cards.get('module', []):
            selected = enrich[enrich.module.eq(module)].sort_values('FDR').head(5)
            module_functions += '<h3>' + e(module) + ' : functional evidence' + '</h3>'
            module_functions += '<p>' + 'Up to five lowest-FDR terms are shown, not necessarily significant. The enrichment background comprises eligible network genes, counted once per gene rather than once per event. Enrichment indicates functional concentration, not activation, inhibition or causality.' + '</p>'
            cols = [c for c in ['term', 'module_genes', 'universe_genes', 'term_genes_in_universe', 'overlap_count', 'overlap_genes', 'FDR'] if c in selected]
            module_functions += table(selected[cols]) if len(selected) else '<p>' + 'No functional test available; this does not establish absence of function.' + '</p>'
        module_functions += '<p><a href="https://doi.org/10.1186/1471-2105-9-559">WGCNA: module–trait association, functional enrichment and intramodular hubs</a></p>'
    add('7. CoSplice: membership, function and confidence', 'Consensus C identifiers are used consistently rather than substituting initial M identifiers. Modules are co-varying events, not condition-exclusive entities. Hubs are representative members, not proven upstream regulators. Functional labels require enrichment FDR support; otherwise identifiers are retained. Mean co-assignment measures resampling consistency and is distinct from association significance. Overall and condition-adjusted module–expression/RBP correlations have separate family-wide FDR control and remain same-cohort exploratory associations.', table(cards) + figure(mp, 'Module scores by sample; signs are arbitrary and color is not interaction significance.'))
    sections[-1] = sections[-1].replace('</section>', module_functions + '</section>')
    usage = base / '20_Presentation_report/manuscript_assets/DTU_DEXSeq_source_data.tsv'
    add('8. Usage evidence and candidate selection', 'DTU tests changes in transcript proportions within a gene; ordinary DEXSeq tests relative exon usage. Neither is the four-group interaction test, and neither requires significant total expression. A supported gene in this table has at least one passing feature; these are not newly calculated gene-level FDR values. Candidates are ordered by supporting dimensions (stable AS, DEG, DTU and DEXSeq), then stable significant event count, then the largest absolute dPSI among those events. These dimensions are not independent replications. Gene cards display selected genes and cached abundance/splicing evidence; they are not GeneCards disease scores or proof of mechanism.', table(pd.read_csv(usage, sep='\t')) if usage.is_file() else '<p>Source summary unavailable, not a zero result.</p>')
    evidence_rows = [('02_AS_reliability', 'Separates measurability, replicate reliability and significance; filtering is not evidence of no biological effect.'), ('07_Cross_comparison', 'Exact-event overlap and direction agreement answer separate questions; shared controls are not independent validation.'), ('09_GO_KEGG', 'Tests selected gene sets against an eligible background; distinct from ranked expression enrichment.'), ('10_RBP', 'Layer 1: original RBP identity/declared target-set annotations remain gene-level prioritization, not event-site binding.'), ('11_Isoform_usage', 'Describes transcript composition, not demonstrated protein isoform change.'), ('13_Isoform_consequences', 'Annotates possible coding, domain or NMD consequences; predictions depend on transcript/reference versions.'), ('14_RBP_motifs', 'Layer 2: boundary motifs/reference CLIP are candidate upstream clues, not automatically mature-RNA consequences. Layer 3 MotifTrace requires actual RNA pairs, loci and usage evidence.')]
    add('9. How the remaining analyses support interpretation', 'The attachment inventory retains all analysis outputs. These layers answer different questions; predictions and missing inputs are not presented as validated findings.', table(pd.DataFrame(evidence_rows, columns=['Result directory', 'Purpose and interpretation boundary'])))
    sections.append(complete_sections(base, output, 'EN'))
    add('10. Integrated interpretation and validation priorities', f'The experiment can be interpreted across related but distinct abundance, event-level splicing and exon-usage responses. The {len(cards)} consensus modules organize events but do not establish mechanism. Validation should prioritize reproducible, localized events with interpretable effects and functional evidence. GeneCard evidence ranking is prioritization, not causal probability. NMD, domain and motif annotations generate hypotheses. Missing BAM prevents Sashimi validation; missing independent cohorts prevents external preservation claims; absent CLIP or perturbation/protein assays prevents direct binding or causal claims.' + ' Expression-interaction themes prioritized for inspection: ' + ('; '.join(dict.fromkeys(pathway_themes)) if pathway_themes else 'No pathways passed the declared FDR.') + ' The new four-group expression model remains distinct from imported pairwise DEG results: effective-length handling, pooled dispersion estimation and eligibility can change statistical results. Differences between pipelines are not independent validation.')
    appendix = '<section><h2>' + 'Attachments and reproducibility' + '</h2><p>' + 'All significant and non-significant outputs are retained. Existing ZIP containers and superseded HTML are not nested; underlying analysis files are included. SHA256 manifests support verification.' + '</p><details><summary>' + 'Open complete file inventory' + '</summary><ul>'
    for row in source_inventory:
        appendix += '<li><a href="' + e(Path(os.path.relpath(base / row['path'], output)).as_posix()) + '">' + e(row['path']) + '</a></li>'
    appendix += '</ul></details><p><a href="SOURCE_MANIFEST_SHA256.tsv">SHA256</a></p><p>Methods: <a href="https://www.nature.com/articles/ncomms14060">DEXSeq factorial precedent</a>; <a href="https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html">DESeq2</a>; <a href="https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html">fgsea</a>.</p></section>'
    toc = '<nav class="toc">' + ''.join(('<a href="#' + anchor + '">' + title + '</a>' for anchor, title in re.findall('<section id="([^"]+)"><h2>(.*?)</h2>', ''.join(sections)))) + '</nav>'
    page = '<!doctype html><html lang="' + 'en' + '"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>' + e(study) + '</title><style>' + STYLE + '.toc{display:flex;flex-wrap:wrap;gap:8px 20px;padding:20px;background:white}.toc a{font:14px sans-serif}figure{margin:24px 0;padding:12px;border:1px solid #e3e8eb}section{scroll-margin-top:20px}details[open]{padding:12px;background:#f7f9fa}h3{font-size:17px;overflow-wrap:anywhere}@media(max-width:700px){main{padding:12px}section{padding:16px}}</style><main><h1>' + e(study) + ' · ' + 'Findings and interpretation' + '</h1><p>' + 'Findings first, followed by evidence; methods, results and limitations are separated.' + '</p>' + toc + ''.join(sections) + appendix + '</main></html>'
    if manifest.get('metadata', {}).get('gene_filter'):
        from katana_next.domain.models import ProjectManifest
        from katana_next.science.gene_filter import gene_filter_footnote
        footer = '<footer id="gene-exclusion-footnote"><p>' + e(gene_filter_footnote(ProjectManifest.model_validate(manifest))) + '</p></footer>'
        page = page.replace('</main>', footer + '</main>')
    (output / ('REPORT_' + 'EN' + '.html')).write_text(page, encoding='utf-8')
    return all_files

def package_report(study, base, output, archive):
    base = Path(base)
    output = Path(output)
    with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED, compresslevel=1, allowZip64=True) as z:
        for p in base.rglob('*'):
            if is_delivery_file(p, base) and p.suffix.lower() != '.zip' and not is_retired_report(p) and (p.suffix.lower() != '.html' or output in p.parents):
                z.write(p, study + '/' + p.relative_to(base).as_posix())
    with zipfile.ZipFile(archive) as z:
        bad = z.testzip()
        if bad:
            raise ValueError('Archive CRC failure: ' + bad)
