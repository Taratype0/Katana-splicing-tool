"""Render existing evidence explicitly; never rerun or relax statistical tests."""
from pathlib import Path
import html
import json
import os
import numpy as np
import pandas as pd
from katana_next.science.heatmap_style import row_scaled, pastel_diverging

def polished_heatmaps(base):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import dendrogram
    base = Path(base)
    for source in (base / '03_AS_landscape/Native_PSI_clustering').glob('*/*/display_PSI_clustered.tsv'):
        frame = pd.read_csv(source, sep='\t', index_col=0)
        if frame.empty:
            continue
        meta = pd.read_csv(source.parent / 'event_annotations.tsv', sep='\t', index_col=0)
        meta = meta.loc[~meta.index.duplicated()]
        cfg = json.loads((source.parent / 'parameters.json').read_text())
        out = base / '24_Narrative_report/figures' / source.parent.parent.name / source.parent.name
        out.mkdir(parents=True, exist_ok=True)
        for mode in ('PSI', 'rowZ'):
            values = frame.to_numpy(float)
            if mode == 'rowZ':
                values = row_scaled(values)
            with plt.rc_context({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False}):
                fig = plt.figure(figsize=(10.5, 12.5), layout='constrained')
                gs = fig.add_gridspec(3, 3, height_ratios=[1, 0.22, 9], width_ratios=[2, 20, 1])
                top = fig.add_subplot(gs[0, 1])
                linkage = source.parent / 'sample_linkage.tsv'
                if linkage.exists():
                    dendrogram(pd.read_csv(linkage, sep='\t').to_numpy(), ax=top, no_labels=True, above_threshold_color='#45535b', color_threshold=0)
                top.axis('off')
                strip = fig.add_subplot(gs[1, 1])
                groups = [str(s).rsplit('_rep', 1)[0] for s in frame.columns]
                codes = {g: i for i, g in enumerate(dict.fromkeys(groups))}
                strip.imshow([[codes[g] for g in groups]], aspect='auto', cmap='Set2', vmin=0, vmax=max(1, len(codes) - 1))
                strip.axis('off')
                for group in codes:
                    positions = [i for i, g in enumerate(groups) if g == group]
                    strip.text(np.mean(positions), 0, group, ha='center', va='center', fontsize=8)
                rowtree = fig.add_subplot(gs[2, 0])
                rowlink = source.parent / 'event_linkage.tsv'
                if rowlink.exists():
                    dendrogram(pd.read_csv(rowlink, sep='\t').to_numpy(), ax=rowtree, orientation='left', no_labels=True, above_threshold_color='#45535b', color_threshold=0)
                    rowtree.invert_yaxis()
                rowtree.axis('off')
                ax = fig.add_subplot(gs[2, 1])
                cax = fig.add_subplot(gs[2, 2])
                cmap = pastel_diverging()
                lim = (0, 1) if mode == 'PSI' else (-2.5, 2.5)
                im = ax.imshow(values, aspect='auto', interpolation='nearest', cmap=cmap, vmin=lim[0], vmax=lim[1])
                labels = []
                for i, event in enumerate(frame.index):
                    row = meta.loc[event] if event in meta.index else {}
                    labels.append(f"{i + 1:03d}  {row.get('gene_symbol', '')} | {row.get('event_type', '')}")
                ax.set_yticks(range(len(frame)), labels, fontsize=6.5)
                ax.set_xticks(range(len(frame.columns)), frame.columns, rotation=90, fontsize=9)
                ax.tick_params(axis='both', length=0, pad=4)
                ax.set_ylabel('Displayed events (full exact IDs in source matrix)')
                fig.colorbar(im, cax=cax, label='PSI (0-1)' if mode == 'PSI' else 'Within-event z score (display clipped at +/-2.5)')
                fig.suptitle(source.parent.parent.name.replace('_vs_', ' vs ') + '\n' + ('Raw inclusion levels' if mode == 'PSI' else 'Relative sample pattern, not raw PSI') + f" | {len(frame)} displayed / {cfg['eligible_events']} eligible\n" + source.parent.name.replace('_', ' ') + '; existing cluster order; gray = missing', fontsize=12)
                for ext in ('png', 'svg', 'pdf'):
                    fig.savefig(out / f'{mode}.{ext}', dpi=220)
                plt.close(fig)
        (out / 'display_policy.json').write_text(json.dumps(dict(cfg, render_only=True, source=str(source), rowZ='sample standard deviation ddof=1; no imputation in display'), indent=2))
    for source in (base / '12_Gene_cards_heatmaps').rglob('selected_gene_psi.tsv'):
        frame = pd.read_csv(source, sep='\t').set_index('event_id')
        display_ids = frame['display_event_id']
        frame = frame.drop(columns=['row_id', 'display_event_id'], errors='ignore')
        config = json.loads((source.parent / 'parameters.json').read_text()).get('parameters', {})
        step = int(config.get('rows_per_page', 40))
        events = pd.read_csv(source.parent / 'selected_gene_events.tsv', sep='\t').drop_duplicates('event_id').set_index('event_id')
        out = base / '24_Narrative_report/figures/selected_genes' / source.parent.name
        out.mkdir(parents=True, exist_ok=True)
        for start in range(0, len(frame), step):
            chunk = frame.iloc[start:start + step]
            fig, ax = plt.subplots(figsize=(11, max(5, 0.22 * len(chunk) + 2.5)), layout='constrained')
            art = ax.imshow(row_scaled(chunk.to_numpy(float)), cmap=pastel_diverging(), vmin=-2.5, vmax=2.5, aspect='auto', interpolation='nearest')
            ax.set_xticks(range(len(chunk.columns)), chunk.columns, rotation=90, fontsize=8)
            labels = [str(events.loc[x, 'gene_symbol']) + ' | ' + str(events.loc[x, 'event_type']) + ' | ' + str(display_ids.loc[x]) for x in chunk.index]
            ax.set_yticks(range(len(chunk)), labels, fontsize=7)
            ax.set_title(f'Selected-gene AS | row-scaled PSI | page {start // step + 1}')
            fig.colorbar(art, ax=ax, label='Row z score; gray = missing; display limits +/-2.5')
            for ext in ('png', 'svg', 'pdf'):
                fig.savefig(out / f'selected_gene_heatmap_{start // step + 1:03d}.{ext}', dpi=220)
            plt.close(fig)
MODULES = [('02_AS_reliability', '输入可信度与筛选损耗', 'Reliability and filtering', '哪些事件可测量、哪些有重复支持？筛掉的事件不是已证明无效。', 'Which events are measurable and replicate-supported? Exclusion is not evidence of no effect.'), ('03_AS_landscape', '全局PCA与样本相似度', 'Global PCA and sample similarity', '观察全部样本的全局关系；它与上方单个comparison筛选后的热图回答不同问题。', 'Global sample relationships differ from the contrast-selected heatmaps above.'), ('04_DAS', '剪接事件分布与差异', 'AS event distribution and differential splicing', '比较SE、RI、A3SS、A5SS、MXE的计数和效应分布；事件数不等于基因数。', 'Compare event-type counts and effect distributions; event counts are not gene counts.'), ('05_DEG', '基因表达响应', 'Gene expression response', '总表达变化与剪接比例变化分别解释；不显著不等于不变。', 'Abundance and splicing proportion are distinct; non-significance is not equivalence.'), ('07_Cross_comparison', '精确事件交叉、包含度与方向', 'Exact-event overlap, containment and direction', '按精确坐标识别同一事件。以下交叉采用原输出的stable-significant集合；重合本身不要求方向一致，方向需单独核对。', 'Exact coordinates define matching events. Existing overlap uses the stable-significant set; membership alone does not require concordant direction.'), ('08_Candidate_genes', '候选基因证据与排名', 'Candidate evidence and ranking', '排名依次使用支持维度数、稳定显著AS数、最大绝对dPSI，不是因果概率。', 'Ranking uses supporting dimensions, stable-significant event count and maximum absolute dPSI; it is not a causal probability.'), ('12_Gene_cards_heatmaps', 'Gene cards：逐基因表达及剪接证据', 'Gene cards: expression and splicing per gene', '每张卡对应原来选中的基因；保留所有已生成的卡片，不只挑阳性结果。完整选中基因审计表记录未能作图的原因。', 'Cards retain the original selected genes, not just positive findings. The selection audit records genes unavailable for plotting.'), ('06_DTU_DEXSeq', '转录本及外显子使用', 'Transcript and exon usage', '这是普通两组使用差异，不是基因型×处理交互；完整未显著结果也保留。', 'These are ordinary usage contrasts, not genotype-by-treatment interactions; non-significant results are retained.'), ('09_GO_KEGG', 'GO与KEGG剪接功能解释', 'GO and KEGG interpretation', '检查原候选集合与指定背景的富集，不能与表达排序富集混作同一个检验。', 'Enrichment uses the original selected set and background, separate from ranked expression enrichment.'), ('10_RBP', '第一层 RBP身份注释', 'Layer 1 RBP identity annotations', '保留原候选基因身份/声明靶集合富集；不是事件位点结合、motif变化或上游因果。', 'Retains original gene identity and declared target-set enrichment; not event-site binding, motif change or upstream causality.'), ('11_Isoform_usage', '异构体组成', 'Isoform composition', '转录本使用改变不自动代表蛋白异构体改变。', 'Transcript usage does not directly establish protein isoform changes.'), ('13_Isoform_consequences', '异构体功能后果', 'Isoform consequences', '编码、结构域和NMD注释形成预测，不替代蛋白验证。', 'Coding, domain and NMD annotations are predictions, not protein validation.'), ('14_RBP_motifs', '第二层 上游motif/CLIP调控线索', 'Layer 2 Upstream motif / CLIP regulatory clues', '邻近motif和匹配CLIP支持候选调控假说；不等于DAS实际改变的成熟RNA序列，也不证明当前样本结合或因果。', 'Neighbourhood motifs and compatible reference CLIP support regulatory hypotheses; not mature-RNA sequence consequences or proof of current-sample binding or causality.'), ('15_CoSplice', 'CoSplice完整结果：模块、条件效应与稳定性', 'CoSplice results: modules, condition effects and stability', '先看模块大小和成员，再看条件关联及效应，最后核对bootstrap可信度与功能支持。初始网络M编号与共识C编号是不同阶段；下列原生图表保留各阶段身份，不把M模块的P值套给C模块。关联显著不等于稳定，也不证明因果。', 'Read membership and size, then condition associations and effects, then bootstrap confidence and functional evidence. Initial M identifiers and consensus C identifiers refer to different stages; original stage identities are preserved, and M-module P-values are not assigned to C modules. Association significance is not stability or causality.'), ('23_Module_explanation', '共识模块功能、hub与表达背景', 'Consensus module function, hubs and expression context', '这些结果使用固定共识C模块；功能命名只在富集FDR支持时给出。总体相关与条件调整后的相关分开报告。', 'These results use fixed consensus C modules. Functional labels require enrichment FDR support. Overall and condition-adjusted correlations are reported separately.')]

def event_distribution(base, output):
    import matplotlib.pyplot as plt
    paths = list((Path(base) / '04_DAS').rglob('differential_as_counts.vl.json'))
    if not paths:
        return None
    source = max(paths, key=lambda p: p.stat().st_mtime)
    frame = pd.DataFrame(json.loads(source.read_text())['data']['values'])
    out = Path(output) / 'figures'
    out.mkdir(exist_ok=True)
    frame.to_csv(out / 'AS_type_distribution_source.tsv', sep='\t', index=False)
    colors = ['#3577a5', '#db6846', '#55a293', '#ba9450', '#8e7897']
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5), layout='constrained')
    for ax, key, title in zip(axes, ['tested_events', 'significant_events', 'stable_significant_events'], ['All tested events', 'Significant events', 'Stable significant events']):
        table = frame.pivot(index='comparison_id', columns='event_type', values=key).fillna(0).reindex(columns=['SE', 'RI', 'A3SS', 'A5SS', 'MXE'], fill_value=0)
        left = np.zeros(len(table))
        for name, color in zip(table.columns, colors):
            ax.barh(range(len(table)), table[name], left=left, color=color, label=name, height=0.65)
            left += table[name].to_numpy()
        ax.set_yticks(range(len(table)), [s.replace('_vs_', ' vs ') for s in table.index], fontsize=8)
        ax.invert_yaxis()
        ax.set_title(title)
        ax.set_xlabel('Number of exact AS events')
        for i, n in enumerate(left):
            ax.text(n, i, f' {int(n):,}', va='center', fontsize=8)
        ax.set_xlim(0, max(left.max(), 1) * 1.2)
        ax.spines[['top', 'right']].set_visible(False)
    axes[-1].legend(ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.14), frameon=False)
    for ext in ('png', 'svg', 'pdf'):
        fig.savefig(out / f'AS_type_distribution.{ext}', dpi=220)
    plt.close(fig)
    return out / 'AS_type_distribution.png'

def complete_sections(base, output, language):
    if str(language).upper() != 'EN':
        raise ValueError('Only English reports are supported')
    from katana_cosplice.branding import display_module_name
    base = Path(base)
    output = Path(output)
    e = lambda value, quote=True: html.escape(display_module_name(str(value)), quote=quote)

    def link(p):
        return e(Path(os.path.relpath(p, output)).as_posix(), quote=True)
    from katana_next.science.rbp_evidence_layers import overview_html
    sections = [overview_html(language=language, inventory=False)]
    audit = []
    distribution = event_distribution(base, output)
    for folder, title_zh, title_en, why_zh, why_en in MODULES:
        root = base / folder
        plots = {}
        for p in sorted(root.rglob('*.png'), key=lambda p: p.stat().st_mtime):
            if 'Native_PSI_clustering' not in p.parts:
                if folder == '15_CoSplice' and 'benchmark' in str(p):
                    continue
                alternative = output / 'figures/selected_genes' / p.parent.parent.name / p.name
                plots[p.name] = alternative if folder == '12_Gene_cards_heatmaps' and p.name.startswith('selected_gene_heatmap_') and alternative.is_file() else p
        tables = {}
        for p in sorted(root.rglob('*.tsv'), key=lambda p: p.stat().st_mtime):
            if folder == '15_CoSplice' and (not any((word in p.name for word in ('association', 'coefficient', 'module_summary', 'module_assignments', 'hub', 'enrichment', 'confidence', 'diagnostics', 'group_summary')))):
                continue
            if 'Native_PSI_clustering' not in p.parts:
                tables[p.name] = p
        body = '<section id="module-' + folder + '"><h2>' + e(title_en) + '</h2><p>' + e(why_en) + '</p>'
        if folder == '04_DAS' and distribution and distribution.is_file():
            body += '<figure><a href="' + link(distribution) + '"><img src="' + link(distribution) + '"></a><figcaption>' + 'Tested, significant and stable-significant sets counted separately; original thresholds unchanged.' + '</figcaption></figure>'
        if folder == '06_DTU_DEXSeq' and 'dtu_features.tsv' in tables:
            frame = pd.read_csv(tables['dtu_features.tsv'], sep='\t')
            for kind, title in [('dexseq_result', 'DEXSeq: exon usage'), ('dtu_result', 'DTU: transcript usage')]:
                part = frame.loc[frame.source_kind.eq(kind)].copy()
                body += '<h3>' + title + '</h3>'
                if part.empty:
                    body += '<p>' + 'No input result of this kind; not zero significance.' + '</p>'
                    continue
                summary = part.groupby('comparison_id').agg(tested_features=('dtu_fdr', lambda x: x.notna().sum()), significant_features=('significant', 'sum')).reset_index()
                body += '<div class="scroll">' + summary.to_html(index=False, escape=True) + '</div>'
                passing = part[part.significant.astype(str).str.lower().isin(['true', '1'])].sort_values('dtu_fdr')
                body += '<p>' + 'Top 12 passing features by the original adjusted P-value; complete significant and non-significant results are linked below. Feature-level significance is not gene-level FDR.' + '</p><div class="scroll">' + passing.drop(columns=['source_path', 'row_id'], errors='ignore').head(12).to_html(index=False, escape=True) + '</div>'
        if not plots:
            body += '<p>' + 'No existing figure for this module; inspect tables or missing-input records below. This is not a zero result.' + '</p>'
        for i, (name, p) in enumerate(plots.items()):
            body += '<figure><a href="' + link(p) + '"><img loading="lazy" src="' + link(p) + '" alt="' + e(name) + '"></a><figcaption>' + e(name.replace('_', ' ')) + ' · <a href="' + link(p) + '">' + 'Open full size' + '</a></figcaption></figure>'
        body += '<details><summary>' + 'Full tables and result previews' + '</summary>'
        for name, p in tables.items():
            body += '<h3><a href="' + link(p) + '">' + e(name) + '</a></h3>'
            try:
                frame = pd.read_csv(p, sep='\t', nrows=8).drop(columns=['row_id'], errors='ignore')
                body += '<p>' + 'Preview: first 8 rows in source order, not a significance filter.' + '</p><div class="scroll">' + frame.to_html(index=False, escape=True) + '</div>'
            except (pd.errors.EmptyDataError, UnicodeError):
                body += '<p>Empty / unavailable table</p>'
        body += '</details></section>'
        sections.append(body)
        audit.append(dict(module=folder, figures=[str(p) for p in plots.values()], tables=[str(p) for p in tables.values()]))
    (output / 'MODULE_COVERAGE.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
    return ''.join(sections)
