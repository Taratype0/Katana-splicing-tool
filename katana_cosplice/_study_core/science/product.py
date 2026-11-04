"""End-to-end analysis and portable report, used by the desktop/API worker."""
import json,sys,time,traceback,html,hashlib,zipfile,shutil,csv,gzip
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from cosplice_research import ROOT,prepare,fit,bootstrap,write_json
from advanced import fit_advanced,nearest_edges,OASNetwork,graph_outputs

def read_matrix(path):
    p=Path(path);sep=',' if p.suffix=='.csv' else '\t'
    opener=gzip.open if p.suffix=='.gz' else open
    with opener(p,'rt',encoding='utf-8-sig',newline='') as handle:header=next(csv.reader(handle,delimiter=sep))
    if len(header)<2 or len(set(header[1:]))!=len(header)-1 or any(not s for s in header[1:]):raise ValueError('Duplicate or empty sample headers')
    raw=pd.read_csv(p,sep=sep,dtype=str)
    if raw.iloc[:,0].isna().any():raise ValueError('Missing event IDs')
    tab=raw.iloc[:,1:].copy();tab.index=pd.Index(raw.iloc[:,0].to_numpy(),dtype=object)
    tab.index=tab.index.astype(str);tab.columns=tab.columns.astype(str)
    if not tab.index.is_unique or not tab.columns.is_unique:raise ValueError('Duplicate event or sample IDs')
    return tab.astype(float)

def raw_inputs(cfg):
    meta=pd.read_csv(cfg['samples'],sep=',' if str(cfg['samples']).endswith('.csv') else '\t',dtype={'sample_id':str,'biological_sample_id':str,'condition':str})
    from katana_cosplice.study_design import validate_study_design
    validate_study_design(meta,paired=cfg.get('paired',False))
    if cfg.get('counts'):
        c=read_matrix(cfg['counts']);d=read_matrix(cfg['totals'])
        if not c.index.equals(d.index) or not c.columns.equals(d.columns):raise ValueError('Counts and totals identities/order differ')
        for tab in [c,d]:
            v=tab.to_numpy();finite=v[np.isfinite(v)]
            if np.isinf(v).any() or (finite<0).any() or (finite%1!=0).any():raise ValueError('Counts must be nonnegative integers or NA')
        if (c>d).any().any():raise ValueError('Numerator exceeds total')
        x=c/d.where(d>0)
    else:x=read_matrix(cfg['psi'])
    if cfg.get('counts') and cfg.get('value_kind')=='response':raise ValueError('Response input requires delta PSI, not counts')
    support=d if cfg.get('counts') else read_matrix(cfg['read_support']) if cfg.get('read_support') else None
    return x,meta,support

def load_inputs(cfg):
    x,meta,support=raw_inputs(cfg)
    measured,complete,audit,means,meta=prepare(x,meta,value_kind=cfg.get('value_kind','usage'))
    if cfg.get('qc_profile','legacy')=='unit_robust':
        if cfg.get('value_kind','usage')!='usage':raise ValueError('Unit influence QC requires individual sample usage')
        from unit_robust import unit_audit
        masked,audit=unit_audit(x,meta,support);measured=masked.loc[audit.eligible];means=measured.mean(axis=1);complete=measured.T.fillna(means).T
        from audit_labels24 import explain_audit
        audit=explain_audit(audit)
    if cfg.get('feature_cap') is not None:
        from feature_selection36 import METHODS,select_events
        if cfg.get('method') not in METHODS:raise ValueError('Event limits currently apply to WGCNA engines')
        measured,audit=select_events(measured,audit,cfg['feature_cap'])
        means=measured.mean(axis=1);complete=measured.T.fillna(means).T
    else:
        audit['qc_eligible']=audit.eligible; audit['network_selected']=audit.eligible
        audit['selection_reason']=np.where(audit.eligible,'selected','failed_qc')
    if len(complete)<40:raise ValueError('Fewer than 40 variable, sufficiently observed events after QC')
    genes=None
    if cfg.get('annotation'):
        a=pd.read_csv(cfg['annotation'],sep='\t',dtype=str)
        if not {'event_key','gene_id'}.issubset(a) or not a.event_key.is_unique:raise ValueError('Annotation requires unique event_key and gene_id')
        genes=a.set_index('event_key').gene_id.reindex(complete.index)
    expression=read_matrix(cfg['expression']) if cfg.get('expression') else None
    return measured,complete,audit,means,meta,genes,expression

def bh(p):
    p=np.asarray(p);order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];return np.minimum(q,1)

def pair_corr(x,left,right,cols):
    a=x[left][:,cols];b=x[right][:,cols]
    # Fixed training imputation used for descriptive network rewiring; logged.
    a=a-a.mean(axis=1,keepdims=True);b=b-b.mean(axis=1,keepdims=True)
    den=np.sqrt((a*a).sum(axis=1)*(b*b).sum(axis=1))
    return np.divide((a*b).sum(axis=1),den,out=np.full(len(left),np.nan),where=den>1e-12)

def differential_edges(complete,meta,edges,out,condition_a,condition_b,permutations=199):
    if condition_a==condition_b:raise ValueError('Choose two distinct conditions')
    selected=meta.condition.isin([condition_a,condition_b]);m=meta[selected].reset_index(drop=True)
    if min(m.groupby('condition').biological_sample_id.nunique())<8:raise ValueError('At least 8 biological units per condition required')
    if m.condition.nunique()!=2:raise ValueError('Contrast condition absent')
    counts=m.groupby('biological_sample_id').condition.nunique()
    paired=counts.gt(1).any()
    if paired and not counts.eq(2).all():raise ValueError('Partially paired contrasts require an explicit model; not treated as independent')
    if paired and m.groupby('biological_sample_id').size().ne(2).any():raise ValueError('Pairing requires exactly one observation per condition per unit')
    keys={k:i for i,k in enumerate(complete.index)}
    tab=edges[edges.source.isin(keys)&edges.target.isin(keys)].copy()
    tab['pair']=tab.apply(lambda r:tuple(sorted([r.source,r.target])),axis=1);tab=tab.drop_duplicates('pair')
    x=complete.loc[:,m.sample_id].to_numpy();left=np.array([keys[k] for k in tab.source]);right=np.array([keys[k] for k in tab.target])
    labels=m.condition.eq(condition_a).to_numpy();rng=np.random.default_rng(20261030)
    permutation_labels=[]
    for _ in range(permutations):
        if paired:
            l=labels.copy()
            for _,ids in m.groupby('biological_sample_id').groups.items():
                if rng.integers(2):l[np.asarray(list(ids))]=~l[np.asarray(list(ids))]
        else:l=rng.permutation(labels)
        permutation_labels.append(l)
    rs=[]
    for start in range(0,len(tab),2000):
        a=left[start:start+2000];b=right[start:start+2000]
        r_a=pair_corr(x,a,b,labels);r_b=pair_corr(x,a,b,~labels);delta=r_a-r_b
        exceed=np.zeros(len(a));valid=np.isfinite(delta)
        for l in permutation_labels:
            v=pair_corr(x,a,b,l)-pair_corr(x,a,b,~l)
            # A nonestimable permutation cannot artificially reduce p-values.
            exceed+=(~np.isfinite(v))|(np.abs(v)>=np.abs(delta))
        p=(1+exceed)/(permutations+1);p[~valid]=np.nan
        rs.extend(zip(r_a,r_b,delta,p))
    for i,name in enumerate(['r_a','r_b','delta_r','p_value']):tab[name]=[r[i] for r in rs]
    tab['q_value']=np.nan;ok=tab.p_value.notna();tab.loc[ok,'q_value']=bh(tab.loc[ok,'p_value'])
    tab.drop(columns='pair').sort_values('p_value').to_csv(out/'differential_edges.tsv',sep='\t',index=False)
    return {'a':condition_a,'b':condition_b,'paired':bool(paired),'units':int(m.biological_sample_id.nunique()),'permutations':permutations,
        'tested_edges':int(ok.sum()),'nonestimable_edges':int((~ok).sum()),'q05_edges':int(tab.q_value.lt(.05).sum()),
        'scope':'pooled-network candidate edges; Pearson rewiring on train-mean imputed usage; exchangeability assumed, no covariate adjustment'}

def write_report(out,summary):
    if summary.get('study'):
        from study_report import build
        build(out,out,summary=summary)
        # The default interpretation is independently portable and has the same
        # final QC/stability context as the main study report.
        interpretation=out/'interpretations/default'
        build(out,interpretation,interpretation,summary)
        with zipfile.ZipFile(interpretation/'results.zip','w',zipfile.ZIP_DEFLATED) as z:
            for p in interpretation.rglob('*'):
                if p.is_file() and p.name!='results.zip':z.write(p,p.relative_to(interpretation))
        with zipfile.ZipFile(out/'results.zip','w',zipfile.ZIP_DEFLATED) as z:
            for p in out.rglob('*'):
                if p.is_file() and p.name not in ['results.zip','status.json','worker.json'] and 'refit_' not in str(p):z.write(p,p.relative_to(out))
        return
    tables=[]
    for name,title in [('modules.tsv','模块概览'),('pathway_enrichment.tsv','模块通路关联（全部模块×基因集共同 BH）'),('model/selection.tsv','训练内交叉验证集成权重'),('model/factor_variance.tsv','各视图因子解释度'),('events.tsv','事件审计与归属'),('model/loadings.tsv','因子载荷'),('differential_edges.tsv','条件间关联变化')]:
        p=out/name
        if p.exists():
            tab=pd.read_csv(p,sep='\t')
            if name=='pathway_enrichment.tsv':tab=tab.sort_values('q_value',kind='stable')
            tables.append(f'<h2>{title}</h2><p>共 {len(tab):,} 行，以下展示前 100 行；完整结果位于下载包。</p>'+tab.head(100).to_html(index=False,escape=True,border=0))
    safe=json.dumps(summary,ensure_ascii=False,indent=2)
    page='''<!doctype html><html lang="zh"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>CoSplice Next 分析报告</title><style>body{font:15px/1.65 system-ui;margin:40px auto;max-width:1150px;padding:0 25px;color:#233a44;background:#f9faf8}h1,h2{color:#12575b}table{border-collapse:collapse;display:block;overflow:auto;font-size:12px}td,th{padding:7px;border-bottom:1px solid #d9e1df;text-align:left}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:white;padding:18px;border:1px solid #d9e1df}aside{padding:16px;background:#fff1d6}</style><h1>CoSplice Next · 分析报告</h1><aside>网络关联、预测重要性和因子载荷各有不同含义。它们不直接证明因果调控，也不自动构成已验证 hub。未归属、缺失与失败保留在输出中。</aside>'''
    def bars(labels,values,title):
        if not len(labels):return ''
        height=40+len(labels)*25;maximum=max(max(values),1e-9);svg=[]
        for i,(label,v) in enumerate(zip(labels,values)):
            svg.append(f'<text x="2" y="{30+i*25}" font-size="12">{html.escape(str(label))}</text><rect x="170" y="{17+i*25}" height="17" width="{450*v/maximum:.2f}" fill="#398475"/><text x="{175+450*v/maximum:.2f}" y="{30+i*25}" font-size="11">{v:.3g}</text>')
        return f'<h2>{title}</h2><svg role="img" aria-label="{title}" viewBox="0 0 720 {height}">'+''.join(svg)+'</svg>'
    if 'eligible_events' in summary:
        page+='<div style="display:flex;gap:35px;margin:25px 0;flex-wrap:wrap">'+''.join(f'<div><strong style="font-size:28px">{summary[k]:,}</strong><br>{label}</div>' for k,label in [('input_events','输入事件'),('eligible_events','质控合格'),('assigned_events','模块归属'),('biological_units','干预目标轮廓' if summary.get('value_kind')=='response' else '生物学单位')])+'</div>'
        page+=f'<p>网络方法：<strong>{html.escape(summary["method"])}</strong>。实际建网 {summary.get("network_events",summary["eligible_events"]):,} 个事件；上限 {summary.get("feature_cap") or "不限量"}。质控排除和上限未选入分别记录。</p>'
        if summary.get('value_kind')=='response':page+='<p>本次输入是扰动 ΔPSI 百分点；列表示干预响应轮廓。不是样本 PSI、计数或独立患者。</p>'
        if summary['method']=='consensus':page+='<p>内部 24 次单位子抽样用于建立共识图；内部共聚类支持不能当作独立稳定性。下方只有外层重新构图、重新聚类才属于重拟合验证。</p>'
        if summary.get('qc_profile')=='unit_robust':
            page+='<h2>单单位影响与读数质控</h2><p>每个生物学单位总权重为 1；排除超过或等于 50% 变异由单个单位贡献的事件。所有排除原因保留在事件审计表，所有条件的逐样本观测共同建网。</p>'
            page+='<p>'+('已按逐样本有效支持读数 ≥10 掩蔽低覆盖测量。' if summary.get('read_support_available') else '未提供读数支持矩阵，未执行读数覆盖质控；仅根据已观测 PSI 检查缺失、变异与单位影响。')+'</p>'
        if summary['method']=='jackknife':page+='<p>每条候选边逐一删除生物学单位的全部观测，取 Spearman 与 Pearson 相关的最小值构图。该步骤控制构建时的单单位依赖，不能替代下方的外层完整重拟合与独立生物学验证。</p>'
        if summary.get('merge_modules'):page+='<p>已启用原生 WGCNA 模块特征向量合并（cutHeight=0.25），保留合并前归属；外层每次重拟合也执行相同合并。这是借鉴标准方法的探索性扩展，尚未证明能提高独立数据上的稳定性。</p>'
        if not summary.get('stability'):page+='<aside>本次尚未附入外层完整重拟合，模块稳定性未验证。</aside>'
        if summary.get('pathway_enrichment'):
            f=summary['pathway_enrichment'];page+=f'<h2>功能解释</h2><p>背景：实际入网事件（含未分配事件）对应的 {f["background_genes"]} 个唯一明确基因，再取数据库可注释部分；{f["tested_modules"]} 个可检验模块，共 {f["tests"]} 次基因集检验。至少一个条目 q≤0.05 的模块为 {f["significant_modules"]} 个。相同基因不因多个剪接事件而重复计数；这不是调控边或剪接机制的实验验证。基因集：{html.escape(f["gene_set_collection"])}。</p>'
    if (out/'modules.tsv').exists():
        m=pd.read_csv(out/'modules.tsv',sep='\t');m=m[m.module.ne(0)].nlargest(20,'events')
        page+=bars(['M'+str(int(x)) for x in m.module],m.events.tolist(),'模块规模（最多显示 20 个）')
    if (out/'model/selection.tsv').exists():
        w=pd.read_csv(out/'model/selection.tsv',sep='\t');page+=bars(w.expert.tolist(),w.weight.tolist(),'训练内学习的集成权重')
        page+='<p>权重仅由训练内个体划分学习。以下训练内误差用于选择权重，不作为独立验证结果。</p>'
    if (out/'stability/stability.tsv').exists():
        b=pd.read_csv(out/'stability/stability.tsv',sep='\t');valid=b.within_module_coassignment.dropna()
        page+=f'<h2>模块重采样</h2><p>请求 {int(b.requested_replicates.iloc[0])} 次，完成 {int(b.successful_replicates.iloc[0])} 次；已归属事件支持度中位数 {valid.median():.3f}。这不是 hub 置信度或 P 值。</p>'
        if len(valid) and valid.median()<.5:page+='<aside>模块归属在重采样中变化较大；分析完成不代表模块稳定。请结合逐事件支持度使用结果。</aside>'
        if summary.get('stability',{}).get('development_replicates'):
            scope=summary['stability'];page+=f'<p>其中前 {scope["development_replicates"]} 次用于开发。未用于调参的后 {scope["confirmation_replicates"]} 次内部确认支持度中位数为 {scope["confirmation_summary"]["event_support_median"]:.3f}；仍属于同一队列，不是外部人群复现。</p>'
        stats=summary.get('stability',{}).get('partition_metrics')
        if stats:page+='<p>覆盖加权支持 '+str(round(stats['coverage_weighted_support'],3))+'；共成员 pair-F1 '+str(round(stats['pair_f1'],3))+'；未归属作为单例的 ARI '+str(round(stats['singleton_grey_ari'],3))+'。同时检查覆盖率与模块大小，避免把大片未归属误算为稳定模块。</p>'
    if summary.get('differential'):
        d=summary['differential'];page+=f'<h2>条件关联变化</h2><p>{html.escape(d["a"])} 对比 {html.escape(d["b"])}：{d["tested_edges"]:,} 条可估计边，{d["q05_edges"]:,} 条 q&lt;0.05；{d["nonestimable_edges"]} 条不可估计边保留为 NA。</p>'
    if summary.get('study'):
        s=summary['study'];a=summary['study_interpretation']
        page+='<h2>研究级共同网络</h2><p>'+html.escape(s['study_id'])+' · '+str(s['samples'])+' 个逐样本观测 · '+str(s['biological_units'])+' 个独立单位；'+str(len(s['conditions']))+' 个条件共同建网。</p><p>固定网络 ID：'+s['network_id']+'</p><p>默认模块—条件分析：'+str(a['tests_estimable'])+'/'+str(a['tests_requested'])+' 个可估计检验；全部对比复用固定模块。</p><p><a href="interpretations/default/report.html">查看逐模块条件效应、配对设计、NA 与方差诊断</a>（结果 ZIP 内可直接打开）。</p><aside>共同网络反映跨条件共同响应，不证明条件内部共调控。模块—条件 P 值是同研究探索性结果，不能替代外部生物学验证。</aside>'
    page+=''.join(tables)+'<details><summary>完整参数、来源与解释边界</summary><pre>'+html.escape(safe)+'</pre></details></html>'
    (out/'report.html').write_text(page,encoding='utf-8')
    with zipfile.ZipFile(out/'results.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in out.rglob('*'):
            if p.is_file() and p.name not in ['results.zip','status.json','worker.json'] and 'refit_' not in str(p):z.write(p,p.relative_to(out))

def run(cfg,out):
    if cfg.get('method')=='factor_adaptive':
        from product32 import run as run_adaptive
        return run_adaptive(cfg,out)
    out=Path(out);status={'status':'running','stage':'input','started':time.time(),'method':cfg['method']};write_json(out/'status.json',status)
    try:
        assembly=None
        if cfg.get('value_kind','usage')=='usage' and cfg['method'] not in ['joint','independent'] and (cfg.get('psi') or cfg.get('psi_sources')):
            from study_inputs import assemble
            cfg,assembly=assemble(cfg,out/'inputs')
        measured,complete,audit,means,meta,genes,expression=load_inputs(cfg)
        study_mode=cfg.get('value_kind','usage')=='usage' and cfg['method'] not in ['joint','independent']
        if study_mode:
            from study import validate_contract
            validate_contract(cfg,meta)
        audit.to_csv(out/'event_audit.tsv',sep='\t',index_label='event_key');measured.to_csv(out/'usage.tsv.gz',sep='\t',index_label='event_key')
        meta.to_csv(out/'samples.tsv',sep='\t',index=False)
        status.update(stage='network',events=len(complete),samples=complete.shape[1]);write_json(out/'status.json',status)
        method=cfg['method'];model=out/'model'
        if cfg.get('revision22_cache'):
            source=Path(cfg['revision22_cache']).resolve();allowed=(ROOT/'evidence/revision22').resolve()
            if not source.is_relative_to(allowed):raise ValueError('Revision cache must belong to the versioned validation directory')
            original=pd.read_csv(source/'input.tsv',sep='\t',index_col=0)
            original_meta=pd.read_csv(source/'samples.tsv',sep='\t')
            if not original.index.equals(complete.index) or not original.columns.equals(complete.columns) or not np.allclose(original,complete,atol=1e-12):raise ValueError('Revision cache input differs')
            fields=['sample_id','biological_sample_id','condition']
            if not original_meta[fields].equals(meta[fields]):raise ValueError('Revision cache sample identities differ')
            result=json.loads((source/'status.json').read_text())
            if result['status']!='completed' or result['method']!=method:raise ValueError('Revision cache engine mismatch')
            shutil.copytree(source,model);result['cache_provenance']={'source':str(source),'matrix_and_units_verified':True,'refitted':False,
                'assignment_sha256':hashlib.file_digest((source/'assignments.tsv').open('rb'),'sha256').hexdigest()}
        elif cfg.get('validated_model'):
            source=Path(cfg['validated_model']).resolve()
            if not source.is_relative_to((ROOT/'evidence/comparison').resolve()):raise ValueError('Only this release verified comparison cache may be imported')
            original=pd.read_csv(source/'wgcna/input.tsv',sep='\t',index_col=0)
            if not original.index.equals(complete.index) or not original.columns.equals(complete.columns) or not np.allclose(original,complete,atol=1e-12):raise ValueError('Cached model training matrix differs')
            original_meta=pd.read_csv(source/'wgcna/samples.tsv',sep='\t').set_index('sample_id')
            current_meta=meta.set_index('sample_id')
            if not original_meta.loc[complete.columns,['biological_sample_id','condition']].equals(current_meta.loc[complete.columns,['biological_sample_id','condition']]):raise ValueError('Cached model biological metadata differs')
            result=json.loads((source/'status.json').read_text())
            if result['status']!='completed' or result['method']!=method:raise ValueError('Cached model does not match requested engine')
            shutil.copytree(source,model)
            result['cache_provenance']={'source':str(source),'matrix_verified':True,'refitted':False}
        elif method=='stacking':
            from stacking import train_stacking
            if genes is None:raise ValueError('Stacking requires event-to-gene annotation for cross-gene prediction')
            result=train_stacking(measured,meta,genes.fillna('UNKNOWN'),model)
        elif method=='global_local':
            if cfg.get('qc_profile')!='unit_robust' or cfg.get('value_kind','usage')!='usage':raise ValueError('Global/local requires individual-sample usage and unit_robust QC')
            from algorithm31 import fit_global_local
            result=fit_global_local(measured,complete,meta,model)
        elif method=='factor_core':
            if cfg.get('qc_profile')!='unit_robust' or cfg.get('value_kind','usage')!='usage':raise ValueError('Factor core requires individual-sample usage and unit_robust QC')
            if meta.biological_sample_id.duplicated().any():raise ValueError('Factor core currently requires independent biological units')
            from algorithm28 import fit_latent
            result=fit_latent(method,measured,complete,meta,model)
        elif method=='jackknife':
            if cfg.get('qc_profile')!='unit_robust':raise ValueError('Jackknife requires unit_robust QC')
            from unit_robust import fit_jackknife
            result=fit_jackknife(measured,complete,meta,model)
        elif method in ['consensus','leiden']:
            from stable_network import fit_stable
            result=fit_stable(measured,complete,meta,model,consensus=method=='consensus')
        elif method in ['factor','megena','wgcna_robust']:
            from validation22 import fit_partition
            result=fit_partition(method,measured,complete,meta,model)
        else:result=fit_advanced(method,measured,complete,meta,model,expression) if method in ['oas','mofa','grnboost2'] else fit(method,complete,meta,model)
        if result['status']!='completed':raise RuntimeError(result.get('error','Model failed; see model/engine.log'))
        if not (model/'assignments.tsv').exists():
            if method in ['factor','pca']:
                data=np.load(model/'factor.npz');w=data['loadings'];sd=np.sqrt(np.sum(w*w,axis=1)+data['noise']);w=w/sd[:,None]
                scorer=lambda a,b:w[a:b]@w.T
            elif method=='megena':
                membership=pd.read_csv(model/'memberships.tsv',sep='\t',dtype=str);position={k:i for i,k in enumerate(complete.index)}
                groups=[[position[k] for k in g.event_key if k in position] for _,g in membership.groupby('module')]
                def scorer(a,b):
                    v=np.zeros((b-a,len(complete)))
                    for ids in groups:
                        here=[i-a for i in ids if a<=i<b]
                        if here and len(ids)>1:v[np.ix_(here,ids)]=np.maximum(v[np.ix_(here,ids)],1/(len(ids)-1))
                    return v
            elif method in ['joint','independent']:
                matrices=[]
                for p in sorted(model.glob('precision_*.tsv')):
                    t=np.loadtxt(p);sd=np.sqrt(np.diag(t));matrices.append(-t/sd[:,None]/sd)
                scorer=lambda a,b:sum(np.abs(t[a:b]) for t in matrices)/len(matrices)
            else:raise ValueError('No graph adapter for engine')
            adapted=nearest_edges(scorer,complete.index)
            adapted.to_csv(model/'edges.tsv',sep='\t',index=False);graph_outputs(adapted,complete,model,mutual=True)
            result['display_partition']='Louvain partition of adapted mutual top20 graph; native method objects retained separately'
        if cfg.get('merge_modules'):
            if method!='jackknife':raise ValueError('Optional eigengene hierarchy currently requires jackknife discovery')
            from module_merge import merge
            merged=merge(model,out/'module_hierarchy')
            if merged['status']!='completed':raise RuntimeError(merged.get('error','Module merge failed'))
            shutil.copyfile(model/'assignments.tsv',model/'native_assignments.tsv')
            shutil.copyfile(out/'module_hierarchy/assignments.tsv',model/'assignments.tsv')
            result['hierarchical_merge']=merged;result['native_modules']=result.get('modules')
            result['modules']=int(pd.read_csv(model/'assignments.tsv',sep='\t').query('module != 0').module.nunique())
        assignment=pd.read_csv(model/'assignments.tsv',sep='\t',dtype={'event_key':str}).set_index('event_key')
        if (model/'edges.tsv').exists():edges=pd.read_csv(model/'edges.tsv',sep='\t',dtype={'source':str,'target':str})
        else:
            labels=assignment.reindex(complete.index).module.to_numpy();z=complete.to_numpy();z=(z-z.mean(axis=1,keepdims=True))/np.maximum(z.std(axis=1,keepdims=True),1e-8)
            def rows(a,b):return (z[a:b]@z.T/z.shape[1])*((labels[a:b,None]==labels[None,:])&(labels[a:b,None]!=0))
            edges=nearest_edges(rows,complete.index)
        edges.to_csv(out/'edges.tsv',sep='\t',index=False)
        event=audit.join(assignment,how='left');event['gene_id']=genes if genes is not None else None
        if cfg.get('annotation'):
            declared_annotation=pd.read_csv(cfg['annotation'],sep='\t',dtype=str).set_index('event_key')
            for field in ('event_type','gene_symbol'):
                if field in declared_annotation:
                    event[field]=declared_annotation[field].reindex(event.index)
        if cfg.get('revision22_stability'):
            source=Path(cfg['revision22_stability']).resolve()
            if not source.is_relative_to((ROOT/'evidence/revision22').resolve()):raise ValueError('Untrusted stability cache')
            tab=pd.read_csv(source/'stability.tsv',sep='\t').set_index('event_key')
            if not tab.index.equals(measured.index) or not np.array_equal(tab.module,assignment.reindex(measured.index).module):raise ValueError('Stability cache partition differs')
            target=out/'stability';target.mkdir()
            for file in source.iterdir():
                if file.is_file():shutil.copy2(file,target/file.name)
            write_json(target/'source.json',{'source':str(source),'partition_verified':True,'full_native_refits_retained_at_source':True})
            event=event.join(tab[['within_module_coassignment']],how='left')
        elif int(cfg.get('bootstrap',0))>0:
            status['stage']='stability';write_json(out/'status.json',status)
            if cfg.get('feature_cap') is not None:
                from feature_selection36 import refit_stability
                raw,rawmeta,support=raw_inputs(cfg)
                refit_stability(cfg,raw,rawmeta,support,assignment,measured.index,out/'stability',int(cfg['bootstrap']))
            elif cfg.get('qc_profile')=='unit_robust':
                from validation24 import product_stability
                raw,rawmeta,support=raw_inputs(cfg)
                if cfg.get('merge_modules'):
                    from product_hierarchy24 import hierarchical_stability
                    hierarchical_stability(raw,rawmeta,support,assignment,measured.index,out/'stability',int(cfg['bootstrap']))
                else:product_stability(method,raw,rawmeta,support,assignment,measured.index,out/'stability',int(cfg['bootstrap']))
            elif method=='wgcna':bootstrap(measured,meta,model,out/'stability',iterations=int(cfg['bootstrap']))
            elif method in ['consensus','leiden','wgcna_enhanced']:
                from validation22 import make_draws,outer_refit,summarize_stability
                draws=make_draws(meta,int(cfg['bootstrap']));stableout=out/'stability';stableout.mkdir();write_json(stableout/'planned_draws.json',{'draws':draws,'seed':22041})
                refits=[]
                for i,units in enumerate(draws):
                    lab,s=outer_refit(method,measured,meta,units,stableout/f'refit_{i:03}');refits.append((i,lab,s));status.update(stability_completed=i+1);write_json(out/'status.json',status)
                summarize_stability(assignment.reindex(measured.index).module.to_numpy(),measured.index,refits,stableout)
            else:raise ValueError('Outer refit stability supports WGCNA, enhanced WGCNA, consensus and Leiden')
            stable=pd.read_csv(out/'stability/stability.tsv',sep='\t').set_index('event_key');event=event.join(stable[['within_module_coassignment']],how='left')
        event.to_csv(out/'events.tsv',sep='\t',index_label='event_key')
        modules=event.dropna(subset=['module']).groupby('module').agg(events=('eligible','size'),mean_observed=('observed_fraction','mean'))
        if genes is not None:
            known=event.dropna(subset=['module','gene_id']);known=known[known.module.ne(0)]
            gene_table=known.reset_index().rename(columns={'index':'event_key'})[['event_key','module','gene_id']]
            gene_table.to_csv(out/'module_genes.tsv',sep='\t',index=False)
            modules=modules.join(known.groupby('module').gene_id.nunique().rename('unique_annotated_genes'))
        if 'within_module_coassignment' in event:modules=modules.join(event.groupby('module').within_module_coassignment.agg(['median','min']).rename(columns={'median':'median_member_support','min':'minimum_member_support'}))
        modules.to_csv(out/'modules.tsv',sep='\t')
        if 'within_module_coassignment' in event:
            event[event.within_module_coassignment.ge(.5)].to_csv(out/'supported_candidates.tsv',sep='\t',index_label='event_key')
        summary={'version':'2.9.0rc1+20261101','method':method,'merge_modules':bool(cfg.get('merge_modules')),'qc_profile':cfg.get('qc_profile','legacy'),'read_support_available':bool(cfg.get('counts') or cfg.get('read_support')),'value_kind':cfg.get('value_kind','usage'),'input_events':len(audit),'eligible_events':len(complete),'excluded_events':int((~audit.eligible).sum()),
            'assigned_events':int(assignment.module.ne(0).sum()),'samples':len(meta),'biological_units':meta.biological_sample_id.nunique(),
            'feature_cap':cfg.get('feature_cap'),'model':result,'bootstrap_requested':int(cfg.get('bootstrap',0)),
            'interpretation':'Predictive/associational network; no causal or validated-hub claim',
            'inputs':{k:{'path':str(v),'sha256':hashlib.file_digest(Path(v).open('rb'),'sha256').hexdigest()} for k,v in cfg.items() if k in ['psi','counts','totals','read_support','samples','expression','annotation'] and v}}
        from feature_selection36 import selection_summary
        summary['feature_selection']=selection_summary(audit,cfg.get('feature_cap'))
        summary['network_events']=len(complete)
        summary['eligible_events']=int(audit.eligible.sum())
        summary['analysis_unit_kind']='perturbation target response profile; shared controls can induce dependence, not independent patients' if cfg.get('value_kind')=='response' else 'declared biological block, retaining all within-block conditions'
        if (out/'stability/stability.tsv').exists():
            b=pd.read_csv(out/'stability/stability.tsv',sep='\t');v=b.within_module_coassignment.dropna()
            summary['stability']={'requested':int(b.requested_replicates.iloc[0]),'successful':int(b.successful_replicates.iloc[0]),'assigned_event_denominator':len(v),'median_support':float(v.median()) if len(v) else None,'interpretation':'Full-QC biological-unit deletion refits; failures NA' if cfg.get('qc_profile')=='unit_robust' else 'Within-original-module coassignment under biological-unit bootstrap; not a P-value or validated hub score'}
            summary['stability'].update(events_support_ge05=int(v.ge(.5).sum()),modules_median_support_ge05=int(modules.median_member_support.ge(.5).sum()),candidate_threshold=.5,
                candidate_scope='Descriptive post-fit view only; no change to full discovery network or benchmark denominator; selected candidates not independently revalidated')
            if (out/'stability/summary.json').exists():summary['stability']['partition_metrics']=json.loads((out/'stability/summary.json').read_text())
            if cfg.get('development_replicates'):
                summary['stability'].update(development_replicates=cfg['development_replicates'],confirmation_replicates=cfg['confirmation_replicates'],confirmation_summary=cfg['confirmation_stability_summary'])
        if study_mode:
            from study import initialize,analyze
            write_json(out/'report_revision.json',{'schema':'study-report-3.1' if method=='global_local' else 'study-report-2.9'})
            status.update(stage='study_interpretation');write_json(out/'status.json',status)
            summary['study']=initialize(out,cfg,measured,meta,assignment)
            default_analysis={'traits':{'condition':'categorical'} if meta.condition.nunique()>1 else {},'contrasts':cfg.get('contrasts',[])}
            summary['study_interpretation']=analyze(out,default_analysis,out/'interpretations/default')
            summary['analysis_scope']='whole_study'
            summary['input_assembly']=assembly
            summary['report_schema']='study-report-2.9'
        if cfg.get('condition_a') and cfg.get('condition_b'):
            status['stage']='differential';write_json(out/'status.json',status)
            summary['differential']=differential_edges(complete,meta,edges,out,cfg['condition_a'],cfg['condition_b'],cfg.get('permutations',199))
        if cfg.get('gene_sets') or cfg.get('gene_set_libraries'):
            if genes is None:raise ValueError('Gene-set interpretation requires event-to-gene annotation')
            from enrichment42 import enrich_libraries
            status.update(stage='functional_interpretation');write_json(out/'status.json',status)
            annotation=pd.read_csv(cfg['annotation'],sep='\t',dtype=str).set_index('event_key')
            # Actual network input, including module 0; do not use QC-only genes
            # that had no chance to enter a module under the event limit.
            libraries=cfg.get('gene_set_libraries') or [{'path':cfg['gene_sets']}]
            summary['pathway_enrichment']=enrich_libraries(assignment.reindex(complete.index),annotation,libraries,out)
            sources=[]
            for item in libraries:
                path=Path(item['path'])
                with path.open('rb') as handle:digest=hashlib.file_digest(handle,'sha256').hexdigest()
                sources.append({**item,'sha256':digest})
            summary['inputs']['gene_set_libraries']=sources
            if cfg.get('gene_sets'):summary['inputs']['gene_sets']=sources[0]
        write_json(out/'summary.json',summary);status.update(stage='report');write_json(out/'status.json',status)
        write_report(out,summary)
        status.update(status='completed',stage='complete',seconds=time.time()-status['started']);write_json(out/'status.json',status)
    except Exception as exc:
        traceback.print_exc();status.update(status='failed',error=f'{type(exc).__name__}: {exc}',seconds=time.time()-status['started']);write_json(out/'status.json',status)
        write_report(out,{'failure':status,'model_outputs_may_be_partial':True})
        raise

if __name__=='__main__':
    folder=Path(sys.argv[1]);run(json.loads((folder/'config.json').read_text(encoding='utf-8')),folder)
