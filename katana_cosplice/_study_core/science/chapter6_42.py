"""Export a reviewable Katana chapter from a completed, frozen CoSplice report."""
import json,re,html,hashlib
from pathlib import Path
import pandas as pd

def build(run,destination):
    run=Path(run);destination=Path(destination)
    report=(destination/'report.html').read_text(encoding='utf8')
    data=json.loads(re.search(r'<script type="application/json" id="reportData">(.*?)</script>',report,re.S)[1])
    report=re.sub(r'<section\b[^>]*id="(?:donor-inference-calibration|donor-validation-5000)".*?</section>','',report,flags=re.S)
    report=re.sub(r'<aside\b[^>]*id="event-selection-summary".*?</aside>','',report,flags=re.S)
    report=report.replace('CO-SPLICING ANALYSIS · STUDY-WIDE NETWORK','CHAPTER 06 · KATANA-SPLICE-NETWORK')
    report=re.sub(r'<h1>.*?</h1>','<h1>CoSplice<br>Shared splicing programmes</h1>',report,count=1,flags=re.S)
    report=re.sub(r'<p class="eyebrow">.*?</p>','<p class="eyebrow">Chapter 06 · '+html.escape(data['manifest']['study_id'])+'</p>',report,count=1,flags=re.S)
    for n in range(1,5):
        report=report.replace(f'<b>0{n}</b>',f'<b>6.{n}</b>').replace(f'<span class="chapter-number">0{n}</span>',f'<span class="chapter-number">6.{n}</span>')
    intro='''<div class="signal chapter6-summary"><strong>One study-wide network; separate module gene lists and downstream interpretations.</strong><p>Each sample remains a separate column. Quality and variance filtering precede WGCNA/TOM module discovery and kME membership filtering; events are not preselected by differential-splicing significance. Frozen module scores are then related to the study design.</p><div class="flow42"><span>Sample-level splicing</span><span>QC + event limit</span><span>WGCNA + kME</span><span>Donor-aware evidence</span><span>Module genes + functions</span></div></div>'''
    report=report.replace('<div class="scopebar">',intro+'<div class="scopebar">',1)
    if (run/'validation37/module_validation.tsv').exists():
        report=report.replace('<div class="facts" id="facts">','<p class="lead">M1 and M2 are two directional modules of approximately one shared differentiation programme. Their memberships remain separate; neither is an independently confirmed mechanism.</p><div class="facts" id="facts">',1)
        tab=pd.read_csv(run/'validation37/module_validation.tsv',sep='\t')
        assoc={r['module']:r for r in data['associations']}
        rows=''
        for r in tab.to_dict('records'):
            a=assoc.get(r['module'],{});q=a.get('q_value');p=a.get('p_value')
            rows+='<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in [r['module'],r['stability_events'],r['genes'],f'{p:.3g}' if p is not None else 'NA',f'{q:.3g}' if q is not None else 'NA',f"{r['median_match_jaccard']:.3f}",f"{r['median_trajectory_r']:.3f}",'Yes' if a.get('qualified_candidate') else 'No'])+'</tr>'
        audit='<details><summary>All four modules: association, stability and held-out response</summary><div class="table-scroll"><table class="metadata"><thead><tr>'+''.join('<th>'+s+'</th>' for s in ['Module','Events','Genes','P','BH FDR','Refit Jaccard','Held-out trajectory r','Qualified'])+'</tr></thead><tbody>'+rows+'</tbody></table></div><p>50 distinct delete-four-donor refits, not a with-replacement bootstrap. Five held-out folds cover 19 donors from this same study. Candidate qualification additionally requires Holm FWER ≤ 0.05 and the unchanged stability gates.</p><p>M1 and M2 have opposite scores (r = −0.995), consistent with approximately one shared time-response programme. M3 and M4 pass structural stability but fail the donor-trajectory test; 90–92% of their score variance lies between donors. These are not four independent biological programmes.</p></details>'
        report=report.replace('<div id="topModules" class="table-scroll"></div>','<div id="topModules" class="table-scroll"></div>'+audit,1)
    methods='''<details><summary>What this chapter can establish</summary><p>Module genes share splicing-usage patterns across this study. kME measures event membership, not regulatory causality. Full donor refits test recovery; held-out donors test transfer within the study. GO/KEGG/Reactome enrichment describes annotation over-representation conditional on the selected network genes. It does not show that a pathway is activated or independently validate a splicing mechanism.</p><p>Repeated visits are not independent donors. Module associations use the comparisons and pairing declared for this project. Optional refits or held-out validation are reported only when actually computed. General superiority over WGCNA or other co-splicing tools has not been established.</p></details>'''
    report=report.replace('<div id="downloads" class="file-list"></div>',methods+'<div id="downloads" class="file-list"></div>',1)
    report=report.replace('</style>',' .flow42{display:flex;gap:10px;flex-wrap:wrap;margin-top:16px}.flow42 span{padding:8px 12px;border:1px solid #c9d5e2;background:white;color:#2863a5}.chapter6-summary{margin-top:24px}.chapter6-summary p{max-width:980px}.preview42{padding:8px 7%;background:#f4f7fa;color:#526171;font-size:12px}</style>',1)
    report=report.replace('</header>','</header><div class="preview42">Whole-study CoSplice · All declared comparable samples · <a href="report.html">Open the complete study report</a></div>',1)
    report=report.replace('<a href="chapter6">Chapter 06</a>','').replace('<a href="cosplice_chapter6.html">Chapter 06</a>','')
    report=report.replace('KATANA · CoSplice Next · report 3.9','KATANA · CoSplice · Chapter 06')
    (destination/'cosplice_chapter6.html').write_text(report,encoding='utf8')
    manifest={'schema':'katana-cosplice-chapter-preview-1','chapter':6,'analysis_id':'programs.cosplice','scope':'whole_study','contrast_id':None,'status':'production','network_id':data['manifest']['network_id'],'assignment_sha256':data['manifest']['assignment_sha256'],'source_run':str(run),'source_report_sha256':hashlib.sha256((destination/'report.html').read_bytes()).hexdigest(),'entrypoint':'cosplice_chapter6.html','report_tables':'report_tables','network_refitted':False}
    (destination/'chapter6_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
    return manifest
