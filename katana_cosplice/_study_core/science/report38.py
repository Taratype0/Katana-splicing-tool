"""Explicit estimand, calibration and results for repeated-donor trajectories."""
from pathlib import Path
import base64,html,json
import pandas as pd
from inference38 import validate_certificate,verify_associations

def panel(run,interpretation,associations):
    root=Path(run)/'inference38';binding=validate_certificate(run)
    valid,reason=verify_associations(run,associations)
    if not valid:raise ValueError(reason)
    cal=json.loads((root/'calibration.json').read_text());scores=json.loads((root/'score_calibration.json').read_text())
    tab=associations.sort_values(['q_value','p_value','module'],kind='stable')
    display=tab[['module','biological_units','p_value','q_value','p_fwer','qualified_candidate']].rename(columns={
        'module':'Module','biological_units':'Independent donors','p_value':'Module P (3 components)',
        'q_value':'BH FDR (all 4 modules)','p_fwer':'Holm FWER','qualified_candidate':'Qualified trajectory candidate'})
    display['Qualified trajectory candidate']=display['Qualified trajectory candidate'].map({True:'Yes',False:'No'})
    parts=pd.read_csv(Path(interpretation)/'trajectory_components.tsv',sep='\t')
    parts=parts[['module','degree','nonzero_donors','positive_donors','negative_donors','median_coefficient','median_ci_low','median_ci_high','p_value']]
    links=[]
    for name in ['analysis_plan.md','analysis_summary.md','calibration.json','score_calibration.tsv','native_scenarios.tsv','native_calibration.tsv']:
        p=root/name
        if p.exists():links.append(f'<a download="{html.escape(name)}" href="data:application/octet-stream;base64,{base64.b64encode(p.read_bytes()).decode()}">{html.escape(name)}</a>')
    sc=pd.DataFrame(scores['scenarios'])[['scenario','experiments','false_positive_families','fwer','upper95','legacy_fwer']]
    native=pd.read_csv(root/'native_scenarios.tsv',sep='\t')
    ledger=pd.read_csv(root/'native_calibration.tsv',sep='\t',keep_default_na=False);null=ledger[ledger.kind.eq('null')]
    fmt=lambda x:f'{x:.4g}'
    return ('<section class="chapter" id="donor-inference-calibration" style="margin:28px clamp(16px,3.2vw,40px)"><div class="chapter-head"><span class="chapter-number">✓</span>'
        '<h2>Repeated-donor trajectory inference / 重复供体推断校准</h2></div>'
        f'<p class="chapter-note"><strong>{int(tab.qualified_candidate.sum())} / 4 modules qualify as within-study trajectory candidates.</strong> '
        'The network, 5,000 selected events, module memberships and stability cutoffs are unchanged. '
        '297 visits are analyzed as 19 independent donor trajectories. '
        'The primary endpoint is a typical-donor cubic time response, replacing the legacy equality-of-all-condition-means endpoint.</p>'
        '<div class="table-scroll">'+display.to_html(index=False,border=0,classes='metadata',float_format=fmt)+'</div>'
        '<p>M1 and M2 have opposite scores (r = −0.995): two directional gene/event lists, approximately one shared time-response program. '
        'M3/M4 remain structurally stable but fail the trajectory test; 90–92% of their variance lies between donors. '
        'Their exclusion is not a stability failure.</p>'
        '<details><summary>Statistical design, completed calibration, costs and limits</summary>'
        '<p><strong>What changed:</strong> the primary endpoint is a typical-donor time response, using linear, quadratic and cubic Legendre coefficients fitted separately in each donor. '
        'Each coefficient uses an exact two-sided sign test; three tests are Bonferroni-corrected within each module, then all four modules receive BH and Holm correction. '
        'Qualification requires both the existing BH FDR ≤ 0.05 and the additional Holm FWER ≤ 0.05 gate. No threshold was relaxed.</p>'
        '<p><strong>This is a different, explicitly defined hypothesis:</strong> balanced donor-coefficient signs about zero, not equality of all 16 condition means. '
        'Legacy categorical tests remain downloadable as diagnostics and cannot qualify repeated-donor candidates. '
        'The reported effect is a median linear coefficient on time scaled to [−1,1], not a per-day slope or total start-to-end change. '
        'A module may also be significant through its planned quadratic or cubic component.</p>'
        '<details><summary>All planned components, donor signs and unadjusted median intervals</summary><div class="table-scroll">'+parts.to_html(index=False,border=0,classes='metadata',float_format=fmt)+
        '</div><p>Degree 1/2/3 = linear/quadratic/cubic. Intervals are distribution-free, at least 95% marginal median intervals; they are not simultaneous confidence intervals.</p>'
        '<p><a download href="report_tables/trajectory_components.tsv">Component tests</a> · <a download href="report_tables/donor_coefficients.tsv">Every donor coefficient</a> · '
        '<a download href="report_tables/legacy_associations.tsv">Legacy categorical tests (diagnostic)</a></p></details>'
        f'<p><strong>Calibration:</strong> {scores["experiments"]:,} fixed-score null experiments across eight scenarios; family-wise false-positive rates '
        f'{min(s["fwer"] for s in scores["scenarios"]):.2%}–{max(s["fwer"] for s in scores["scenarios"]):.2%}. '
        f'Full QC → variance selection → native WGCNA → inference: {cal["native_false_positive_families"]}/{cal["native_null_experiments"]} null datasets had any false call '
        f'(pooled one-sided 95% upper bound {cal["native_upper95"]:.2%}). Every native dataset starts from 97,385 synthetic events and selects 5,000. '
        f'{cal["positive_controls"]}/{cal["positive_control_total"]} positive controls detected a signal; '
        f'{cal["mixed_control_any_null"]}/{cal["mixed_control_total"]} mixed controls called a module containing no signal events.</p>'
        '<details><summary>Every calibration scenario, uncertainty and legacy comparison</summary>'
        '<p>Fixed-score calibration. Legacy family-wise error is reported only where both mean-zero and median-zero nulls hold; the asymmetric median-null case is not comparable. '
        'The legacy 15-restriction cluster-Wald test is poorly calibrated with 19 donors in these simulations; these rates do not describe all uses of cluster-robust inference.</p>'
        '<div class="table-scroll">'+sc.to_html(index=False,border=0,classes='metadata',float_format=fmt,na_rep='Not comparable')+'</div>'
        '<p>Full native calibration scenarios: 0 = independent noise; 1 = autocorrelation + donor offsets; 2 = heterogeneous heavy-tailed noise; '
        '3 = strong autocorrelation + offsets + missing values. There are 32 datasets per scenario; their individual confidence bounds are wider than the pooled bound. '
        f'Synthetic null networks contained {int(null.tested_modules.min())}–{int(null.tested_modules.max())} tested modules '
        f'and {int(null.assigned_events.min())}–{int(null.assigned_events.max())} assigned events; they are stress-test generators, not exact replicas of the observed four-module partition.</p>'
        '<div class="table-scroll">'+native.to_html(index=False,border=0,classes='metadata',float_format=fmt)+'</div></details>'
        '<p><strong>Scope:</strong> exact sign-test calculations apply to fixed scores. Network discovery on the same data is assessed empirically under the documented generators; '
        'this is not a universal selective-inference guarantee or external replication. Native simulations start from usage values with synthetic missingness; read-count sampling is not simulated. Independent donors, comparable time measurements and noninformative missingness are required. '
        'Temporal association does not establish causation. Other traits, covariates, contrasts or networks do not inherit this calibration.</p>'
        f'<p>Calibration cost: {cal["elapsed_seconds"]/60:.1f} minutes for 144 native fits with four workers; summed native fit time '
        f'{cal["native_seconds_sum"]/60:.1f} minutes; largest native process tree {cal["native_peak_rss_bytes"]/2**30:.2f} GiB. '
        'This is a one-time validation cost; ordinary reinterpretation reuses the frozen network.</p></details>'
        '<p>'+ ' · '.join(links)+'</p></section>')

def adapt_template(rendered):
    # Only this explicitly calibrated interpretation receives these labels.
    rendered=rendered.replace('Permutation P values have resolution 0.0001; ties are retained and ordered by module ID.',
        'P values are exact binomial tails with planned multiplicity corrections; no Monte Carlo resolution floor applies. Ties are retained and ordered by module ID.')
    rendered=rendered.replace("{key:'absolute_loading',label:'|Factor loading|'},{key:'loading_margin',label:'Loading margin'},",'')
    rendered=rendered.replace('they did not define the factor cores','they are a display summary, not the native WGCNA adjacency/TOM used for clustering')
    rendered=rendered.replace('No attached negative control for this run.','Donor-inference null, positive and mixed controls are documented in the calibration panel above. No separate external causal-edge validation is attached.')
    rendered=rendered.replace("'Slope: '+fmt(r.effect)+' SD per trait unit. '","'Median donor linear Legendre coefficient: '+fmt(r.effect)+' (time scaled to [−1,1]); marginal 95% median interval ['+fmt(r.ci_low)+', '+fmt(r.ci_high)+']. '")
    rendered=rendered.replace("'All conditions / '+r.factor+' (overall)'","'Typical-donor cubic time trajectory / '+r.factor")
    rendered=rendered.replace("+' · '+E(m.latent_factor)+' ('+E(m.loading_direction)+')</strong><span>'", "+' · WGCNA module</strong><span>'")
    rendered=rendered.replace("'Member support'},{key:'median_match_jaccard'","'Member support'},{key:'p_fwer',label:'Holm FWER'},{key:'median_match_jaccard'")
    rendered=rendered.replace("let audit=D.statistical_audit;","let audit=null;")
    rendered=rendered.replace("conditions=[...new Set(D.samples.map(s=>s.condition))]","conditions=[...new Set(D.samples.map(s=>s.condition))].sort((a,b)=>Number(D.samples.find(s=>s.condition===a).time_day)-Number(D.samples.find(s=>s.condition===b).time_day))")
    rendered=rendered.replace("E(c)+'</text><text x=\"'", "E('d'+D.samples.find(s=>s.condition===c).time_day)+'</text><text x=\"'")
    rendered=rendered.replace('Dots: samples. Lines: condition means.', 'Days are chronological (d0–d15). Dots: visits, grouped by day; the same donor contributes repeated visits. Lines: descriptive condition means. Inference uses donor coefficients.')
    rendered=rendered.replace("inference_not_calibrated_for_candidate_gate:'Inference not qualified for this design',","inference_not_calibrated_for_candidate_gate:'Inference not qualified for this design',donor_calibration_invalid_or_mismatched:'Donor calibration does not match this evidence',donor_Holm_FWER_above_0_05:'Holm FWER > 0.05',")
    rendered=rendered.replace("'Stability evidence: '+D.qualification.status+'. FDR is computed before the stability gate; no selected-only FDR recalculation.'","'Stability evidence: '+D.qualification.status+'. Donor trajectory candidates additionally require matching empirical calibration and Holm FWER ≤ 0.05. BH and Holm use all four modules before filtering.'")
    return rendered
