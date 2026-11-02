"""Explain an unmeasured stability endpoint without treating it as failure."""
def adapt(rendered, result):
    if int(result.get('bootstrap_requested',0)) != 0:
        return rendered
    rendered=rendered.replace('Stability-qualified module associations','Module associations')
    start='<p class="chapter-note">Default candidates must pass '
    if start in rendered:
        a=rendered.index(start);b=rendered.index('</p>',a)+4
        rendered=rendered[:a]+('<p class="chapter-note">One whole-study WGCNA/TOM network with kME membership filtering. '
            'Modules are shown by association FDR, with separate gene and event lists. Stability resampling is off; '
            'these are exploratory associations, not stability-validated candidates. kME is a membership measure, not a stability test.</p>')+rendered[b:]
    rendered=rendered.replace('<option value="qualified">Stability-qualified candidates</option>',
        '<option value="qualified" disabled>Stability not assessed</option>')
    rendered=rendered.replace('<option value="exploratory">','<option value="exploratory" selected>')
    rendered=rendered.replace('href="report_tables/qualified_module_lists.zip">Qualified module lists for this interpretation (separate files in ZIP)',
        'href="report_tables/all_module_lists.zip">Complete gene/event lists for all modules (separate files in ZIP)')
    hook="""
const originalScopeWithoutRefits=updateScope;
updateScope=function(){
 originalScopeWithoutRefits();
 const sig=scopeRows().filter(r=>r.status==='estimated'&&num(r.q_value)&&r.q_value<=.05);
 $('headline').textContent=sig.length+' FDR-significant associations from '+D.modules.length+' shared modules. One whole-study WGCNA + kME fit; stability not assessed.';
 $('facts').innerHTML=[['Stability','Not assessed'],['FDR-significant associations',sig.length],['Discovered modules',D.modules.length],['Individual samples / conditions',D.manifest.samples+' / '+conditions.length]].map(([l,v])=>'<div class="fact"><strong>'+E(v)+'</strong><span>'+E(l)+'</span></div>').join('');
 $('scopeSummary').textContent=sig.length+' FDR-significant · stability not assessed';
 $('gateSummary').textContent='Stability refits are off (0 requested). Modules are shown using the unchanged association FDR <= 0.05 criterion. Missing validation is not a failed stability test. Neither kME nor association significance establishes reproducibility; repeated-donor inference retains its recorded calibration limits.';
};
$('scopeSelect').onchange=updateScope;$('evidenceView').onchange=updateScope;
"""
    needle="$('downloads').innerHTML=D.downloads.map"
    assert needle in rendered
    return rendered.replace(needle,hook+'\n'+needle,1)
