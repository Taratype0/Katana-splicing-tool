"""English presentation of native MotifTrace candidate-gene annotations."""
import html
import re

def candidate_html(card):
    esc=lambda value:html.escape(str(value or 'Not provided'),quote=True)
    rows=[]
    for gene in card.get('candidate_genes',[]):
        symbol=esc(gene['symbol']);stable_id=gene.get('gene_id','')
        if stable_id.startswith(('ENSG','ENSMUSG')) and re.fullmatch('[A-Za-z0-9.]+',stable_id):
            stable='<a href="https://www.ensembl.org/id/'+stable_id+'">'+esc(stable_id)+'</a>'
        else:stable=esc(stable_id)
        evidence=esc(gene['library_sources'])+'; model '+esc(gene['native_model_ids'])
        if gene.get('identity_url'):evidence+='; <a href="'+esc(gene['identity_url'])+'">gene mapping evidence</a>'
        rows.append('<tr><td><strong>'+symbol+'</strong><br>'+stable+'</td><td>'+esc(gene['species'].replace('_',' '))+'</td><td>'+esc(gene['mapping'])+
            '<br>Source: '+esc(gene['source_symbol'])+' / '+esc(gene['source_gene_id'])+' / '+esc(gene['source_species'].replace('_',' '))+'</td><td>'+evidence+'</td></tr>')
    table=('<div class="table-scroll"><table class="metadata"><tr><th>Candidate RBP gene / stable ID</th><th>Candidate species</th><th>Mapping to the original motif</th><th>Evidence source</th></tr>'+''.join(rows)+'</table></div>') if rows else '<p>No RBP gene assignment was provided for this motif model.</p>'
    out='<section class="candidate-rbp-genes"><h4>Candidate RBP genes</h4>'+table
    if len(rows)>1:out+='<p class="fig-notes">Shared motif: retain all listed candidates; this sequence does not uniquely identify one RBP.</p>'
    out+='<p class="fig-notes">Gene mapping annotates the original motif model. Cross-species candidates require orthology and are not direct binding evidence in the study species. Motif FWER/FDR values above are not RBP gene P/q values.</p>'
    target=card.get('target_summary',{})
    out+='<h4>Target RNA genes and events with changed motif opportunities</h4>'
    if target.get('available'):
        out+='<p>'+str(target['genes'])+' target genes; '+str(target['events'])+' selected events ('+str(target['matched_events'])+' matched events). '+str(target['increased_sites'])+' increased and '+str(target['decreased_sites'])+' decreased candidate-site opportunities.</p>'
        if target['preview']:out+='<p>Examples (alphabetical): '+', '.join(esc(r['symbol'])+' ('+esc(r['gene_id'])+')' for r in target['preview'])+'.</p>'
        out+='<p class="fig-notes">Descriptive counts across all selected events for this motif, including explicitly labelled unmatched events. These are not individually significant targets or proof that a specific target drives a significant region. The complete event table is downloadable below.</p>'
    else:out+='<p>Target-level sequence-change data are unavailable for this report; no targets are inferred.</p>'
    return out+'</section>'
