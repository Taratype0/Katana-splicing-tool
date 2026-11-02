"""Visible, source-traceable RBP gene annotations; no scans or statistical tests."""
from collections import defaultdict
from pathlib import Path
import csv
import gzip
import html
import json
import re
import io
import os
import hashlib


def species(value):
    value=str(value or '').strip().lower().replace(' ', '_')
    return {'human':'homo_sapiens','mouse':'mus_musculus','9606':'homo_sapiens','10090':'mus_musculus'}.get(value,value)


def project_species(project):
    annotation=project.get('annotation',{})
    value=annotation.get('species') or annotation.get('organism') or project.get('metadata',{}).get('species')
    if value:return species(value), 'declared project species'
    build=str(annotation.get('genome_build','')).lower()
    inferred={'hg38':'homo_sapiens','grch38':'homo_sapiens','hg19':'homo_sapiens','grch37':'homo_sapiens',
              'mm10':'mus_musculus','grcm38':'mus_musculus','mm39':'mus_musculus','grcm39':'mus_musculus'}.get(build,'')
    return inferred, ('project reference '+build if inferred else 'project species not declared')


def references():
    path=Path(__file__).with_name('motif_gene_reference.json')
    return json.loads(path.read_text('utf-8')) if path.is_file() else {'records':{}}


def candidate_genes(model, target_species='', reference=None):
    """Retain every declared RBP; never treat a motif/native ID as a gene ID."""
    reference=reference if reference is not None else references()
    records=reference.get('records',{})
    assignments=list(model.get('assignments') or [])
    if not assignments:
        assignments=[{'RBP':name.strip()} for name in str(model.get('RBPs','')).split(';') if name.strip()]
    collected={}
    for assignment in assignments:
        names=[s.strip() for s in str(assignment.get('RBP') or assignment.get('RBPs') or '').split(';') if s.strip()]
        for name in names:
            origin=species(assignment.get('species'))
            known=records.get(origin+'|'+name.upper(),{})
            source_id=assignment.get('gene_id') or assignment.get('ensembl_gene_id') or known.get('gene_id','')
            source_symbol=known.get('symbol',name)
            candidates=[]
            if target_species and origin and target_species != origin:
                for ortholog in known.get('orthologs',[]):
                    if species(ortholog.get('species'))==target_species:
                        candidates.append(dict(symbol=ortholog['symbol'],gene_id=ortholog['gene_id'],species=target_species,
                            mapping='Cross-species candidate: '+ortholog['relationship'].replace('_',' '),
                            identity_evidence=ortholog.get('evidence',''),identity_url=ortholog.get('url',''),hgnc_id=ortholog.get('hgnc_id','')))
                if not candidates:
                    candidates=[dict(symbol=source_symbol,gene_id=source_id,species=origin,
                        mapping='Source-species annotation; project-species mapping unavailable',identity_evidence=known.get('identity_evidence',''),identity_url=known.get('identity_url',''))]
            else:
                candidates=[dict(symbol=source_symbol,gene_id=source_id,species=origin,
                    mapping='Same-species library assignment' if origin and origin==target_species else 'Source annotation; project or source species unavailable',
                    identity_evidence=known.get('identity_evidence',''),identity_url=known.get('identity_url',''))]
            for candidate in candidates:
                key=(candidate['symbol'],candidate['gene_id'],candidate['species'],origin,source_id,candidate['mapping'])
                row=collected.setdefault(key,{**candidate,'source_symbol':source_symbol,'source_gene_id':source_id,'source_species':origin,
                    'library_sources':set(),'native_model_ids':set(),'library_evidence':set()})
                for field,source in [('library_sources','source'),('native_model_ids','native_id'),('library_evidence','evidence')]:
                    if assignment.get(source):row[field].add(str(assignment[source]))
    rows=[]
    for row in collected.values():
        rows.append({k:'; '.join(sorted(v)) if isinstance(v,set) else v for k,v in row.items()})
    return sorted(rows,key=lambda r:(r['symbol'],r['gene_id'],r['source_species']))


def _read_rows(path):
    with (gzip.open(path,'rt',encoding='utf-8-sig',newline='') if str(path).endswith('.gz') else Path(path).open(encoding='utf-8-sig',newline='')) as stream:
        yield from csv.DictReader(stream,delimiter='\t')


def enrich(summary, models, directory, comparison_id, *, target_species='', species_basis='', events_path=None, matching_path=None, changes_path=None):
    directory=Path(directory); directory.mkdir(parents=True,exist_ok=True)
    reference=references(); all_candidates=[]
    for token,model in models.items():
        rows=candidate_genes(model,target_species,reference)
        all_candidates.extend({'comparison_id':comparison_id,'model_id':token,**r} for r in rows)
    columns=['comparison_id','model_id','symbol','gene_id','species','mapping','source_symbol','source_gene_id','source_species',
             'library_sources','native_model_ids','library_evidence','identity_evidence','identity_url','hgnc_id']
    candidates_path=directory/'MOTIF_REPORT_CANDIDATE_RBP_GENES.tsv'
    _tsv(candidates_path,all_candidates,columns)
    reference_path=directory/'MOTIF_REPORT_GENE_MAPPING_REFERENCE.json'
    atomic(reference_path,json.dumps(reference,indent=2,ensure_ascii=False))
    files=[candidates_path,reference_path]
    cards={c['model_id']:c for c in summary['cards']}
    events={}
    if events_path:
        for row in _read_rows(events_path):
            if row.get('comparison_id',comparison_id)!=comparison_id:raise ValueError('Target events belong to a different comparison')
            if row.get('selected','true').lower() not in ('true','1','1.0'):continue
            events[row['event_id']]=row
    matching={r['foreground_event']:r.get('status','not recorded') for r in _read_rows(matching_path)} if matching_path else {}
    target_rows={token:{} for token in models}
    if changes_path and events_path:
        for site in _read_rows(changes_path):
            token=site.get('model_id'); event=site.get('event_id')
            if token not in models or event not in events:continue
            if site.get('comparison_id',comparison_id)!=comparison_id:raise ValueError('Target-site data belongs to a different comparison')
            direction=site.get('consequence_direction')
            if direction not in ('increased','decreased'):continue
            source=events[event]; gene_id=site.get('gene_id') or source.get('gene_id','')
            key=(event,gene_id)
            row=target_rows[token].setdefault(key,dict(comparison_id=comparison_id,model_id=token,event_id=event,
                target_gene_id=gene_id,target_gene_symbol=source.get('gene_symbol') or source.get('geneSymbol',''),
                event_type=site.get('event_type',''),dPSI=source.get('dpsi',source.get('IncLevelDifference','')),
                matching_status=matching.get(event,'not recorded'),increased_candidate_sites=0,decreased_candidate_sites=0,
                interpretation='Descriptive sequence-opportunity change; no individual target significance inferred'))
            row[direction+'_candidate_sites']+=1
    full_targets=directory/'MOTIF_REPORT_ALL_TARGET_EVENTS.tsv'
    target_columns=['comparison_id','model_id','target_gene_id','target_gene_symbol','event_id','event_type','dPSI',
        'matching_status','increased_candidate_sites','decreased_candidate_sites','interpretation']
    _tsv(full_targets,[r for token in sorted(target_rows) for r in target_rows[token].values()],target_columns)
    files.append(full_targets)
    for token,card in cards.items():
        card['candidate_genes']=[r for r in all_candidates if r['model_id']==token]
        card['target_species']=target_species;card['target_species_basis']=species_basis
        rows=sorted(target_rows[token].values(),key=lambda r:(r['target_gene_symbol'],r['target_gene_id'],r['event_id']))
        target_path=directory/('MOTIF_REPORT_'+hashlib.sha256(token.encode()).hexdigest()[:16]+'-target-events.tsv')
        _tsv(target_path,rows,target_columns)
        files.append(target_path)
        genes=sorted({(r['target_gene_symbol'],r['target_gene_id']) for r in rows})
        card['target_events_path']=str(target_path);card['candidate_genes_path']=str(candidates_path)
        card['gene_mapping_reference_path']=str(reference_path)
        card['target_summary']=dict(available=bool(changes_path and events_path),genes=len({r['target_gene_id'] for r in rows}),
            events=len({r['event_id'] for r in rows}),matched_events=len({r['event_id'] for r in rows if r['matching_status']=='matched'}),
            increased_sites=sum(r['increased_candidate_sites'] for r in rows),decreased_sites=sum(r['decreased_candidate_sites'] for r in rows),
            preview=[dict(symbol=s or g,gene_id=g) for s,g in genes[:8]],ordering='Alphabetical display only; full event list retained')
    summary['candidate_gene_annotation']=dict(version='motif-rbp-genes-v1',target_species=target_species,species_basis=species_basis,
        qualifying_candidate_genes=len({(r['species'],r['gene_id'] or r['symbol']) for c in cards.values() for r in c['candidate_genes']}),
        displayed_motifs_with_gene_annotations=sum(bool(c['candidate_genes']) for c in cards.values()),
        inference='Candidate binding partners; not RBP differential expression or gene-level significance')
    return files




def atomic(path,content):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_name(path.name+'.annotation.tmp')
    temporary.write_bytes(content.encode('utf-8') if isinstance(content,str) else content)
    os.replace(temporary,path)


def _tsv(path,rows,columns):
    stream=io.StringIO(newline='');writer=csv.DictWriter(stream,fieldnames=columns,delimiter='\t',extrasaction='ignore')
    writer.writeheader();writer.writerows(rows);atomic(path,stream.getvalue())


def annotate_completed_results(source, destination, *, project=None):
    """Native annotation-only reuse API: validate a completed run, never refit it.

    Writes complete motif-to-gene and descriptive target-event tables in a new
    directory. Statistics and original results remain byte-identical.
    """
    from katana_next.reports.motif_sequences import evidence, _best_windows
    source=Path(source).resolve();destination=Path(destination).resolve()
    if destination.is_relative_to(source):raise ValueError('Annotations must use a new directory outside the completed source')
    complete=json.loads((source/'COMPLETE.json').read_text('utf-8-sig'))
    if not complete.get('complete') or not complete.get('members'):raise ValueError('A verified complete native run is required')
    def digest(path):
        with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
    originals={}
    for name,record in complete['members'].items():
        path=(source/name).resolve()
        if not path.is_relative_to(source) or not path.is_file():raise ValueError('Invalid source member: '+name)
        expected=record['sha256'] if isinstance(record,dict) else record
        if digest(path)!=expected:raise ValueError('Completed source hash mismatch: '+name)
        originals[str(path)]=expected
    plan=json.loads((source/'ANALYSIS_PLAN_BEFORE_SCORING.json').read_text('utf-8-sig'))
    result=json.loads((source/'RESULTS.json').read_text('utf-8-sig'))
    comparison_id=str(plan['comparison_id'])
    if result['comparison_id']!=comparison_id:raise ValueError('Mismatched completed comparisons')
    models=plan['models'];regions=list(_read_rows(source/'NEW_CLUSTER_FWER.tsv'))
    windows=_best_windows(source/'NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz') if result.get('new_window_q_lt_0_05',0) else []
    summary=evidence(models,regions,windows)
    provenance=plan.get('provenance',{})
    project=project or {'annotation':{'genome_build':provenance.get('genome_build')},'metadata':{'species':provenance.get('species')}}
    target,basis=project_species(project)
    destination.mkdir(parents=True,exist_ok=False)
    def existing(name):
        path=source/name
        return path if path.is_file() else None
    files=enrich(summary,models,destination,comparison_id,target_species=target,species_basis=basis,
        events_path=existing('FROZEN_SELECTED_EVENTS.tsv'),matching_path=existing('MATCHING_ALL_SELECTED_WITH_LOSSES.tsv'),changes_path=existing('PAIRED_FORM_INSTANCE_CHANGES.tsv.gz'))
    if not all(digest(Path(p))==h for p,h in originals.items()):raise RuntimeError('Completed source changed during annotation')
    receipt=dict(schema='katana.motiftrace_rbp_gene_annotation/v1',complete=True,comparison_id=comparison_id,
        source_complete_sha256=digest(source/'COMPLETE.json'),source=str(source),originals_unchanged=True,
        new_scans=0,new_statistical_tests=0,annotation=summary['candidate_gene_annotation'],
        members={p.name:dict(sha256=digest(p),bytes=p.stat().st_size) for p in files})
    atomic(destination/'MOTIF_RBP_ANNOTATION.json',json.dumps(receipt,indent=2))
    return receipt


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description='Annotate completed MotifTrace motifs with candidate RBP genes without refitting.')
    parser.add_argument('--source',required=True,type=Path);parser.add_argument('--outdir',required=True,type=Path)
    parser.add_argument('--project',type=Path)
    args=parser.parse_args()
    project=json.loads(args.project.read_text('utf-8-sig')) if args.project else None
    print(json.dumps(annotate_completed_results(args.source,args.outdir,project=project),indent=2))
