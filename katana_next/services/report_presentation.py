"""Read delivered figures and network data without executing bundled HTML."""
from functools import lru_cache
from pathlib import Path
import json
import re
from .report_import import member


def report_root(context):
    return Path(context.manifest.metadata['source_report_directory']).resolve()


@lru_cache(maxsize=4)
def _manifest(root, stamp):
    return json.loads(member(Path(root), 'report-manifest.json').read_text('utf-8-sig'))


def manifest(context):
    root=report_root(context)
    return _manifest(str(root), (root/'report-manifest.json').stat().st_mtime_ns)


def presentation(context):
    data=manifest(context)
    has_splicenet = any(r.get('analysis_id') == 'programs.splicenet' for r in data.get('complete_results', []))
    has_splicetrace = any(str(r.get('analysis_id', '')).startswith('prioritize.splicetrace') for r in data.get('complete_results', []))
    return {'panels': [{k:p.get(k) for k in ('id','chapter','title','image','svg','contrast_id','caption','metadata','data_links')}
                       for p in data.get('panels', []) if not p.get('output_only')],
            'network_available': bool(data.get('network_browsers')) or has_splicenet,
            'splicetrace_available': has_splicetrace,
            'motif_summary': data.get('motif_summaries', []),
            'complete_results': len(data.get('complete_results', []))}


def _comparisons(value):
    found=set()
    if isinstance(value, dict):
        for key,item in value.items():
            if key in ('comparison_id','contrast_id') and isinstance(item,str): found.add(item)
            elif isinstance(item,(dict,list)): found |= _comparisons(item)
    elif isinstance(value,list):
        for item in value: found |= _comparisons(item)
    return found


@lru_cache(maxsize=8)
def _figure(root, manifest_stamp, panel_id):
    data=_manifest(root, manifest_stamp); root=Path(root)
    panel=next((p for p in data.get('panels',[]) if p.get('id')==panel_id),None)
    if panel is None: raise ValueError('Unknown delivered figure')
    filename=re.sub(r'^[a-f0-9]{12}-','',Path(panel.get('image','')).name)
    filename=re.sub(r'\.(png|svg)$','.json',filename)
    candidates=[row for row in data.get('source_files',[]) if row.get('source_name')==filename and str(row.get('report_path','')).endswith('.json')]
    bindings = data.get('figure_data_links', {}).get('bindings', [])
    direct = [r for r in bindings if r.get('panel_id') == panel_id and r.get('role') == 'methods' and str(r.get('path','')).endswith('.vl.json')]
    if direct:
        candidates = [{'report_path': r['path']} for r in direct]
    for row in candidates:
        spec=json.loads(member(root,row['report_path']).read_text('utf-8-sig'))
        comparisons=_comparisons(spec)
        cid=panel.get('contrast_id')
        if direct or comparisons=={cid} or (len(candidates)==1 and (not comparisons or cid in (None,'all'))):
            # External data URLs are never fetched by the local presentation.
            def check(value):
                if isinstance(value,dict):
                    if 'data' in value and isinstance(value['data'],dict) and 'url' in value['data']:raise ValueError('Figure requires external data; use its delivered image')
                    for item in value.values():check(item)
                elif isinstance(value,list):
                    for item in value:check(item)
            check(spec)
            return {'spec':spec,'source':row['report_path'],'panel':panel_id}
    return {'spec':None,'panel':panel_id,'image':panel.get('image'),'reason':'Only the delivered static figure can be bound unambiguously'}


def figure(context,panel_id):
    root=report_root(context)
    return _figure(str(root),(root/'report-manifest.json').stat().st_mtime_ns,panel_id)


@lru_cache(maxsize=4)
def _network(root, stamp):
    data=_manifest(root,stamp); views=data.get('network_browsers',[])
    if not views:raise ValueError('No completed CoSplice result delivered')
    path=member(Path(root),views[0]['source'])
    match=re.search(r'<script\s+type="application/json"\s+id="reportData">(.*?)</script>',path.read_text('utf-8'),re.S)
    if not match:raise ValueError('Network report has no structured result payload')
    payload=json.loads(match.group(1))
    payload['source_directory']=path.parent.relative_to(Path(root)).as_posix()
    return payload


def network(context):
    root=report_root(context)
    return _network(str(root),(root/'report-manifest.json').stat().st_mtime_ns)
