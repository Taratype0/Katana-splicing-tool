import json
from types import SimpleNamespace
from pathlib import Path
import pytest
from katana_next.services.report_presentation import presentation, figure, network

def context(tmp_path, data):
    (tmp_path/'report-manifest.json').write_text(json.dumps(data),encoding='utf-8')
    return SimpleNamespace(manifest=SimpleNamespace(metadata={'source_report_directory':str(tmp_path)}))

def test_figure_binding_preserves_exact_comparison_and_data(tmp_path):
    for cid in ('A','B'):
        (tmp_path/f'{cid}.vl.json').write_text(json.dumps({'$schema':'https://vega.github.io/schema/vega-lite/v5.json','data':{'values':[{'comparison_id':cid,'value':0.0}]},'mark':'point'}))
    data={'panels':[{'id':'A','title':'Motif','chapter':7,'contrast_id':'A','image':'same.png'}],
          'source_files':[{'source_name':'same.json','report_path':'B.vl.json'},{'source_name':'same.json','report_path':'A.vl.json'}],
          'figure_data_links':{'bindings':[{'panel_id':'A','role':'methods','path':'A.vl.json'}]}}
    ctx=context(tmp_path,data)
    assert figure(ctx,'A')['spec']['data']['values']==[{'comparison_id':'A','value':0.0}]
    assert presentation(ctx)['panels'][0]['chapter']==7
    with pytest.raises(ValueError,match='Unknown'):figure(ctx,'unknown')

def test_external_data_is_not_loaded_and_ambiguous_comparison_is_not_guessed(tmp_path):
    (tmp_path/'external.vl.json').write_text(json.dumps({'data':{'url':'https://invalid.example/data'},'mark':'point'}))
    ctx=context(tmp_path,{'panels':[{'id':'p','image':'external.png'}],'source_files':[{'source_name':'external.json','report_path':'external.vl.json'}]})
    with pytest.raises(ValueError,match='external'):figure(ctx,'p')

def test_network_reads_data_without_executing_html_and_keeps_all_members(tmp_path):
    payload={'modules':[{'module':'M1'}],'members':[{'event_key':f'e{i}'} for i in range(301)]}
    (tmp_path/'network.html').write_text('<script>throw Error("must not execute")</script><script type="application/json" id="reportData">'+json.dumps(payload)+'</script>')
    ctx=context(tmp_path,{'network_browsers':[{'source':'network.html'}]})
    result=network(ctx)
    assert len(result['members'])==301
    assert result['source_directory']=='.'

def test_bound_paths_cannot_leave_report_root(tmp_path):
    ctx=context(tmp_path,{'network_browsers':[{'source':'../outside.html'}]})
    with pytest.raises(ValueError):network(ctx)
