"""Product outer refits must repeat the exact optional hierarchy operation."""
from pathlib import Path
from validation24 import deletion_draws,fit_common,summarize
from validation22 import labels_at
from module_merge import merge
from cosplice_research import write_json

def hierarchical_stability(raw,meta,support,assignment,keys,out,iterations):
    out=Path(out);out.mkdir(exist_ok=False);draws=deletion_draws(meta,iterations)
    write_json(out/'planned_draws.json',{'training_sample_draws':draws,'seed':24031,'mergeCloseModules_cutHeight':.25});records=[]
    for i,cols in enumerate(draws):
        target=out/f'refit_{i:03}';result,_,_=fit_common('jackknife',raw.loc[:,cols],meta,support.loc[:,cols] if support is not None else None,target)
        if result['status']=='completed':
            merged=merge(target/'model',target/'merged_model');result['hierarchy_status']=merged
            if merged['status']!='completed':result.update(status='failed',error=merged.get('error'))
        labels=labels_at(target/'merged_model',keys) if result['status']=='completed' else None;records.append((i,labels,result))
    return summarize(assignment.reindex(keys).module.to_numpy(),keys,records,out)
