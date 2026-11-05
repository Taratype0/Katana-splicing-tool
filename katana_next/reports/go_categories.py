"""Offline namespace annotation for report copies of completed GO statistics."""
import gzip
import json
from functools import lru_cache
from pathlib import Path


@lru_cache(maxsize=1)
def category_index():
    path = Path(__file__).parent/'assets/go-categories.json.gz'
    data = json.loads(gzip.decompress(path.read_bytes()))
    if data.get('schema') != 'katana.report_go_categories/v1': raise ValueError('Unknown GO category index')
    return data


def annotate_completed_categories(frame):
    """Preserve explicitly recorded categories; fill only missing GO namespaces."""
    data = frame.copy(deep=True)
    data['category'] = data.ontology_category.replace({'GO BP':'BP','GO CC':'CC','GO MF':'MF'})
    missing = ~data.category.isin(['BP','CC','MF','KEGG'])
    identifiers = data.term.astype(str).str.extract(r'(?:^|/\s*)(GO:\d{7})(?:\s|$)',expand=False)
    needed = missing & identifiers.notna()
    if not needed.any(): return data, None
    index = category_index()
    assigned = identifiers.loc[needed].map(index['categories'])
    resolved = assigned.dropna().index
    data.loc[resolved, 'category'] = assigned.loc[resolved]
    audit = {k:v for k,v in index.items() if k!='categories'}
    audit.update(rows_resolved=len(resolved),unresolved_identifiers=sorted(set(identifiers.loc[needed & ~data.index.isin(resolved)])),
                 statistics_recomputed=False, original_categories_changed=False)
    return data, audit
