"""Build a reproducible, display-only GO ID to namespace index from official OBO."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def build(obo, source_url, out):
    content = Path(obo).read_bytes()
    categories, aliases, current = {}, {}, {}
    namespaces = {'biological_process': 'BP', 'cellular_component': 'CC', 'molecular_function': 'MF'}
    version = None
    for line in content.decode('utf8').splitlines() + ['[end]']:
        if line.startswith('data-version: '): version = line.partition(': ')[2]
        if line.startswith('['):
            if current.get('id', '').startswith('GO:') and current.get('namespace') in namespaces:
                category = namespaces[current['namespace']]
                categories[current['id']] = category
                for alias in current.get('alt_id', []):
                    if alias in aliases and aliases[alias] != category: raise ValueError('Conflicting GO namespace alias')
                    aliases[alias] = category
            current = {}
        elif line.startswith('alt_id: '): current.setdefault('alt_id', []).append(line[8:])
        elif line.startswith(('id: ', 'namespace: ')):
            key, value = line.split(': ', 1); current[key] = value
    for alias, category in aliases.items():
        if alias in categories and categories[alias] != category: raise ValueError('Conflicting GO namespace')
        categories[alias] = category
    if not version or len(categories) < 10000: raise ValueError('Incomplete GO ontology')
    data = dict(schema='katana.report_go_categories/v1', ontology_version=version,
        source_url=source_url, source_sha256=hashlib.sha256(content).hexdigest(),
        attribution='Gene Ontology Consortium', license='CC BY 4.0',
        license_url='https://geneontology.org/docs/go-citation-policy/',
        purpose='Display categories only; no gene membership, enrichment, propagation or P/FDR calculation.',
        alias_policy='Official alt_id namespace lookup; original result identifiers remain unchanged.',
        categories=dict(sorted(categories.items())))
    encoded=json.dumps(data,sort_keys=True,separators=(',', ':')).encode()
    Path(out).write_bytes(gzip.compress(encoded,mtime=0))
    print(json.dumps({k:v for k,v in data.items() if k!='categories'}|{'terms':len(categories),'bytes':Path(out).stat().st_size}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--obo',required=True);p.add_argument('--source-url',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();build(a.obo,a.source_url,a.out)
