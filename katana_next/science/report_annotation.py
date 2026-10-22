"""Non-destructive gene annotation for report tables."""
import re
import pandas as pd


def annotate_exons(frame, annotation):
    result = frame.copy()
    symbol = 'gene_symbol' if 'gene_symbol' in annotation else 'gene_name'
    mapping = {}
    for _, row in annotation.dropna(subset=['gene_id', symbol]).iterrows():
        key = str(row.gene_id).split('.')[0]
        mapping.setdefault(key, set()).add(str(row[symbol]))
    def resolve(value):
        ids = re.findall(r'ENS[A-Z]*G\d+(?:\.\d+)?', str(value))
        return '; '.join('/'.join(sorted(mapping.get(g.split('.')[0], {g}))) for g in dict.fromkeys(ids))
    source = 'groupID' if 'groupID' in result else 'bin_id'
    result['gene_name'] = result[source].map(resolve)
    return result
