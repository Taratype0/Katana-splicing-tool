"""Traceable display names from the report's own unambiguous identifier maps.

Original result columns and files are never rewritten. Ambiguous aliases are
not guessed, and unknown IDs remain visible as IDs.
"""
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import re
import pandas as pd


def gene_token(value):
    return re.sub(r'\.\d+$', '', str(value).strip().strip('"'))


@lru_cache(maxsize=8)
def _mapping(receipt_path, modified):
    receipt = Path(receipt_path)
    records = json.loads(receipt.read_text('utf-8'))['tables']
    candidates = {}
    for entry in records:
        source = entry.get('source_relative_path', '')
        if not source.endswith('gene_identifier_mapping.parquet'):
            continue
        # Locate the immutable imported copy by its independently recorded hash.
        # Its filename is not part of the public receipt contract.
        candidates[entry['imported_table_sha256']] = source
    pairs = {}
    if not candidates:
        return pd.DataFrame({'_name_id': pd.Series(dtype=str), '_mapped_name': pd.Series(dtype=str)})
    for path in (receipt.parent / 'results' / 'imported').glob('*.parquet'):
        import pyarrow.parquet as pq
        if not {'alias', 'canonical_gene_ids', 'mapping_status'} <= set(pq.read_schema(path).names):
            continue
        digest = hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()
        if digest not in candidates:
            continue
        frame = pd.read_parquet(path, columns=['alias', 'canonical_gene_ids', 'mapping_status'])
        for alias, ids, status in frame.itertuples(index=False, name=None):
            symbol, identifier = gene_token(alias), gene_token(ids)
            if status != 'mapped' or not identifier or ';' in identifier or not symbol or symbol.startswith('ENS'):
                continue
            pairs.setdefault(identifier, set()).add(symbol)
    return pd.DataFrame([(identifier, next(iter(names))) for identifier, names in sorted(pairs.items()) if len(names) == 1],
                        columns=['_name_id', '_mapped_name'])


def display_names(artifact):
    receipt = artifact.path.parent.parent.parent / 'report-import.json'
    if artifact.analysis_id != 'imported.report' or not receipt.is_file():
        return pd.DataFrame({'_name_id': pd.Series(dtype=str), '_mapped_name': pd.Series(dtype=str)})
    return _mapping(str(receipt), receipt.stat().st_mtime_ns)
