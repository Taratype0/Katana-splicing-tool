"""Lossless storage for newly generated portable reports.

This pass changes storage and portable links only. Native result tables are
retained byte for byte; repeated copies are linked to a canonical copy. Complete
table-browser JSON can be gzip encoded without discarding a row or scalar.
Call before creating the report ZIP, in an isolated new output directory.
"""
from __future__ import annotations

import base64
import ast
from collections import defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
import html
import json
import os
from pathlib import Path
import re
import tempfile
from urllib.parse import quote, unquote, urlsplit, urlunsplit


_DATA_SUFFIXES = {'.tsv', '.csv', '.parquet', '.gz', '.png', '.fa', '.fasta', '.npy'}
_TEXT_SUFFIXES = {'.html', '.json', '.js', '.css'}
_JSON_STRING = re.compile(r'"(?:[^"\\]|\\.)*"')
_QUOTED_STRING = re.compile(r'''"(?:[^"\\]|\\.)*"|'(?:[^'\\]|\\.)*' ''', re.X)
_ATTR = re.compile(r'(?P<head>\b(?:href|src|data-src|poster)\s*=\s*)(?P<q>["\'])(?P<url>.*?)(?P=q)', re.I | re.S)
_SCRIPT = re.compile(r'(<script\b[^>]*>)(.*?)(</script\s*>)', re.I | re.S)
_READ_TABLE = re.compile(r"JSON\.parse\(\s*document\.getElementById\(\s*(['\"])table-data\1\s*\)\.textContent\s*\)")


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(4 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def _atomic_write(path, data):
    """Replace a directory entry, never modify an old hardlinked inode."""
    path = Path(path)
    fd, temporary = tempfile.mkstemp(prefix='.' + path.name + '.', suffix='.tmp', dir=path.parent)
    try:
        with os.fdopen(fd, 'wb') as handle:
            handle.write(data)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _inside(root, path):
    resolved = path.resolve()
    if not resolved.is_relative_to(root):
        raise ValueError('Portable report path escapes the report directory: ' + str(path))
    return resolved


def _canonical_key(path):
    preferred = ('outputs/results/', 'outputs/figure-data/', 'outputs/figures/', 'figures/')
    return next((i for i, p in enumerate(preferred) if path.startswith(p)), len(preferred)), len(path), path


def _format_suffix(path):
    suffixes = Path(path).suffixes
    return ''.join(suffixes[-2:]).lower() if suffixes and suffixes[-1].lower() == '.gz' else Path(path).suffix.lower()


def _aliases(root, files, supplied):
    """Only byte payloads are deduplicated; positional HTML/JSON are excluded."""
    if supplied is not None:
        if isinstance(supplied, (str, Path)):
            supplied = json.loads(Path(supplied).read_text(encoding='utf-8'))
        supplied = supplied.get('aliases', supplied)
        result = {}
        for alias, target in supplied.items():
            source, canonical = _inside(root, root / alias), _inside(root, root / target)
            if source == canonical or source.suffix.lower() not in _DATA_SUFFIXES:
                continue
            if _format_suffix(source) != _format_suffix(canonical):
                # An extension can determine a TSV/CSV/gzip decoder. Equal
                # bytes alone do not justify changing that reader contract.
                continue
            if not source.is_file() or not canonical.is_file():
                raise ValueError('Missing duplicate payload: ' + str(alias))
            if source.stat().st_size != canonical.stat().st_size or _sha(source) != _sha(canonical):
                raise ValueError('Duplicate payloads differ: ' + str(alias))
            result[source.relative_to(root).as_posix()] = canonical.relative_to(root).as_posix()
        # Reject chains/cycles; a canonical file must remain present.
        if any(target in result for target in result.values()):
            raise ValueError('Canonical payload aliases must not form chains or cycles.')
        return result
    candidates = defaultdict(list)
    for path in files:
        if path.suffix.lower() in _DATA_SUFFIXES:
            candidates[(path.stat().st_size, _format_suffix(path))].append(path)
    result = {}
    for paths in candidates.values():
        if len(paths) < 2:
            continue
        by_hash = defaultdict(list)
        inode_hash = {}
        for path in paths:
            st = path.stat()
            identity = st.st_dev, st.st_ino
            if identity not in inode_hash:
                inode_hash[identity] = _sha(path)
            by_hash[inode_hash[identity]].append(path.relative_to(root).as_posix())
        for matches in by_hash.values():
            ordered = sorted(matches, key=_canonical_key)
            for alias in ordered[1:]:
                result[alias] = ordered[0]
    return result


def _rewrite_link(value, document, root, aliases, alias_names=None):
    # Source-machine identities and URLs are provenance, not portable links.
    if not value or value.startswith(('#', '/', '\\')) or re.match(r'^[A-Za-z][A-Za-z0-9+.-]*:', value):
        return value
    try:
        parsed = urlsplit(value)
    except ValueError:
        return value
    if not parsed.path:
        return value
    decoded = unquote(parsed.path)
    if alias_names is not None and decoded.replace('\\', '/').rsplit('/', 1)[-1] not in alias_names:
        return value
    candidates = [(document.parent / decoded).resolve()]
    if decoded.startswith(('outputs/', 'assets/', 'figures/')):
        candidates.append((root / decoded).resolve())
    for path in candidates:
        if not path.is_relative_to(root):
            continue
        alias = path.relative_to(root).as_posix()
        if alias in aliases:
            relative = Path(os.path.relpath(root / aliases[alias], document.parent)).as_posix()
            if '%' in parsed.path:
                relative = quote(relative, safe='/-._~')
            return urlunsplit(('', '', relative, parsed.query, parsed.fragment))
    return value


def _portable_json_document(document, root, text):
    """Only generated navigational JSON can change portable path strings."""
    relative = document.relative_to(root)
    if relative.as_posix() in {'report-manifest.json', 'provenance.json', 'gene-set-inputs.json'}:
        return True
    # This is an explicitly generated portable index, distinct from the
    # author's native MAPP catalog and scientific JSON result payloads.
    return bool(re.search(r'"schema"\s*:\s*"katana\.portable_native_MAPP_fit_index/v1"', text))


def _opaque_data_script(opening_tag):
    """JSON/octet-stream script bodies are immutable data, not navigation."""
    content_type = re.search(r'(?<![\w-])type\s*=\s*(?:(["\'])(.*?)\1|([^\s>]+))', opening_tag, re.I | re.S)
    if content_type:
        value = (content_type[2] or content_type[3]).split(';', 1)[0].strip().lower()
        if value in {'application/json', 'application/octet-stream'} or value.endswith('+json'):
            return True
    return bool(re.search(r'\bid\s*=\s*["\']table-data["\']', opening_tag, re.I))


def _rewrite_text(text, document, root, aliases):
    if not aliases:
        return text
    alias_names = {Path(path).name for path in aliases}
    if document.suffix.lower() == '.json' and not _portable_json_document(document, root, text):
        # A producer's original file-index/path identities are evidence. They
        # remain verbatim; portable view/catalog URLs are rewritten elsewhere.
        return text
    def attribute(match):
        old = html.unescape(match['url'])
        new = _rewrite_link(old, document, root, aliases, alias_names)
        if old == new:
            return match[0]
        return match['head'] + match['q'] + html.escape(new, quote=True) + match['q']
    # Replacing only JSON string tokens preserves every original numeric token,
    # including tiny P/FDR values, and all scientific list/table structure.
    def token(match):
        # JSON keys (including original file-index keys) are identities.
        position = match.end()
        while position < len(text) and text[position].isspace():
            position += 1
        if position < len(text) and text[position] == ':':
            return match[0]
        try:
            old = json.loads(match[0])
        except (ValueError, TypeError):
            return match[0]
        new = _rewrite_link(old, document, root, aliases, alias_names)
        return json.dumps(new, ensure_ascii=False) if new != old else match[0]
    def quoted_token(match):
        try:
            old = json.loads(match[0]) if match[0][0] == '"' else ast.literal_eval(match[0])
        except (ValueError, SyntaxError):
            return match[0]
        new = _rewrite_link(old, document, root, aliases, alias_names)
        return json.dumps(new, ensure_ascii=False) if new != old else match[0]
    def chunk(value):
        if document.suffix.lower() in {'.html', '.css', '.js'}:
            value = _ATTR.sub(attribute, value)
            return _QUOTED_STRING.sub(quoted_token, value)
        return _JSON_STRING.sub(token, value)
    if document.suffix.lower() == '.html':
        # Scientific JSON payloads are evidence, not navigational links. Native
        # network browsers can contain millions of member/metric strings here.
        # Preserve every payload byte and do not tokenize its cells as code.
        parts, end = [], 0
        for script in _SCRIPT.finditer(text):
            if _opaque_data_script(script[1]):
                parts.extend((chunk(text[end:script.start()]), script[0]))
                end = script.end()
        parts.append(chunk(text[end:]))
        return ''.join(parts)
    return chunk(text)


_LOADER = """async function katanaReadTableData(){
const e=document.getElementById('table-data');
if(e.dataset.encoding!=='gzip-base64')return JSON.parse(e.textContent);
if(typeof DecompressionStream==='undefined')throw new Error('This complete table requires a current browser with gzip DecompressionStream. The original source-table download remains available.');
const binary=atob(e.textContent.trim()),bytes=new Uint8Array(binary.length);
for(let i=0;i<binary.length;i++)bytes[i]=binary.charCodeAt(i);
const stream=new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
const raw=await new Response(stream).arrayBuffer();
if(raw.byteLength!==Number(e.dataset.jsonBytes))throw new Error('Complete table decompressed size mismatch.');
if(globalThis.crypto&&crypto.subtle){const actual=Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',raw)),x=>x.toString(16).padStart(2,'0')).join('');if(actual!==e.dataset.jsonSha256)throw new Error('Complete table payload checksum mismatch.');}
return JSON.parse(new TextDecoder('utf-8',{fatal:true}).decode(raw));
}
"""


def _compress_table(text, minimum):
    scripts = list(_SCRIPT.finditer(text))
    payload = next((m for m in scripts if re.search(r'\bid\s*=\s*["\']table-data["\']', m[1], re.I)), None)
    if payload is None or 'data-encoding=' in payload[1]:
        return text, None
    consumers = [m for m in scripts if m != payload and _READ_TABLE.search(m[2])]
    raw = payload[2].encode('utf-8')
    if len(raw) < minimum or not consumers:
        return text, None
    # Validate JSON but retain its original UTF-8 bytes, never serialize rows.
    parsed = json.loads(payload[2])
    if not isinstance(parsed, dict) or not isinstance(parsed.get('rows'), list):
        return text, None
    compressed = gzip.compress(raw, compresslevel=6, mtime=0)
    encoded = base64.b64encode(compressed).decode('ascii')
    if len(encoded) >= len(raw):
        return text, None
    if gzip.decompress(base64.b64decode(encoded)) != raw:
        raise ValueError('Complete-table lossless encoding failed.')
    sha = hashlib.sha256(raw).hexdigest()
    new_payload = '<script id="table-data" type="application/octet-stream" data-encoding="gzip-base64" data-json-bytes="' + str(len(raw)) + '" data-json-sha256="' + sha + '">' + encoded + '</script>'
    replacements = [(payload.start(), payload.end(), new_payload)]
    for consumer in consumers:
        body = _READ_TABLE.sub('(await katanaReadTableData())', consumer[2])
        wrapped = '<script>' + _LOADER + '(async()=>{try{' + body + '}catch(error){window.tableLoadError=String(error);console.error(error);const e=document.getElementById("count");if(e)e.textContent=String(error);}})();</script>'
        replacements.append((consumer.start(), consumer.end(), wrapped))
    for start, end, replacement in sorted(replacements, reverse=True):
        text = text[:start] + replacement + text[end:]
    return text, dict(json_bytes=len(raw), json_sha256=sha, gzip_bytes=len(compressed), rows=len(parsed['rows']), columns=len(parsed.get('columns', [])))


def _refresh_portable_hashes(root, manifest, *, changed_paths=None, hash_cache=None):
    """Refresh portable bindings only; native/source SHA identities stay intact."""
    updates = []
    hash_cache = {} if hash_cache is None else hash_cache
    def actual_sha(path):
        if not isinstance(path, str):
            return None
        candidate = (root / path).resolve()
        if not candidate.is_relative_to(root) or not candidate.is_file():
            return None
        if candidate not in hash_cache:
            hash_cache[candidate] = _sha(candidate)
        return hash_cache[candidate]
    def metadata_bindings(value):
        metadata = value.get('metadata')
        if not isinstance(metadata, dict):
            return
        if 'portable_sha256' in metadata:
            path = value.get('source', value.get('path', value.get('view')))
            actual = actual_sha(path)
            previous = metadata['portable_sha256']
            if actual and actual != previous:
                metadata.setdefault('pre_storage_portable_sha256', previous)
                metadata['portable_sha256'] = actual
                updates.append(dict(path=path, binding='metadata.portable_sha256', previous_portable_sha256=previous,
                                    portable_sha256=actual, source_sha256=metadata.get('source_sha256')))
        if 'portable_fit_index_sha256' in metadata:
            candidates = {link['path'] for link in value.get('data_links', [])
                          if isinstance(link.get('path'), str) and link['path'].endswith('complete_native_fits_portable.json')}
            if len(candidates) != 1:
                raise ValueError('MAPP portable fit-index SHA needs one explicit portable index path.')
            path = candidates.pop()
            actual = actual_sha(path)
            previous = metadata['portable_fit_index_sha256']
            if actual and actual != previous:
                metadata.setdefault('pre_storage_portable_fit_index_sha256', previous)
                metadata['portable_fit_index_sha256'] = actual
                updates.append(dict(path=path, binding='metadata.portable_fit_index_sha256', previous_portable_sha256=previous,
                                    portable_sha256=actual, source_sha256=metadata.get('original_native_fit_index_sha256')))
    def walk(value):
        if isinstance(value, list):
            for item in value:
                walk(item)
        elif isinstance(value, dict):
            metadata_bindings(value)
            path = value.get('path', value.get('report_path'))
            if isinstance(path, str) and any(k in value for k in ('portable_sha256', 'source_sha256')):
                candidate = (root / path).resolve()
                if candidate.is_relative_to(root) and candidate.is_file():
                    relative = candidate.relative_to(root).as_posix()
                    if changed_paths is not None and relative not in changed_paths and re.fullmatch(r'[a-fA-F0-9]{64}', str(value.get('portable_sha256', ''))) and value.get('bytes') == candidate.stat().st_size:
                        # Retain an immediately preceding completed compaction
                        # verification, rather than re-reading immutable native
                        # GB tables for each metadata-only reconciliation.
                        for item in value.values():
                            walk(item)
                        return
                    actual = actual_sha(path)
                    previous = value.get('portable_sha256')
                    if 'report_path' in value and 'bytes' in value:
                        value.setdefault('source_bytes', value['bytes'])
                    value['portable_sha256'] = actual
                    value['bytes'] = candidate.stat().st_size
                    if actual != previous:
                        updates.append(dict(path=path, previous_portable_sha256=previous, portable_sha256=actual, source_sha256=value.get('source_sha256')))
            for item in value.values():
                walk(item)
    walk(manifest)
    return updates


def compact_report(report_dir, *, alias_map=None, deduplicate=True, compress_html=True,
                   html_min_bytes=32768, require_unarchived=True):
    """Compact a NEW completed report before archiving; return its storage audit.

    Existing native tables, statistics, figures and row identities are retained.
    Atomic replacements protect any original file hardlinked into this output.
    The caller owns this output directory; an existing ZIP is rejected by
    default to prevent accidentally compacting a previously delivered report.
    ``alias_map`` may be a verified mapping, audit dict or JSON path; omitted
    mappings are derived by SHA-256. Position-dependent HTML/JSON files are
    deliberately never deduplicated, even when their bytes happen to match.
    """
    root = Path(report_dir).expanduser().resolve()
    if not (root / 'REPORT.html').is_file() or not (root / 'report-manifest.json').is_file():
        raise ValueError('A completed report HTML and manifest are required.')
    if require_unarchived and (root / 'REPORT.zip').exists():
        raise ValueError('Compact only a new unarchived report; REPORT.zip already exists.')
    files = [p for p in root.rglob('*') if p.is_file() and p.name not in {'REPORT.zip', 'storage-audit.json'}]
    for path in files:
        _inside(root, path)
    before = sum(p.stat().st_size for p in files)
    aliases = _aliases(root, files, alias_map) if deduplicate else {}
    rewritten, tables = [], []
    for path in files:
        if path.suffix.lower() not in _TEXT_SUFFIXES:
            continue
        original = path.read_bytes()
        try:
            text = original.decode('utf-8')
        except UnicodeDecodeError:
            continue
        text = _rewrite_text(text, path, root, aliases)
        table = None
        if compress_html and path.suffix.lower() == '.html':
            text, table = _compress_table(text, html_min_bytes)
        data = text.encode('utf-8')
        if data != original:
            _atomic_write(path, data)
            rewritten.append(dict(path=path.relative_to(root).as_posix(), previous_sha256=hashlib.sha256(original).hexdigest(), portable_sha256=hashlib.sha256(data).hexdigest(), before_bytes=len(original), after_bytes=len(data)))
        if table:
            tables.append(dict(path=path.relative_to(root).as_posix(), **table))
    # A lost reference is worse than a large report. Check supported portable
    # text references after rewriting, before deleting any aliases.
    for path in files:
        if path.suffix.lower() not in _TEXT_SUFFIXES:
            continue
        try:
            text = path.read_text(encoding='utf-8')
        except UnicodeDecodeError:
            continue
        if _rewrite_text(text, path, root, aliases) != text:
            raise ValueError('Unresolved duplicate link in ' + str(path))
    removed = []
    for alias, canonical in aliases.items():
        path = root / alias
        size = path.stat().st_size
        # Content identity was checked before any aliases could be removed.
        path.unlink()
        removed.append(dict(path=alias, canonical=canonical, bytes=size))
    # Generated report provenance contains portable bindings in addition to
    # immutable source identities; update those bindings as well. The native
    # producer proofs under outputs/provenance were left byte-identical above.
    provenance_path = root / 'provenance.json'
    bindings = []
    if provenance_path.is_file():
        provenance = json.loads(provenance_path.read_text(encoding='utf-8'))
        bindings.extend(dict(document='provenance.json', **item) for item in _refresh_portable_hashes(root, provenance))
        provenance['portable_storage'] = dict(schema='katana.report_storage/v1', audit='storage-audit.json', lossless=True)
        _atomic_write(provenance_path, json.dumps(provenance, ensure_ascii=False, indent=2).encode('utf-8'))
    manifest_path = root / 'report-manifest.json'
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    bindings.extend(dict(document='report-manifest.json', **item) for item in _refresh_portable_hashes(root, manifest))
    # Source identities and P/FDR stay untouched. This metadata is explicitly
    # a storage transformation, not a recomputation or a new source receipt.
    manifest['storage'] = dict(schema='katana.report_storage/v1', audit='storage-audit.json', lossless=True,
                               native_result_tables='unchanged bytes; canonical copies retained',
                               complete_table_browsers='all original rows; gzip-base64 embedded JSON')
    _atomic_write(manifest_path, json.dumps(manifest, ensure_ascii=False, indent=2).encode('utf-8'))
    remaining = [p for p in root.rglob('*') if p.is_file() and p.name not in {'REPORT.zip', 'storage-audit.json'}]
    after = sum(p.stat().st_size for p in remaining)
    audit = dict(schema='katana.report_storage_audit/v1', created_utc=datetime.now(timezone.utc).isoformat(),
                 lossless=True, native_statistics_recomputed=False, before_bytes=before, after_bytes=after,
                 saved_bytes=before-after, deleted_duplicate_aliases=removed, rewritten_files=rewritten,
                 complete_table_payloads=tables, portable_binding_updates=bindings,
                 size_scope='all report files except REPORT.zip and storage-audit.json',
                 note='Logical portable payload sizes; hardlinked storage may be shared with source directories.')
    _atomic_write(root / 'storage-audit.json', json.dumps(audit, ensure_ascii=False, indent=2).encode('utf-8'))
    return audit


def reconcile_compact_metadata(report_dir, storage_audit=None, *, source_report_dir=None,
                               rehash_unchanged=True):
    """Reconcile generated portable metadata after a completed compaction.

    Alias payloads may already have been removed. Their recorded canonical
    targets must exist; native/source identities are not reinterpreted. An
    optional immutable source report can restore a scientific JSON modified by
    an older compactor, but only when its original SHA matches that compactor's
    recorded pre-change SHA. No scientific process is run or source written.
    """
    root = Path(report_dir).expanduser().resolve()
    if storage_audit is None:
        storage_audit = root / 'storage-audit.json'
    if isinstance(storage_audit, (str, Path)):
        audit = json.loads(Path(storage_audit).read_text(encoding='utf-8'))
    else:
        audit = dict(storage_audit)
    aliases = {item['path']: item['canonical'] for item in audit.get('deleted_duplicate_aliases', [])}
    for alias, canonical in aliases.items():
        _inside(root, root / alias)
        target = _inside(root, root / canonical)
        if not target.is_file():
            raise ValueError('Recorded canonical payload is missing: ' + canonical)
    restored = []
    if source_report_dir is not None:
        original_root = Path(source_report_dir).expanduser().resolve()
        if original_root == root:
            raise ValueError('Immutable source report must differ from output report.')
        for item in audit.get('rewritten_files', []):
            path = _inside(root, root / item['path'])
            if path.suffix.lower() != '.json' or not path.is_file():
                continue
            text = path.read_text(encoding='utf-8')
            if _portable_json_document(path, root, text):
                continue
            original = _inside(original_root, original_root / item['path'])
            if not original.is_file():
                continue
            sha = _sha(original)
            if sha != item['previous_sha256']:
                raise ValueError('Native JSON source does not match recorded pre-compaction SHA: ' + item['path'])
            current = _sha(path)
            if current != sha:
                _atomic_write(path, original.read_bytes())
                restored.append(dict(path=item['path'], previous_portable_sha256=current,
                                     restored_source_sha256=sha, source_report=str(original_root)))
    generated = [root / name for name in ('gene-set-inputs.json', 'provenance.json', 'report-manifest.json')]
    generated.extend((root / 'outputs' / 'figure-data').glob('*complete_native_fits_portable.json'))
    rewritten = []
    for path in generated:
        if not path.is_file():
            continue
        original = path.read_bytes()
        text = _rewrite_text(original.decode('utf-8'), path, root, aliases)
        if text.encode('utf-8') != original:
            _atomic_write(path, text.encode('utf-8'))
            rewritten.append(dict(path=path.relative_to(root).as_posix(), previous_sha256=hashlib.sha256(original).hexdigest(), portable_sha256=_sha(path)))
    portable_updates = []
    changed_paths = None if rehash_unchanged else {item['path'] for item in rewritten + restored}
    shared_hashes = {}
    # Provenance first so a manifest binding to it sees its final bytes.
    for name in ('provenance.json', 'report-manifest.json'):
        path = root / name
        if not path.is_file():
            continue
        value = json.loads(path.read_text(encoding='utf-8'))
        portable_updates.extend(dict(document=name, **item) for item in _refresh_portable_hashes(root, value, changed_paths=changed_paths, hash_cache=shared_hashes))
        value['portable_storage'] = dict(schema='katana.report_storage/v1', audit='storage-audit.json', lossless=True)
        _atomic_write(path, json.dumps(value, ensure_ascii=False, indent=2).encode('utf-8'))
        if changed_paths is not None:
            changed_paths.add(name)
    receipt = dict(schema='katana.report_compact_metadata_reconciliation/v1', created_utc=datetime.now(timezone.utc).isoformat(),
                   rewritten_generated_metadata=rewritten, portable_binding_updates=portable_updates,
                   restored_native_json=restored, native_statistics_recomputed=False,
                   unchanged_native_hashes='reverified' if rehash_unchanged else 'retained from immediately preceding completed compaction proof; matching portable SHA and byte count required')
    audit.setdefault('metadata_reconciliation', []).append(receipt)
    audit['after_bytes'] = sum(path.stat().st_size for path in root.rglob('*') if path.is_file() and path.name not in {'REPORT.zip', 'storage-audit.json'})
    audit['saved_bytes'] = audit['before_bytes'] - audit['after_bytes']
    _atomic_write(root / 'storage-audit.json', json.dumps(audit, ensure_ascii=False, indent=2).encode('utf-8'))
    return receipt


def _inspect_encoded_table(path, source_root, root):
    text = path.read_text(encoding='utf-8')
    if 'data-encoding="gzip-base64"' not in text:
        return None
    block = next((m for m in _SCRIPT.finditer(text) if re.search(r'\bid\s*=\s*["\']table-data["\']', m[1], re.I)), None)
    if block is None or 'data-encoding="gzip-base64"' not in block[1]:
        return None
    sha_match = re.search(r'data-json-sha256="([a-f0-9]{64})"', block[1])
    bytes_match = re.search(r'data-json-bytes="(\d+)"', block[1])
    if not sha_match or not bytes_match:
        raise ValueError('Encoded complete table lacks its lossless proof: ' + str(path))
    packed = base64.b64decode(block[2].strip(), validate=True)
    raw = gzip.decompress(packed)
    sha = hashlib.sha256(raw).hexdigest()
    if sha != sha_match[1] or len(raw) != int(bytes_match[1]):
        raise ValueError('Interrupted table payload proof failed: ' + str(path))
    original = source_root / path.relative_to(root)
    original_verified = False
    if original.is_file():
        before_text = original.read_text(encoding='utf-8')
        before = next((m for m in _SCRIPT.finditer(before_text) if re.search(r'\bid\s*=\s*["\']table-data["\']', m[1], re.I)), None)
        if before is not None and 'data-encoding=' not in before[1]:
            if hashlib.sha256(before[2].encode('utf-8')).hexdigest() != sha:
                raise ValueError('Complete-table rows differ from the immutable source report: ' + str(path))
            original_verified = True
    payload = json.loads(raw)
    return dict(path=path.relative_to(root).as_posix(), json_bytes=len(raw), json_sha256=sha,
                gzip_bytes=len(packed), rows=len(payload['rows']), columns=len(payload.get('columns', [])),
                preserved_from_interrupted_pass=True, immutable_source_payload_verified=original_verified)


def recover_compact_report(report_dir, source_report_dir, alias_map, *, progress=None):
    """Finish storage-only compaction after an interrupted presentation pass.

    The source is read-only. Previously encoded complete tables are verified
    and indexed; scientific JSON can be restored from matching source paths.
    Partial starting size and the older source-report size are kept separate.
    No original archive is deleted and no analysis or ZIP build is performed.
    """
    root = Path(report_dir).expanduser().resolve()
    original_root = Path(source_report_dir).expanduser().resolve()
    if root == original_root or not (original_root / 'REPORT.html').is_file():
        raise ValueError('Recovery needs a distinct immutable source report.')
    if (root / 'REPORT.zip').exists():
        raise ValueError('Recover only the new, unarchived output report.')
    emit = progress or (lambda value: None)
    def size(folder):
        return sum(p.stat().st_size for p in folder.rglob('*') if p.is_file() and p.name not in {'REPORT.zip', 'storage-audit.json'})
    partial_start = size(root)
    original_bytes = size(original_root)
    emit(dict(phase='started', partial_start_bytes=partial_start, immutable_source_payload_bytes=original_bytes))
    restored = []
    for path in root.rglob('*.json'):
        _inside(root, path)
        if _portable_json_document(path, root, path.read_text(encoding='utf-8')):
            continue
        relative = path.relative_to(root)
        original = _inside(original_root, original_root / relative)
        if not original.is_file():
            continue
        before_sha, original_sha = _sha(path), _sha(original)
        if before_sha != original_sha:
            _atomic_write(path, original.read_bytes())
            restored.append(dict(path=relative.as_posix(), interrupted_portable_sha256=before_sha,
                                 restored_source_sha256=original_sha, basis='same native relative source path; immutable source bytes verified'))
    emit(dict(phase='native_json_preserved', restored=len(restored)))
    encoded = []
    for path in root.rglob('*.html'):
        _inside(root, path)
        record = _inspect_encoded_table(path, original_root, root)
        if record:
            encoded.append(record)
    emit(dict(phase='encoded_complete_tables_verified', tables=len(encoded), rows=sum(x['rows'] for x in encoded)))
    supplied = json.loads(Path(alias_map).read_text(encoding='utf-8')) if isinstance(alias_map, (str, Path)) else alias_map
    supplied = supplied.get('aliases', supplied)
    present, already_removed = {}, []
    for alias, canonical in supplied.items():
        path, target = _inside(root, root / alias), _inside(root, root / canonical)
        if path.suffix.lower() not in _DATA_SUFFIXES or _format_suffix(path) != _format_suffix(target):
            continue
        if path.is_file():
            present[alias] = canonical
        else:
            source_alias = _inside(original_root, original_root / alias)
            if not target.is_file() or not source_alias.is_file() or _sha(target) != _sha(source_alias):
                raise ValueError('Previously removed alias lacks immutable byte proof: ' + alias)
            already_removed.append(dict(path=alias, canonical=canonical, bytes=source_alias.stat().st_size))
    audit = compact_report(root, alias_map=present)
    emit(dict(phase='fast_compaction_finished', deleted_aliases=len(audit['deleted_duplicate_aliases']), newly_encoded_tables=len(audit['complete_table_payloads'])))
    audit['deleted_duplicate_aliases'] = already_removed + audit['deleted_duplicate_aliases']
    complete = {item['path']: item for item in encoded}
    complete.update({item['path']: item for item in audit['complete_table_payloads']})
    audit['complete_table_payloads'] = list(complete.values())
    audit['recovery'] = dict(schema='katana.interrupted_report_storage_recovery/v1', created_utc=datetime.now(timezone.utc).isoformat(),
                             partial_start_bytes=partial_start, immutable_source_payload_bytes=original_bytes,
                             source_report=str(original_root), restored_native_json=restored,
                             already_encoded_tables=len(encoded), source_report_untouched=True,
                             scientific_results_rerun=0, zip_built=False)
    _atomic_write(root / 'storage-audit.json', json.dumps(audit, ensure_ascii=False, indent=2).encode('utf-8'))
    reconcile_compact_metadata(root, audit, rehash_unchanged=False)
    final = json.loads((root / 'storage-audit.json').read_text(encoding='utf-8'))
    emit(dict(phase='completed', after_bytes=final['after_bytes'], tables=len(final['complete_table_payloads']),
              rows=sum(x['rows'] for x in final['complete_table_payloads']), deleted_aliases=len(final['deleted_duplicate_aliases'])))
    return final
