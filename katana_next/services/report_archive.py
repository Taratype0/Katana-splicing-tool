"""Open directories directly; safely unpack portable reports into an owned cache."""
import json
import shutil
import stat
import tempfile
import os
from pathlib import Path, PurePosixPath, PureWindowsPath
from zipfile import ZipFile, BadZipFile

from katana_next.services.report_import import sha, member


def compact_copy(source, destination):
    """Storage-only export: retain the original presentation and every result byte."""
    from katana_next.reports.compact import compact_report
    from katana_next.reports.refresh import write_archive
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if destination.exists() or destination.is_relative_to(source) or source.is_relative_to(destination):
        raise ValueError('Choose a new report directory outside the original report')
    before = sha(source / 'report-manifest.json')
    manifest = json.loads((source / 'report-manifest.json').read_text('utf-8-sig'))
    if manifest.get('schema') != 'katana.approved_report/v1':
        raise ValueError('Unsupported report schema')
    files = [p for p in source.rglob('*') if p.is_file() and p.relative_to(source).as_posix() != 'REPORT.zip']
    for path in files:
        member(source, path.relative_to(source).as_posix())
    destination.mkdir(parents=True)
    for path in files:
        target = destination / path.relative_to(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        try:
            os.link(path, target)
        except OSError:
            shutil.copy2(path, target)
    audit = compact_report(destination)
    write_archive(destination)
    if sha(source / 'report-manifest.json') != before:
        raise ValueError('Original report changed while exporting')
    return audit


def report_source(location, cache_root):
    source = Path(location).expanduser().resolve()
    if source.is_dir():
        return source
    if source.name in {'report-manifest.json', 'REPORT.html'} and source.is_file():
        return source.parent
    if not source.is_file() or source.suffix.lower() != '.zip':
        raise ValueError('Choose an extracted report folder, report-manifest.json, or REPORT.zip')
    digest = sha(source)
    cache = Path(cache_root).resolve() / 'report-archives'
    cache.mkdir(parents=True, exist_ok=True)
    destination = cache / digest
    if destination.exists():
        receipt = json.loads((destination / 'extraction.json').read_text('utf-8'))
        if receipt['archive_sha256'] != digest:
            raise ValueError('Report archive cache identity changed')
        for name, expected in receipt['members'].items():
            if sha(member(destination / 'payload', name)) != expected:
                raise ValueError('Extracted archive changed; use a fresh workspace cache')
        return destination / 'payload' / receipt['report_root']
    with tempfile.TemporaryDirectory(prefix='extract-', dir=cache) as temporary:
        staging = Path(temporary)
        payload = staging / 'payload'
        payload.mkdir()
        try:
            with ZipFile(source) as archive:
                entries = archive.infolist()
                if len(entries) > 100000 or sum(i.file_size for i in entries) > 100 * 1024**3:
                    raise ValueError('Report exceeds local import limits (100,000 members / 100 GiB); extract and inspect it separately')
                required = sum(i.file_size for i in entries) + 256 * 1024**2
                if shutil.disk_usage(cache).free < required:
                    raise ValueError(f'Not enough space to extract this report: need at least {required / 1024**3:.1f} GiB in {cache}. Set KATANA_WORKSPACE_HOME to a larger drive or use an extracted report folder.')
                seen, checked = set(), []
                for info in entries:
                    name = info.filename.rstrip('/')
                    parts = name.split('/')
                    reserved = {'CON', 'PRN', 'AUX', 'NUL', *(f'COM{i}' for i in range(1, 10)), *(f'LPT{i}' for i in range(1, 10))}
                    if (not name or '\\' in name or ':' in name or '\x00' in name
                            or PurePosixPath(name).is_absolute() or PureWindowsPath(name).drive
                            or any(p in {'', '.', '..'} or p.endswith(('.', ' ')) or p.split('.')[0].upper() in reserved for p in parts)
                            or stat.S_ISLNK(info.external_attr >> 16) or info.flag_bits & 1):
                        raise ValueError('Unsafe or encrypted report archive member')
                    if name.casefold() in seen:
                        raise ValueError('Duplicate archive member')
                    seen.add(name.casefold())
                    checked.append((info, name))
                roots = [name for info, name in checked if not info.is_dir() and
                         (name == 'report-manifest.json' or (name.count('/') == 1 and name.endswith('/report-manifest.json')))]
                if len(roots) != 1:
                    raise ValueError('ZIP must contain one report-manifest.json at its root or in one enclosing folder')
                manifest = json.loads(archive.read(roots[0]))
                if manifest.get('schema') != 'katana.approved_report/v1':
                    raise ValueError('This is not a supported study report ZIP; target results use Import results')
                ledger = {}
                for info, name in checked:
                    target = payload.joinpath(*name.split('/'))
                    if info.is_dir():
                        target.mkdir(parents=True, exist_ok=True)
                        continue
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with archive.open(info) as incoming, target.open('xb') as outgoing:
                        shutil.copyfileobj(incoming, outgoing, 1024 * 1024)
                    ledger[name] = sha(target)
        except BadZipFile as exc:
            raise ValueError('Invalid or damaged report ZIP') from exc
        if sha(source) != digest:
            raise ValueError('Report ZIP changed during import')
        relative_root = str(PurePosixPath(roots[0]).parent)
        (staging / 'extraction.json').write_text(json.dumps({'archive_sha256': digest,
            'report_root': relative_root, 'members': ledger}), encoding='utf-8')
        # Only promote a fully extracted and checked cache. TemporaryDirectory cleans failures.
        try:
            staging.rename(destination)
        except FileExistsError:
            return report_source(source, cache_root)
    return destination / 'payload' / relative_root
