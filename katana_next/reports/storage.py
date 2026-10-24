"""Optional immutable payload store shared by independent portable reports."""
from pathlib import Path
import os
import shutil
import tempfile
from .sources import digest


def copy_payload(source, destination, *, store=None, expected_sha=None):
    source, destination = Path(source), Path(destination)
    identity = expected_sha or digest(source)
    before = source.stat()
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if not destination.is_file() or digest(destination) != identity:
            raise ValueError('Refusing to replace an existing report payload: '+str(destination))
        return 'existing_verified'
    mode = 'copied'
    if store:
        root = Path(store).resolve()
        if destination.resolve().is_relative_to(root) or root.is_relative_to(destination.parent.resolve()):
            raise ValueError('Result store must be separate from the report')
        root.mkdir(parents=True,exist_ok=True)
        canonical = root / (identity + source.suffix.lower())
        if not canonical.exists():
            fd, temporary = tempfile.mkstemp(prefix='.payload-',dir=root)
            os.close(fd)
            temporary = Path(temporary)
            try:
                shutil.copyfile(source,temporary)
                if digest(temporary) != identity: raise ValueError('Source changed while populating result store')
                os.utime(temporary, ns=(before.st_atime_ns, before.st_mtime_ns))
                try: os.link(temporary,canonical)
                except FileExistsError: pass
            finally:
                temporary.unlink(missing_ok=True)
        if canonical.is_symlink() or digest(canonical) != identity:
            raise ValueError('Immutable result-store payload identity differs')
        try:
            if canonical.stat().st_mtime_ns != before.st_mtime_ns:
                raise OSError('A separate copy preserves this source timestamp without retiming shared files')
            os.link(canonical,destination)
            mode = 'store_hardlink'
        except OSError:
            shutil.copyfile(canonical,destination)
            os.utime(destination, ns=(before.st_atime_ns, before.st_mtime_ns))
            mode = 'store_copy'
    else:
        shutil.copyfile(source,destination)
        os.utime(destination, ns=(before.st_atime_ns, before.st_mtime_ns))
    if (source.stat().st_size,source.stat().st_mtime_ns)!=(before.st_size,before.st_mtime_ns) or digest(destination)!=identity:
        raise RuntimeError('Source changed during report payload collection')
    return mode
