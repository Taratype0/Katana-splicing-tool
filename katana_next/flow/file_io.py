"""Avoid re-reading unchanged files and copying immutable task products twice."""
from __future__ import annotations

from collections import OrderedDict
import hashlib
import os
from pathlib import Path
import shutil
from threading import RLock

_digests = OrderedDict()
_lock = RLock()


def identity(path):
    stat = Path(path).stat()
    changed = stat.st_ctime_ns
    if os.name == 'nt':
        # Python's Windows st_ctime is creation time, not metadata change time.
        # FILE_BASIC_INFO.ChangeTime catches same-size writes even if a caller
        # restores LastWriteTime (mtime) afterwards.
        import ctypes
        import msvcrt
        class BasicInfo(ctypes.Structure):
            _fields_ = [(name, ctypes.c_longlong) for name in
                        ('creation', 'access', 'write', 'change')] + [('attributes', ctypes.c_uint32)]
        info = BasicInfo()
        api = ctypes.windll.kernel32.GetFileInformationByHandleEx
        api.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p, ctypes.c_uint32]
        api.restype = ctypes.c_int
        with Path(path).open('rb') as handle:
            if not api(msvcrt.get_osfhandle(handle.fileno()), 0, ctypes.byref(info), ctypes.sizeof(info)):
                raise OSError('Cannot read file change identity: ' + str(path))
        changed = info.change
    return [stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, changed]


def digest(path):
    path = Path(path).resolve(strict=True)
    before = identity(path)
    key = (("inode" if before[1] else str(path)), *before)
    with _lock:
        cached = _digests.get(key)
        if cached is not None and os.environ.get("KATANA_VERIFY_CONTENT") != "1":
            _digests.move_to_end(key)
            return cached
    with path.open('rb') as handle:
        value = hashlib.file_digest(handle, 'sha256').hexdigest()
    if identity(path) != before:
        raise ValueError(f'File changed during content verification: {path}')
    with _lock:
        _digests[key] = value
        while len(_digests) > 8192:
            _digests.popitem(last=False)
    return value


def retain_owned_file(source, destination):
    """Retain a completed, owned immutable product; never use for user inputs.

    A hard link survives removal of work. Across filesystems use an ordinary
    copy. Neither route leaves a symbolic link back to the execution directory.
    """
    source, destination = Path(source), Path(destination)
    if source.is_symlink() or not source.is_file():
        raise ValueError(f'Expected a regular completed product: {source}')
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    try:
        os.link(source, destination)
    except OSError:
        # copy2 overwrites by default; this destination belongs to a newly
        # created private tree and is checked above.
        shutil.copy2(source, destination)
    return str(destination)
