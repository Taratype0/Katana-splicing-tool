"""Content verification regression; generated files, never scientific reruns."""
import hashlib
import os
from pathlib import Path

import pytest

from katana_next.flow import file_io, runtime


def test_record_accepts_digest_verified_metadata_transition(tmp_path, monkeypatch):
    path = tmp_path / 'dtu_features.tsv'
    path.write_bytes(b'transcript_id\tfdr\nT1\t0.02\n')
    original = file_io.hashlib.file_digest
    calls = []
    def read_and_touch(handle, algorithm):
        result = original(handle, algorithm)
        calls.append(1)
        if len(calls) == 1:
            stamp = path.stat()
            os.utime(path, ns=(stamp.st_atime_ns, stamp.st_mtime_ns + 1000000000))
        return result
    monkeypatch.setenv('KATANA_VERIFY_CONTENT', '1')
    monkeypatch.setattr(file_io.hashlib, 'file_digest', read_and_touch)
    record = runtime.file_record(path)
    assert record['sha256'] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert record['file_identity'] == file_io.identity(path)
    runtime.verify_files([record])


@pytest.mark.parametrize('change', ['bytes', 'inode', 'size', 'unstable_metadata'])
def test_record_rejects_unverified_transition(tmp_path, monkeypatch, change):
    path = tmp_path / 'dtu_features.tsv'; path.write_bytes(b'abcd')
    original = runtime.file_hash
    count = 0
    def race(target):
        nonlocal count
        result = original(target)
        count += 1
        if count == 1 or change == 'unstable_metadata':
            if change == 'bytes':
                path.write_bytes(b'wxyz')
            elif change == 'size':
                path.write_bytes(b'larger')
            elif change == 'inode':
                other = path.with_suffix('.new'); other.write_bytes(b'abcd'); other.replace(path)
            stamp = path.stat()
            os.utime(path, ns=(stamp.st_atime_ns, stamp.st_mtime_ns + 1000000000))
        return result
    monkeypatch.setenv('KATANA_VERIFY_CONTENT', '1')
    monkeypatch.setattr(runtime, 'file_hash', race)
    with pytest.raises(ValueError, match='changed'):
        runtime.file_record(path)


def test_unchanged_record_does_not_add_another_digest(tmp_path, monkeypatch):
    path = tmp_path / 'stable.tsv'; path.write_bytes(b'stable')
    original = runtime.file_hash
    calls = []
    def counted(target):
        calls.append(1)
        return original(target)
    monkeypatch.setattr(runtime, 'file_hash', counted)
    runtime.file_record(path)
    assert len(calls) == 1
