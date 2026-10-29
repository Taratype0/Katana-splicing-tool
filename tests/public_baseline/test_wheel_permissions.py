"""Regression for Windows-built wheels losing Unix pipeline execute bits."""
import hashlib
import os
from pathlib import Path
import stat
from zipfile import ZipFile, ZipInfo

import pytest


def test_windows_wheel_metadata_preserves_payload_and_record(tmp_path):
    from katana_build import normalize_pipeline_permissions
    wheel = tmp_path / 'test.whl'
    payloads = {
        'katana_next/flow/_pipeline/bin/katana_differential.R': b'#!/usr/bin/env Rscript\nprint(1)\n',
        'katana_next/flow/_pipeline/bin/tximport.R': b'#!/usr/bin/env Rscript\nprint(2)\n',
        'katana_next/flow/_pipeline/bin/shared.R': b'# sourced helper\n',
        'example.dist-info/RECORD': b'original RECORD is preserved\n',
    }
    with ZipFile(wheel, 'w') as archive:
        for name, data in payloads.items():
            info = ZipInfo(name)
            info.create_system = 0  # Simulate the Windows release builder.
            info.external_attr = (stat.S_IFREG | 0o644) << 16
            archive.writestr(info, data)
    changed = normalize_pipeline_permissions(wheel)
    assert len(changed) == 2
    with ZipFile(wheel) as archive:
        for name, data in payloads.items():
            assert archive.read(name) == data
            mode = archive.getinfo(name).external_attr >> 16
            assert bool(mode & 0o111) == data.startswith(b'#!')
    first = hashlib.sha256(wheel.read_bytes()).hexdigest()
    normalize_pipeline_permissions(wheel)
    assert hashlib.sha256(wheel.read_bytes()).hexdigest() == first


@pytest.mark.skipif(os.name != 'posix', reason='POSIX executable permission semantics')
def test_readonly_package_is_checked_without_rewriting_it(tmp_path):
    from katana_next.flow.distribution import validate_pipeline_executables
    root = tmp_path / 'installed'
    (root / 'bin').mkdir(parents=True)
    script = root / 'bin/tximport.R'
    script.write_text('#!/usr/bin/env Rscript\nprint(1)\n')
    script.chmod(0o444)
    before = script.stat()
    with pytest.raises(PermissionError, match='No jobs were submitted'):
        validate_pipeline_executables(root)
    assert script.stat().st_mtime_ns == before.st_mtime_ns
    assert stat.S_IMODE(script.stat().st_mode) == 0o444
    script.chmod(0o555)
    validate_pipeline_executables(root)


@pytest.mark.skipif(os.name != 'posix', reason='POSIX executable permission semantics')
def test_managed_source_rejects_broken_install_before_nextflow(tmp_path):
    from katana_next.flow.managed import _pipeline_source
    (tmp_path / 'bin').mkdir()
    script = tmp_path / 'bin/katana_differential.R'
    script.write_text('#!/usr/bin/env Rscript\n')
    script.chmod(0o644)
    main = tmp_path / 'main.nf'
    main.write_text('workflow {}\n')
    with pytest.raises(PermissionError, match='katana_differential.R'):
        _pipeline_source(main)
