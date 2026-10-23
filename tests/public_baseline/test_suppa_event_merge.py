"""Header-only SUPPA classes must not hang or drop adjacent event rows."""
from pathlib import Path
import shutil
import subprocess

import pytest


@pytest.mark.skipif(shutil.which('awk') is None, reason='Native awk check runs on Linux')
@pytest.mark.parametrize('rows', [([], [], []), ([], ['a\tb'], []), (['a\tb'], [], ['c\td'])])
def test_real_module_merge_handles_empty_classes(tmp_path, rows):
    module = Path(__file__).resolve().parents[2] / 'pipeline/modules/local/suppa_generateevents.nf'
    command = next(line.strip() for line in module.read_text().splitlines() if line.strip().startswith('awk '))
    sources = []
    for i, values in enumerate(rows):
        path = tmp_path / f'events_{i}_strict.ioe'
        path.write_text('seqname\tevent\n' + ''.join(row + '\n' for row in values))
        sources.append((path, path.read_bytes(), path.stat().st_mtime_ns))
    subprocess.run(['sh', '-c', command], cwd=tmp_path, timeout=5, check=True)
    expected = 'seqname\tevent\n' + ''.join(row + '\n' for values in rows for row in values)
    assert (tmp_path / 'events.ioe').read_text() == expected
    assert all(path.read_bytes() == content and path.stat().st_mtime_ns == mtime for path, content, mtime in sources)
