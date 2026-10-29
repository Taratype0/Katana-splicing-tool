"""Build-only wheel metadata normalization, including Windows build hosts."""
from pathlib import Path, PurePosixPath
import os
import stat
import tempfile
from zipfile import ZipFile


def normalize_pipeline_permissions(wheel):
    """Preserve payload/RECORD bytes while marking bundled shebang scripts executable."""
    wheel = Path(wheel)
    fd, temporary = tempfile.mkstemp(prefix=wheel.name + '.', suffix='.tmp', dir=wheel.parent)
    os.close(fd)
    changed = []
    try:
        with ZipFile(wheel) as source, ZipFile(temporary, 'w') as destination:
            destination.comment = source.comment
            for info in source.infolist():
                data = source.read(info.filename)
                path = PurePosixPath(info.filename)
                is_bin = path.parent.as_posix() == 'katana_next/flow/_pipeline/bin'
                is_fixture = path.as_posix() == 'katana_next/flow/_pipeline/tests/protocol_fixture.py'
                if (is_bin or is_fixture) and data.startswith(b'#!'):
                    info.create_system = 3  # Unix permissions, regardless of build OS.
                    info.external_attr = ((stat.S_IFREG | 0o755) << 16) | (info.external_attr & 0xFFFF)
                    changed.append(info.filename)
                destination.writestr(info, data)
        if not changed:
            raise ValueError('Wheel contains no executable pipeline scripts')
        os.replace(temporary, wheel)
    finally:
        Path(temporary).unlink(missing_ok=True)
    return changed
