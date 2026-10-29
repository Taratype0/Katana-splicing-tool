"""Export the exact packaged Nextflow sources to a new execution directory."""
from pathlib import Path
import os
import shutil


def validate_pipeline_executables(pipeline):
    """Fail before submitting jobs if an installation/copy lost executable bits."""
    root = Path(pipeline)
    if os.name != 'posix' or not (root / 'bin').is_dir():
        return
    missing = []
    for script in sorted((root / 'bin').iterdir()):
        if script.is_file():
            with script.open('rb') as handle:
                shebang = handle.read(2) == b'#!'
            if shebang and not os.access(script, os.X_OK):
                missing.append(script.name)
    if missing:
        raise PermissionError('Installed pipeline scripts are not executable: ' + ', '.join(missing) +
                              '. Install the corrected release or use katana-flow export-pipeline '
                              'and pass its main.nf with --pipeline. No jobs were submitted.')


def export_pipeline(outdir):
    package = Path(__file__).resolve().parent
    source = package / "_pipeline"
    if not (source / "main.nf").is_file():
        source = package.parents[1] / "pipeline"
    if not (source / "main.nf").is_file():
        raise FileNotFoundError("This installation does not contain the integrated Nextflow workflow")
    destination = Path(outdir).expanduser().resolve()
    if destination.exists():
        raise FileExistsError("Choose a new pipeline export directory; existing versions are preserved")
    shutil.copytree(source, destination)
    for script in (destination / "bin").iterdir():
        if script.is_file() and script.suffix in {".py", ".R", ".awk"}:
            script.chmod(script.stat().st_mode | 0o111)
    fixture = destination / "tests" / "protocol_fixture.py"
    if fixture.is_file():
        fixture.chmod(fixture.stat().st_mode | 0o111)
    validate_pipeline_executables(destination)
    return destination
