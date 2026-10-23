"""Export the exact packaged Nextflow sources to a new execution directory."""
from pathlib import Path
import shutil


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
    return destination
