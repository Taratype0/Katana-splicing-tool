"""The same portable native-runtime check for CLI and Nextflow."""
import os
import shutil
import subprocess
import sys


def check_runtime():
    from katana_next.flow.container_runtime import sanitize_container_environment, require_image_path
    sanitize_container_environment()
    from katana_next.flow.dependency_checks import require_python_imports
    require_python_imports(sys.executable, ('numpy','pandas','scipy','statsmodels','sklearn','psutil'))
    rscript = os.environ.get('KATANA_RSCRIPT') or shutil.which('Rscript')
    if not rscript:
        raise ValueError('Rscript is missing from the selected environment; use the Katana container or standalone Conda environment')
    rscript = require_image_path(rscript, 'Rscript')
    try:
        check = subprocess.run([rscript,'--vanilla','-e','stopifnot(requireNamespace("WGCNA",quietly=TRUE))'],
                               capture_output=True,text=True,timeout=60)
    except subprocess.TimeoutExpired as exc:
        raise ValueError('The R/WGCNA environment probe did not respond within 60 seconds; no network calculation was started') from exc
    if check.returncode:
        raise ValueError('WGCNA failed to load in the selected R runtime: '+check.stderr[-1000:])
    return rscript
