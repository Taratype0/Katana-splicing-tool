"""Import-only checks in the interpreters used by enabled scientific stages.

No analysis, dependency installation or input-file access is performed here.
Successful checks are cached only for this process and executable/environment.
"""
from __future__ import annotations

from functools import lru_cache
import json
import os
from pathlib import Path
import re
import subprocess
import sys


_MARKER = "KATANA_PYTHON_IMPORTS="
_IMPORT_SCRIPT = r'''
import importlib, json, sys
failures = []
for name in json.loads(sys.argv[1]):
    try:
        importlib.import_module(name)
    except Exception as exc:
        failures.append({"module": name, "error": type(exc).__name__ + ": " + str(exc)})
print("KATANA_PYTHON_IMPORTS=" + json.dumps({"executable": sys.executable, "failures": failures}))
sys.exit(1 if failures else 0)
'''

# These are dependencies actually imported by the adapter's Python code. Native
# R tools and reference resources retain their own readiness checks in stages.py.
_HOST_IMPORTS = {
    "as.input_projection": ("numpy", "pandas", "pyarrow"),
    "mapp.prepare": ("numpy", "pandas", "Bio.Seq", "Bio.SeqIO"),
    "mapp.maei": ("numpy", "pandas"),
    "network.fit": ("numpy", "pandas", "scipy", "sklearn", "networkx", "igraph", "dynamicTreeCut"),
    "network.bootstrap": ("numpy", "pandas", "scipy", "sklearn", "networkx", "igraph", "dynamicTreeCut"),
    "network.associate": ("numpy", "pandas", "scipy", "statsmodels"),
    "network.annotate": ("numpy", "pandas", "scipy"),
}
# Official MAEI Python postprocessors import numpy/pandas. The kmer counter
# itself uses only the standard library; Bio is required in the Katana host by
# mapp_inputs -> isoform_consequences, not by that counter.
_WORKER_IMPORTS = {"mapp.maei": ("numpy", "pandas")}


@lru_cache(maxsize=64)
def _import_probe(executable, executable_identity, modules, environment):
    del executable_identity  # Included in the cache identity, not the command.
    try:
        result = subprocess.run(
            [executable, "-B", "-c", _IMPORT_SCRIPT, json.dumps(modules)],
            env=dict(environment), capture_output=True, text=True,
            timeout=45, check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise ValueError(f"Python dependency probe could not run in {executable}: {exc}") from exc
    records = [line[len(_MARKER):] for line in result.stdout.splitlines() if line.startswith(_MARKER)]
    try:
        if len(records) != 1:
            raise ValueError("missing or repeated probe result")
        report = json.loads(records[0])
        failures = report["failures"]
        if not isinstance(failures, list):
            raise ValueError("invalid import results")
    except (KeyError, TypeError, ValueError) as exc:
        detail = (result.stderr.strip() or result.stdout.strip())[-1500:]
        raise ValueError(f"Python dependency probe returned no valid result in {executable}: {detail}") from exc
    if failures:
        detail = "; ".join(f"{item['module']}: {item['error']}" for item in failures)
        raise ValueError(f"Required Python imports failed in {executable}: {detail}")
    if result.returncode:
        raise ValueError(f"Python dependency probe exited {result.returncode} in {executable}: {result.stderr.strip()[-1500:]}")
    return report["executable"]


def require_python_imports(executable, modules):
    """Test real imports, including binary/transitive failures, without fitting."""
    from .container_runtime import sanitize_container_environment, require_image_path
    sanitize_container_environment()
    require_image_path(executable, "Python dependency worker")
    modules = tuple(dict.fromkeys(modules))
    if not modules:
        return str(executable)
    if any(not isinstance(name, str) or not re.fullmatch(r"[A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*", name)
           for name in modules):
        raise ValueError("Dependency probes accept Python module names only")
    # A POSIX virtualenv's bin/python is often a symlink. Resolving that link
    # invokes the base interpreter and silently discards the selected venv.
    path = Path(os.path.abspath(Path(executable).expanduser()))
    stat = path.stat()
    environment = dict(os.environ)
    environment.update(PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1", MPLBACKEND="Agg")
    return _import_probe(str(path), (stat.st_size, stat.st_mtime_ns), modules,
                         tuple(sorted(environment.items())))


def python_dependency_issues(spec, *, resolve_executable):
    """Check only enabled compute dependencies; historical adoption has none.

    Parameters selecting an external Python interpreter are resolved separately
    from the current Katana process. Finding an executable never proves imports.
    """
    adapter = str(spec.get("adapter") or spec.get("workflow_id", "").removeprefix("flow."))
    if adapter not in _HOST_IMPORTS:
        return []
    probes = [(sys.executable, _HOST_IMPORTS[adapter], "Katana runtime")]
    parameters = spec.get("parameters") or {}
    if adapter in _WORKER_IMPORTS:
        try:
            probes.append((resolve_executable(parameters.get("python"), sys.executable),
                           _WORKER_IMPORTS[adapter], "configured worker runtime"))
        except (OSError, ValueError) as exc:
            return [f"{adapter}: {exc}"]
    issues = []
    for executable, modules, role in probes:
        try:
            require_python_imports(executable, modules)
        except (OSError, ValueError) as exc:
            package_hint = " Bio is supplied by biopython." if "Bio." in str(exc) else ""
            issues.append(
                f"{adapter} ({role}): {exc}. Repair the selected pipeline image/module environment; "
                "for standalone use, repair the active Katana environment and run preflight again. "
                f"{package_hint} No analysis was started."
            )
    return issues
