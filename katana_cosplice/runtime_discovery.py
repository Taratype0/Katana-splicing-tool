"""Portable, standard-library-only discovery of user-installed R executables."""
from __future__ import annotations

import os
from pathlib import Path
import shutil
import sys


def rscript_candidates() -> list[Path]:
    """Return existing candidates in configured/PATH/environment priority order.

    This only locates executables. Callers still verify the R packages required
    by their selected method; discovering R never authorizes a fallback method.
    """
    candidates: list[Path] = []
    configured = os.environ.get("KATANA_RSCRIPT", "").strip()
    if configured:
        candidates.append(Path(shutil.which(configured) or configured).expanduser())
    on_path = shutil.which("Rscript")
    if on_path:
        candidates.append(Path(on_path))
    r_home = os.environ.get("R_HOME", "").strip()
    if r_home:
        candidates.extend((Path(r_home) / "bin" / "Rscript", Path(r_home) / "bin" / "Rscript.exe"))
    prefixes = [Path(value) for value in (os.environ.get("CONDA_PREFIX"), sys.prefix) if value]
    conda_exe = os.environ.get("CONDA_EXE") or shutil.which("conda") or shutil.which("conda.exe")
    if conda_exe:
        prefixes.append(Path(conda_exe).resolve().parent.parent)
    prefixes.extend(Path.home() / name for name in ("anaconda3", "miniconda3", "mambaforge"))
    for prefix in dict.fromkeys(path.resolve() for path in prefixes):
        for suffix in ("bin/Rscript", "Scripts/Rscript.exe", "Lib/R/bin/Rscript.exe"):
            candidates.append(prefix / suffix)
        envs = prefix / "envs"
        if envs.is_dir():
            for pattern in ("*/bin/Rscript", "*/Scripts/Rscript.exe", "*/Lib/R/bin/Rscript.exe"):
                candidates.extend(sorted(envs.glob(pattern)))
    if sys.platform == "win32":
        for variable in ("ProgramFiles", "ProgramFiles(x86)"):
            installed = os.environ.get(variable)
            if installed:
                candidates.extend(sorted((Path(installed) / "R").glob("R-*/bin/Rscript.exe"), reverse=True))
    result: list[Path] = []
    seen: set[str] = set()
    for candidate in candidates:
        resolved = candidate.resolve()
        key = str(resolved).casefold() if os.name == "nt" else str(resolved)
        if key not in seen and resolved.is_file():
            result.append(resolved)
            seen.add(key)
    return result
