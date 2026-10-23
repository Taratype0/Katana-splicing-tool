"""Keep pipeline containers independent of mounted host software environments."""
from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import sys

_IMAGE_PREFIX = Path("/opt/conda")


def in_container():
    return os.environ.get("KATANA_CONTAINER") == "1"


def sanitize_container_environment():
    """Clear inherited interpreter/library overrides only in pipeline images."""
    if not in_container():
        return []
    removed = []
    exact = {"PYTHONPATH", "PYTHONHOME", "PYTHONUSERBASE", "R_HOME", "R_ENVIRON",
             "R_ENVIRON_USER", "R_PROFILE", "R_PROFILE_USER", "VIRTUAL_ENV",
             "LD_LIBRARY_PATH", "LD_PRELOAD", "LIBRARY_PATH", "CPATH", "C_INCLUDE_PATH",
             "CPLUS_INCLUDE_PATH", "PKG_CONFIG_PATH", "CMAKE_PREFIX_PATH"}
    for name in list(os.environ):
        if name in exact or name.startswith("R_LIBS"):
            removed.append(name)
            os.environ.pop(name)
    os.environ["PYTHONNOUSERSITE"] = "1"
    os.environ["PYTHONSAFEPATH"] = "1"
    # The image prefix owns scientific tools; retain only image system tools
    # for shell/R compilation support, never mounted host Conda PATH entries.
    if os.name != "nt":
        os.environ["PATH"] = os.pathsep.join([str(_IMAGE_PREFIX / "bin"),
            "/usr/local/sbin", "/usr/local/bin", "/usr/sbin", "/usr/bin", "/sbin", "/bin"])
    return sorted(removed)


def require_image_path(value, role):
    """Reject mounted executables/libraries, including escaping symlinks."""
    if not in_container():
        return str(value)
    prefix = _IMAGE_PREFIX.absolute()
    if Path(sys.prefix).resolve() != prefix.resolve():
        raise ValueError(f"Container Python prefix {sys.prefix} does not match the fixed image runtime {prefix}; mounted host launchers are not allowed")
    candidate = Path(value).expanduser().absolute()
    if not candidate.is_relative_to(prefix) or not candidate.resolve().is_relative_to(prefix.resolve()):
        raise ValueError(f"Container {role} must belong to the image runtime {prefix}; external host path is not allowed: {value}")
    return str(candidate)


def require_image_executable(value, role="executable"):
    if not in_container():
        return str(value)
    resolved = shutil.which(str(value))
    if not resolved:
        raise ValueError(f"Container {role} is unavailable in its image: {value}")
    return require_image_path(resolved, role)


def validate_runtime_settings(settings):
    """Validate selected native interpreters and R libraries before invocation."""
    if not in_container():
        return
    if isinstance(settings, dict):
        if settings.get("wsl_distribution"):
            raise ValueError("Container tasks require native image tools, not a WSL runtime")
        for name, value in settings.items():
            if name in {"python", "rscript", "rscript_path", "wgcna_rscript", "rseqc_executable"} and value:
                require_image_executable(value, name)
            elif name in {"r_library", "r_libraries"} and value:
                paths = value if isinstance(value, list) else [value]
                for path in paths:
                    require_image_path(path, name)
            elif isinstance(value, (dict, list)):
                validate_runtime_settings(value)
    elif isinstance(settings, list):
        for item in settings:
            validate_runtime_settings(item)


def project_asset_path(asset, project, project_path):
    """Resolve assets against input_root, as the project manifest reader does."""
    path = Path(asset["path"])
    if path.is_absolute():
        return path
    root = Path(project.get("input_root", "."))
    if not root.is_absolute():
        root = Path(project_path).resolve().parent / root
    return root / path


def validate_task_runtime_boundaries(task, *, project=None, project_path=None):
    """Shared probe/direct-task gate; it does not create any output directory."""
    if not in_container():
        return
    sanitize_container_environment()
    require_image_path(sys.executable, "Python")
    parameters = task.get("parameters") or {}
    validate_runtime_settings(parameters)
    project = project or {}
    validate_runtime_settings(project.get("metadata", {}))
    adapter = task.get("adapter", "")
    config_path = task.get("inputs", {}).get("config") if adapter.startswith("network.") else None
    if config_path and Path(config_path).is_file():
        from katana_cosplice.config import load_config
        validate_runtime_settings({"wgcna_rscript": load_config(config_path).network.wgcna_rscript})
    config_path = parameters.get("config_path")
    if config_path and Path(config_path).is_file():
        validate_runtime_settings(json.loads(Path(config_path).read_text(encoding="utf-8-sig")))
    if task.get("workflow_id") in {"prioritize.event_motifs", "prioritize.rbp_mechanism", "prioritize.se_rbp_impact"}:
        for asset in project.get("assets", []):
            if asset.get("kind") != "fimo_runtime":
                continue
            path = project_asset_path(asset, project, project_path)
            validate_runtime_settings(json.loads(path.read_text(encoding="utf-8-sig")))
