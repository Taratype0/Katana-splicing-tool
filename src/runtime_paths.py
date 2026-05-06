from __future__ import annotations

import sys
from pathlib import Path


def app_root() -> Path:
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
    return Path(__file__).resolve().parents[1]


def resource_path(*parts: str) -> Path:
    return app_root().joinpath(*parts)


def bundled_software_root() -> Path:
    return resource_path("software")


def bundled_configs_root() -> Path:
    return resource_path("configs")


def bundled_assets_root() -> Path:
    return resource_path("assets")
