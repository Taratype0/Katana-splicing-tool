"""Locate the built browser interface in a source checkout or installed wheel."""
from pathlib import Path


def frontend_directory() -> Path | None:
    package = Path(__file__).resolve().parent
    for candidate in (package.parent / "apps" / "web" / "dist", package / "_web"):
        if (candidate / "index.html").is_file():
            return candidate
    return None
