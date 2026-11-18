"""Directory metadata for the local application's click-to-select browser."""
import os
from pathlib import Path
from urllib.parse import urlsplit

from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel


def workspace_home() -> Path:
    if os.environ.get('KATANA_WORKSPACE_HOME'):
        return Path(os.environ['KATANA_WORKSPACE_HOME']).expanduser().resolve()
    return Path(os.environ.get("LOCALAPPDATA", str(Path.home() / ".local" / "share"))) / "Katana" / "workspaces"


class BrowseLocalFiles(BaseModel):
    path: str = ""


def create_local_files_router() -> APIRouter:
    router = APIRouter()

    @router.post("/api/local-files/browse")
    def browse(payload: BrowseLocalFiles, request: Request):
        # Do not expose the local filesystem browser to remote origins/clients.
        allowed = {"127.0.0.1", "localhost", "::1", "testclient"}
        if not request.client or request.client.host not in allowed or request.url.hostname not in allowed:
            raise HTTPException(403, "File selection is available on this computer only")
        origin = request.headers.get("origin")
        if origin and urlsplit(origin).hostname not in allowed:
            raise HTTPException(403, "File selection requires a local application origin")
        home = Path.home()
        try:
            path = Path(payload.path).expanduser().resolve() if payload.path else home
            if path.is_file():
                path = path.parent
            if not path.is_dir():
                raise ValueError("Folder not found")
            entries = []
            with os.scandir(path) as items:
                for item in items:
                    try:
                        directory = item.is_dir()
                        if directory or item.is_file():
                            entries.append({"name": item.name, "path": str(path / item.name), "directory": directory})
                    except OSError:
                        continue
            entries.sort(key=lambda item: (not item["directory"], item["name"].casefold()))
        except (OSError, ValueError) as exc:
            raise HTTPException(422, "Cannot open this folder: " + str(exc)) from exc
        locations = [("Home", home), ("Desktop", home / "Desktop"), ("Downloads", home / "Downloads"),
                     ("Workspaces", workspace_home())]
        shortcuts = [{"name": name, "path": str(p)} for name, p in locations if p.is_dir()]
        if os.name == "nt":
            # Listing drive letters must not probe disconnected network/removable drives.
            import ctypes
            mask = ctypes.windll.kernel32.GetLogicalDrives()
            shortcuts.extend({"name": f"{chr(65+i)}:", "path": f"{chr(65+i)}:/"} for i in range(26) if mask & (1 << i))
        else:
            shortcuts.append({"name": "Computer", "path": "/"})
        return {"path": str(path), "parent": str(path.parent) if path.parent != path else None,
                "entries": entries, "locations": shortcuts}

    return router
