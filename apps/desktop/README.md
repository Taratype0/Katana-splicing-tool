# Katana Next desktop host

Katana Next uses the same tested React workspace in browsers and in the native
Windows host. The desktop entry point starts an isolated FastAPI server on
loopback, waits for its health endpoint, and renders it with PySide6
QtWebEngine. Closing the window shuts the local server down.

```powershell
pnpm --dir apps/web build
python -m katana_next.desktop --check
python -m katana_next.desktop
```

Heavy analysis outputs remain in the project-selected `output_root`; the
desktop host does not write result caches into the installation directory.
