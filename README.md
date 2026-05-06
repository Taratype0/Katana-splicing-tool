# Katana Splicing Tool

Local desktop application for RNA splicing, transcriptome integration, and visualization.

## Current scope

- Project scanning and input validation
- Comparison selection with direction reversal
- `JC` / `JCEC` selection
- Modular analysis pipeline skeleton
- Desktop GUI shell using `PySide6`
- Hooks for `rmats2sashimiplot`, `Jutils`, and `IsoformSwitchAnalyzeR`

## Quick start

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -e .
python -m app.main
```

Or on Windows:

```powershell
.\scripts\setup_windows_dev.ps1
.\launch_katana.bat
```

For Windows source-mode development, the bundled `software/` directory is used for
`Jutils`, `SUPPA`, and `rmats2sashimiplot`, so they do not need to be installed
through `conda`/`pip` during initial setup.

On Windows, `pysam` is not included in the default environment bootstrap yet.
This means the core app can be tested first, while `rmats2sashimiplot` may need
an extra dependency step later.

## Distribution targets

- Portable Windows bundle via `PyInstaller`
- Git checkout + local environment setup

Current release focus:
- Windows first
- Linux compatibility later

See `docs/PACKAGING.md` and `docs/WINDOWS_RELEASE_CHECKLIST.md` for packaging details.

## Planned integrations

- rMATS landscape and program comparison
- DEG / SUPPA / DEXSeq / DTU support matrix
- transcriptome-splicing integration
- cards, sashimi, and isoform switch views
