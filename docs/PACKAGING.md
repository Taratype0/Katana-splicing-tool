# Packaging

## Supported delivery modes

### 1. Portable EXE bundle

Recommended for end users on Windows.

Contents:
- desktop app
- bundled Python runtime
- bundled `software/Jutils`
- bundled `software/rmats2sashimiplot`
- bundled `software/SUPPA`

Still required on user machine:
- `Rscript` for IsoformSwitchAnalyzeR unless you later bundle an R runtime
- `samtools` if not bundled separately
- any missing Python-side runtime dependency not already frozen into the bundle

### 2. Git checkout + one-step setup

Recommended for lab/internal users.

```powershell
conda env create -f configs/environment.yml
conda activate katana-splicing-tool
pip install -e .
python -m app.main
```

## Build EXE

Install PyInstaller in the build environment:

```powershell
pip install pyinstaller
pyinstaller build.spec
```

The generated folder will contain:
- `KatanaSplicingTool.exe`
- `configs/`
- `software/`

The current spec intentionally bundles only:
- `software/Jutils`
- `software/SUPPA`
- `software/rmats2sashimiplot/src`

This avoids shipping the old Python 2 `build/lib` tree from `rmats2sashimiplot`.

## Current packaging status

Implemented:
- application entrypoint
- portable folder layout
- bundled tool directories in spec
- runtime dependency reporting in app

Not yet fully bundled:
- R runtime / Bioconductor stack for IsoformSwitchAnalyzeR
- automatic `samtools` bundling
- final installer wrapper

## Recommendation

Short term:
- ship a portable Windows folder produced by PyInstaller
- include `software/` inside the bundle
- let users configure `Rscript` once in Settings
- use `launch_katana.bat` for source-mode startup

Long term:
- add a true installer
- optionally bundle R runtime or provide a guided installer
