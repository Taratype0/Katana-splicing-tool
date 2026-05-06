# Windows Release Checklist

## Goal

Ship a Windows-first desktop package that users can unzip and run with minimal setup.

## Before packaging

- Confirm `software/Jutils` exists
- Confirm `software/SUPPA` exists
- Confirm `software/rmats2sashimiplot/src` exists
- Confirm the app starts with `python -m app.main`
- Confirm `PySide6` imports successfully

## Required user-configurable items

These are still expected to be configured by the user inside the app:

- `Rscript.exe`
- annotation `GTF`
- optional transcript `FASTA`

## Recommended optional dependencies

- `pysam` for `rmats2sashimiplot`
- R packages for `IsoformSwitchAnalyzeR`

## Build

```powershell
.\scripts\build_windows.ps1
```

## Smoke test after build

1. Open `dist\KatanaSplicingTool\KatanaSplicingTool.exe`
2. Load a small project
3. Confirm Settings shows bundled `Jutils`, `rmats2sashimiplot`, `SUPPA`
4. Confirm output directory can be changed
5. Run main analysis
6. Generate shortlist
7. Generate cards
8. Open Sashimi page
9. Open Isoform page and verify `Rscript` prompt path works

## Current known limitations

- `IsoformSwitchAnalyzeR` still depends on external R/Bioconductor
- `rmats2sashimiplot` still depends on runtime compatibility and `pysam`
- Linux packaging is not the current target
