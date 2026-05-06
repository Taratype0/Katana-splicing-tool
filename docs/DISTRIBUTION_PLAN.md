# Distribution Plan

## Goal

Allow users to either:
- download one Windows bundle and double-click an `.exe`
- or clone the repo and run after one environment setup

## Packaging layers

### Layer 1: Python app
- bundle with PyInstaller

### Layer 2: bundled tools
- `Jutils`
- `rmats2sashimiplot`
- `SUPPA`

### Layer 3: external runtimes
- `Rscript` for IsoformSwitchAnalyzeR
- `samtools` for sashimi support

## Current blocker summary

- `PyInstaller` is not yet installed in the current environment
- `Rscript` is not available on PATH
- `pysam` is required for `rmats2sashimiplot`

## Shipping recommendation

First public release:
- provide a portable bundle
- provide a release note listing required external runtimes
- keep repo setup instructions in parallel

