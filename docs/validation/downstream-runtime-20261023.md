# Downstream runtime validation — 23 October 2026

Katana `2.0.0b5.dev20261023` passed the Linux standalone runtime checks and the
Windows Python/browser checks below. These are installation and native-tool
checks using small generated inputs, not a real-study end-to-end validation.

## Distribution identity

The wheel was rebuilt from its source distribution, then checked for browser
assets, Nextflow resources, R scripts, bundled MAPP/DJExpress author source and
license notices. PySide6 is an optional desktop dependency; these installations
did not include it.

| Artifact | SHA256 |
| --- | --- |
| Final wheel, 3,836,184 bytes | `bc6ed2b54d189a8e3b4684309d3147cf53ef4e25d53d0c2153af3d2cd4c5fa0c` |
| Source distribution | `a769b0e10a2771d5c7acda90404a5f5f179a95978aeb5003cd77975bf8594cd5` |
| Linux native dependency lock | `cc6b9d9c17b98fdcbf59b0dfffe77c20356b786830c7dc503ac2dc5f4b434b99` |
| Supplementary Linux pip lock | `4bbfdd6b3510821fce032d52d42a52bef1313238fad82b7c2d831adda5cb87f1` |

## Completed checks

| Environment | Observed result | Scope |
| --- | --- | --- |
| Linux x86-64, Python 3.12.15, R 4.4.3 | **51 passed; 0 failed** | Imports, seven CLI entry points, package resources, `pip check`, headless graphics and native checks listed below |
| Windows x86-64, Python 3.13.9 | **41 passed; 0 failed** | Portable Python, CLI, browser dependencies, resources, `pip check` and headless graphics |
| Windows installed browser service | **6 HTTP checks passed** | Health/version, workspace inventory, empty job list, homepage, JavaScript and CSS; all returned HTTP 200 |
| Windows installed-file comparison | **605 matched; 0 mismatches** | Installed package files compared byte-for-byte with the final wheel |
| Linux x86-64 container, disconnected, non-root | **51 passed; 0 failed** | The same installed native runtime checks, without a host environment mount |
| Nextflow 25.10.4 with the actual Docker image | **3 completed; 3 cached on resume** | Independent network, MAPP and MotifTrace dependency probes; output hashes unchanged |

Both environments were created separately from existing installations. Initial
candidate results were retained; the corrected final application wheel was
reinstalled with `--no-deps` into those isolated environments before the final
checks. Tests ran outside the checkout with `PYTHONPATH` cleared. The browser
test started and stopped only its own temporary service.

Linux native checks exercised:

- WGCNA adjacency/TOM and the installed network interface on 60 generated events
  across six samples, producing 60 assignments and 59 linkage rows.
- A small DESeq2 fit, DEXSeq construction/size-factor estimation and fgsea
  permutations; GO.db and GenomeInfoDbData were loaded explicitly.
- Private Rcpp compilation of the unchanged bundled MAPP model.
- The actual FIMO worker, a synthetic PyHMMER search, and generated BAM reads
  through samtools and RSeQC gene-body coverage.

The native lock contains 369 distributions and the supplementary pip lock 59.
Its versions include WGCNA 1.73, DESeq2 1.46.0, DEXSeq 1.52.0, fgsea 1.32.2
and PyMEME Suite 0.1.0a4. Explicit-package linking order and locked Python build
tools were checked. An ordinal-rank integer-index compatibility correction
preserves the existing dynamicTreeCut ranking values and analysis parameters.

## Container and Nextflow checks

The accepted Linux container has image ID
`sha256:6df830ac7cc41cf4726abf352646965ac0da7d13120a103b500207c77821f4dd`.
It contains the application wheel identified above and the locked native/Python
dependencies. Three separately hash-locked Debian packages provide `ps` for
Nextflow task metrics. An earlier candidate passed native checks but failed
actual Nextflow execution because `ps` was missing; that failure was retained
and the corrected image was tested again. No host system upgrade was needed.

The final image passed all 51 native checks with networking disabled and UID
1000. The reusable Nextflow check then executed three independent installed
CLI/native-runtime probes. All three completed with exit code zero, and all
three were cached on `-resume`, preserving every output hash. An existing,
mounted host interpreter path was rejected in each task. These are actual
container executions, not a mocked launcher or a claim of scientific completion.
The image ID is local and has not been published to a registry.

Reproduce the Nextflow check on Linux with an already built or loaded image:

```sh
python pipeline/tools/check_container_runtime.py \
  --image katana-downstream:2.0.0b5-dev20261023 \
  --outdir /tmp/katana-image-check
```

## Reproduction and limits

Use the exact-lock installation commands in
[Runtime environments](../RUNTIME_ENVIRONMENTS.md). Run the installed checker
outside the checkout and choose a new output directory:

```sh
python -m scripts.check_downstream_environment --out runtime-check
# Windows portable scope:
python -m scripts.check_downstream_environment --portable-only --out portable-check
```

Machine-readable receipts retain command outcomes, tool versions and output
hashes. No study inputs or existing study results were accessed. These checks
do not certify a complete scientific workflow, biological conclusions, all
optional resource combinations, native Windows R/FIMO execution, macOS, or
remote scheduler integration.
