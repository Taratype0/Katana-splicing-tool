# Runtime environments

The integrated Nextflow workflow must supply the tools needed by its downstream
tasks through versioned module environments or containers, just as the upstream
modules do. Users of that workflow should not need to activate or repair a shared
host Conda environment. A standalone Conda installation is for users who invoke
Katana's Python API, CLI or browser service directly.

The current development repository does **not** yet distribute a validated
Katana analysis image. `--katana_container` remains an explicit image supplied
by the operator; no published image URL or digest is assumed. Choosing the
Docker, Apptainer or Singularity profile alone does not install Katana or its
dependencies. Container builds, solved platform locks and clean-container
acceptance remain release work.

## Enabled-stage checks

`katana_next/flow/dependency_checks.py` now performs actual imports in a separate
Python process. It catches missing packages, broken binary extensions and failed
transitive imports. An available Python executable is insufficient. Checks use
the executable that will run the stage; an overridden MAEI worker is checked
separately from the Katana host. The probe imports libraries only: it does not
open study inputs, fit a model, install software or submit a job.

| Enabled adapter | Python imports checked | Other existing checks |
| --- | --- | --- |
| MAPP input preparation | NumPy, pandas, `Bio.Seq`, `Bio.SeqIO` in the Katana host | Native counter source identity and declared inputs |
| MAPP MAEI fit | NumPy and pandas in the host and selected worker | Rcpp, optparse, Rscript, seven native source identities and input tables |
| Network discovery/bootstrap | NumPy, pandas, SciPy, scikit-learn, NetworkX, igraph, dynamicTreeCut | The discovery stage checks Rscript/WGCNA and its configuration |
| Network association | NumPy, pandas, SciPy, statsmodels | Project, sample and declared comparison validation |
| Network enrichment | NumPy, pandas, SciPy | Declared matched resources and module inputs |
| Existing-AS input projection | NumPy, pandas, PyArrow | Explicit source and project identities |

Biopython supplies the `Bio` import. In the current MAPP bridge it is required
by Katana's `isoform_consequences` module, imported by the preparation adapter.
The official kmer counter uses the standard library; the official MAEI Python
postprocessors use NumPy/pandas. These runtime roles are checked independently.
Successful checks are cached within the current preflight process only. Missing
dependencies block execution before a stage output directory is created.

This table describes implemented Python checks, not complete installation
acceptance. RSeQC/FIMO workers, optional native tools, package versions, reference
compatibility and every distributed container still require their own checks.
Historical-result adoption remains scoped to transport validation; it does not
require or execute unrelated analysis backends. An authenticated zero-eligible
MAEI result remains `not_applicable`, rather than forcing an unused fit runtime.

## Standalone installation

The development bootstrap is `environments/standalone-python.yml`. Its Python
3.13 version follows the existing Windows/Linux CI target. From a checkout:

```sh
conda env create -f environments/standalone-python.yml
conda activate katana-standalone
python -m pip install -e '.[flow,next]'
python -m pip check
```

Python dependency constraints come from `pyproject.toml`; the YAML does not
duplicate them. This bootstrap has **not** been solved or installed as part of
the import-check repair. It is not a locked scientific environment. R, RSeQC,
the selected FIMO runtime and optional native methods must also be provided for
the workflows that use them. Native versions that require different R or Python
stacks should have separate module environments rather than silently replacing
one another. The current package still includes desktop dependencies; removing
those from the headless/browser installation is separate packaging work.

## Container packaging still to finish

The remaining deliverable is a tested container specification for each required
runtime family, built from the same package constraints and pinned native tools,
with resolved package locks and image digests recorded after actual builds.
Nextflow must select those images automatically for enabled tasks. No candidate
recipe should be presented as a working image before it has been built and
checked. Until that work is complete the pipeline is a development workflow,
not a self-contained release.

FASTA/GTF, motif libraries, CLIP peaks and pathway databases are versioned input
resources. Container/Conda installation supplies their reader and analysis
software; it does not substitute reference data or silently download a different
genome build. Completed study results are reused when their source, scope,
parameters and artifacts match; repairing an environment does not authorize a
new analysis of those studies.
