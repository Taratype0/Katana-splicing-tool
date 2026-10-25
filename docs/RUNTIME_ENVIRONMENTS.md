# Runtime environments

The integrated Nextflow pipeline uses a versioned Katana image for its downstream
tasks. Existing upstream nf-core tool modules keep their own environments and
images. Its users do not need to activate or repair a shared host Conda environment.
Standalone users can install the same
analysis tools with Conda and run the Python API, CLI or browser directly.

Qt is an optional `desktop` extra. The `standalone` extra installs the analysis,
report and browser dependencies without PySide6. Neither route supplies study
references: FASTA/GTF, motif libraries, CLIP peaks and pathway databases remain
explicit, versioned inputs.

The current Linux x86-64 development image, wheel and installed Conda binding
are identified in the [entry/runtime record](validation/entrypoints-runtime-20261025.md).
The [23 October container record](validation/current-container-20261023.md)
and older archives remain historical evidence. Build a matching wheel and image
when changing the software revision. No registry image or PyPI release is claimed.
See [container build/load instructions](../containers/downstream/README.md). A complete runtime is a
prerequisite for scientific workflow acceptance, not a substitute for it.
The [generic workflow record](validation/generic-pipeline-20261023.md) tracks
existing-result execution separately from those installation checks.

## Standalone Linux environment

`environments/standalone-linux.yml` describes the complete native runtime family:

| Component | Used by |
| --- | --- |
| Python 3.12 and the `standalone` extra | Analysis, reports, browser service, target cards and portable result handling |
| R 4.4, WGCNA 1.73, dynamicTreeCut, NetRep | Network discovery, module statistics and external preservation |
| DESeq2, DEXSeq, fgsea, limma | Declared expression/exon interaction, expression-response and ranking analyses |
| Rcpp, optparse, C/C++ compilers and make | Original MAPP MAEI native model; private compilation keeps packaged source unchanged |
| RSeQC, samtools and pysam | BAM/reference validation and gene-body coverage |
| PyMEME Suite 0.1.0a4 | The existing FIMO Python worker, including its fixed terminal-site handling |
| PyHMMER | Optional protein-domain worker; Pfam models are separate inputs |

MAPP's unchanged MAEI model and counter are bundled at the recorded author
commit, with per-file hashes and original license notices. This does not include
or claim execution of all MAPP subworkflows. DJExpress's existing author-source
adapter and R scripts are included in the wheel too.

Use the exact accepted dependency locks to reproduce an installation. Run from
this checkout, with the wheel located in `dist/`, and choose a new environment:

```sh
micromamba create -n katana --file environments/locks/linux-64.explicit.txt
micromamba run -n katana python -m pip install --no-deps --require-hashes --no-build-isolation \
  -r environments/locks/pip-linux-64.lock
micromamba run -n katana python -m pip install --no-deps dist/katana_splicing_tool-*.whl
micromamba run -n katana python -m pip check
micromamba run -n katana katana-flow --help
micromamba run -n katana katana-open --no-open
```

Conda can consume the same explicit native lock using `conda create --file`.
The pip lock supplements that base; it is not a standalone requirements file.
Install native tools first, then the additional pip distributions, then the
Katana wheel. Do not update only one side of this tested combination.
The explicit native lock retains the successful solver's linking order, with
Bioconductor annotation post-link packages last. Alphabetically sorting it can
run their installers before Python, R, yq and curl are available. Package-manager
exit status alone is insufficient: the native acceptance checks load GO.db and
GenomeInfoDbData through their dependent R packages. Source-built Python wheels
use the locked setuptools/wheel tools; build isolation is disabled to avoid an
unrecorded temporary build environment.

The `dist/` glob must select one matching Katana wheel. A source checkout does
not contain that built artifact. To build from the current development branch,
use the [wheel build instructions](INSTALL.md#development-wheels), or use a
supplied wheel with its recorded checksum. Do not combine an older image or
wheel with newer task code without rebuilding and checking that combination.

For a deliberate new solve, use the YAML and then install the wheel with its
`standalone` extra. This may select newer dependencies and requires another
clean acceptance before release:

```sh
micromamba create --strict-channel-priority --override-channels \
  -c conda-forge -c bioconda -n katana-new -f environments/standalone-linux.yml
micromamba run -n katana-new python -m pip install 'dist/katana_splicing_tool-<version>-py3-none-any.whl[standalone]'
```

## What is checked

`python -m scripts.check_downstream_environment --out <new-directory>` performs
real imports, CLI invocations, wheel-resource checks and headless figure exports.
On Linux it also exercises native WGCNA/Rcpp, a tiny DESeq2 fit, DEXSeq
construction/normalization, fgsea, the installed six-sample network bridge,
private compilation of the original MAPP model, the FIMO worker, PyHMMER and
synthetic BAM-to-RSeQC coverage. It writes machine-readable results and logs.
All inputs are generated locally. It never accesses a study or downloads genome
resources. Run outside the checkout with `PYTHONPATH` unset to ensure imports
come from the installed wheel.

This verifies installation and native linkage, not a real-data end-to-end
scientific run. Full pipeline acceptance, statistical validity of a particular
experiment and study completion remain separate checks. A successful import
alone is not reported as a native fit.

The historical `2.0.0b5.dev20261023` Linux x86-64 standalone build was installed in a new
isolated environment and passed all 51 installation checks. The wheel was
rebuilt from its source distribution and checked for browser assets, R scripts,
author source and license notices. This environment uses Python 3.12.15,
R 4.4.3, WGCNA 1.73, DESeq2 1.46.0, DEXSeq 1.52.0 and fgsea 1.32.2. The
recorded native and supplementary pip locks contain 369 and 59 distributions,
respectively. No study datasets were processed during this acceptance.

The matching Linux image also passed all 51 native checks while disconnected
and running as a non-root user. Actual Nextflow image tasks passed for network,
MAPP and MotifTrace dependencies, then reused all three tasks on resume without
changing output files. Windows passed 41 portable checks and actual browser
HTTP checks. See the [validation record](validation/downstream-runtime-20261023.md)
for artifact identities, scope and reproducible commands, and the
[container instructions](../containers/downstream/README.md) for image builds.
The current `2.0.0b5.dev20261025` artifacts retain the verified native layers and
locks, check every installed package/resource against the new wheel, run
`pip check`, and exercise public target restore/plan/cache reuse. The original
51 native checks remain linked to their actual receipts and unchanged workers;
they are not reported as 51 new executions for every application overlay.
Complete-mode and report acceptance remain separate in the
[current record](validation/entrypoints-runtime-20261025.md).

The existing dynamicTreeCut 0.1.1 adapter converts finite, integral ordinal
ranks to integer indices for recent SciPy versions. Ranking values, ordering,
non-ordinal methods and analysis thresholds are unchanged; invalid indices
are rejected rather than truncated. The installed six-sample native fit
exercises this compatibility path.

`--portable-only` checks the Python/browser/CLI side without claiming native R,
RSeQC, FIMO or HMMER readiness. It is the Windows portable acceptance scope.
Native computation can run on Linux directly or in a Linux container; WSL is an
optional local Linux host. The Linux lock must not be used as a Windows or macOS
lock. `standalone-python.yml` remains a Python-only bootstrap and does not
install the native analysis toolchain.

## Enabled-stage preflight

Preflight checks actual Python imports in the runtime that will execute a stage,
including the `Bio` imports used by MAPP preparation and NumPy/pandas in an
explicitly selected MAEI worker. It also checks selected native tools/R packages,
packaged author-source identity and declared resources. A missing module or
broken binary extension blocks the stage before creating analysis outputs.
Historical-result adoption only validates the transport it uses; it does not
run unrelated methods. Authenticated zero-eligible MAPP inputs remain
`not_applicable` rather than forcing a fit.

Completed results are reused only when their input identity, comparison,
parameters and output evidence match. Installing or repairing an environment
does not authorize rerunning a completed study.
