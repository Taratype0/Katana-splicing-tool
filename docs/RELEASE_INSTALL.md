# Install the downloadable Katana release

[Katana 2.0.0 beta 5, 25 October 2026](https://github.com/Taratype0/Katana-splicing-tool/releases/tag/v2.0.0-beta.5.dev20261025)
provides the accepted Python wheel, Linux Docker images, exported Nextflow
workflow and native dependency locks. This is a beta release. The wheel and
images retain their tested identities; publishing them does not rebuild the
software or recompute studies.

## Current main versus the downloadable release

Current `main` includes automatic strandedness, the simplified parameter entry,
stage cleanup and matched-subject comparisons. The frozen 25 October downloads
do not include all of these updates. For these features, use a current installed
package, its workflow and a matching downstream image together. Updating a Git
checkout does not update an installed wheel, exported workflow or container.

The paired implementation and matching HPC image have passed the checks in
[paired-entry validation](validation/paired-entry-20261028.md). That image has
not yet been published as a new downloadable release; the release linked above
still provides the older artifacts.

Install the current controller in a fresh environment on the Linux execution
host (Python 3.11 or newer, plus native Nextflow/Java and the container engine):

```sh
git clone --branch main --single-branch https://github.com/Taratype0/Katana-splicing-tool.git tools/katana-source
python3 -m venv tools/katana-controller
. tools/katana-controller/bin/activate
export PYTHONNOUSERSITE=1
python -m pip install 'tools/katana-source[flow]'
python -m pip check
katana-flow export-pipeline --outdir tools/katana-pipeline
katana-flow run --help
```

The controller launches Nextflow and manages completed output handoffs; task
containers supply the scientific Python/R tools. This installation does not
depend on an older Katana study directory. Avoid `--system-site-packages` when
the goal is an independent installation: it makes the new environment depend
on the base environment. Keep reference and database paths separate from both.

Build a matching image using the [container build instructions](../containers/downstream/README.md)
or, when the native and pip dependency locks are unchanged, the
[verified-release SIF update](validation/paired-entry-20261028.md#updating-an-hpc-application-image).
Pass the new image using `--katana_container` and the matching site configuration
with `-c`. The standard `--input`, `--contrasts`, reference and `--outdir`
arguments are unchanged. Matched animals share `subject` in the samplesheet;
their contrast selects `paired=true` as described in the
[contrast guide](NEXTFLOW_USAGE.md#4-contrastsheet-input-and-direction).

The installation sections below describe the frozen downloadable release.

Use the [Nextflow guide](NEXTFLOW_USAGE.md), [complete parameter reference](NEXTFLOW_PARAMETERS.md)
and [Quick start](QUICKSTART.md) for sample/contrast formats,
entry points and analysis choices. Copyable sheets are in
[examples/nextflow](../examples/nextflow/README.md).

## Choose downloads

| Use | Release assets |
| --- | --- |
| Windows/Linux browser, completed reports and target cards | Python wheel, `RELEASE_MANIFEST.json`, `SHA256SUMS` |
| Linux bulk FASTQ/BAM complete pipeline | Both Docker archives, workflow ZIP, manifest and checksums |
| Linux existing-results downstream pipeline | Katana downstream Docker archive, workflow ZIP, manifest and checksums |
| Standalone native Linux CLI/Python/browser | Python wheel, runtime-locks ZIP, manifest and checksums |

References, motif/pathway databases and study data are separate inputs. The
Python-only installation does not supply R and all external analysis tools.
Docker archives are loaded locally; this release does not publish a registry
image or a PyPI package. The release page also provides source ZIP/tar archives.

Download the files for your route into a new directory. On Linux, check the
downloaded assets against `SHA256SUMS`:

```sh
sha256sum --ignore-missing --check SHA256SUMS
```

Every selected file must be present and report `OK`. On Windows, compare each
download's `Get-FileHash -Algorithm SHA256` value with `SHA256SUMS` or the
manifest. The manifest records the wheel version, image IDs and original
software commits separately from the release's documentation revision.

## Browser, reports and selected targets

Use Python 3.11 or newer in a fresh environment. On Windows PowerShell:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install ".\katana_splicing_tool-2.0.0b5.dev20261025-py3-none-any.whl[next,flow,contracts,targets]"
python -m pip check
katana-open --runtime-root .\katana-runtime
```

On Linux, activate with `. .venv/bin/activate` and use the same wheel filename
using a normal `./` path prefix. `katana-open --no-open` starts the
service on a server without launching a local browser. The browser talks to
that Python service; an offline report can display saved results independently.

Generate a report from explicitly registered completed sources:

```sh
katana-report --project results.json --outdir exports/study-report --compact --no-archive
```

For targets exported from the browser, bind resources on the execution machine:

```sh
katana-targets inspect --request targets.zip
katana-targets plan --request targets.zip --bindings bindings.json
katana-targets run --request targets.zip --bindings bindings.json --outdir target-results
katana-targets results --outdir target-results
```

Saved cards can use completed evidence. New Sashimi requires accessible matching
BAMs/indexes and references; the request does not embed them. See
[target inputs, return import and restore](TARGET_COMMANDS.md).

## Linux Nextflow with Docker

The execution host needs compatible Java, Nextflow and Docker. Native acceptance
used Nextflow 25.10.4 on Linux x86-64. Unzip the workflow archive into a new
directory and load the required images:

```sh
unzip katana-pipeline-2.0.0b5.dev20261025.zip -d tools/katana-pipeline
docker load -i katana-downstream-2.0.0b5.dev20261025-linux-amd64.docker.tar.gz
docker load -i kanata-dexseq-tools-20261024-linux-amd64.docker.tar.gz

export KATANA_IMAGE=katana-downstream:acceptance-20261025-098ceae-presentation09
export DEXSEQ_IMAGE=kanata-dexseq-tools:dev-20261024
docker image inspect "$KATANA_IMAGE" --format '{{.Id}} {{.Os}}/{{.Architecture}}'
docker image inspect "$DEXSEQ_IMAGE" --format '{{.Id}} {{.Os}}/{{.Architecture}}'
```

Compare those IDs with `RELEASE_MANIFEST.json`. The separate DEXSeq image is
needed for new raw-input processing; completed-result downstream-only execution
does not require it. Upstream nf-core modules retain their declared public tool
images, which the execution host must also be able to obtain.

FASTQ and eligible genome-BAM samples use the same complete entry:

```sh
nextflow run tools/katana-pipeline/main.nf -profile docker \
  --mode complete --input inputs/samples.csv \
  --contrasts inputs/contrasts.csv \
  --fasta references/genome.fa --gtf references/annotation.gtf \
  --project_template inputs/project-template.json \
  --flow_config inputs/flow.json \
  --katana_container "$KATANA_IMAGE" \
  --dexseq_tools_container "$DEXSEQ_IMAGE" \
  --outdir results -work-dir work
```

Supply actual sample/contrast sheets, matching references and configured
resources. For existing upstream results, map the files and experimental design
into a valid Katana project, then use `--mode downstream-only`. See
[entry schemas and commands](ENTRY_POINTS.md) and
[complete execution, missing-input gates and resume](INTEGRATED_PIPELINE.md).
An arbitrary nf-core output directory is not automatically a complete project.

Apptainer/Singularity and scheduler profiles require site qualification. Docker
acceptance does not certify every Slurm site, ARM64/macOS host or raw scRNA/snRNA
workflow.

## Standalone Linux native environment

Extract the runtime-locks ZIP; it contains the versioned `environments/` files:

```sh
unzip katana-runtime-locks-2.0.0b5.dev20261025.zip
micromamba create -n katana --file environments/locks/linux-64.explicit.txt
micromamba run -n katana python -m pip install --no-deps --require-hashes --no-build-isolation \
  -r environments/locks/pip-linux-64.lock
micromamba run -n katana python -m pip install --no-deps katana_splicing_tool-2.0.0b5.dev20261025-py3-none-any.whl
micromamba run -n katana python -m pip check
micromamba run -n katana katana-open --no-open
```

Conda can consume the same native lock with `conda create --file`. Preserve lock
ordering and install the supplementary pip lock before the Katana wheel. See
[runtime details](RUNTIME_ENVIRONMENTS.md). Nothing in this installation starts
a study analysis automatically.
