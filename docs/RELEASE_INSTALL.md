# Install the downloadable Katana release

[Katana 2.0.0 beta 5, 25 October 2026](https://github.com/Taratype0/Katana-splicing-tool/releases/tag/v2.0.0-beta.5.dev20261025)
provides the accepted Python wheel, Linux Docker images, exported Nextflow
workflow and native dependency locks. This is a beta release. The wheel and
images retain their tested identities; publishing them does not rebuild the
software or recompute studies.

The current `main` source also supports `strandedness=auto`. That addition is
not in these frozen beta 5 wheel/workflow downloads. For auto inference, keep
the matching downloaded downstream images but run the current source workflow:

```sh
git clone --branch main --single-branch https://github.com/Taratype0/Katana-splicing-tool.git tools/katana-source
nextflow run tools/katana-source/pipeline/main.nf --help
```

Use that `pipeline/main.nf` path in the complete FASTQ/BAM command below.
RSeQC/gffread use their own declared upstream environments. Exporting a workflow
from an older installed wheel does not update it. See
[auto-strandedness requirements and validation](validation/auto-strandedness-20261025.md).

Use the [Nextflow guide](NEXTFLOW_USAGE.md), [complete parameter reference](NEXTFLOW_PARAMETERS.md)
and [Chinese quick start](NEXTFLOW_QUICKSTART_ZH.md) for sample/contrast formats,
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
