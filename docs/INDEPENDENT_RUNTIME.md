# Install the independent Katana runtime

Use the launcher and wheel from the [29 October permission-fix release](https://github.com/Taratype0/Katana-splicing-tool/releases/tag/v2.0.0-beta.5.dev20261029).
This beta includes paired comparisons, automatic strandedness, direct Nextflow
parameters and managed cleanup. The older downloads remain available for
reproducing earlier installations.

## Linux cluster

The controller archive includes Python, installed dependencies, Nextflow
25.10.4 and the exported pipeline. It does not use another Katana installation,
a Conda environment or your Python user-site packages. The host supplies
**Java 17+**, **Singularity/Apptainer** and scheduler clients. Load your site's
modules before running. Scientific tasks use containers; references and
databases are separate inputs.

Download the controller from the new release. Reuse the unchanged scientific
images from the [28 October release](https://github.com/Taratype0/Katana-splicing-tool/releases/tag/v2.0.0-beta.5.dev20261028); their exact hashes are also listed in the new manifest. This patch changes installation and launch validation, not scientific methods or image dependencies.

Required assets:

| Asset | Purpose |
| --- | --- |
| `katana-controller-2.0.0b5.dev20261029-linux-x86_64.tar.gz` | Portable launcher, Python and Nextflow |
| `katana-downstream-2.0.0b5.dev20261028-linux-amd64.sif` | Matching Katana scientific runtime |
| `kanata-dexseq-tools-20261024-linux-amd64.sif` | Unchanged DEXSeq counting runtime |
| `SHA256SUMS`, `RELEASE_MANIFEST.json` | Checksums and tested versions |

```bash
sha256sum --ignore-missing --check SHA256SUMS
mkdir -p tools
tar -xzf katana-controller-2.0.0b5.dev20261029-linux-x86_64.tar.gz -C tools
export PATH="$PWD/tools/katana-controller/bin:$PATH"
katana-flow run --help
nextflow -version
```

Every selected download must be listed in the checksums and report `OK`.
The archive can be extracted to a different directory on another compatible
Linux x86-64 host. Its launchers locate files relative to their installation.
The bundled Python serves the controller; scientific tools run in task images.

Put caches on the allocation intended for the study:

```bash
export PROJECT=/path/to/your/project
export NXF_HOME="$PROJECT/cache/katana/nextflow"
export NXF_SINGULARITY_CACHEDIR="$PROJECT/cache/katana/images"
export SINGULARITY_CACHEDIR="$PROJECT/cache/katana/singularity"
mkdir -p "$NXF_HOME" "$NXF_SINGULARITY_CACHEDIR" "$SINGULARITY_CACHEDIR"
```

For Apptainer use `NXF_APPTAINER_CACHEDIR` and `APPTAINER_CACHEDIR` instead.
Pass a writable `NUMBA_CACHE_DIR` through the container environment when the
site disallows writes to a default cache location.

## FASTQ or BAM to the report

The scientific arguments are unchanged. A samplesheet selects FASTQ or
genome-aligned BAM; the contrast sheet selects comparisons and their direction.
Katana generates its internal project/settings JSON.

```bash
katana-flow run \
  -profile slurm,singularity -c /path/to/site.config \
  --slurm_account YOUR_ALLOCATION \
  --input /data/study/samples.csv \
  --contrasts /data/study/contrasts.csv \
  --fasta /data/reference/genome.fa \
  --gtf /data/reference/annotation.gtf \
  --reference_build YOUR_BUILD --annotation_release YOUR_RELEASE \
  --species 'Mus musculus' \
  --gene_sets /data/resources/pathways.gmt \
  --rbp_gene_sets /data/resources/rbp_gene_sets.gmt \
  --motif_database /data/resources/motifs.json \
  --clip_peaks /data/resources/matched_clip_peaks.tsv \
  --katana_container /data/images/katana-downstream-2.0.0b5.dev20261028-linux-amd64.sif \
  --dexseq_tools_container /data/images/kanata-dexseq-tools-20261024-linux-amd64.sif \
  --run-dir "$PROJECT/cache/katana/study-run" \
  --outdir "$PROJECT/results/study"
```

Replace the paths and site settings. The launcher retains published results and
target-analysis inputs, then cleans validated completed stage work. Failed work
stays available; repeat with `-resume`. `--keep-work` retains completed work too.
Use separate run/output directories for different studies. Queue selection is
site configuration; Katana does not automatically select a memory partition.

The direct Nextflow entry is `nextflow run tools/katana-controller/pipeline/main.nf`
with the same scientific parameters. Use `-work-dir` for that route, which keeps
native Nextflow work. Other upstream tool images come from the workflow or the
site image cache, not the controller's Python environment.

For paired animals, share `subject` between matching samples and set
`paired=true` for that contrast. See [samples and contrasts](NEXTFLOW_USAGE.md),
the [parameter reference](NEXTFLOW_PARAMETERS.md) and
[editable template](../examples/nextflow/params.complete.yaml). Optional analyses
have switches; supply resources for the analyses you enable.

## Docker workstation

Load the unchanged 28 October downstream archive and DEXSeq tools archive:

```bash
docker load -i katana-downstream-2.0.0b5.dev20261028-linux-amd64.docker.tar.gz
docker load -i kanata-dexseq-tools-20261024-linux-amd64.docker.tar.gz
```

Use `-profile local,docker`, omit the cluster config/account and select
`--katana_container katana-downstream:2.0.0b5.dev20261028` plus
`--dexseq_tools_container kanata-dexseq-tools:dev-20261024`.
The controller runs on Linux. These are downloadable local images; the release
does not publish registry tags or a PyPI package.

## Python, browser and existing results

Install the release wheel in a fresh Python 3.11+ environment:

```bash
python -m venv .venv
# Linux: source .venv/bin/activate
# Windows PowerShell: .\.venv\Scripts\Activate.ps1
python -m pip install './katana_splicing_tool-2.0.0b5.dev20261029-py3-none-any.whl[standalone]'
python -m pip check
katana-open --help
katana-report --help
katana-targets --help
```

Completed results and reports can be read without a cluster. New computations
also need their documented R/external tools and resources; Nextflow supplies
those through containers. Standalone native Linux users can use the locked
[Conda environment](RUNTIME_ENVIRONMENTS.md).

Existing bundles enter downstream analysis without repeating read processing;
completed collections can generate reports directly. Gene/event deep dives use
`katana-targets`. See the [entry-point guide](NEXTFLOW_USAGE.md#1-choose-an-entry).

## Validation and updates

The release manifest records hashes, CI and native runtime checks. Wheels built
on both Windows and Linux are installed on Linux and tested by directly running
the packaged gene-import and transcript-import R modules through Nextflow,
including a cached second run. This catches executable-bit loss that an exported
pipeline or `--help` check would miss. The launcher checks packaged script
permissions before submitting work, without modifying the installation. Checks use
installed commands, generated paired counts and generated Nextflow tasks;
they do not imply every study or cluster has been run. The approved seven-section
report format is unchanged.

Install updates into new software directories and use matching images. Keep
data/results outside software directories. Remove an old runtime after launchers
and configs reference the new install. Updating Git alone does not update an
installed package or image.
