# Katana downstream image

This image supplies the Python, R and native tools used by Katana downstream
tasks. Existing upstream modules retain their own nf-core tool environments.
The image contains no research data, genome references, motif libraries or
study-specific configuration. Qt is not installed.

Use `codex/package-development` for the commands below. The accepted runtime is
Linux x86-64; native ARM64 and macOS execution have not been certified. Docker
or a compatible Apptainer/Singularity setup must already be available on the
execution host. Nextflow and Java run on that host; Katana's native dependencies
run inside the image.

There is no published registry image yet. The current local development image
is `katana-downstream:dev-d05c76c`; its matching wheel, image ID, archive checksum
and observed checks are in the [current acceptance record](../../docs/validation/pipeline-closure-20261024.md).
The [earlier runtime record](../../docs/validation/downstream-runtime-20261023.md)
refers to an older image. For later source revisions, build a matching wheel and
image and retain their receipt. [Generic workflow evidence](../../docs/validation/generic-pipeline-20261023.md)
is recorded separately from image installation checks.

## Build locally

Start with a browser-enabled wheel and the checked-in Linux dependency locks.
A source checkout includes the built browser assets; in a separate build
environment, install `build` and build the wheel from the repository root:

```sh
python -m pip install build
KATANA_REQUIRE_WEB=1 python -m build --wheel
```

Set `wheel` to the single file just built. From a clean checkout, use a commit
tag to keep this build separate from earlier runtime-check images:

```sh
wheel=dist/katana_splicing_tool-2.0.0b5.dev20261023-py3-none-any.whl
image="katana-downstream:dev-$(git rev-parse --short HEAD)"
python scripts/build_downstream_container.py \
  --wheel "$wheel" \
  --context /tmp/katana-container-build \
  --tag "$image"
```

Use a new context directory outside the checkout. The helper pins the base by
digest and records the wheel, dependency locks and resulting image ID. It does
not publish to a registry. A local Conda archive cache may be supplied through
`--conda-cache`; only archives named in the lock are copied, then verified by
the installer. The final image excludes the installation cache.

Nextflow's task wrapper also needs the system `ps` utility. The separate
`system-packages-linux-amd64.json` lock records official Debian package URLs,
versions, sizes and SHA256 values for `procps` and its two missing libraries.
The helper verifies those downloads against the lock; the final image installs
only these `.deb` files and runs `ps` as a build check. It does not run an apt
upgrade or resolve additional packages. This system lock is tied to the pinned
Debian 13 base image and does not change the Python wheel or Conda dependency
locks. A different base image requires its own validated system package lock.

## Load a supplied image archive

If you have received an image archive, set `archive` and `image` to the filenames
and tag in its delivery record, then verify its SHA256 and load it on the Docker
execution host:

```sh
sha256sum "$archive"
docker load -i "$archive"
docker image inspect "$image" --format '{{.Id}} {{.Os}}/{{.Architecture}}'
```

Confirm the image ID against the accompanying validation record. A source ZIP
or Git clone is not a container archive. These commands do not download a
Katana image from an unpublished registry. Match the archive's wheel/source
identity to your pipeline revision; the earlier installation-check image is not
an archive of current development code.

## Check before analysis

Check the actual image with generated installation inputs, without a network
connection or mounted host environment:

```sh
docker run --name katana-runtime-check --network none --user 1000:1000 \
  "$image" \
  python /opt/katana/check_downstream_environment.py --out /tmp/acceptance
docker cp katana-runtime-check:/tmp/acceptance ./container-acceptance
docker rm katana-runtime-check
```

`acceptance.json` must report `passed`. This verifies installed native tool
execution; it is not a real-data scientific workflow acceptance.

To exercise the actual image under Nextflow and verify cache reuse, use a new
output directory outside the checkout:

```sh
python pipeline/tools/check_container_runtime.py \
  --image "$image" --outdir /tmp/katana-nextflow-runtime-check
```

This check requires a non-root Linux user and Nextflow 25.10.4 or newer on
`PATH`. It executes three independent native dependency probes, verifies their
outputs, then requires all three to be cached on resume. It does not run a
study or download sequencing data.

## Use in the integrated workflow

Run existing-result analysis with the built image:

```sh
nextflow run pipeline/main.nf --mode downstream-only \
  --project inputs/project.json --flow_config config/flow.json \
  --katana_container "$image" \
  -profile docker --outdir results -resume
```

Docker workers on another machine need the same image loaded there or an
explicitly published immutable registry reference. Slurm users can select
`slurm,apptainer` or `slurm,singularity` and provide an image/SIF available to
their compute nodes. Site account, partition and resources are user settings.
No Dardel-specific path or account is embedded.

For new FASTQ/BAM processing, select `--mode complete` and supply the sample
sheet, contrast sheet, references and downstream resource configuration from
the [integrated workflow guide](../../docs/INTEGRATED_PIPELINE.md). The host does
not need a standalone Katana Conda environment for containerized tasks.
Upstream tool images are separate and must also be available to the executor.
The DEXSeq annotation/count steps require the separately
[built or loaded DEXSeq tools image](../dexseq-tools/README.md), selected with
`--dexseq_tools_container`. This is required for containerized `complete` and
`upstream-only` modes; no published DEXSeq registry image is assumed.

Container tasks use `/opt/conda`, clear inherited Python/R overrides, reject
host interpreter/library paths and probe required dependencies before analysis.
Use portable tool names such as `Rscript` in project metadata. Explicit
external runtimes remain available to standalone users outside the image.
FASTA, GTF and other analysis resources are separate validated inputs.
