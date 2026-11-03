# Install and run Katana

## Nextflow users

Use the [Nextflow quick start](NEXTFLOW_SIMPLE.md). If your server already has
Nextflow >=25.10.4, Java >=17 and a supported container engine, no Katana Python
installation is required. Nextflow downloads the versioned containers and
Katana fetches selected databases once into your resource cache. No server-specific
Katana installation path is used.

## Standalone downstream and browser users

Install the wheel from the matching [release](https://github.com/Taratype0/Katana-splicing-tool/releases/tag/v2.0.0-beta.5.dev2026110401)
into an isolated Python environment. The standalone scientific dependencies and
Conda instructions are in [installation](INSTALL.md). The Python API, `katana-flow`
and target-analysis commands remain available independently of Nextflow.

## Offline or portable controller installation

The optional controller archive includes Python, installed controller dependencies,
Nextflow and an exported workflow. Scientific tasks still use containers. Download
the matching controller and downstream SIF plus the unchanged DEXSeq SIF linked
in the release manifest. Check `SHA256SUMS` once when installing, then extract:

```bash
mkdir -p tools
tar -xzf katana-controller-2.0.0b5.dev2026110401-linux-x86_64.tar.gz -C tools
export PATH="$PWD/tools/katana-controller/bin:$PATH"
nextflow run tools/katana-controller/pipeline/main.nf --help
```

Offline sites can pass `--katana_container /images/katana.sif` and
`--dexseq_tools_container /images/dexseq.sif`, and use a pre-populated resource
cache or explicit database paths. These overrides are optional for online use.

## Cache locations

Put reusable caches on your project's storage, separately from task work:

```bash
export NXF_HOME=/data/cache/katana/nextflow
export NXF_SINGULARITY_CACHEDIR=/data/cache/katana/images
export SINGULARITY_CACHEDIR=/data/cache/katana/singularity
mkdir -p "$NXF_HOME" "$NXF_SINGULARITY_CACHEDIR" "$SINGULARITY_CACHEDIR"
```

For Apptainer use `NXF_APPTAINER_CACHEDIR` and `APPTAINER_CACHEDIR` instead.
Pass `--resource_cache /data/cache/katana/resources` and
`--digest_cache /data/cache/katana/checksums` to the pipeline. Download checksums
and changed-file checks remain enabled. Repeated reads of unchanged local files
reuse a small checksum database; `--verify_content true` forces full rereads.

Keep a running study on its original release and image. Updating a release does
not silently replace an active run's code or invalidate its successful work.
