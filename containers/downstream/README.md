# Katana downstream image

This image supplies the Python, R and native tools used by Katana downstream
tasks. Existing upstream modules retain their own nf-core tool environments.
The image contains no research data, genome references, motif libraries or
study-specific configuration. Qt is not installed.

Build from the accepted wheel and checked-in Linux dependency locks:

```sh
python scripts/build_downstream_container.py \
  --wheel dist/katana_splicing_tool-<version>-py3-none-any.whl \
  --context /tmp/katana-container-build \
  --tag katana-downstream:<version>
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

Check the actual image with generated installation inputs, without a network
connection or mounted host environment:

```sh
docker run --name katana-runtime-check --network none \
  katana-downstream:<version> \
  python /opt/katana/check_downstream_environment.py --out /tmp/acceptance
docker cp katana-runtime-check:/tmp/acceptance ./container-acceptance
docker rm katana-runtime-check
```

`acceptance.json` must report `passed`. This verifies installed native tool
execution; it is not a real-data scientific workflow acceptance.

Run existing-result analysis with the built image:

```sh
nextflow run pipeline/main.nf --mode downstream-only \
  --project inputs/project.json --flow_config config/flow.json \
  --katana_container katana-downstream:<version> \
  -profile docker --outdir results -resume
```

Docker workers on another machine need the same image loaded there or an
explicitly published immutable registry reference. Slurm users can select
`slurm,apptainer` or `slurm,singularity` and provide an image/SIF available to
their compute nodes. Site account, partition and resources are user settings.
No Dardel-specific path or account is embedded.

Container tasks use `/opt/conda`, clear inherited Python/R overrides, reject
host interpreter/library paths and probe required dependencies before analysis.
Use portable tool names such as `Rscript` in project metadata. Explicit
external runtimes remain available to standalone users outside the image.
FASTA, GTF and other analysis resources are separate validated inputs.
