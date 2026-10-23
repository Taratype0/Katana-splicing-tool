# Actual container installation and transport acceptance

This fixture starts three independent Docker tasks through the standard
Nextflow launcher: network, MAPP and motif runtime probes. Each invokes the
installed `katana-flow` CLI and its native runtime checks. It queries software
imports/R namespaces and verifies bundled author source identity; it does not
fit scientific models or load a study.

Build the runtime image first, then run on Linux/WSL with Docker and Nextflow:

```bash
python pipeline/tools/check_container_runtime.py \
  --image your-locally-built-katana-image \
  --nextflow nextflow \
  --outdir /absolute/new/path/katana-container-acceptance
```

The image must already exist locally. The checker resolves its immutable image
ID, disables container networking, and uses the caller's UID/GID. A root caller
must provide a non-root `--uid` and `--gid`; arrange write access to the new
output directory for that identity. `--nextflow` accepts a launcher path, not a
shell command. Configure `JAVA_HOME`/Docker connection normally in the caller.

Acceptance requires exactly three `COMPLETED` tasks with exit 0, then exactly
three `CACHED` tasks on `-resume`, with unchanged output hashes. Each probe must
reject an existing executable in mounted task storage outside `/opt/conda`
before it runs. This tests the boundary rather than a missing executable.

The evidence directory retains the image identity, copied fixtures, both
Nextflow logs/traces, task work, per-probe CLI/runtime records and final receipt.
Existing directories are never overwritten or cleaned. A failed attempt is not
automatically retried. These are installation/transport checks, not scientific
FASTQ/BAM end-to-end acceptance.
