# Installation

## Choose an installation

Use the **`codex/package-development`** branch for the packaged workflow described here. `main` can lag behind it. Choose one of these routes:

| Route | What it supplies | Host requirements |
| --- | --- | --- |
| [Nextflow with the downstream image](../containers/downstream/README.md) | Katana Python, R and native analysis tools inside each downstream task; upstream modules use their own tool images | Linux execution host, Nextflow 25.10.4 or newer, compatible Java, Docker or supported Apptainer/Singularity setup |
| [Standalone Linux Conda](RUNTIME_ENVIRONMENTS.md#standalone-linux-environment) | The full native downstream runtime for direct CLI, Python and browser use | Linux x86-64 and Conda or micromamba |
| [Python source installation](#source-checkout) | Browser, CLI, reports and Python analyses | Python 3.11 or newer; selected native computations need their additional runtime |

For the current development image and its matching wheel, use the
[current container acceptance record](validation/pipeline-closure-20261024.md).
The [earlier runtime record](validation/downstream-runtime-20261023.md) documents
the standalone Linux and Windows Python/browser installations; its older image
does not include subsequent workflow/report fixes. Match the wheel, image and
source revision instead of selecting by version string alone. No Katana image
has been published to a registry and no PyPI release is advertised.
[Generic workflow acceptance](validation/generic-pipeline-20261023.md) is tracked
separately from installation checks.

### Load the current development container

On a Linux x86-64 Docker host, place the supplied archive and checksum file in
the same directory, then run:

```sh
sha256sum --check SHA256SUMS
docker load -i katana-downstream-dev-d05c76c-linux-amd64.tar.gz
docker image inspect katana-downstream:dev-d05c76c --format '{{.Id}} {{.Os}}/{{.Architecture}}'
docker run --rm --network none katana-downstream:dev-d05c76c katana-flow --help
```

For **new FASTQ/BAM processing**, also load the supplied DEXSeq tools archive:

```sh
docker load -i kanata-dexseq-tools-dev-20261024-linux-amd64.tar.gz
docker image inspect kanata-dexseq-tools:dev-20261024 --format '{{.Id}}'
```

Set `--dexseq_tools_container kanata-dexseq-tools:dev-20261024` when starting
`complete` or `upstream-only` mode with Docker. Existing-result downstream-only
analysis does not need this upstream image. If no archive was supplied,
[build it locally](../containers/dexseq-tools/README.md); do not try to pull the
previously assumed, unpublished registry tag. The runtime record lists both
archive hashes and image IDs. Upstream modules also use their pinned public
tool images, which the execution host must be able to obtain.

Compare the image ID with the current acceptance record. The archive is a local
development deliverable, not a registry download; users without it can follow
the [build instructions](../containers/downstream/README.md#build-locally).
The host needs Docker and, for the integrated pipeline, Java/Nextflow. Katana's
Python/R/native dependencies are inside the image; a host Conda environment is
not needed. Standalone CLI/browser users can use the locked Conda route below.

Installation does not require recreating an existing analysis plan. To resume
unchanged downstream work, use the original `--flow_plan` path with `-resume`
and retain its inputs and Nextflow cache; no `--flow_config` is needed for that
route. See [plan validation and recovery commands](INTEGRATED_PIPELINE.md#resume-an-unchanged-downstream-plan).
The same guide documents the explicit, presentation-only recovery of a verified
MotifTrace display failure; installing a new wheel does not automatically
rerun a completed computation.

## Source checkout

Use a fresh environment so an older installation cannot supply hidden dependencies:

```sh
git clone --branch codex/package-development https://github.com/Taratype0/Katana-splicing-tool.git
cd Katana-splicing-tool
python -m venv .venv
```

Activate it on Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Or on Linux/macOS:

```sh
. .venv/bin/activate
```

Then install and start the browser interface:

```sh
python -m pip install -e ".[next]"
python -m pip check
katana-open --port 8765 --runtime-root ../katana-runtime
```

The `next` extra supplies the application service, artifact storage and report dependencies. Developers can add test tools with `python -m pip install -e '.[next,next-dev]'`. Qt is optional: install the `desktop` extra only for the compatibility desktop mode. Browser and CLI installations do not require PySide6. Editable installation resolves the declared dependency ranges; use the standalone locks or accepted image when reproducing the complete native runtime.

The launcher opens a browser by default and prints its actual URL. If the requested port is occupied, it can select another free port. Launch logs and service metadata go into the external runtime directory. `python -m katana_next.launcher` provides the same entry point when invoking a specific Python environment.

The browser is the current interface and release target. Explicit `--mode desktop` and `--mode choose` retain the existing compatibility modes where the `desktop` extra, Qt WebEngine and a graphical session are available.

## Full standalone downstream runtime

For reproducible Linux downstream computation, install the native and pip locks,
then the matching wheel as shown in [Runtime environments](RUNTIME_ENVIRONMENTS.md#standalone-linux-environment).
This supplies the browser/API/CLI Python
dependencies alongside R WGCNA, DESeq2, DEXSeq, fgsea, NetRep, Rcpp, RSeQC and
samtools. The exact accepted Conda builds and additional pip distributions are
kept under `environments/locks`. The YAML plus the wheel's `standalone` extra is
an alternative when deliberately solving a new environment; it can choose
different dependency builds and must be checked again.

The Python-only route above is useful for reading results, reports and the
browser, but does not install R or native analysis tools. Windows supports this
portable route; a Linux container or a Linux installation runs the native
toolchain. WSL is one optional local Linux host, not a requirement for a remote
server or container host.

## Development wheels

A browser-enabled development wheel includes the built interface; users do not
need Node.js or the original source checkout. Install its service and selected
target dependencies with:

```sh
python -m pip install 'katana_splicing_tool-<version>-py3-none-any.whl[next,targets]'
katana-open
katana-targets --help
```

Replace the placeholder with the actual wheel filename. Wheels are development
artifacts until a formal release is accepted. To build a browser wheel from
source, first run `pnpm build` in `apps/web`, then `python -m build --wheel` in
the repository root. Set `KATANA_REQUIRE_WEB=1` in the build environment to fail
if browser assets have not been built. A CLI-only build can omit those assets.
The same wheel supports local and server execution; BAM/reference paths are
configured on the machine running Python. See [target workflows](TARGET_COMMANDS.md).

### Verify the installed wheel outside the checkout

Use the newly created environment after installing the wheel. Choose a writable
directory outside any source checkout and clear implicit source imports with
Python's isolated mode:

```sh
mkdir -p ~/katana-install-check
cd ~/katana-install-check
python -I -m pip check
python -I -m katana_next.flow.cli --help
python -I -m katana_motiftrace.cli --help
python -I -m scripts.check_downstream_environment --portable-only --out portable-check
python -I -m katana_next.flow.cli export-pipeline --outdir exported-pipeline
```

On a fully installed Linux native environment, omit `--portable-only` and use a
new output directory for the full generated-runtime checks. These checks do not
open a research project. An existing output/export directory is preserved;
choose a new one for another check. The wheel contains the built browser,
flow/report modules and Nextflow sources, so these commands must work without
the original checkout on `PYTHONPATH`.

## Browser interface

`katana-open` stores launcher state in your user data directory: `%LOCALAPPDATA%/Katana/runtime` on Windows, `$XDG_DATA_HOME/katana/runtime` (or `~/.local/share/katana/runtime`) on Linux, and `~/Library/Application Support/Katana/runtime` on macOS. It does not write into the installed Python package. Set `--runtime-root /writable/path` to select another location; existing explicit locations remain supported.

The browser communicates with a Katana Python service on the machine that holds the project files and runs its calculations. Closing a browser tab does not turn an active job into a completed result. Offline report HTML can display saved results but does not provide this execution service.

The checkout includes `apps/web/dist`, so a frontend build is not required for ordinary startup. To modify and rebuild it, install a compatible Node.js and pnpm environment and run:

```sh
cd apps/web
pnpm install --frozen-lockfile
pnpm build
cd ../..
```

Do not copy `node_modules` into a release. Preserve bundled dependency and font licenses when redistributing a built interface.

### Server access

On the analysis server, start or reuse the service without opening a graphical window:

```sh
katana-open --no-open --port 8765 --runtime-root ~/katana-runtime
```

The service listens on the server's loopback address. From your local machine, forward its printed port over your existing SSH connection, for example:

```sh
ssh -L 8765:127.0.0.1:8765 user@server
```

Open `http://127.0.0.1:8765/` in your local browser. Replace `user@server` with your login and server; if the launcher selected another port, use that port on the server side of the forwarding rule. Analysis inputs and BAMs stay on the machine running the service. Cluster job submission still uses the configured execution route; opening the browser does not submit a Slurm job or start an analysis.

## Scientific tools and separate data resources

The accepted downstream image and standalone Linux lock include the supported
native toolchain below. Integrated container users do not need to repair a host
Conda environment. These additional installations apply to a Python-only or
custom standalone installation. FASTA/GTF, motif libraries, CLIP and pathway
data are always explicit project resources.

| Analysis | Additional requirements |
| --- | --- |
| Original WGCNA network route | R with WGCNA and its dependencies, supplied by the native runtime. For a custom standalone environment, configure or expose `Rscript`; `scripts/install_wgcna.R` is a helper, not an exact dependency lock. |
| DJExpress author core | Its declared R dependencies and suitable junction-count inputs. The vendored author source retains its own notices and source hashes. |
| Motif and sequence analyses | Matched reference FASTA/GTF, the appropriate motif library and species/build mappings. The included FIMO worker requires `pymemesuite` in its selected Python runtime; reference CLIP is a separate optional input. |
| Protein/domain annotations | Matched transcript/protein references, a suitable Pfam HMM resource and `pyhmmer` in the selected worker runtime. |
| Sashimi and BAM-derived routes | Eligible aligned BAMs/indexes. Depending on the selected route, use `pysam`, configured samtools or the supported BAM reader; these choices are not interchangeable for every workflow. |
| Native MAPP MAEI | The bundled, source-verified author model, R/Rcpp and compilers are included in the native runtime. Matching transcript quantities, reference FASTA/GTF and motif-count inputs remain required; other MAPP subworkflows are outside this route. |
| Other external methods | The original method's executable/runtime and its supported input schema. Importing previous results does not install or execute that method. |

Readiness checks identify the resources needed for each workflow. Missing data or tools do not represent zero significant findings. Do not substitute a different genome build, statistical method or experimental unit merely to pass a check.

Install optional worker dependencies only in the runtime used for the selected route. The Linux `standalone` extra includes the fixed PyMEME FIMO worker and PyHMMER. Some adapters can use WSL when explicitly configured; WSL is not a global Katana installation requirement. The report reader uses Beautiful Soup, and the included resource-download helper uses Requests; both are supplied by the application's `next` extra.

## Distribution boundary

This repository contains software and synthetic fixtures, not a research dataset or preconfigured workstation. Large references, raw reads, BAMs, personal project manifests, private runtime configuration and generated results belong outside the checkout. Installation checks verify the software runtime; they do not establish completion or scientific validity of an individual study.
