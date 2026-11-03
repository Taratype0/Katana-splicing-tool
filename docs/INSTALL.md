# Installation

For downloadable wheels and matching Docker images, use the
[release installation guide](RELEASE_INSTALL.md). Source installation and local
builds remain available below. The published assets retain their original
accepted wheel/image identities.

## Choose an installation

Use **`main`** for the accepted workflow described here. Further development
continues on `codex/` development branches; use a recorded commit or backup tag
when reproducing a specific version. Choose one of these routes:

First choose the [entry point for your inputs](ENTRY_POINTS.md). Use the
[independent runtime guide](INDEPENDENT_RUNTIME.md) and the matching release's
`RELEASE_MANIFEST.json`, `RUNTIME_VALIDATION.json` and `SHA256SUMS` for current
software, image and validation identities. The November integration uses
whole-study CoSplice and EventGraph MotifTrace; older wheels/images do not gain
these engines from an updated source checkout.

| Route | What it supplies | Host requirements |
| --- | --- | --- |
| [Nextflow with the downstream image](../containers/downstream/README.md) | Katana Python, R and native analysis tools inside each downstream task; upstream modules use their own tool images | Linux execution host, Nextflow 25.10.4 or newer and compatible Java; actual native acceptance uses Linux x86-64 Docker. Apptainer/Singularity profiles need site qualification. |
| [Standalone Linux Conda](RUNTIME_ENVIRONMENTS.md#standalone-linux-environment) | The full native downstream runtime for direct CLI, Python and browser use | Linux x86-64 and Conda or micromamba |
| [Python source installation](#source-checkout) | Browser, CLI, reports and Python analyses | Python 3.11 or newer; selected native computations need their additional runtime |

Earlier checks are retained in the dated [October entry/runtime record](validation/entrypoints-runtime-20261025.md)
and [completed-result record](validation/completed-results-runtime-20261025.md).
They describe their recorded versions, not the current engines. New engine
integration checks use generated inputs; they do not imply a new real-study run.
Python wheels are distributed through GitHub Releases. Nextflow profiles select
the matching Docker registry image or downloadable Apptainer/Singularity SIF.
No PyPI publication is advertised.

## Versions and rollback

`main` contains the accepted software. `codex/` branches are used for
subsequent development. Record the exact Git commit with each installation;
existing wheel/image acceptance records keep their original source hashes.
The downloadable assets preserve the accepted builds. Match the release tag
in your command to its wheel and container versions; PyPI is not required.

The previous `main` is preserved by the annotated tag
[`archive/pre-main-package-20261025`](https://github.com/Taratype0/Katana-splicing-tool/tree/archive/pre-main-package-20261025),
at commit `c873ee8e97c0ef626234f4f9698280d79c94c57e`.
To inspect that version, use a separate checkout and environment:

```sh
git clone --branch archive/pre-main-package-20261025 https://github.com/Taratype0/Katana-splicing-tool.git Katana-previous
```

Do not switch a running installation or existing study directory in place.
The tag preserves the old code; it does not rewrite projects or completed
results into an older format.

## Installed-package report commands

Use one matched software revision across the wheel, image and exported
Nextflow sources. A shared package version string alone is not sufficient:
verify the supplied wheel SHA256, source commit, image ID and dependency-lock
identities in its delivery/acceptance record. Existing accepted artifacts keep
their original identities; a newer source change does not update an old image.
For the current public image addresses, see the
[independent runtime guide](INDEPENDENT_RUNTIME.md).

For report assembly only, install the supplied wheel with its report extra in
a fresh Python 3.11+ environment, replacing the example filename with the exact
wheel you received:

```sh
python -m pip install "/path/to/katana_splicing_tool-<version>-py3-none-any.whl[report,next]"
python -m pip check
katana-report --help
```

Run from the study workspace, outside the source checkout:

```sh
katana-report --project results.json --outdir exports/study-report \
  --compact --no-archive --result-store payload-store
```

`results.json` can be an existing project/descriptor or an explicitly bound
multi-source completed-result set. See the [two-source example, path mappings
and identity checks](EXISTING_RESULTS.md#generate-one-report-from-several-completed-sources).
Use `--archive downloads/study-report.zip` instead of `--no-archive` to produce a
ZIP at a separate new location. The same command is available as
`python -m katana_next.reports.cli`; it uses the installed package and preserves
the approved seven-section format. Native differential/network/MAPP fits are
not rerun to assemble a report.

For new downstream calculations, use the [locked Linux Conda runtime](RUNTIME_ENVIRONMENTS.md#standalone-linux-environment)
plus the matching wheel, or a matching downstream image. A report-only wheel
installation does not supply R and the external analysis toolchain. Integrated
Nextflow users obtain that toolchain from task containers and do not need a
second host Conda installation.

An installed image can run the same report command. Set `KATANA_IMAGE` to the
locally built/loaded tag or immutable image reference in its delivery record.
From a workspace containing the completed inputs, this Linux Docker example
mounts those inputs read-only and writes a new report separately:

```sh
mkdir -p exports
docker run --rm --network none --user "$(id -u):$(id -g)" \
  --mount "type=bind,src=$PWD,dst=/data,readonly" \
  --mount "type=bind,src=$PWD/exports,dst=/exports" \
  --workdir /data "$KATANA_IMAGE" \
  katana-report --project results.json --outdir /exports/study-report \
  --compact --no-archive
```

Source descriptors must use relative paths within the mounted inputs or an
explicit `--path-map` for the container's locations. Standalone Conda users run
the same command with `micromamba run -n katana katana-report ...` after
installing the native/pip locks and matching wheel. Genome references and
analysis databases remain project resources rather than installer contents.

### Automatic container downloads

The normal Nextflow command needs no manual image download or host Katana
Python installation. `-profile docker` selects versioned registry images;
`-profile singularity` or `-profile apptainer` selects the matching SIF files.
This includes the DEXSeq tools image for FASTQ/BAM analysis. Other upstream
modules retain their own pinned nf-core tool images.

For offline sites, preload the images and pass `--katana_container` and
`--dexseq_tools_container` with local paths. Keep reusable image and database
caches outside task work. The [runtime guide](INDEPENDENT_RUNTIME.md) documents
cache locations and the optional portable controller; the
[database guide](RESOURCE_DATABASES.md) lists automatic and custom resources.

Local builds are optional: follow the [container instructions](../containers/downstream/README.md#build-locally)
when developing or changing dependencies. An existing running study should
continue with its original release and image so `-resume` can reuse its work.

## Source checkout

Use a fresh environment so an older installation cannot supply hidden dependencies:

```sh
git clone --branch main https://github.com/Taratype0/Katana-splicing-tool.git
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
python -m pip install -e ".[next,flow,contracts]"
python -m pip check
katana-open --port 8765 --runtime-root ../katana-runtime
```

The `next` extra supplies the application service, artifact storage and report dependencies; `flow` and `contracts` include the workflow interface and schema validation. Developers can add test tools with `python -m pip install -e '.[next,flow,contracts,next-dev]'`. Qt is optional: install the `desktop` extra only for the compatibility desktop mode. Browser and CLI installations do not require PySide6. Editable installation resolves the declared dependency ranges; use the standalone locks or accepted image when reproducing the complete native runtime.

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
source, use the built browser assets already included in the checkout. If the
frontend has changed, rebuild it with `pnpm build` in `apps/web` first, then run
`python -m build --wheel` in the repository root. Set `KATANA_REQUIRE_WEB=1` in the build environment to fail
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
