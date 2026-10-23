# Installation

## Source checkout

For the distinction between integrated Nextflow task environments and a direct
Python/Conda installation, see [Runtime environments](RUNTIME_ENVIRONMENTS.md).
That page records the tested runtime, exact dependency locks and remaining
platform limits. A container profile alone does not install a missing image.

The supported installation route for this baseline is an editable source checkout with Python 3.11 or newer. Use a fresh environment so an older installation cannot supply hidden dependencies.

```sh
git clone https://github.com/Taratype0/Katana-splicing-tool.git
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
python -m pip install -e ".[next,next-dev]"
katana-open --port 8765 --runtime-root ../katana-runtime
```

The `next` extra supplies the application service, artifact storage and report dependencies. `next-dev` supplies test dependencies. Qt is optional: install the `desktop` extra only for the compatibility desktop mode. Browser and CLI installations do not require PySide6.

The launcher opens a browser by default and prints its actual URL. If the requested port is occupied, it can select another free port. Launch logs and service metadata go into the external runtime directory. `python -m katana_next.launcher` provides the same entry point when invoking a specific Python environment.

The browser is the current interface and release target. Explicit `--mode desktop` and `--mode choose` retain the existing compatibility modes where the `desktop` extra, Qt WebEngine and a graphical session are available.

## Full standalone downstream runtime

For Linux downstream computation, use `environments/standalone-linux.yml` and
install the wheel's `standalone` extra. This supplies the browser/API/CLI Python
dependencies alongside R WGCNA, DESeq2, DEXSeq, fgsea, NetRep, Rcpp, RSeQC and
samtools. The exact accepted Conda builds and additional pip distributions are
kept under `environments/locks`; see [Runtime environments](RUNTIME_ENVIRONMENTS.md)
for reproducible installation and verification commands.

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

## Optional scientific tools and resources

| Analysis | Additional requirements |
| --- | --- |
| Original WGCNA network route | An accessible R installation with WGCNA and its dependencies. Configure or expose `Rscript`; use `scripts/install_wgcna.R` as an installation helper. |
| DJExpress author core | Its declared R dependencies and suitable junction-count inputs. The vendored author source retains its own notices and source hashes. |
| Motif and sequence analyses | Matched reference FASTA/GTF, the appropriate motif library and species/build mappings. The included FIMO worker requires `pymemesuite` in its selected Python runtime; reference CLIP is a separate optional input. |
| Protein/domain annotations | Matched transcript/protein references, a suitable Pfam HMM resource and `pyhmmer` in the selected worker runtime. |
| Sashimi and BAM-derived routes | Eligible aligned BAMs/indexes. Depending on the selected route, use `pysam`, configured samtools or the supported BAM reader; these choices are not interchangeable for every workflow. |
| Native external methods | The original method's executable/runtime and its supported input schema. Importing previous results does not install or execute that method. |

Readiness checks identify the resources needed for each workflow. Missing data or tools do not represent zero significant findings. Do not substitute a different genome build, statistical method or experimental unit merely to pass a check.

Install optional worker dependencies only in the runtime used for the selected route. The Linux `standalone` extra includes the fixed PyMEME FIMO worker and PyHMMER. Some adapters can use WSL when explicitly configured; WSL is not a global Katana installation requirement. The report reader uses Beautiful Soup, and the included resource-download helper uses Requests; both are supplied by the application's `next` extra.

## Distribution boundary

This repository contains software and synthetic fixtures, not a research dataset or preconfigured workstation. Large references, raw reads, BAMs, personal project manifests, private runtime configuration and generated results belong outside the checkout. Installation checks verify the software runtime; they do not establish completion or scientific validity of an individual study.
