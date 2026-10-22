# Installation

## Source checkout

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

Then install and start the application:

```sh
python -m pip install -e ".[next,next-dev]"
python -m katana_next.launcher --mode browser --port 8765 --runtime-root ../katana-runtime
```

The `next` extra supplies the application service, artifact storage and report dependencies. `next-dev` supplies test dependencies. The current package still declares desktop/GUI dependencies even in browser mode; this is not a minimal headless distribution.

The launcher prints its actual URL. If the requested port is occupied, it can select another free port. Launch logs and service metadata go into the external runtime directory. Desktop mode is available with `--mode desktop` where Qt WebEngine and a graphical session are supported.

## Browser interface

The checkout includes `apps/web/dist`, so a frontend build is not required for ordinary startup. To modify and rebuild it, install a compatible Node.js and pnpm environment and run:

```sh
cd apps/web
pnpm install --frozen-lockfile
pnpm build
cd ../..
```

Do not copy `node_modules` into a release. Preserve bundled dependency and font licenses when redistributing a built interface.

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

Install optional worker dependencies only in the runtime used for the selected route. Some adapters can use WSL when explicitly configured; WSL is not a global Katana installation requirement. The English approved-reader report route uses Beautiful Soup, and the included resource-download helper uses Requests; both are supplied by the application's `next` extra.

## Distribution boundary

This repository contains software and synthetic fixtures, not a research dataset or preconfigured workstation. Large references, raw reads, BAMs, personal project manifests, private runtime configuration and generated results belong outside the checkout. A self-contained wheel/container and the future GUI-independent package require separate packaging and acceptance work.
