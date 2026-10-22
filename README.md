# Katana

**An interactive workspace for RNA splicing analysis, regulatory evidence and co-splicing networks.**

Katana brings differential expression, splicing, transcript/exon usage and candidate evidence into one project. Explore individual events and genes, trace sequence and RBP evidence with **MotifTrace**, investigate modules and hubs with **Katana-Splice-Network**, and share an **English one-click report** with its supporting results.

This repository contains the complete application source baseline, its browser interface, scientific modules and command-line interfaces. The current registry contains **54 workflow entries**. Each entry requires its declared inputs and resources; a registered workflow is not a promise that every analysis runs on every dataset.

## What you can do

| Area | Capabilities |
| --- | --- |
| Bulk analysis | Import supported upstream results; inspect DAS, DEG, DTU and exon usage; compare selected contrasts; review replicate patterns and candidate evidence. |
| MotifTrace | Reconstruct supported event-level RNA sequence consequences, inspect motif changes and RBP evidence, and trace results to source loci, models and inputs. |
| Networks | Build qualified co-splicing networks, inspect modules and hubs, relate module scores to sample metadata, and examine preservation and functional evidence. |
| Single-cell results | Use supported processed single-cell inputs and dedicated result adapters. Input format, experimental unit and method eligibility remain explicit. |
| Reports and exports | Generate an English reproducibility evidence package or the existing interpreted study report, retaining tables, figures, parameters and provenance. |

## Install and open

Use Python 3.11 or newer in an isolated environment. From a source checkout:

```sh
git clone https://github.com/Taratype0/Katana-splicing-tool.git
cd Katana-splicing-tool
python -m pip install -e ".[next,next-dev]"
python -m katana_next.launcher --mode browser --port 8765 --runtime-root ../katana-runtime
```

Open the URL printed by the launcher. The source checkout includes the built browser interface; Node.js is needed only when rebuilding it. Keep project data, analysis outputs and runtime files outside the checkout. See [installation](docs/INSTALL.md) for environment setup, rebuilding and optional tools.

## Start an analysis

1. Open a supported Katana project manifest, or migrate an existing compatible Katana project.
2. Check sample identities, biological units, references and contrast direction.
3. Select the analysis and the contrasts you need; inspect readiness before running.
4. Browse completed results and generate the appropriate report. Reporting reads completed evidence and does not silently run missing analyses.

See [usage](docs/USAGE.md) for CLI examples, project migration and the two report routes.

## Scope and requirements

- Katana consumes supported analysis inputs and results. This application baseline is **not an official nf-core pipeline release** and does not itself provide a general raw FASTQ-to-results workflow.
- R/WGCNA, author-native tools, genome annotations, sequence files, motif libraries and CLIP resources are required only by the workflows that use them. They are not bundled as experimental data.
- Motif/CLIP evidence suggests regulatory hypotheses; hubs and associations are not proof of causality. Imported results retain their original method and statistical meaning.
- The current installation uses an editable source checkout. A fully self-contained wheel, a GUI-independent analysis package and the expanded Nextflow integration remain development work.

[Methods and limitations](docs/METHODS.md) describe the scientific boundaries. [Development architecture](docs/development/ARCHITECTURE.md) and the [roadmap](docs/development/ROADMAP.md) are separate from the software guide. Ongoing package development continues on the development branch.

## License and attribution

Katana's own code is distributed under the [MIT License](LICENSE). Bundled third-party components retain their original licenses, including the original DJExpress author core and browser fonts. See [third-party notices](THIRD_PARTY_NOTICES.md).
