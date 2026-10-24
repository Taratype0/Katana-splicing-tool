# Development and reproducibility audit — 24 October 2026

This audit separates an accepted study deliverable from acceptance of the
software that another user installs. It supersedes earlier statements that
the Primary culture report was still paused or incomplete; the historical
runtime records retain the identities of the artifacts they actually tested.

## Source and CI

At this audit, the public `codex/package-development` software revision was
`9bd4877b3a29aac38ad88cb476d64915183d0b23`. The default `main` branch was still
`c873ee8e97c0ef626234f4f9698280d79c94c57e` and does not contain the subsequent
packaging and network fixes. Follow the development-branch installation
instructions when evaluating those changes.

The development branch contains the full-event network and bounded edge-export
repair (`6105ec2`), bootstrap-off default and continuous kME ranking (`1885a44`),
tree-cut and thread-budget performance repair (`5ca6533`), exact report kME
ordering (`ad0ab60`), and ZIP timestamp compatibility (`9bd4877`).

[Software CI at 9bd4877](https://github.com/Taratype0/Katana-splicing-tool/actions/runs/40130647595)
completed successfully, covering Linux/Windows application checks, the browser
build, and installed-wheel tests outside the checkout.
[Native Nextflow CI at f69d953](https://github.com/Taratype0/Katana-splicing-tool/actions/runs/40112127794)
also completed successfully. Neither CI run establishes a full native-toolchain
release test for every optional analysis.

## What has already been exercised

- Generated FASTQ and BAM entries completed 101/82 upstream tasks respectively,
  followed by 55 downstream tasks per entry and seven-section report checks.
  See the [historical integration record](pipeline-closure-20261024.md).
- A generated two-comparison existing-result workflow exercised separate
  networks, multiple motif models, significant enrichment with hit genes,
  collection and the seven-section report. See the
  [existing-result acceptance](generic-pipeline-20261023.md).
- The target-analysis workflow passed browser selection, portable request
  export, installed-package Gene/Event cards and indexed-BAM Sashimi checks.
  Those checks used generated fixtures and retained completed-source bytes
  and caches; they are not biological validation.
- The real Primary culture deliverable now includes 12 accepted independent
  networks and 168 completed MAPP fit families, alongside its retained
  MotifTrace and differential results. Its seven-section report, figure-data
  links, offline browser views and complete ZIP were accepted. This report
  integration did not rerun completed scientific analyses or modify the
  original Katana application.

## Remaining before a current one-command release claim

The public package already provides `katana-report`, `katana-flow report` and
the Nextflow report task. However, the final study above also used local
presentation/integration adapters that are not yet public generic entry points:

1. Explicit identity-checked assembly of separately completed historical
   collections into one report, without assuming a prior report snapshot or
   workstation paths.
2. Portable relative MAPP fit-index links for the complete multi-collection
   report, direct main-report links to complete module browsers, and the final
   network-enrichment plot-padding corrections.
3. Regression acceptance of these adapters through the installed public
   package and its documented report command, using completed results only.
4. A matching final wheel and native container built from the final public
   software revision, followed by installed-runtime, report, and resume checks.
   The earlier accepted image `dev-9546632` predates the network repairs listed
   above and must not be presented as containing them.
5. Updated installation/artifact identities and release documentation after
   those checks pass. No current PyPI or container-registry release is implied.

Passing source CI and delivering one complete real report do not close these
remaining distribution checks. In particular, the exact final historical study
assembly is not yet reproducible through a single documented public command.
New users still supply their own samples, selected contrasts, matching
references and databases; study data and local project configurations are not
distributed with the software.

Accepted platform scopes remain Linux x86-64 native/container execution and
Windows Python/browser/report use. Native macOS/ARM and arbitrary scheduler or
site combinations have not been certified by these checks.
