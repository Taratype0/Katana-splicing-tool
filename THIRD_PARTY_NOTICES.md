# Third-party notices

Katana's own source is covered by the repository MIT license. Third-party source, browser code and fonts retain their original licenses; that MIT notice does not replace theirs.

## nf-core upstream modules

The Nextflow workflow includes and adapts nf-core modules under their retained
[MIT license](pipeline/UPSTREAM_LICENSE). The auto-strandedness interface records
its nf-core revision, RSeQC/gffread versions and adaptations in
[method provenance](pipeline/modules/local/infer_strandedness.provenance.json).
RSeQC and gffread execute from their declared external tool environments; their
source is not included in this adapter and retains its own licensing.

## DJExpress author core

The original DJExpress core is included under `katana_next/vendor/DJExpress` at the recorded upstream revision `a0396272fd0522565190d97f4ecd2e9ba8a98366`. Copyright is held by Lina Gallego-Paez. Its MIT notices are preserved in [LICENSE](katana_next/vendor/DJExpress/LICENSE) and [LICENSE.md](katana_next/vendor/DJExpress/LICENSE.md). The original DESCRIPTION, NAMESPACE and [source-hash manifest](katana_next/vendor/DJExpress_SOURCE_SHA256.json) remain with that source. Katana does not claim authorship of this method.

## MAPP author MAEI source

The original MAEI scripts and kmer counter are included under
`katana_next/vendor/MAPP` from upstream commit
`6ad6e12275e6d7aa1254d1a8933cfbce6766b16c`. Their contents are unchanged;
`SOURCE.json` records the source URLs, sizes and SHA256 values. This supplies the
existing standalone MAEI adapter, not the complete upstream MAPP pipeline.

Licenses differ by file: the original `LICENSE` and the counter/postprocessing
scripts specify Apache-2.0; `model.cpp` and
`optimize-MARA-model-parameters.R` specify GPL-3.0; the original R Core Team
`optim.c` and `zeroin.c` specify GPL-2.0-or-later. Their original copyright and
license headers and full `COPYING-GPL-2.0` / `COPYING-GPL-3.0` texts are retained.
Katana's MIT notice does not replace any of these licenses. The native model is
compiled from private copies of this supplied source at execution time.

## Browser interface and fonts

The source checkout includes a built browser interface. Original license texts for a conservative production dependency closure of the committed pnpm lockfile are retained below. This includes optional/transitive dependencies and can be larger than the set of modules retained by the JavaScript bundler; inclusion in this inventory is not a claim that a package executes in the browser.

IBM Plex Sans and Newsreader fonts are licensed under SIL Open Font License 1.1. React/React DOM and TanStack packages use MIT; Lucide uses ISC; Vega family packages use BSD-3-Clause. The exact package/version notices below govern each component. License files were copied from installed npm packages whose package names and versions match `apps/web/pnpm-lock.yaml`; package metadata and copied-file hashes are recorded in [inventory.json](THIRD_PARTY_LICENSES/inventory.json). No workstation paths or installed package directories are included in that provenance.

| Package | Version | Declared license | Original text |
| --- | --- | --- | --- |
| `@fontsource/ibm-plex-sans` | 5.3.0 | OFL-1.1 | [LICENSE](THIRD_PARTY_LICENSES/npm-fontsource__ibm-plex-sans-5.3.0/LICENSE) |
| `@fontsource/newsreader` | 5.3.0 | OFL-1.1 | [LICENSE](THIRD_PARTY_LICENSES/npm-fontsource__newsreader-5.3.0/LICENSE) |
| `@tanstack/query-core` | 5.102.8 | MIT | [LICENSE](THIRD_PARTY_LICENSES/npm-tanstack__query-core-5.102.8/LICENSE) |
| `@tanstack/react-query` | 5.102.8 | MIT | [LICENSE](THIRD_PARTY_LICENSES/npm-tanstack__react-query-5.102.8/LICENSE) |
| `@tanstack/react-table` | 8.21.3 | MIT | [LICENSE](THIRD_PARTY_LICENSES/npm-tanstack__react-table-8.21.3/LICENSE) |
| `@tanstack/react-virtual` | 3.14.11 | MIT | [LICENSE](THIRD_PARTY_LICENSES/npm-tanstack__react-virtual-3.14.11/LICENSE) |
| `@tanstack/table-core` | 8.21.3 | MIT | [LICENSE](THIRD_PARTY_LICENSES/npm-tanstack__table-core-8.21.3/LICENSE) |
| `@tanstack/virtual-core` | 3.17.9 | MIT | [LICENSE](THIRD_PARTY_LICENSES/npm-tanstack__virtual-core-3.17.9/LICENSE) |
| `@types/estree` | 1.0.9 | MIT | [LICENSE](THIRD_PARTY_LICENSES/npm-types__estree-1.0.9/LICENSE) |
| `@types/geojson` | 7946.0.16 | MIT | [LICENSE](THIRD_PARTY_LICENSES/npm-types__geojson-7946.0.16/LICENSE) |
| `ansi-regex` | 6.3.0 | MIT | [license](THIRD_PARTY_LICENSES/npm-ansi-regex-6.3.0/license) |
| `ansi-styles` | 6.2.3 | MIT | [license](THIRD_PARTY_LICENSES/npm-ansi-styles-6.2.3/license) |
| `cliui` | 9.0.1 | ISC | [LICENSE.txt](THIRD_PARTY_LICENSES/npm-cliui-9.0.1/LICENSE.txt) |
| `commander` | 2.20.3 | MIT | [LICENSE](THIRD_PARTY_LICENSES/npm-commander-2.20.3/LICENSE) |
| `commander` | 7.2.0 | MIT | [LICENSE](THIRD_PARTY_LICENSES/npm-commander-7.2.0/LICENSE) |
| `d3-array` | 3.2.4 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-array-3.2.4/LICENSE) |
| `d3-color` | 3.1.0 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-color-3.1.0/LICENSE) |
| `d3-delaunay` | 6.0.4 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-delaunay-6.0.4/LICENSE) |
| `d3-dispatch` | 3.0.1 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-dispatch-3.0.1/LICENSE) |
| `d3-dsv` | 3.0.1 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-dsv-3.0.1/LICENSE) |
| `d3-ease` | 3.0.1 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-ease-3.0.1/LICENSE) |
| `d3-force` | 3.0.0 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-force-3.0.0/LICENSE) |
| `d3-format` | 3.1.2 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-format-3.1.2/LICENSE) |
| `d3-geo-projection` | 4.0.0 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-geo-projection-4.0.0/LICENSE) |
| `d3-geo` | 3.1.1 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-geo-3.1.1/LICENSE) |
| `d3-hierarchy` | 3.1.2 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-hierarchy-3.1.2/LICENSE) |
| `d3-interpolate` | 3.0.1 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-interpolate-3.0.1/LICENSE) |
| `d3-path` | 3.1.0 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-path-3.1.0/LICENSE) |
| `d3-quadtree` | 3.0.1 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-quadtree-3.0.1/LICENSE) |
| `d3-scale-chromatic` | 3.1.0 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-scale-chromatic-3.1.0/LICENSE) |
| `d3-scale` | 4.0.2 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-scale-4.0.2/LICENSE) |
| `d3-shape` | 3.2.0 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-shape-3.2.0/LICENSE) |
| `d3-time-format` | 4.1.0 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-time-format-4.1.0/LICENSE) |
| `d3-time` | 3.1.0 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-time-3.1.0/LICENSE) |
| `d3-timer` | 3.0.1 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-d3-timer-3.0.1/LICENSE) |
| `delaunator` | 5.1.0 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-delaunator-5.1.0/LICENSE) |
| `emoji-regex` | 10.6.0 | MIT | [LICENSE-MIT.txt](THIRD_PARTY_LICENSES/npm-emoji-regex-10.6.0/LICENSE-MIT.txt) |
| `escalade` | 3.2.0 | MIT | [license](THIRD_PARTY_LICENSES/npm-escalade-3.2.0/license) |
| `fast-json-patch` | 3.1.1 | MIT | [LICENSE.txt](THIRD_PARTY_LICENSES/npm-fast-json-patch-3.1.1/LICENSE.txt) |
| `get-caller-file` | 2.0.5 | ISC | [LICENSE.md](THIRD_PARTY_LICENSES/npm-get-caller-file-2.0.5/LICENSE.md) |
| `get-east-asian-width` | 1.6.0 | MIT | [license](THIRD_PARTY_LICENSES/npm-get-east-asian-width-1.6.0/license) |
| `iconv-lite` | 0.6.3 | MIT | [LICENSE](THIRD_PARTY_LICENSES/npm-iconv-lite-0.6.3/LICENSE) |
| `internmap` | 2.0.3 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-internmap-2.0.3/LICENSE) |
| `json-stringify-pretty-compact` | 4.0.0 | MIT | [LICENSE](THIRD_PARTY_LICENSES/npm-json-stringify-pretty-compact-4.0.0/LICENSE) |
| `lucide-react` | 1.42.0 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-lucide-react-1.42.0/LICENSE) |
| `react-dom` | 19.2.8 | MIT | [LICENSE](THIRD_PARTY_LICENSES/npm-react-dom-19.2.8/LICENSE) |
| `react` | 19.2.8 | MIT | [LICENSE](THIRD_PARTY_LICENSES/npm-react-19.2.8/LICENSE) |
| `robust-predicates` | 3.0.3 | Unlicense | [LICENSE](THIRD_PARTY_LICENSES/npm-robust-predicates-3.0.3/LICENSE) |
| `rw` | 1.3.3 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-rw-1.3.3/LICENSE) |
| `safer-buffer` | 2.1.2 | MIT | [LICENSE](THIRD_PARTY_LICENSES/npm-safer-buffer-2.1.2/LICENSE) |
| `scheduler` | 0.27.0 | MIT | [LICENSE](THIRD_PARTY_LICENSES/npm-scheduler-0.27.0/LICENSE) |
| `semver` | 7.8.5 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-semver-7.8.5/LICENSE) |
| `string-width` | 7.2.0 | MIT | [license](THIRD_PARTY_LICENSES/npm-string-width-7.2.0/license) |
| `strip-ansi` | 7.2.0 | MIT | [license](THIRD_PARTY_LICENSES/npm-strip-ansi-7.2.0/license) |
| `topojson-client` | 3.1.0 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-topojson-client-3.1.0/LICENSE) |
| `tslib` | 2.8.1 | 0BSD | [LICENSE.txt](THIRD_PARTY_LICENSES/npm-tslib-2.8.1/LICENSE.txt) |
| `vega-canvas` | 2.0.0 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-canvas-2.0.0/LICENSE) |
| `vega-crossfilter` | 5.1.3 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-crossfilter-5.1.3/LICENSE) |
| `vega-dataflow` | 6.1.3 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-dataflow-6.1.3/LICENSE) |
| `vega-embed` | 7.2.0 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-embed-7.2.0/LICENSE) |
| `vega-encode` | 5.2.2 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-encode-5.2.2/LICENSE) |
| `vega-event-selector` | 4.0.0 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-event-selector-4.0.0/LICENSE) |
| `vega-expression` | 6.1.0 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-expression-6.1.0/LICENSE) |
| `vega-expression` | 6.2.2 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-expression-6.2.2/LICENSE) |
| `vega-force` | 5.1.3 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-force-5.1.3/LICENSE) |
| `vega-format` | 2.1.3 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-format-2.1.3/LICENSE) |
| `vega-functions` | 6.2.0 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-functions-6.2.0/LICENSE) |
| `vega-geo` | 5.1.3 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-geo-5.1.3/LICENSE) |
| `vega-hierarchy` | 5.1.3 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-hierarchy-5.1.3/LICENSE) |
| `vega-interpreter` | 2.3.2 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-interpreter-2.3.2/LICENSE) |
| `vega-label` | 2.1.3 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-label-2.1.3/LICENSE) |
| `vega-lite` | 6.4.3 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-lite-6.4.3/LICENSE) |
| `vega-loader` | 5.1.3 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-loader-5.1.3/LICENSE) |
| `vega-parser` | 7.1.3 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-parser-7.1.3/LICENSE) |
| `vega-projection` | 2.1.3 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-projection-2.1.3/LICENSE) |
| `vega-regression` | 2.1.3 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-regression-2.1.3/LICENSE) |
| `vega-runtime` | 7.1.3 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-runtime-7.1.3/LICENSE) |
| `vega-scale` | 8.1.3 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-scale-8.1.3/LICENSE) |
| `vega-scenegraph` | 5.3.0 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-scenegraph-5.3.0/LICENSE) |
| `vega-schema-url-parser` | 3.0.2 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-schema-url-parser-3.0.2/LICENSE) |
| `vega-selections` | 6.1.5 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-selections-6.1.5/LICENSE) |
| `vega-statistics` | 2.0.0 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-statistics-2.0.0/LICENSE) |
| `vega-themes` | 3.0.0 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-themes-3.0.0/LICENSE) |
| `vega-time` | 3.3.0 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-time-3.3.0/LICENSE) |
| `vega-tooltip` | 1.1.0 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-tooltip-1.1.0/LICENSE) |
| `vega-transforms` | 5.2.2 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-transforms-5.2.2/LICENSE) |
| `vega-typings` | 2.3.0 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-typings-2.3.0/LICENSE) |
| `vega-util` | 2.1.3 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-util-2.1.3/LICENSE) |
| `vega-view-transforms` | 5.2.2 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-view-transforms-5.2.2/LICENSE) |
| `vega-view` | 6.2.0 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-view-6.2.0/LICENSE) |
| `vega-voronoi` | 5.1.3 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-voronoi-5.1.3/LICENSE) |
| `vega-wordcloud` | 5.1.3 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-wordcloud-5.1.3/LICENSE) |
| `vega` | 6.4.0 | BSD-3-Clause | [LICENSE](THIRD_PARTY_LICENSES/npm-vega-6.4.0/LICENSE) |
| `wrap-ansi` | 9.0.2 | MIT | [license](THIRD_PARTY_LICENSES/npm-wrap-ansi-9.0.2/license) |
| `y18n` | 5.0.8 | ISC | [LICENSE](THIRD_PARTY_LICENSES/npm-y18n-5.0.8/LICENSE) |
| `yargs-parser` | 22.0.0 | ISC | [LICENSE.txt](THIRD_PARTY_LICENSES/npm-yargs-parser-22.0.0/LICENSE.txt) |
| `yargs` | 18.0.0 | MIT | [LICENSE](THIRD_PARTY_LICENSES/npm-yargs-18.0.0/LICENSE) |

## Separately installed scientific dependencies

Python/R runtimes, Qt/PySide6, numerical libraries, WGCNA, FIMO and other external native methods are installed separately and retain their respective upstream licenses. Installing them does not relicense their code under Katana's MIT license. The repository does not bundle experimental datasets, large biological databases, raw reads, BAMs or complete external runtime installations. Resource acquisition and use must preserve the source resource's terms and attribution.

When redistributing a rebuilt browser bundle or a packaged application, retain the third-party notices corresponding to the actual included components. A future container or executable requires a fresh dependency and license inventory for what it distributes.
