# Completed-result reporting and runtime validation — 25 October 2026

**Completed: Windows completed-result assembly and its download archive, current
Windows repeated-refresh commands, and full Linux Conda/container portable
reporting passed their declared acceptance checks.** Installation, source
protection, native results, complete tables, links, browser behaviour and the
approved seven-section format were verified. The environments exercise
different declared entry points: Windows
assembles completed source collections, while standalone Linux Conda and the
Linux container consume the resulting portable report. This is not a claim
that three independent fresh source assemblies have completed.

This run reuses an existing study's completed analyses through the public
`katana-report` entry point. No sequencing dataset, network discovery,
MotifTrace model or MAPP fit was rerun. Report rendering and lossless packaging
are separate from scientific computation. Private study inputs, results and
machine-specific paths are not included in this repository.

The final summary is `katana.cross_environment_installed_report06_acceptance/v1`,
recorded at 02:40:17 UTC on 25 October 2026, with SHA256
`c03827ecd91561307309c8e31e6426a90960c4aea34c9200dfaace53a46b0bdf`.
The matching-runtime closure record is
`katana.final_matching_runtime_portable_report_acceptance/v1`, SHA256
`5fde7dfcf82e78380c3ebb413fffd3020549827384d9e9bf8be44582b4484de5`.

## Exact software baseline

| Item | Identity |
| --- | --- |
| Software commit | `d025eafb9591f859911e102c5d9bfe4ac3cd4778` |
| Python distribution version | `2.0.0b5.dev20261025` |
| Final wheel | `katana_splicing_tool-2.0.0b5.dev20261025-py3-none-any.whl` |
| Final wheel SHA256 | `9a86e98dad72f6bb163de76971afdb149f694ec91a43cf2f913c498c4bc53c41` |
| Wheel size | 3,955,889 bytes |
| Approved report specification | `1.1`, retaining the seven-section structure |
| Format approval SHA256 | `0f7bcae52fe53a9f4375556210803a1b44483bab58325ab162123a236f406ecb` |
| Matching final Linux image ID | `sha256:45286083debdabac5caad10ac2c8637cb2153f21120d6fbd7a494e4904b53725` |
| Passing CI revision | `d025eafb9591f859911e102c5d9bfe4ac3cd4778` |

This wheel was built from a clean committed Git archive. Uncommitted
documentation drafts were excluded. Of the 628 package/resource files, only
`katana_next/reports/refresh.py` differs from the preceding wheel; distribution
metadata also differs and is recorded separately. The version remains
`2.0.0b5.dev20261025`, so select artifacts by their full SHA256.

The accepted preceding Windows report and archive retain their original
software identity: revision `bf8732dc71b6208cdc2bce431630651f276d8094`, wheel
SHA256 `d387c1af245f57b25fedd0d323066f3c78ee6917e7201842db8efeac03dcf437`.
Its matched preceding image was
`sha256:75dd901c73e9bd91712be81b7cb5cbf8ef17f8b2e801ecde912e92a1084c4e11`.
Their receipts are not silently relabeled as executions of the new wheel. A
local image ID does not imply a published registry image or a PyPI release.

## GitHub CI

[Run 40150538654](https://github.com/Taratype0/Katana-splicing-tool/actions/runs/40150538654)
completed successfully for revision
`e1ebb827113ed2a20ef3f137eee8f7d285b3e263`. All four jobs passed: the
Ubuntu application checks, Windows application checks, browser interface, and
installed package-wheel checks. The final job finished at 01:36:09 UTC on
25 October 2026.

The only file changed between the accepted software revision `bf8732d` and
this CI revision is a synthetic test-fixture file. Packaged application,
scientific worker, dependency and browser-resource sources are unchanged.
That run belongs to the preceding software baseline. The current filesystem
correction at revision `d025eaf` also passed all four jobs in
[CI run 40151643701](https://github.com/Taratype0/Katana-splicing-tool/actions/runs/40151643701).
Its Windows application job finished last at 01:54:17 UTC on 25 October 2026.
A CI wheel is not silently substituted for an accepted local artifact. These CI checks complement the
runtime and real completed-result report checks below; they do not execute
every scientific workflow or HPC configuration.

## Observed status

| Environment or boundary | Installation/runtime result | Real completed-result report |
| --- | --- | --- |
| New Windows x86-64 venv, Python 3.13.9 | Passed: current 628-file exact match, isolated origins, `pip check`, three copied source-protection tests and four real public refresh commands with compaction on/off | Preceding baseline passed fresh source assembly, public portable refresh, native-data/links/format checks, inherited browser acceptance and complete ZIP validation; current wheel is not claimed as a fresh source assembly |
| New standalone Linux x86-64 Conda prefix, Python 3.12.15 | Passed: current install, `pip check`, exact 628 files, browser, five reader/mapping and three source-protection checks; original 51 native passes explicitly bound unchanged | Passed: full public CLI, complete source guard, independent native-data/table/link checks and browser acceptance inherited from exact offline-file equality with the actual container browser run |
| Matching Linux x86-64 Docker image, Python 3.12.15 | Passed: current image, exact 628 files, browser, five reader/mapping and three source-protection checks; original 51 native passes explicitly bound unchanged | Passed: full public CLI, complete source guard, independent native-data/table/link checks, actual 24-scope browser run and visual review |
| Actual Nextflow 25.10.4 → Docker image | Previous actual three completed and three cached probes retained as unchanged runtime/transport evidence; no new scientific or transport tasks executed for the application-only overlay | Runtime transport proof only |
| Read-only assembly of completed study sources | Passed: all 14 source identities, biological design and oriented contrasts checked before explicit target binding | Passed: Windows fresh assembly and final native-data/archive checks |

The retained six Windows CLI checks cover `katana-report`, `katana-flow`,
`katana-motiftrace`, `katana-targets`, `katana-next` and `katana-open`.
Windows tests used isolated Python mode outside the checkout, with no editable
install, source `PYTHONPATH`, system-site-packages or imports from the old
environment. Tests were copied to an unrelated directory. The 15 resource
checks cover installed report templates, browser assets, MotifTrace static
renderers and exported workflow sources. Three project/report/target execution
tests were deliberately excluded from this resource-only check. The old
environment's package set was unchanged.

The retained bounded-reader/path-map overlay passed ten copied tests; four
native Windows symlink cases were skipped for unavailable symlink privileges.
Linux symlink enforcement has separate evidence. These files are unchanged in
the current wheel. The reader retains native values, types, row order and
projections while limiting parallel readahead; mapping retains strict actual
file bounds without repeatedly resolving metadata paths.

An earlier Windows overlay added report-status reconciliation and upgrading
of copied completed network browsers; all 36 copied targeted tests passed. The
preceding display overlay changed only the network display-reuse predicate. Four copied
installed tests passed: current complete displays are reused, a missing required
display still needs native completion authentication, zero assigned modules do
not trigger a new drawing, and the historical bootstrap display keeps its
native edges. Six unrelated cases were deselected for this focused overlay.
All 628 installed package/resource files matched that preceding wheel; `Requires-Dist`
and frozen scientific/vendor bytes match the preceding wheel. The earlier 39
dependency and 15 resource passes remain linked to their original receipts,
with unchanged artifacts checked, rather than being relabeled as fresh runs.

One real completed comparison was upgraded with the preceding installed public copy
hook. Its native numerical payload and every module drawing matched the
approved reference display, including the full module gene table. The original
source and approved report bytes/timestamps were unchanged. Network discovery,
TOM, linkage, association and enrichment execution counts were zero. This
bounded copy check is retained by exact receipt SHA because its payload and
drawing code are unchanged in the filesystem correction. That correction performs
no network discovery or drawing export. This bounded check preceded the
separate complete report acceptance recorded below.

The accepted preceding Windows public `--from-report` command exited 0. All 11,693 offline
files, including 3,094 HTML files, were byte-identical to its immediately
preceding browser-reviewed report. The retained browser check covered all 24
network/RBP comparison scopes with screenshots. Browser acceptance is inherited
through complete file equality; the final report is not described as having
repeated those browser sessions. Separate final native-data, relative-link,
format and ZIP checks also passed. The archive contains 11,697 entries; every
entry passed CRC and SHA validation with no corruption.

| Accepted preceding Windows archive | Identity |
| --- | --- |
| Size | 4,557,562,038 bytes |
| SHA256 | `0aef94a63091d23fa4f2b7c679d08e024d1edf8882b5e4ad2f9dca7edc03722b` |
| Entries fully checked | 11,697 |
| Original sources | Bytes and timestamps preserved |

This archive identifies the checked private acceptance output, not a public
example dataset or an attachment stored in the repository.

## Repeated-refresh filesystem protection

A portable refresh may hardlink retained files when input and output share a
filesystem. The previous implementation copied an old
`presentation-refresh.json` this way, then overwrote that shared inode when
recording the new lifecycle. The correction writes the new receipt to a
temporary file and atomically replaces the destination directory entry.
Native result files and report rendering are unchanged.

Existing source control files, including a previous refresh receipt when
present, are bound by SHA256, size and modification timestamp before copying.
They are checked before recording completion and again after every output
write, including archive creation. Twenty-six targeted refresh, comparison,
compaction and RBP presentation tests passed. The new regressions use real
same-filesystem hardlinks with compaction enabled and disabled, check all
source bytes and timestamps, cover a first refresh without an older receipt,
and reject a source-receipt change occurring after the final archive write.
No scientific test, model, default or native result format changed.

The current installed Windows wheel passed three copied source-protection
regressions and four actual public CLI invocations. Each compaction mode ran a
first and repeated refresh, with real native-file hardlinks, distinct new
receipt inodes, unchanged source SHA256/size/mtime inventories and valid small
fixture ZIP CRCs. The current standalone Conda environment and image each
passed the same three source-protection boundaries against installed code.
These are focused runtime checks, not a second fresh assembly of the private
study or a rerun of its analyses.

An earlier report-rendering process exhausted host memory. Its stderr remains
part of the retained execution history. This was a rendering/resource failure,
not a failure of the completed scientific analyses. Earlier attempts are
retained separately from successful later commands.

A redundant standalone Conda fresh assembly was explicitly stopped before
completion. Its elapsed time was dominated by cross-filesystem metadata I/O;
the owned process identity, partial checkpoint, source read-only mounts and
confirmed stop are recorded. After the complete replacement passed validation,
the incomplete owned output was inventoried with all file hashes and retired;
its external logs and checkpoint remain. That attempt is not counted as a
successful report. A separate preceding-wheel portable attempt used a read-only
mount and fell back to redundant copies. That process was also explicitly
stopped, with exit code -15 and its partial output retained. It is distinct from
the retired fresh-assembly directory and is not counted as successful.

The two current-wheel public `--from-report --compact --no-archive` commands
started at 01:52:55 UTC. They consume a staged copy of the complete portable
report. Input and output share one owned writable staging mount so that
retained payloads can use hardlinks. The accepted Windows source and original
study directories are not mounted into the container. This is not a strict
read-only-input layout: the input file hashes, sizes and timestamps are bound
before execution, and their complete final preservation check remains a gate.
Neither replacement fits or repeats scientific analyses. Native-data, scope,
format and file checks all still apply; generating output files alone is not
recorded as full-report success.

The current container command finished at 02:33:41.755 UTC and the standalone
Conda command at 02:35:31.314 UTC, both with actual exit code 0. Their separate
full-source guards passed at 02:34:08 and 02:36:11 UTC respectively. A native
Windows verifier freshly hashed all 11,697 portable
payload files (12,209,662,836 logical bytes) in the source/staging input and
checked their complete file set, SHA256, size and modification timestamps.
Every payload passed. The previously accepted download ZIP was neither copied
nor rehashed in this portable run: its size and modification timestamp were
checked, and its earlier complete CRC/SHA acceptance remains the ZIP proof.
Inode change-time and link-count stability are not claimed because creating
hardlinks legitimately changes those fields.

Each output contains 11,697 files, with 11,687 output paths sharing retained
source inodes. Unshared output inode payloads total approximately 58.9 million
logical bytes per report; logical sizes are not a claim about physical disk allocation.
No new ZIP or scientific model was created. Independent native-data, full-table
and relative-link checks passed for both outputs. Actual container browser and
visual review passed; the standalone Conda report passed its own static checks
and exact offline-file binding to that browser-reviewed container output.

## Cross-platform references and browser acceptance

Strict Windows-to-container browser inheritance did not pass, and its failed
receipt is retained. The main `REPORT.html` contains 60 Top 10 result links
whose Windows backslash separators are serialized as forward slashes on Linux.
The other 11,692 offline files have strictly identical SHA256 hashes. Comparing
the main page after normalizing only `href`/`src` URL separators gives exact
HTML equality; this is a diagnostic comparison, not a rewrite of either report.
No native value, PNG, network payload or compressed table content differs.

The container therefore received an actual new browser run across all 24
comparison/DAS-selection scopes and visual review. Both passed; the Windows
main-page browser proof was not substituted for this run. The container browser
receipt has SHA256
`1802e8336696ade55c16fbf2df3b80d390660c9a7999bd40774e9f3bbdebef3f`.

The standalone Conda output has all 11,693 offline files, including 3,094 HTML
files, strictly byte-identical to the container output. Its successful browser
inheritance receipt, SHA256
`048309ab73495bdae576a12efebadad56bde14e3538cdd74adb2f89b005a1a90`,
binds the actual container browser run rather than the Windows run. The seven
sections, styling, content and analysis meaning are unchanged by URL
serialization. No new wheel, model fit or report-data calculation was needed
for this comparison.

## Verified report coverage

| Retained output | Accepted coverage |
| --- | --- |
| Selected comparisons | 12 |
| Top 10 records | 60 |
| Main figures | 554 |
| Independent network browsers | 12 |
| Native MAPP fit families | 168 |
| Lossless complete tables | 1,648 |
| Complete table rows | 14,142,605 |
| Original Katana files protected | 1,356 |
| Scientific analyses rerun | 0 |

These are acceptance-output counts, not published study findings. Native
statistics, gene/event identities, figures and complete result payloads remain
bound to their original evidence and comparison scopes.

The retained preceding Linux native checks exercise the installed WGCNA bridge, R namespaces,
unchanged bundled MAPP compilation, FIMO, PyHMMER, BAM/RSeQC handling and headless
graphics using small generated runtime inputs. They do not fit the real study.
The preceding Nextflow probes use the recorded image's interpreter and native tools with
networking disabled. Each rejects an existing mounted host interpreter rather
than using it as an unrecorded fallback. Their exact prior image identity is
`sha256:c56ce85b5aa01bca513d3f20f450fffc164ed1f3d875cdd8ee9a5690bcb7daea`;
its wheel was `d053dae0651b0ea7246f2aec0c15945df2a4ce4e924d8aaef303ae4668e81fe4`.
The current image has installed-browser and five-check bounded-reader/mapping
receipts, and all 628 package/resource files match the current wheel. The
preceding nine presentation checks remain explicit historical evidence; three
current installed source-protection checks passed separately. The current
`katana.unchanged_native_runtime_binding/v1` receipt, SHA256
`89c14e0f976afe49cc84b9eb8592e5c11342a19760c521590cb3b8a39478c354`,
binds both earlier 51-check proofs by their exact hashes to the current wheel.
Only nine reviewed report/presentation modules changed: native imports,
workers, model sources, checker and dependency bytes are unchanged. The new
filesystem-only wheel has its own byte binding and image identity. Neither
image is presented as having rerun those 51 checks or the three Nextflow
probes. Its portable-report execution is tracked separately from native runtime
installation and Windows source assembly.

## Locked Linux dependencies

The new Conda prefix was constructed from the existing exact dependency locks,
not copied from an installed prefix and not resolved against newer packages.
It contains 369 native distributions and the 59 locked supplementary pip
packages. Verified package archives were reused; the previous environment and
its cache were not modified. The current filesystem-only wheel was installed
and checked in this new prefix; its real portable-report acceptance remains
separate. Actual installed distribution metadata confirms that version and
`Requires-Dist` are unchanged, rather than assuming this from their filenames.
The explicit locks are unchanged by these report-code fixes.

| Tracked lock | SHA256 |
| --- | --- |
| `environments/locks/linux-64.explicit.txt` | `cc6b9d9c17b98fdcbf59b0dfffe77c20356b786830c7dc503ac2dc5f4b434b99` |
| `environments/locks/linux-64.sha256.explicit.txt` | `400b5800d191138128a27f6e2c7fdcec3b251eeb05dd6ff1445145c7793f62a2` |
| `environments/locks/pip-linux-64.lock` | `4bbfdd6b3510821fce032d52d42a52bef1313238fad82b7c2d831adda5cb87f1` |

These are Linux x86-64 locks. Windows used its own native-platform dependency
installation and verification; the Linux package URLs were not used to claim
Windows compatibility. Windows report generation from completed files does
not certify native Windows R/FIMO execution.

## Completed-result assembly and preservation

`katana.completed_result_set/v1` declares each input source explicitly:
completed artifact descriptors, verified collections or saved view descriptors.
Each source is checked against its original project manifest before any
presentation-only binding to the target report project. A changed manifest
revision requires an exact target SHA; samples, biological units, oriented
contrast arms and genome build must match. A saved view from a different
revision must supply its native `project_manifest`; its identity is not simply
overwritten to make it load.

Source files and native P values, FDR, model identities and comparison scopes
remain unchanged. Native MAPP Bonferroni values remain distinct from MotifTrace
FDR. A separate portable MAPP catalogue links every retained completed fit by
relative URL and SHA256 while preserving the original native catalogue bytes.
Missing or altered fits, scope mismatches and paths escaping their declared
root are rejected.

This round inherits recorded statistical outputs and filtering history. A
historical post-fit report filter remains described as post-fit; this report
acceptance does not claim a new pre-fit DEG/DAS exclusion or refitted model.
Biotype-based exclusions use the recorded annotation, not gene-name suffixes.
An annotated protein-coding retrogene is retained when its recorded biotype is
protein coding, even if its name resembles a pseudogene. No name-based filter
was added to manufacture a different result.

Fresh reports use the approved comparison selector and seven-section format.
Per-comparison network browsing exposes its retained network edges and complete
module/event/gene exports. MotifTrace sequence displays use native significance
and distinct-model ranking; a comparison without significant motifs receives a
short notice, while complete results remain in the output files. Lossless
compaction changes storage references, not numerical results. Its audit must
continue to validate every resulting reference and payload hash.

## Public CLI usage

Obtain or build the wheel identified above, then install it in a separate
environment using the [installation guide](../INSTALL.md). With that environment
activated, this command works in both a Linux shell and Windows PowerShell:

```text
python -m pip install "./katana_splicing_tool-2.0.0b5.dev20261025-py3-none-any.whl[report,motiftrace,targets,next]"
python -m pip check
katana-report --help
```

Create a `completed-results.json` descriptor for the completed study. For the
simple case of one existing completed artifact descriptor:

```json
{
  "schema": "katana.completed_result_set/v1",
  "project_manifest": "study/project.json",
  "sources": [
    {"source_id": "analysis", "kind": "descriptor", "path": "study/completed/report-project.json"}
  ]
}
```

When moving files between machines, declare the original path prefix and its
actual destination explicitly in `path-map.json`, for example:

```json
[
  {"source_prefix": "/source/study", "target_prefix": "./study"}
]
```

The target prefix is resolved relative to the mapping file. Use the prefix
actually recorded by the producer; a Windows source such as `X:/study` can be
mapped explicitly on either platform. Prefix boundaries and traversal checks
are enforced. Mapping is an in-memory consumer projection and does not edit the
producer manifests or their native hashes.

Generate into a new directory on either platform:

```text
katana-report --project completed-results.json --path-map path-map.json --outdir report-new --compact --no-archive
```

`--no-archive` keeps the complete portable report directory without creating a
second ZIP copy. To produce a download archive, use `--archive report-new.zip`
instead of `--no-archive`. A result-store cache can be explicitly selected with
`--result-store`; it does not replace the files needed to move the report.
Existing report directories are not overwritten.

For a complete portable report copied from another machine, use the saved
report entry point instead of repeating collection assembly:

```text
katana-report --from-report report-moved --outdir report-refreshed --compact --no-archive
```

This reads the report's relative assets and retained native records. It does
not require the source machine's original absolute paths and does not rerun
the analyses. The resulting report must still pass its own output validation.

## Evidence and completed acceptance

The retained local evidence is private where it contains study identifiers or
machine paths. This public summary records only software identities, check
outcomes and the evidence schemas:

| Evidence | Schema and required fields |
| --- | --- |
| Software build | `katana.installed_report_software_build/v1`: commit, wheel SHA/size, format version/approval SHA, packaged entry points/resources |
| Windows wheel | `katana.clean_windows_wheel_acceptance/v1`: exact wheel SHA, independent prefix, installed origins, CLI outcomes, no editable/source imports |
| Windows dependency and resource checks | `katana.clean_windows_requested_extras/v1`, `katana.clean_windows_installed_resources/v1`: applicable requirements, `pip check`, copied test SHAs and JUnit outcomes |
| Preceding Windows reader/path-map overlay | `katana.windows_wheel_overlay_acceptance/v1`: the preceding wheel SHA, 626 installed-file hashes, changed-file set, copied reader/path-map tests and inherited receipt hashes |
| Current Windows filesystem overlay | `katana.windows_installed_wheel06_refresh_acceptance/v1`: current 628 installed-file hashes, exact wheel identity, isolated origins, unchanged native/dependency metadata, three copied regressions and four real public repeated-refresh commands |
| Preceding network-reuse overlay | `katana.windows_wheel_overlay_acceptance/v1`: 628 preceding installed-file hashes, unchanged scientific/vendor bytes, four copied predicate tests, exact prior 36-test overlay, single-comparison and dependency receipt SHAs |
| Installed single-comparison copy | `katana.installed_single_comparison_network_display_acceptance/v1`: preceding wheel and installed module origin, approved payload/drawing equivalence, complete gene table, preserved original bytes/timestamps and zero discovery/test runs; later artifact bindings remain explicit |
| Current Linux application overlays | `katana.current_report_application_acceptance/v1`: current wheel/image identities, installed-browser and five reader/mapping check outcomes |
| Installed container presentation checks | `katana.installed_presentation_acceptance/v1`: nine exact display/status/preservation checks, probe SHA, no scientific models or original study data used |
| Current installed source protection | `katana.installed_refresh_source_protection_acceptance/v1`: real same-inode preservation with compaction on/off and expected rejection of a synthetic post-archive source mutation |
| Unchanged native proof binding | `katana.unchanged_native_runtime_binding/v1`: old/new wheel SHA, exact changed application files, byte-unchanged native code/dependencies and hashes of the two earlier 51-check proofs |
| New native prefix | `katana.new_locked_conda_native_prefix/v1`: native lock/archive SHA, package count, source environment unchanged |
| Final matching runtimes and portable reports | `katana.final_matching_runtime_portable_report_acceptance/v1`: current wheel/image, installed application/native bindings, separate actual public CLI/source guards and complete report evidence |
| Native checks | `katana.runtime_acceptance/v1`: 51 checks, zero failures, installed executables, generated output hashes |
| Actual Nextflow/image transport | `katana.actual_nextflow_container_acceptance/v1`: requested/actual image ID, completed/cached trace rows, native probes, preserved output hashes |
| Source assembly | `katana.completed_report_source_assembly/v1`: original/target identities, biological-design checks, preserved status history, zero scientific submissions |
| Installed report execution | `katana.installed_report_execution/v1`: exact public command, installed software identity, actual exit code, distinct entry point and separate final acceptance state |
| Full shared-staging execution | `katana.actual_shared_staging_portable_execution/v1`: installed wheel/image, actual public CLI and exit code, declared writable owned staging, original-source mounts excluded and archive/scientific execution flags |
| Complete portable source guard | `katana.actual_portable_full_source_guard/v1`: fresh payload SHA/size/mtime inventory, complete file set, independent input/original-source checks, retained archive identity and separate logical/inode storage accounting |
| Final read-only report validation | `katana.installed_cli_report_readonly_acceptance/v1`: format, native data, scoped views, relative links, source preservation and complete archive CRC/SHA checks with no report or scientific-result writes |
| Byte-equivalent browser inheritance | `katana.byte_equivalent_report_browser_inheritance/v1`: complete offline file hashes, source browser receipt SHA, matching HTML/assets and explicit statement that browser sessions were inherited rather than rerun |
| Cross-platform reference diagnostic | `katana.cross_os_report_html_link_serialization/v1`: retained strict-inheritance failure, exact differing link count, main-page native hashes, equality after URL-separator-only comparison and actual new-browser requirement |
| Final cross-environment acceptance | `katana.cross_environment_installed_report06_acceptance/v1`: exact software/runtime identities, separate Windows/source and Linux/portable entry points, full static/browser evidence, complete source guards and retained archive proof |
| Stopped redundant renderer | `katana.owned_local_report_supersession/v1`: exact owned process, preserved checkpoint, verified stop, read-only source mounts, incomplete attempt and replacement entry point |
| Retired incomplete output | `katana.owned_incomplete_report_retirement/v1`: bounded owned root, complete file-hash inventory, accepted replacement identity, original sources excluded and external history retained |

The declared acceptance checks are complete:

- [x] The filesystem-only application wheel is installed and checked in all
  three environments; the new container image and unchanged-native proof
  binding have their own exact identities.
- [x] Windows final native/archive checks pass after its successful public
  refresh. Browser acceptance is bound through complete offline-file equality.
- [x] Both current Linux full public CLIs exit 0 and their complete source
  payload guards pass; this does not substitute for final report-content checks.
- [x] Both Linux environments complete and validate the public portable entry
  point. No independent Linux fresh assembly is claimed.
- [x] Final manifest, provenance, figure data, relative links and complete
  result attachments pass validation on each platform.
- [x] Native statistics and original source identities remain unchanged after
  rendering and compaction; no source timestamp or byte mutation is accepted.
- [x] Seven-section specification `1.1`, comparison selection, network browser
  and MotifTrace significance display pass content/browser review.
- [x] Cross-environment report acceptance distinguishes expected rendering
  differences from numerical, identity or missing-file discrepancies.
- [x] Final portable directory/archive integrity and matched software/runtime
  identities are recorded with the Windows archive and current Linux output
  identities kept separate.

This record does not certify new FASTQ/BAM scientific runs, every optional data
mode or external database, other HPC scheduler configurations, ARM/macOS native
execution, an image-registry release or PyPI publication.
