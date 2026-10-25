# Entry-point and runtime validation — 25 October 2026

**The recorded bulk FASTQ and BAM complete-mode executions, unchanged cache
repeats, installed target restore, and matching runtime checks passed.** These
are generated software cases. The separately accepted real-study report and
retained BAM recovery keep their original evidence; they are not newly
recomputed here. See [entry points and commands](../ENTRY_POINTS.md).

## Software and runtime

| Item | Accepted identity |
| --- | --- |
| Current installed/display revision | `098ceaea3ed56e3927318d579dc4d6e725494b92` |
| Distribution | `katana_splicing_tool-2.0.0b5.dev20261025-py3-none-any.whl` |
| Wheel SHA256 | `b9eda865201363edf77209bf1469eb6bbf7ad7fc8b9e390530c1232a54c46890` |
| Wheel bytes | 3,966,687 |
| Matching Linux x86-64 image | `sha256:619b6ace665b686c2c5703098a8ccc1989f424060901cb2c9bfe2d7b1f7d6222` |
| Local image tag | `katana-downstream:acceptance-20261025-098ceae-presentation09` |
| Native dependency lock SHA256 | `cc6b9d9c17b98fdcbf59b0dfffe77c20356b786830c7dc503ac2dc5f4b434b99` |
| Supplementary pip lock SHA256 | `4bbfdd6b3510821fce032d52d42a52bef1313238fad82b7c2d831adda5cb87f1` |
| Installed package/resource identities checked per environment | 630 |
| Installed/exported workflow files checked | 117 |
| Preserved upstream source files / explicit integration records checked | 89 / 36 |

The current wheel was installed and checked on Windows, standalone Linux Conda
and the Linux Docker image. `pip check`, dependency-version preservation and
public target restore/plan/run reuse passed in all three. The installed workflow
export matched the wheel and committed source. The current overlay changes only
the three report modules described below; scientific kernels, workflow sources,
target workers, CSS/JavaScript and dependencies are unchanged. Original native
first-use checks remain bound by exact bytes, image layers and locks; this is
not a claim that all 51 native probes were rerun for each overlay.

The actual whole-workflow runs retain this preceding execution identity:

| Whole-workflow baseline | Original identity, retained |
| --- | --- |
| Source | `c448007faa838586086d5f482ab7e5a632de8ae8` |
| Wheel SHA256 | `36996a2dd5ce4b19fd97bee383ad755c03f5b4f50d5586c1a309d2628f02bc3e` |
| Image | `sha256:9dfe1d3396f77319bb9b377c5f6136df876e7963cb84468fabf967da71d36796` |

Both [application CI](https://github.com/Taratype0/Katana-splicing-tool/actions/runs/40181808030)
(Ubuntu, Windows, browser and installed wheel) and
[native Nextflow CI](https://github.com/Taratype0/Katana-splicing-tool/actions/runs/40181808070)
passed at the `c448007` whole-workflow baseline. CI transport cases and the
actual native executions below are distinct evidence.

The current `098ceae` revision also passed all four
[application CI jobs](https://github.com/Taratype0/Katana-splicing-tool/actions/runs/40184505728).
Its 117 exported workflow files are unchanged from `c448007`; the native
Nextflow CI record remains the actual `c448007` run, not a newly triggered
`098ceae` native run.

## Complete bulk workflow

The actual main `--mode complete` entry at `c448007` performs resource preflight, upstream
processing, native bundle import, the downstream task DAG, collection and
report acceptance under Nextflow 25.10.4 and Docker. The test study has six
generated paired-end bulk samples, declared read length 75, and one explicit
comparison. No public sequencing dataset or existing study analysis was rerun.

| Entry | Observed successful execution | Unchanged repeat | Protected products |
| --- | --- | --- | --- |
| FASTQ | Main complete-mode recovery: 158 completed tasks reused from the preserved attempt, two new terminal tasks (`COLLECT`, `REPORT`); exit 0 | 160/160 `CACHED`, exit 0, no new task | 6,885 product entries retained their SHA256, size and mtime |
| BAM | Fresh main complete-mode run from the six generated BAMs: 141 `COMPLETED`, zero failed, exit 0 | 141/141 `CACHED`, exit 0, no new task | 6,185 product entries unchanged; 257 source FASTQ-bundle files also unchanged |

The FASTQ result is a verified recovery of the whole command, followed by an
entirely cached repeat; it is not advertised as a cold run of all tasks under
the final revision. Earlier failed collection and interrupted diagnostic
attempts remain recorded. Recovery corrected the display check for the native
network writer's exact six-significant-digit export and deterministic DEXSeq
version-file selection. It did not edit scientific outputs, cache records or
receipts, change algorithms, or repeat the 158 successful tasks.

Each final collection verifies 54 completed downstream tasks, including all 14
native MAPP fit families. The original `c448007` fixed seven-section reports
passed the following structural, numerical/catalogue and offline browser checks;
the later content correction is recorded separately below:

| Report check | FASTQ | BAM |
| --- | --- | --- |
| Figures, including output-only figures | 75 | 73 |
| Figures in the main report manifest | 34 | 32 |
| Complete-result entries | 270 | 230 |
| HTML files | 187 | 177 |
| Relative links checked | 1,578 | 1,464 |
| ZIP members checked by CRC and SHA256 | 652 | 608 |
| Required-analysis gaps | 0 | 0 |

Browser checks exercised ordinary DAS/RS-DAS switching, the module selector,
network drawing, complete gene downloads, and offline searchable tables. The
network contains 128 events across six samples, two assigned modules and the
unassigned set. A display-limited graph is not presented as the complete
network topology; full member/gene tables remain available.

All 45 shared native FASTQ/BAM result tables were byte-identical in this case,
including quantification, DEG, DTU/DEXSeq, PSI/read-count matrices and the ten
rMATS event-type × JC/JCEC exports. Geometry, sample vectors, effect direction,
P/FDR and references were checked. The generated case contains nonempty SE
events; the other four event types have valid empty exports, so this is not a
nonempty five-event-type biological benchmark.

### Generated-case scientific limits

The small acceptance configuration explicitly uses 99 motif permutations and
MAPP exon-set fraction 0.5 for this software case. These are not new production
defaults. Its twelve MotifTrace models are authentically **not estimable**;
the report records that state and shows no significant sequence cards. It
does not call missing/invalid FDR "nonsignificant" or invent ten motifs.
Multiple-model Top 10/significance behaviour is covered by separate regression
and previously accepted real-result report checks. A single comparison also
has no fabricated cross-comparison overlap or interaction fit.

This checks execution, data transport, completion, numerical consistency and
report behaviour. It does not establish statistical power or biological
validity for an arbitrary experiment.

## Presentation-only follow-up

A subsequent manual content check found two presentation defects despite the
successful structural checks above: zero-count pie labels could overlap, and
the completed-source loader could merge different enrichment foregrounds or a
statistics table with its same-name plot specification. A custom GMT family
without recognized GO/KEGG labels could also receive an inaccurate generic
"not computed" message. The original statistical files were present; this was
an indexing/display defect, not an absent fit.

Revision `098ceaea3ed56e3927318d579dc4d6e725494b92` corrects those three report
modules. It retains the four native foreground-specific statistics tables and
their four plot specifications, hides only zero-count slice labels, and
describes the actual empty-view state. The same deduplication correction also
restores the saved PSI measurement table and three already approved chapter-2
views: PCA, within-condition consistency and outliers. Those views derive from
completed PSI; this is result-index restoration, not merely a caption edit.
The seven-section layout, CSS/JavaScript, native values, scientific kernels,
task pipeline and target workers are unchanged. Source checks passed 56 focused
tests and verified the original collection files without mutation.

The installed public `katana-flow report` command then rendered the existing
FASTQ collection in Conda and the existing BAM collection in Docker into new
directories. Both exited 0 (14.3 and 10.2 seconds respectively), retaining all
1,370/1,324 source collection files' SHA256, size and mtime. No whole-workflow,
network, differential or MAPP task was resubmitted. Independent static and
offline browser checks passed both reports, including all 54 task identities,
14 MAPP families, ordinary/RS-DAS switching, network members/genes, the restored
PSI views and the four original enrichment statistics tables. The four custom
GMT plots remain output-only and their native statistics are byte-identical.

| Corrected report check | FASTQ collection / Conda | BAM collection / Docker |
| --- | --- | --- |
| Figures / main-report figures | 82 / 37 | 80 / 35 |
| Complete-result entries | 300 | 260 |
| ZIP members checked | 726 | 682 |

The final installed/display acceptance record,
`katana.installed_presentation_accuracy_release09_acceptance/v1`, has SHA256
`62775dca160e2a909daf041eaad320899e994169e24c4e68368bf532a725dd33`.
The independent report-content/browser record has SHA256
`ab93c8c058674de8bb65c6f583f5b1d8e1680526b75b86db85dee1193363bca6`.

The updated reports have their own rendering/acceptance records. They do not
relabel the `c448007` whole-workflow and 160/141 cached-repeat evidence as
executions of a later version. No original study report was regenerated.

## Gene/Event targets and returned-result restoration

The [target commands](../TARGET_COMMANDS.md) support browser selection, manual
request inputs, execution, returned-product import and restoration into a new
workspace. Windows native-product checks and a new Linux first-use case
verified Gene/Event cards and five actual indexed-BAM Sashimi products. The
Linux first-use case used four tiny generated indexed BAMs, two targets,
21 HTML files, 35 PNG images and 217 relative links. No study data was used.

Restore preserves all 202 recorded product members, including the five Sashimi
products. The matching Windows, Conda and Docker installed commands restored
the products, validated the plan and reused both targets. Native statistics,
resources, source archives and cached products retained SHA256, size and mtime.
An occupied destination was rejected. Windows-produced results were readable
on Linux; differing method/dependency fingerprints were not forced into cache
hits. Exact current-runtime reuse and cross-environment readability are
different claims.

Requests do not embed BAMs or assume a particular cluster. New Sashimi needs
verified BAM/index bindings on the execution machine. Restore itself performs
no counting, fitting or rendering.

## Retained completed-result and real-study evidence

- The [completed-result report record](completed-results-runtime-20261025.md)
  retains real-study assembly, full archive and Windows/Linux portable-report
  checks under their original artifact identities and fixed report format.
- Current production terminal handoffs separately reused each historical
  sealed 55-task FASTQ/BAM plan: two terminal tasks completed, then both cached
  on resume, with no scientific task submitted and protected inputs unchanged.
  These terminal checks are not substituted for the main complete-mode runs.
- The [preserved raw-entry record](raw-entry-preserved-evidence-20261025.md)
  retains the earlier generated 101/82 upstream and 55 downstream task scopes.
  It also records a real BAM export recovery: native controller/export jobs
  `26791376` and `26791383` completed with exit `0:0`, with 74 cached tasks and
  one new export. Fresh bounded identity checks matched 397 native result
  files to the accepted FASTQ parent. This is not a new full real-study BAM fit.

## Deployment and release boundary

The accepted native platform is Linux x86-64 with Docker or the locked Conda
runtime; Windows supplies the Python/browser/report route and verified target
products. Apptainer/Singularity and Slurm profiles are available, but each site's
mounts, scheduler settings and image conversion still need qualification.
Native macOS/ARM and arbitrary read-layout or external BAM/reference combinations
are not certified by these checks. Raw scRNA/snRNA processing is outside the bulk
workflow. External nf-core output directories require a matched adapter and
manifest, not merely a directory name.

The pipeline is an independent Katana DAG reusing nf-core modules; it does not
run three whole unmodified nf-core pipelines. No approved nf-core membership,
PyPI release or public container-registry image is claimed. Original Katana,
approved study reports, references and completed study statistics were not
modified by this acceptance.
