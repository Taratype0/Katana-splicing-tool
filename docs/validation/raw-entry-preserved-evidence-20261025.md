# FASTQ/BAM entry: preserved evidence and current export binding

This audit verifies the installed workflow sources and retained results. It
does not turn an older scientific run into a new successful production run.
No completed upstream analysis, BAM audit or statistical model was repeated.

## Current installed workflow

The public `katana-flow export-pipeline` command was executed outside the source
checkout using the installed wheel for software revision `d025eaf`:

`9a86e98dad72f6bb163de76971afdb149f694ec91a43cf2f913c498c4bc53c41`.

All 116 tracked workflow files matched the wheel and the audited `a775c9d`
checkout byte for byte. This covers the main entry, upstream and downstream
workflows, process modules, configuration and helper scripts. The exported
source verifier separately passed all 89 preserved upstream files and 35
explicit integration changes. Twenty focused distribution, native-bundle
import and raw-entry report regressions passed.

This proves the installed/exported code identity and those tested contracts.
It is not an execution receipt for new alignments or differential models.

## Historical raw-entry checks

The [24 October integration record](pipeline-closure-20261024.md) retains actual
generated-input execution: 101 FASTQ upstream tasks and 82 BAM upstream tasks,
followed by 55 completed downstream tasks per entry. Each entry includes the
native 14-window MAPP family. These were generated software inputs, not study
data or a demonstration of arbitrary external BAM compatibility.

The original complete-mode controllers and their repaired downstream resumes
did not finish successfully. Once all 55 task receipts existed, collection was
recovered through a separate native Nextflow finalizer. Reporting was then
corrected and rendered separately. Their verified calculations and final
reports remain valid within those scopes, but must not be described as a
terminal-zero execution of the full main `--mode complete` command.

The missing collection resource mounting and non-finite report-value handling
were subsequently fixed. A new complete-mode first/resume execution needs its
own actual exit status, trace, source/runtime bindings and report acceptance;
plan generation or the sum of older successful pieces cannot substitute for it.

## Retained real BAM route

A retained study route originally completed its BAM identities, headers,
primary-record measurements and reference audits, then failed while exporting
the reused result bundle. A staged manifest was incorrectly treated as if its
parent directory contained the full source bundle. This was an export-adapter
failure, not a failed differential fit.

A later independent recovery preserved that failure and reused the successful
tasks and frozen payload. Both its controller and export job were verified as
native `COMPLETED/0:0`; the recovery trace records 74 cached tasks and one newly
completed export task. The resulting full bundle remains present, with all
1,498 ledger members checked for existence and declared byte count. Large BAMs
were neither downloaded nor rehashed for this audit. Historical full-payload
verification and the recovery's bounded adoption evidence remain distinct from
this fresh metadata check.

The recovered route's sample identities, contrast directions, reference
identities and rMATS vector ordering match its accepted FASTQ parent. Fresh
local SHA256 checks of 397 quantification/statistical/matrix files matched both
the recovered BAM manifest and the already-local parent results. The retained
original-reader projection had independently read 297 tables across all its
declared comparisons. Its receipt explicitly records a results projection,
without local BAM/BAI files.

These checks establish retained native-result identity and reader-projection
compatibility. They do **not** claim a newly executed full-BAM import through
the current package, a fresh current BAM scientific run, or support for BAMs
with incompatible references or missing required records. The public full-bundle
importer still checks the complete declared artifact ledger; this audit did
not weaken that requirement or replace absent BAMs with placeholders.

## Evidence separation

Keep separate records for source/export identity, generated transport, native
scientific completion, recovery/adoption, result import and portable reporting.
The [current installed-report acceptance](completed-results-runtime-20261025.md)
certifies existing-result rendering and its runtimes; it does not supply a raw
FASTQ/BAM execution receipt. The exact research manifests, scheduler records and
private data locations remain in the local acceptance dossier, not this public
software repository.
