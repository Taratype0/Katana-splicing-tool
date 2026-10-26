# Performance candidate

This candidate is developed on `codex/performance-test-20261026`. It is not a
replacement for `main` until its validation is complete. It builds on development
commit `04941d6e9805f48134955275ac9a20ba06b09af6`.

## Changes

1. **MotifTrace scan reuse.** Reuse scores for the same sequence and complete
   model definition within one analysis. A 32 MiB memoization cache evicts old
   entries when full; it never removes events, models, windows or tests.
2. **MotifTrace matrix streaming.** Generate and process one paired matrix at a
   time. All comparisons retain their original within-comparison statistical
   families: the same pair swaps, random seed, global maximum-cluster null and
   final BH adjustment. Final result tables and the permutation signs still
   require memory; this is not a claim of constant-memory execution.
3. **Network memory overlap.** Execute the existing R WGCNA module detection
   before allocating Python correlation, adjacency and TOM matrices. Both
   calculations and their parameters are retained. This reduces overlapping
   allocations; it does not replace WGCNA or change module membership.
4. **MAPP preparation.** Compute transcript membership, measured-event
   eligibility and TPM sums once per comparison, then reuse them in its 14
   windows. Window-specific sequence eligibility and native fitting remain
   independent.
5. **Independent molecular branches.** Read-support processing, AS, DEG and
   DTU processing are separate tasks per comparison. Evidence membership,
   landscape and cross-comparison summaries join the completed branches where
   they genuinely need several comparisons. They retain pairwise results.
6. **Task-specific cache identity.** Workers receive their own task and ancestor
   definitions instead of the complete plan. Unrelated module selections do not
   invalidate those workers. Changed task parameters invalidate that task and
   its descendants. Planning explicitly binds configuration and project-manifest
   contents, rather than relying on their directory timestamps.
7. **Report delivery.** Verified, pipeline-owned payloads can be hardlinked on
   the same filesystem; external user inputs are still copied. Report ZIP
   creation happens after the calculation and figure checks. Already-compressed
   payloads are stored without another compression pass. The approved report
   layout and file contents are unchanged.

## Validation and limits

`tests/public_baseline/test_performance_regression.py` compares generated inputs
with byte-preserved sources from the development baseline. It checks all
positions, missing values, raw P, BH q, cluster FWER, permutation nulls, network
outputs, all 14 MAPP preparation windows and the complete molecular tables after
joining comparison branches. Historical reference sources remain unchanged.

Observed instrumentation: one consensus sequence is scanned once instead of
eight times across the tested methods/windows. Fourteen MAPP windows reuse one
event-quantity preparation. These are counts of eliminated duplicate work, not
claims of an eightfold or fourteenfold complete-pipeline speedup.

`pipeline/tools/check_native_cache_scope.py` uses actual Nextflow execution with
generated transport fixtures. An unrelated configuration change preserves all
12 worker cache hits. Changing one task parameter reruns exactly that task and
its two descendants; the other nine workers remain cached. The check is included
in GitHub CI. The existing native transport checks cover task fan-out, dependency
ordering, failed-task propagation and unchanged-run resume.

No Primary culture or public biological analysis is rerun by these tests. No
new large-cohort wall-time or peak-RSS benchmark is claimed. Matching container
execution and installed-wheel checks remain distinct from source-code tests.

## Trying and resuming the candidate

Use the candidate branch (preferably its exact tested commit) and a matching
package/container. Do not mix a new Nextflow checkout with an older installed
Katana image. The normal commands and analysis configuration remain applicable.

For a native downstream run with retained work, use `nextflow run pipeline/main.nf
--mode downstream-only --project project.json --flow_config flow.json -resume`
with the same launch directory, work directory and input paths. Choose a new
output directory when creating a revised report. If work was already cleaned,
use the retained result project; Nextflow cannot recover deleted intermediates.

The managed `katana-flow run --resume` command still resumes the same scientific
request, while allowing execution-resource changes. A changed scientific request
uses a new managed run directory. This candidate does not silently reinterpret an
existing managed run as a different study or configuration.

Switching from a legacy aggregate plan to this task layout may create a new cache
identity. It is not permission to rerun completed studies: existing results can
continue through the completed-result entry point. Rollback uses the previous
commit and its existing result directories; no original code or results are
modified in place.
