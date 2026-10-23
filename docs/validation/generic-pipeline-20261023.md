# Generic existing-result workflow validation — 23 October 2026

**Status: passed for the generated existing-result entry, native downstream
workflow, collection and fixed seven-section report.** This is a development
acceptance record, not a new FASTQ/BAM upstream or release-image certification.

This check runs the actual Nextflow existing-result entry point and registered
Katana tasks on Linux. It uses generated software fixtures rather than a public
or research dataset. The source inputs include nine generated biological-sample
identities, two comparisons with six samples each, twelve distinct motif models,
generated reference sequences, BAMs and upstream-format result tables. Existing
research inputs and completed study results are not used or recomputed.

## Execution scope

| Check | Scope and completion rule |
| --- | --- |
| Native downstream graph | A 98-task planned analysis graph under Nextflow 25.10.4, followed by collection and the complete seven-section report gate; final native trace and exit status must confirm completion |
| Per-comparison networks | Separate selected sample sets, fits, bootstrap, memberships, activities, associations and report figures; no pooled network substituted for independent fits |
| MotifTrace | Independent tasks for both comparisons; multiple model identities, sequence displays, original statistics, Top 10 ranking and complete output tables |
| QC and imported differential layers | Real generated BAM files support native coverage; upstream-format QC/DEG/DTU/DEXSeq observations are marked synthetic imports, not newly executed upstream analyses |
| MAPP | The fixture has no eligible paired transcript inputs: 28 preparation tasks and 28 dependent fits are verified `not_applicable`; no native fit success is claimed |
| Reporting | All required computed evidence and per-comparison figures must pass collection and report acceptance before a complete archive is delivered |

Fixture-only settings use 10 network bootstrap iterations and 99 motif-cluster
permutations to bound the software check. They do not replace production
defaults or establish empirical biological stability. Generated CLIP-shaped
records and input P/FDR values test input/transport behavior and are not evidence
of experimental binding or differential expression.

## Observed results

- The sealed 98-task plan contains 42 completed calculations/adapters and the
  56 verified MAPP not-applicable tasks above. The final native Nextflow run
  exited zero with 100 cached processes and one completed report process.
  None of the completed calculations was rerun for report corrections.
- Both comparisons independently used six samples and 128 events. Their eight
  fit/bootstrap/association/annotation tasks passed. The exported network
  browsers retain all members and native edges, PSI and module activity.
- Each comparison has 12 completed MotifTrace models and 10 sequence cards.
  Complete families contain 179,184 and 113,280 tests. Each displayed P, FDR
  and effect was checked against the same native test row. Model rankings do
  not relabel test-family FDR as a model-level FDR.
- The complete report and ZIP passed their required-analysis, seven-section,
  figure-data, link and archive checks. There are no main figures missing
  exact data bindings. The 114 generated input files remain byte-identical.
- RBP identity and boundary-motif enrichment use compact per-comparison
  summaries. Only native FDR ≤ 0.05 rows appear in these summaries; hit genes
  expand in place. All tested sets and motif hit-gene/event joins are exported.
  Repeated support graphics remain in outputs. Native bootstrap heatmaps and
  size plots link their actual plotted matrix/summary, not an unrelated table.
- The Windows application/report regression passed 895 tests, with three
  skips; browser compatibility passed 15 tests. Follow-up report/data-binding,
  plan-reuse and failure-handling checks passed after the final corrections.

## Completed failure-handling checks

An independent shell-only check imported the actual production Nextflow config.
One task exited with code 23 while two independent tasks were running. Nextflow
returned nonzero; both independent tasks completed with exit code zero and the
dependent report remained blocked. An explicit repaired resume reused both
successful tasks from cache, preserving output hashes, modification times and
work directories. Only the failed branch was rerun before the dependent report
started. There were no automatic retries or ignored failures.

The final report wrapper has also been tested with actual rendered ZIP files.
Required-analysis gaps, missing chapter figures, missing final artifacts and
acceptance errors retain a diagnostic ZIP and failure receipt, with no successful
report receipt. Existing failed or successful report directories are preserved.

## Distribution and limits

The [earlier runtime validation](downstream-runtime-20261023.md) identifies a
specific installed wheel and image. New workflow corrections are not present in
that old archive. A current wheel candidate was built from an isolated source
snapshot and installed outside the checkout. It passed 41 portable checks,
seven newly added module identity checks, `validate-plan` and presentation
recovery CLI checks, exact pipeline export comparison and real browser HTTP
requests. Served JavaScript/CSS matched the wheel bytes. Its SHA256 is
`b5b4facd446d4e827079c982ee79ccd6434861a16e7a7cd90a2c82c233ac3d3f`.
This check used a new environment with the previously accepted native
dependencies available through `system-site-packages` and installed the wheel
with `--no-deps`. It validates current wheel distribution and installed entry
points; it is not another independent clean dependency installation or a
scientific rerun. The existing runtime environment was not modified.
The installed wheel also generated and passed the complete report from the
verified collection outside the source checkout, using Python isolated mode.
This final check rendered existing results only and executed no scientific task.

The later [current-container acceptance](current-container-20261023.md) records
the matching development image, rebuilt wheel and disconnected installation
checks. That image also rendered this existing collection through the complete
report acceptance gate without rerunning calculations. The older wheel hash
above remains the identity of this original acceptance step.

This check does not add FASTQ/BAM upstream scientific acceptance, execute public
data, resume a paused research study, validate every optional resource/model
combination, or establish remote-scheduler/platform coverage. Native MAPP fit
coverage is not newly certified by a fixture with zero eligible inputs. The
paused research study remains separate from this software acceptance.
