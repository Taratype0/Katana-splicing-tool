# Pipeline integration and usability review

Reviewed against source `ab213c6` on 26 October 2026, with the documentation
and configurable network-minimum changes accompanying this review.
This is a code-path audit and targeted software check, not a new study run.

## Is the intended one-command workflow available?

The native complete FASTQ/BAM route and downstream project route exist.
Nextflow schedules their dependent tasks and report generation automatically.
The managed `katana-flow run` command adds staged cleanup and a retained project
for browser analysis. However, the intended turnkey experience is **not yet
complete**: new users still prepare project/resource configuration, external
nf-core directories need explicit adapters, full-report requirements are
rigid, and the latest source is ahead of downloadable runtime artifacts.

| Entry | Implemented behavior | Remaining gap |
| --- | --- | --- |
| Bulk FASTQ | Native upstream, import, downstream, report | One merged read set per sample; independent `~ condition` model; resource-template setup; current matching distribution still needed. |
| Genome BAM | Uses supplied alignments; quantification/differential/downstream/report | Same design/resource restrictions; genome/reference compatibility remains necessary. |
| Katana upstream bundle/project | Import and downstream without read-processing rerun | Setup is not yet a guided single-file experience. |
| External nf-core products | Mapped compatible results can be reused | No universal directory-to-project importer with automatic missing-step planning. |
| Completed downstream results | Report assembly and browsing | Report assembly alone does not compute missing analyses. |
| Selected genes/events | Public target commands, cards, Sashimi, portable results | BAM/reference locations are needed for new read-backed calculations. |

## Findings and priorities

| Priority | Finding and source | Consequence / next action |
| --- | --- | --- |
| P0 | `flow/network_inputs.py` hard-coded six contrast samples even though the standalone engine supports a configurable minimum. | Corrected in this change: `network.parameters.minimum_network_samples`, integer >=3, defaults preserved. Four samples can pass an explicitly configured network gate; this is not full two-plus-two statistical acceptance. |
| P0 | `flow/planner.py` always schedules network, coverage, boundary/CLIP and MAPP in the complete profile; preflight requires matching assets. | Define explicitly selected analysis modules and distinguish unavailable/not selected from failed. Keep the approved seven-section format. This policy is not changed by this review. |
| P0 | Current source is ahead of beta wheel/images. | Build and qualify one matching distribution before presenting current-source features as available from older downloads. |
| P1 | `pipeline/nextflow.config` falls back to one CPU, 8 GB and 12 h for downstream tasks; the planner does not populate task resource estimates. | Real studies need native site overrides. A copyable site example is added. Establish measured per-module defaults; no unmeasured larger defaults were silently substituted. |
| P1 | The public Slurm profile uses `params.slurm_partition`; no automatic `memory` partition routing was found in the generic workflow. | Keep queue/account policy in user site configuration. The prior study's actual partition decision cannot be reconstructed from generic source alone. |
| P1 | `katana_cosplice/network.py` and the documented WGCNA route use full event-by-event TOM. Python diagnostics and R fitting retain distinct calculations. | Memory grows as O(events squared). Avoidable allocations were already reduced, but replacing the method with true blockwise fitting needs separate numerical/algorithm validation. No new fit or scientific-method change was performed. |
| P1 | `flow/runtime.py` rehashes input/output members during receipt validation; collection and retained delivery validate again. | Repeated SHA256 reads and filesystem metadata operations can become I/O bottlenecks. Profile bytes/read counts; reuse trustworthy per-run identities at appropriate boundaries without accepting changed inputs. Not yet optimized in this change. |
| P1 | Collection/report materialization and archive generation retain additional copies; managed cleanup runs at upstream/downstream boundaries. | Current cleanup reduces retained work, not all peak disk use. Measure collection/archive peak storage and add lossless sharing where supported. Do not delete per-process work while downstream consumers still need it. |
| P1 | Managed resume binds source and configuration identity. | Changing a resource/site configuration cannot currently resume that same managed scope. Separate execution-only configuration from scientific identity with explicit compatibility checks. Not implemented here. |
| P1 | Retained projects have rebased browser/target resources, while the archived native collection preserves its original project provenance paths. | Existing offline reports and target reuse are verified separately. Qualify report regeneration from retained projects after old work deletion before documenting that as a universal refresh route. |
| P1 | Downstream `errorStrategy = 'finish'`, with no automatic retry. | Independent submitted tasks may finish, but a failed task stops complete delivery. Add narrowly scoped transient/resource retry policy with measurable bounds; do not retry invalid data indefinitely. |
| P1 | No `nextflow_schema.json` or nf-schema integration was found in the public pipeline, and its CI is custom integration testing rather than nf-core template/lint compliance. | Add a single discoverable parameter schema and standard help/validation, then run nf-core lint before a community submission. Existing native DAG execution is separate from this conformity gap. |
| P2 | Multi-lane raw input and shared paired/batch designs are not implemented in the complete raw route. | Add supported nf-core-compatible preprocessing/design adapters with independent numerical tests. Keep current limitations visible until done. |
| P2 | README and navigation mixed introductory usage with internal verification terminology and a Chinese quick start. | Rewritten English README and short English quick start; detailed methods/verification remain in linked reference pages. |

## Automation and nf-core

`pipeline/workflows/katana_downstream.nf` connects actual dependency channels.
Network fits, network followups, per-comparison motif tasks and MAPP window
tasks run as separate Nextflow processes. Collection waits for their results;
the report waits for collection. There is no four-task application cap in the
public Slurm profile, and routine progress does not depend on human polling.

Official nf-core membership is separate from using Nextflow correctly.
Reusing nf-core modules does not make the downstream code equally mature, nor
does it remove project-specific database setup, resource sizing or native
method constraints. Nextflow resource requests and partition selection remain
distinct configuration choices. See the official
[Nextflow process documentation](https://www.nextflow.io/docs/latest/process.html)
and [nf-core installation guide](https://nf-co.re/docs/usage/getting_started/configuration).

## Storage policy already implemented

The managed runner verifies a complete upstream handoff before deleting
upstream work; verifies the final report and independent reusable project
before deleting downstream work; and retains failures for recovery. Original
input files, shared images/environments, result tables, models and resources
for later target analysis are preserved. Native direct runs retain work unless
the operator selects a cleanup route. See [storage management](../STORAGE_MANAGEMENT.md).

Cleaning after success cannot prevent peak space exhaustion during a large
running stage. Trace disk use as well as final retained size. A scientific
matrix or statistical model is not disposable solely because it is large.

## Changes made during this review

- Configurable integrated network sample minimum, with validation and tests
  for a four-sample scope, unchanged defaults, plan/followup propagation and
  generated engine configuration.
- English-only public Markdown documentation, shorter README and a quick start
  for the three main input routes, browsing and target analysis.
- Native site-resource example with process selectors; no hard-coded real
  cluster, memory queue or application concurrency cap.
- No upstream study rerun, original-application edit, report-format edit,
  dense-network algorithm replacement or prior-results rewrite.

The code audit identifies optimization opportunities; it does not claim a
measured speedup, a new container build, or full real-data acceptance of every
small-sample configuration.

Validation for this change: 28 network-scope tests and 43 planning, reuse,
bootstrap-policy and report-policy tests passed. The site configuration parsed
with the fixed Nextflow 25.10.4 runtime. The 89-file upstream source check
passed, public documentation links resolved, and public documentation/example
text contained no Chinese. The approved real-study report's SHA256 was unchanged.
