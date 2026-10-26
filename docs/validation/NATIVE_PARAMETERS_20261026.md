# Direct Nextflow parameter validation — 26 October 2026

Scope: `codex/performance-test-20261026`; not merged into `main` by this change.

## Implemented

- FASTQ/BAM complete mode generates internal project/settings JSON from public
  parameters. It still uses the existing upstream workflow and importer.
- Existing-project downstream mode can generate settings without `flow.json`.
- Resource files are Nextflow path inputs; local file identities, optional
  provider metadata and runtime descriptors are recorded by the generated task.
- Ten analysis skip switches, opt-in gene-body coverage/bootstrap, five
  molecular thresholds, event selection and two network module parameters.
- Layer-specific expression and usage cutoffs also reach the candidate evidence
  summary. Its version is 2.2.1; inherited defaults and default numerical output
  remain unchanged. Historical inventory JSON is not rewritten.
- Explicit advanced configuration and original-plan reuse remain supported.

## Observed checks

- 106 Python checks passed: direct parameter generation, invalid values,
  resource provenance conflicts, unchanged source inputs, genuine bundle
  import, generated FIMO descriptor transport, real planner parameter binding,
  network module override binding, optional analysis selection, complete YAML
  option coverage, usable application baseline and legacy candidate parity.
- Real Nextflow 25.10.4: seven transport tasks completed; equivalent values via
  `-params-file` plus `-resume` reused all seven. Changing the DEG cutoff rebuilt
  configuration and relevant transport, retaining three cached fixture tasks.
- The actual raw entry stopped at the configuration task for missing required
  resources, before scheduling read-processing jobs.
- The upstream source verifier passed: 89 source files and 37 previously
  declared integration changes. No upstream scientific module was changed.
- New documentation links resolve; the full native reference lists 114 options.

Native tests use explicit protocol fixtures, not biological calculations.
They test actual Nextflow execution and handoff, not full FASTQ/BAM scientific
acceptance, Slurm execution, image builds or performance on a real cohort.
The matching container build/installation qualification remains separate.

The approved seven-chapter report renderer and existing biological results
were not rerun or changed. New native tests are part of the GitHub workflow;
the remote CI result must be checked on the pushed commit.

Follow-up: an earlier commit's Windows CI exposed a same-clock-tick file
identity collision. Fresh Windows change timestamps now trigger content
verification rather than metadata-only reuse; no production sleep is added.
Only an already verified stable identity can use the shortcut. Older receipts
without this marker still load and receive content verification. Tests cover
stable reuse, rapid same-size edits and a fresh receipt becoming older.
The follow-up runtime, reuse, storage and execution-option suite passed 65
checks on Windows. The entrypoint commit's remote native Nextflow CI also
passed; the follow-up commit receives its own CI run.
