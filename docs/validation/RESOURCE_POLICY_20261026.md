# Portable resource policy: performance test branch

This change separates tool resource requests from site-specific queue names.
It is part of `codex/performance-test-20261026`, not a production/main release.

## Changes

- Gene-body RSeQC: one CPU, 4 GB and 12 hours initially. Sequencing-QC
  collection, reference BED preparation and project audit: one CPU, 4 GB and
  two hours. Allocations can be overridden for a study/site.
- Downstream retry preserves CPU count; default memory/time may increase once.
  No queue change follows from a retry or a memory request.
- DEXSeq counting/annotation, reference indexing/preparation and FASTQ staging
  request one CPU for their current serial invocations. Existing memory/time
  requests remain in place. Imported upstream tool sources are unchanged.
- Generic examples no longer choose `main` or override tool resources through
  the site-config example. An optional Dardel resource example selects `shared`.
- Execution-only JSON supports exact process labels as well as process aliases.
  Process aliases take precedence. Empty overrides produce a valid no-op config.
- Scientific inputs, parameters, measurements, report layout and existing results
  are unchanged. No biological data were rerun for this validation.

## Validation

Local Python regressions: **45 passed, one skipped** across portable resources,
execution options and container-runtime boundary tests. The skipped test retains
its existing platform/dependency condition; it is not counted as passed.

Real Nextflow **25.10.4** on Linux executed ten generated shell tasks using the
actual pipeline configuration: five without a site queue and five with explicit
site mappings. All completed. The test confirms the RSeQC/default requests,
serial CPU override, qualified process-alias precedence and separate queues.
The task bodies record allocations; they do not run RSeQC, DEXSeq or network fits.

An independent native retry/resume test retained the simulated exit-137 failure,
then completed the retry with **one CPU and 32 MB**. Both completed tasks were
cached after a memory/time-only override. Upstream source verification passed
for all **89 source files and 37 declared integration changes**.

These checks establish configuration and recovery behavior. They do not certify
Slurm admission, performance on a new dataset, arbitrary scheduler deployments,
or the matching candidate container. Queue limits and allocation policy belong
to the site. No Dardel jobs were submitted by these checks.

The two native checks are included in GitHub's native integration workflow:

```sh
python pipeline/tools/check_resource_policy.py \
  --nextflow-command-file nextflow-command.json --outdir resource-policy-check
python pipeline/tools/check_resource_resume.py \
  --nextflow-command-file nextflow-command.json --outdir resource-resume-check
```

`nextflow-command.json` contains the verified runtime's command argv. Each output
directory must be new and receives logs and an `acceptance.receipt.json`.

See [execution configuration](../ANALYSIS_AND_EXECUTION_OPTIONS.md) for usage.
