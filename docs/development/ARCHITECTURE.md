# Development architecture

## Current software baseline

The repository contains the working Katana application source, browser interface, CLI/API consumers and scientific modules. MotifTrace, Katana-Splice-Network and the two report routes are retained. Its current registry contains 54 workflow entries, with input-dependent readiness and optional external tools.

The application still has shared GUI/application dependencies and source-checkout resource assumptions. Publication of this baseline does not mean that every module has already been extracted into a separately installable headless package.

## Target structure

```text
Bulk FASTQ / eligible genome BAM
              |
     Nextflow upstream integration
     (reuse mature nf-core methods)
              |
     versioned inputs and results
              |
       Katana analysis package
   core evidence | MotifTrace | networks
              |
       shared results and provenance
              |
      CLI / API / interactive UI
              |
        English reports / exports
```

Upstream generation and custom downstream analysis remain separate layers connected by explicit contracts. Independent legal task partitions can run in parallel; a joint statistical model and its full test family remain intact. Users select relevant contrasts rather than being required to compute every possible pair.

## Fixed decisions

- Adopted functionality is formal scope, with default versus optional activation recorded separately from migration or acceptance status.
- Historical launchers and experimental branches are not silently promoted into the default workflow.
- Project, sample, lane, donor, contrast, method, reference and result identities must be traceable.
- A report reads completed evidence. Import, analysis and reporting are separate intents.
- Missing inputs and unavailable methods yield explicit blocked/omitted states, not fabricated zero findings.
- Original installations, frozen algorithms and completed results are preserved. New work uses versioned development changes and fresh output identities.
- New report output is English. Existing report structure and substantive functionality are retained.

Planning documents belong in this development directory; the homepage describes the software and how to use it.
