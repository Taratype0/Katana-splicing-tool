# Gene exclusions and detected-gene enrichment backgrounds

Mitochondrial-encoded genes and pseudogenes are excluded by default for new
analyses. Identification uses chromosome and gene-biotype annotation, not gene
name guesses. Nuclear genes involved in mitochondrial biology are retained.
Missing required annotation is an actionable input error.

## Where filtering occurs

The upstream pipeline must apply the policy to the statistical inputs **before**
DEG, DAS, DTU and DEXSeq model fitting. Full alignment references and BAM files
remain available for mapping and read inspection. The independent upstream
implementation publishes `kanata.preanalysis_gene_exclusion/v1` provenance;
completed native jobs and manifest hashes are still required to verify a run.

The Python package also supports importing completed results. Filtering those
rows cannot retroactively change a fitted model or its multiple-testing family.
Such imports retain original effects, P values and FDR and explicitly record
`upstream_refit_performed: false`. They are never relabelled as models fitted
after pre-analysis exclusion.

`ProjectService.create_bulk_project(manifest, new_manifest_path)` writes explicit
ON defaults and validates annotation before creating a fresh analysis. Existing
manifests are opened with `open_project`; absent or explicitly saved historical
settings keep their prior meaning. The creation method cannot overwrite a
manifest or reuse a nonempty output directory. The gene-filter API and browser
panel recommend the same ON settings for a newly requested isolated scope.

## Foreground and background

GO/KEGG and RBP identity enrichment expose `foreground_source`:

- `candidate_union`: explicitly labelled supported candidates across assays.
- `das`, `deg`, `dtu`, `dexseq`: the corresponding native supported flags.

`background_source=match_foreground` is the default. Per-assay foregrounds use
the same assay's detected/testable universe; candidate union uses the explicitly
labelled union of measurable assays. A caller can explicitly choose another
listed background source; the actual choice is recorded in every result.

Eligible background rows require native detected/tested/available evidence.
Native measurement eligibility and invalid/failed model status are respected.
No extra count cutoff or significance cutoff is invented by the enrichment
adapter. **Nonsignificant detected/testable genes are retained.** A reference GTF,
genome-wide gene list, or list containing only significant hits is not a valid
detection universe. Missing assay evidence is not filled from another assay.

Numeric P values and `model_status` alone are not detection evidence. Gene
backgrounds also trace to complete local native result artifacts and their
registered input files, with actual content hashes, row counts and parent
lineage. An arbitrary hit table cannot certify itself by adding `detected=True`
or declaring `certified=True`. Explicit `result_scope=significant_only`, subset,
or foreground-only inputs are refused as backgrounds; they can be selected as
foregrounds against a separately available complete native Candidate source.

This verifies the files and scope read locally, **not an external guarantee that
a user never truncated an imported file before registering it**. The audit records
`delivery_completeness=not_independently_verified` and preserves the declared
source scope. Standard complete rMATS/DEG/DTU/DEXSeq imports remain supported via
the Candidate workflow; no extra raw-count catalog or refitting is required.
Missing historical lineage is an actionable request to import the complete native
results, rather than a silent switch to all reference genes or a hit-only list.

Rows remain separated by contrast, cell type and assay when those identities
exist. Event/transcript/exon identities map to unique canonical host genes; a
gene with multiple events counts once. The same exclusion and ID-mapping policy
applies to foreground and background. GO/KEGG then intersects the detected
background with genes mapped to the declared library. RBP identity tests retain
measured **RBP and non-RBP genes** in the background.

The hypergeometric test and BH adjustment are unchanged. Versions are bumped so
new source-contract runs cannot masquerade as old completed enrichment runs.
The historical helper used by frozen motif calculations is unchanged.

## Output and report policy

The output directory contains complete gene inputs, eligibility decisions,
mapping/resource audits and full enrichment tables. Rows include the effective
foreground/background source; the eligibility audit lists its evidence columns.
The report provides links to these inputs and outputs. Long gene lists, full
heatmap pages and complete matrices do not occupy the report body.

No completed study is recomputed or relabelled by this change. Missing native
pre-statistics provenance remains `not_verified`; annotation exclusion counts
are never substituted for detected or tested gene counts.
