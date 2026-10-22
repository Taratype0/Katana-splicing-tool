# Methods and interpretation boundaries

Katana combines methods and evidence within a common project. It does not make different endpoints interchangeable, combine unrelated P values, or claim authorship of external statistical methods.

See [bulk threshold values and evidence](THRESHOLDS.md) for the active `bulk-all-samples-counts-v1` policy and separate significance choices. The AS read gate requires finite `IJC + SJC >= 10` in every expected sample of both arms of the selected contrast, with unambiguous sample mapping and no pooling. Its per-sample JC precedent is distinguished from the explicitly declared JCEC extension; neither is a universal assay-independent validation.

| Evidence | Interpretation and limits |
| --- | --- |
| Differential expression | Gene abundance with the imported or executed method's native effect and significance definitions. |
| DAS, DTU and exon usage | Event, transcript and exon-bin endpoints remain distinct. AS read support and native significance are separate states; no replicate-consistency or approximate-SE gate is added. A feature-level FDR is not automatically a gene-level FDR. JC/JCEC and reference identities remain explicit. |
| Candidate evidence | Support across eligible endpoints and contrasts is prioritization, not independent replication or causal proof. Unavailable evidence is not a negative finding. |
| MotifTrace | Sequence consequences, motif changes, usage support and RBP/CLIP evidence are separate layers. Predictions are not measured occupancy; unmatched public binding data do not establish binding in the user's experiment. |
| Co-splicing networks | Count-bearing input uses the declared whole-cohort `network-all-samples-counts-v1` universe: every participating sample has valid PSI and IJC+SJC >=10, constant events are excluded, and DAS significance is not an input prefilter. Original discovery/merging algorithms remain separate from this eligibility change. PSI-only input cannot claim count-verified support. Modules group co-varying features; membership, association, stability, preservation and hubs have separate meanings. |
| DJExpress and other author-native routes | Preserve author implementations, their input requirements and method identity. Importing output is distinguished from executing the native method. |
| Single-cell routes | Supported processed inputs/results are qualified by their adapter. Cells, donors, pseudobulk samples and features are not interchangeable inferential units. Raw single-cell/nucleus FASTQ processing is outside this application's upstream promise. |
| Reports | Summarize observed completed evidence, parameters and limitations. They do not silently run absent analyses or transform exploratory results into validation. |

## Preserved algorithms

The source baseline retains MotifTrace, network and one-click report modules. Publication cleanup removes private paths and experiment-specific material from the public copy without replacing the scientific kernels. Workflow versions and relevant source identities remain separate from packaging and presentation versions.

R/WGCNA, DJExpress, rMATS, DEXSeq, DESeq2, MAPP, FIMO and other tools should be credited when actually used. Availability of an import adapter, an execution helper or a historical benchmark is not evidence that a new fit has been run for the current project. Optional integrations can require additional resources and environment validation.

## Scope

The planned upstream workflow targets bulk short-read FASTQ and eligible aligned BAM inputs by reusing mature nf-core modules and workflows. This repository's application baseline consumes supported inputs and results; the expanded upstream-to-package integration is separate engineering work and has not become an official nf-core release through publication of this source tree.

Single-cell and single-nucleus claims remain input- and method-specific. Processed result import, supported new downstream computation and raw-data preprocessing are recorded separately. Long-read result import does not imply a general long-read raw-data workflow.

## Reproducibility

Keep sample and contrast identities, biological units, reference versions, native result columns and completed-run provenance. Preserve complete statistical test families when parallelizing work. Record omitted, blocked, failed and unassessable tasks explicitly. Reusing completed output requires compatible inputs and methods; reopening a report alone does not authorize scientific recomputation.

## Current and historical AS support policies

The current all-sample read-count gate does not refit native statistical models or rewrite their P-values, FDR or effect estimates. In its compatibility columns, `stable_event` means current read support and `stable_significant` adds native significance/effect selection. The names are not claims of experimentally validated stability. A fixed sample-count gate does not establish adequate replication or power for a statistical design.

The prior presets are retained only in the [versioned historical archive](thresholds/history/8c368b7/README.md). New runs must record the active policy and cannot silently reuse a previous-policy Quality result as if it had passed the current gate. Selected-contrast Landscape matrices require support in every selected sample and eligibility of every included selected-contrast event row; an OR across contrasts cannot rescue a failing measurement.

## Whole-cohort network eligibility

The independently versioned [network input policy](thresholds/network-input.v1.json) applies to the entire sample set declared for each count-bearing network fit. It retains complete, count-supported, nonconstant PSI events without selecting significant DAS events. Every participating sample needs actual finite, nonnegative IJC/SJC summing to at least ten and a valid finite PSI; pooled or fraction-supported substitutes do not pass. The per-sample JC literature precedent and JCEC extension have the same limits described in the [threshold reference](THRESHOLDS.md#whole-cohort-network-input).

This gate acts on the resolved sample-level matrices. The existing rMATS adapter selects one canonical source for each condition/event type by largest file size, then lexical path/side order, and records that choice. It does not take median/maximum counts across comparisons or assemble partial counts from several tables. It also does not certify agreement with unselected candidate tables; a passing matrix observation cannot be cited as proof of cross-comparison source consistency. Source selection and raw measurement definitions are unchanged by this input-policy update.

Old sample-fraction, condition-coverage, minimum-two-replicate, positive variance-floor, repeatability and strong condition-response/dynamic-range gates no longer determine input eligibility under this policy. Diagnostic values may remain available. Count support is distinct from sufficient biological replication, reliable correlation, module membership confidence or significance of a module-condition association. External PSI-only matrices retain an explicit `read_support_not_assessed` state; their numerical values never imply observed read counts.

This policy changes the event universe used by new eligible runs; it does not claim that old network outputs already used that universe. New output/cache identities must distinguish it from historical defaults. Raw BAM-junction measurement parameters, original network kernels, module size/merging choices, and completed research fits remain separate and unchanged by the input approval.

CoSplice engine `1.3.0` and workflow `programs.cosplice` `3.2.0` implement the input policy and record it in resolved configuration and stage provenance. Historical outputs without that policy cannot be overwritten as current-policy outputs, including under force. Local implementation checks are distinct from a new native R fit or new biological validation; see the [versioned source bindings](thresholds/network-input.v1.json).
