# Methods and interpretation boundaries

Katana combines methods and evidence within a common project. It does not make different endpoints interchangeable, combine unrelated P values, or claim authorship of external statistical methods.

| Evidence | Interpretation and limits |
| --- | --- |
| Differential expression | Gene abundance with the imported or executed method's native effect and significance definitions. |
| DAS, DTU and exon usage | Event, transcript and exon-bin endpoints remain distinct. A feature-level FDR is not automatically a gene-level FDR. JC/JCEC and reference identities must remain explicit. |
| Candidate evidence | Support across eligible endpoints and contrasts is prioritization, not independent replication or causal proof. Unavailable evidence is not a negative finding. |
| MotifTrace | Sequence consequences, motif changes, usage support and RBP/CLIP evidence are separate layers. Predictions are not measured occupancy; unmatched public binding data do not establish binding in the user's experiment. |
| Co-splicing networks | The original WGCNA/TOM backbone and the declared cohort/event universe are retained. Modules group co-varying features; membership, association, stability, preservation and hubs have separate meanings. |
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
