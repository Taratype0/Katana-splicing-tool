# Selected contrasts and experimental designs

Katana follows the upstream nf-core contrast-sheet convention:

```csv
contrast,treatment,control
treated_vs_control,Treated,Control
```

Only declared rows are checked. A third condition does not create another comparison.
This layer supplements mature upstream workflows with experimental identities and
request-specific checks. It does not rewrite read processing or statistical methods,
fit a model, submit jobs, or grant execution authorization.

## Inputs

The experiment metadata sidecar has one row per canonical sample, with columns
`sample,condition,biological_unit_id,donor_id`. Optional columns include `batch`,
`cell_type`, `technical_replicate_of` and explicit `covariate_*` metadata. Lane rows
remain in the upstream samplesheet and are not independent biological replicates.
Blank biological identities stay unknown; no identity is inferred from a name.

Selected-contrast rows use `contrast,treatment,control`, with optional `design`,
`paired_by`, `paired` and `cell_type`. Treatment is the declared numerator and control
the denominator. Supplied `paired` values must be explicit boolean text:
`true`, `1`, `yes`, `false`, `0` or `no`. Pairing requires an explicit identity mapping and complete matching
within the selected comparison; a repeated donor does not authorize automatic pairing.
Technical replicas must be explicitly aggregated by the appropriate upstream route
before they can be used as one model observation. No counts or reads are merged here.

## Current support

| Route | Metadata that can be checked | Execution boundary |
| --- | --- | --- |
| `bulk_complete` | Independent selected two-condition comparisons, `~ condition` | Current complete upstream wrapper; paired, batch-adjusted and interaction requests are blocked |
| `deseq2_adapter` | Additive categorical `condition`, optional `batch` and/or `subject`, explicit complete pairs | DEG adapter exists, but extended designs are not wired into the complete route; successful preflight does not enable dispatch |
| `import_only` | Identity, selected comparisons and source design provenance | Preserves opaque source model text and native results as unverified; does not run or certify that source model |

At least two independent biological units per arm are required for compute preflight.
This is an implementation minimum, not proof of statistical power or independence.
Unknown units, technical duplication, incomplete pair mappings, confounded or
saturated designs, unsupported formulas and missing selected levels produce explicit
diagnostics. No simpler formula is silently substituted. Rank and residual degrees
of freedom are metadata diagnostics; the native R model and event-specific eligibility
must still pass their own checks at execution time.

Numeric/arbitrary covariates, nonlinear terms and interactions need a named compatible
adapter. Katana's adopted `differential.replicate_interaction` remains a **formal,
opt-in** 2x2 DEXSeq feature; its dedicated implementation and acceptance status are
separate from this thin contrast validator, which requires a dedicated adapter for
that request. A parsed or estimable design is not a
fitted result. Existing single-cell/nucleus results retain donor and cell-type context;
this validator does not turn a bulk wrapper into a single-cell raw-data pipeline.

## Preserve native direction and statistics

- DESeq2 reports native `log2(treatment/control)`.
- rMATS uses treatment as group 1 and control as group 2. JC/JCEC and sample-vector
  order remain explicit.
- SUPPA2 reports second minus first; the existing adapter places control first and
  treatment second, preserving that input-order provenance.
- **DEXSeq fold-change headers remain authoritative.** Current export records the
  control-over-treatment exception. Gene and stageR q-values have no effect direction.

Selecting a comparison never rewrites a native effect sign, P-value, adjusted P-value,
tested universe or testing family. Reversing an imported comparison requires an
explicit downstream display/transform contract, not silent mutation of its source.
For imports, the original declared row, including native pairing/model fields, is
retained as provenance. Method input-order assignments describe the pinned wrapper's
expected convention and remain unverified for an imported native run.
`selected_sample_ids` identifies the requested comparison;
it does not identify every sample used in the original statistical fit. Imported
`fit_sample_ids` stays null and `fit_scope` is `native_fit_unverified`: this validator
does not reconstruct or certify a native fit universe from comparison arms.

## Synthetic examples and CLI

The files under `examples/contracts/selected-contrasts-v1/` are synthetic metadata:

- `selected-subset.*.csv`: three available conditions, one requested contrast.
- `paired-adapter.*.csv`: three explicit donor pairs; checks the DEG adapter boundary.
- `confounded-batch.*.csv`: batch is identical to condition; the requested adjusted
  design must be rejected rather than run without adjustment.

Run from the development `source` folder with that folder on `PYTHONPATH`; these
examples do not require installing a package or console entry point. For the public
repository snapshot, the folder is `development/source`.

```powershell
$env:PYTHONPATH = (Get-Location).Path
python -B -m katana_next.contracts contrasts --metadata examples/contracts/selected-contrasts-v1/selected-subset.metadata.csv --contrasts examples/contracts/selected-contrasts-v1/selected-subset.contrasts.csv --route bulk_complete
python -B -m katana_next.contracts contrasts --metadata examples/contracts/selected-contrasts-v1/paired-adapter.metadata.csv --contrasts examples/contracts/selected-contrasts-v1/paired-adapter.contrasts.csv --route deseq2_adapter --design "~ subject + condition" --paired-by donor_id
python -B -m katana_next.contracts contrasts --metadata examples/contracts/selected-contrasts-v1/confounded-batch.metadata.csv --contrasts examples/contracts/selected-contrasts-v1/confounded-batch.contrasts.csv --route deseq2_adapter --design "~ batch + condition"
```

API entry points are `validate_contrast_rows` and `validate_contrastsheets` in
`katana_next.contracts`. Their route/design/pairing/intent/assay/cell-type inputs
are explicit. Validation produces diagnostics and preserves input rows; it does not
generate omitted comparisons. Passing metadata checks does not prove source-science
completion or package/runtime acceptance. The method support resource is versioned
as `katana.method_design_support/v1`, with relative source identities and SHA256
evidence for the current wrapper boundaries.
