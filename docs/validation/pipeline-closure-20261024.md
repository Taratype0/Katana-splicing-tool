# Pipeline integration and installation — 24 October 2026

This is a development acceptance record, not a production release or an
nf-core approval. The original application and completed study results were
preserved. Primary culture remains paused until the outstanding entry checks
and deployment checks pass.

## Accepted downstream behavior

The [generated two-comparison workflow](generic-pipeline-20261023.md) verifies:

- Separate network discovery, bootstrap, association and annotation tasks for
  each selected comparison, using only its biological samples.
- Multiple MotifTrace models per comparison, sequence Top 10 with the matching
  native P/FDR and complete model/test tables in outputs.
- Compact significant RBP and boundary-motif enrichment tables with complete
  recorded hit genes and the full tested results available separately.
- The approved seven-section report, independent comparison identities, plot
  data links and a verified downloadable result package.

These are executed software-fixture checks. They do not certify a real-study
report, six-sample network precision, or an absent analysis resource.

## Current accepted images

Software revision: `086338842b241ba6bdef137b05a783b18eb6ca19` on
`codex/package-development`. Later documentation-only commits retain these
identities. Images have not been published to a registry.

| Item | Identity |
| --- | --- |
| Wheel | `katana_splicing_tool-2.0.0b5.dev20261023-py3-none-any.whl` |
| Wheel SHA256 | `e12d7c493ef5ac83b9af0263fd9394952b4f486c16fa5aa6ec5132b49fd7725e` |
| Downstream image | `katana-downstream:dev-0863388` |
| Downstream image ID | `sha256:32cd5a67bb67311d0c665d3df3b4db3ce7a88f7637d6216c589af8fd2ad6e57c` |
| Downstream archive | `katana-downstream-dev-0863388-linux-amd64.tar.gz`, 1,463,502,232 bytes |
| Downstream archive SHA256 | `089a8a940a60360a7488b8339b600503d895a9b69d2840b98dce0d44d09558ae` |
| DEXSeq tools image | `kanata-dexseq-tools:dev-20261024` |
| DEXSeq tools image ID | `sha256:9da6802a176f77a81c732cbee5a5a9f8ca4a64ee38391f35903e1fa57881c209` |
| DEXSeq tools archive | `kanata-dexseq-tools-dev-20261024-linux-amd64.tar.gz`, 122,383,650 bytes |
| DEXSeq tools archive SHA256 | `9eca391e9009fb7ff4ad7479e882d5f73ea170532ed79d4f8f29de3ad27f2fe7` |
| Platform | Linux x86-64 / Docker |

The downstream native/Python locks are unchanged from the
[earlier runtime acceptance](downstream-runtime-20261023.md). The current wheel
was installed into a new image; this was not a new dependency solver run.

Observed checks on the current downstream image:

- 51 native runtime checks passed offline as UID 1000, without host Python/R.
- 523 installed package/resource files matched the wheel.
- Browser health and six HTTP routes/assets passed.
- Three native Nextflow container probes completed, then all three were cached
  on resume with unchanged outputs.
- Report rendering verified all 98 records in the already accepted collection
  without repeating its calculations. This includes 42 completed calculations
  or adapters and 56 explicitly not-applicable MAPP tasks, not 98 successful fits.
- Both exported image archives were hashed, loaded and matched their image IDs.

The separate DEXSeq image downloaded hash-pinned, unchanged GPL tool scripts
at build time. Actual generated annotation and single/paired BAM count checks
passed offline as UID 1000. Source license and attribution are retained in that
image. See the [DEXSeq image instructions](../../containers/dexseq-tools/README.md).

## Additional integration checks and fixes

Reference declarations, annotation identity and gene-exclusion provenance now
survive bundle export/import. Older bundles can recover these fields only from
their hash-verified source map; conflicting declarations fail validation.
Unavailable metadata remains unavailable. Saved FASTQ/BAM bundles (212 and 199
files respectively) passed export/import checks without rerunning their native
calculations or changing source files.

Upstream BioContainers image names now explicitly use Quay. The formerly
assumed DEXSeq registry image was not published: containerized raw-data entry
now requires `--dexseq_tools_container` and fails before task submission if it
is missing. Versions and scientific scripts were not changed by these fixes.

A separate generated positive MAPP check executed the official sequence
counter and native MAEI fit/significance chain for 120 exact paired SE events,
six samples and one 3-prime window. The counter produced 5,440 kmer columns;
this does not mean every model was fitted or significant. The fixture used an
explicit minimum motif fraction of 0.5 without changing the production default.
This positive branch check is additional to the not-applicable branches above;
it is not acceptance of a new study's full 14-window family.

[Native Nextflow integration](https://github.com/Taratype0/Katana-splicing-tool/actions/runs/40038202713)
and [Linux/Windows, interface and installed-wheel CI](https://github.com/Taratype0/Katana-splicing-tool/actions/runs/40038202761).
Both runs passed for the image revision above. Subsequent fixes are separately
validated before updating this image acceptance record.

## Outstanding acceptance

The raw-entry check exposed an input/output basename collision in reference
preparation: an input named `annotation.gtf` could be truncated through its
staged symlink, and a chained empty-file check did not terminate the task.
Inputs are now staged in separate directories and each required output has an
independent fatal check. A real Nextflow/gffread container regression passed
for same-named plain inputs, compressed inputs and rejected empty transcript
output, with source bytes and modification times preserved in all three cases.
The failed generated fixture and logs were retained. Its annotation was
reconstructed in a new directory and matched the original pre-failure SHA256.
No study inputs were involved. The current image includes this correction.

Both generated raw-entry upstream runs have completed: FASTQ finished
101 tasks (94 cached and seven newly completed in its final upstream resume);
BAM finished 82 tasks, including six BAM-record FastQC tasks. The 256-file
FASTQ bundle and 207-file BAM bundle passed their full SHA256 inventories.
These are upstream results; the complete downstream/report runs remain
under acceptance.

Native execution also exposed and corrected a SUPPA merge loop on
header-only event classes. Completed source IOE files were preserved. A
single selected comparison now explicitly marks cross-comparison intersections
not applicable while retaining chapter 4; it does not fabricate a self-comparison.

Native DEXSeq exon CSV files now enter the exon layer through a string-preserving
TSV projection. Per-gene screen tables remain a separate asset. Completed
stageR confirmation tables take priority for DTU; the underlying transcript
test is retained separately. Actual native bundle checks verified 384 exon bins
and 256 distinct transcript rows with no duplicate units or changed P/FDR.
All source statistical outputs remain unchanged.

- Complete the current raw FASTQ entry, then the BAM entry, through real native
  tools and the integrated downstream/report handoff. A generated local input
  fixture is being used; no public or study dataset is being rerun. The initial
  Docker credential-helper failure is retained, and an isolated Docker client
  configuration is used for recovery. A running task or image download is not
  a completed analysis.
- Check the actual target site's container engine, storage headroom, resource
  paths and retained task/result inventory before resuming the paused study.
  Local Docker acceptance does not certify a Dardel Apptainer deployment.
- Primary culture's missing analyses and its complete final report remain
  unaccepted. They are not covered by software-fixture success.

Native ARM/macOS and every scheduler/site combination are outside this accepted
platform scope. Database files and references remain separately declared inputs.
No PyPI or registry publication is implied by these checks.
