# Matched-subject comparisons

The raw FASTQ/BAM entry now carries an explicit paired design from the
contrastsheet through native differential analysis and downstream project
metadata. `paired=true` requires one observation per subject in each arm.
The independent `~ condition` route retains its existing models.

| Analysis | Paired behavior |
| --- | --- |
| DEG | Native DESeq2 `~ subject + condition`, treatment minus control |
| rMATS | Native PAIRADISE via `--paired-stats`; both BAM lists aligned by subject |
| SUPPA | Native `-pa`; PSI and TPM columns follow the same subject order |
| DEXSeq exon usage / DTU | Full `~ sample + exon + subject:exon + condition:exon`; reduced `~ sample + exon + subject:exon` |
| Network | Observed-PSI discovery and continuous kME unchanged; subject-blocked module association with complete pairs |
| Report | Existing format; model identities and actual parameters travel with result metadata |

DEXSeq's native `estimateExonFoldChanges` remains a condition-only descriptive
summary. It is not a new subject-adjusted effect estimator. The paired P/FDR
come from the subject-adjusted likelihood-ratio test. Both models and this
distinction are recorded in the fitted object metadata.

A six-observation paired network has three independent subjects. The minimum
for the paired module-association adapter is three complete pairs; it does not
treat six observations as six independent animals. A batch column alone does
not request batch adjustment. Incomplete pairs, duplicated within-arm subjects,
contradictory design flags and undeclared shared subjects are rejected.

## Validation

- Python regression checks cover subject ordering, incomplete/duplicated pairs,
  flag conflicts, downstream scope identity and module-association P values
  against an independent paired t-test.
- `pipeline/tools/check_paired_transport.py` executes the actual rMATS and SUPPA
  Nextflow modules with generated command doubles: four native tasks verify
  paired/unpaired flag dispatch. This is transport evidence, not model fitting.
- `pipeline/tests/paired_models.R` fits generated counts using the production
  DEG, DTU and exon adapters and checks sample-order invariance against a native
  reference fit. Dardel job 27023715 completed with exit 0 on 28 October 2026,
  using DEXSeq 1.52.0 / DESeq2 1.46.0. No study data were used.
- CI additionally runs that fixture in the pinned upstream DEXSeq 1.36.0 image.
  Its result must be checked on the exact candidate commit; adding the test is
  not evidence that this second runtime has passed.
- Matching application-image deployment and complete study execution are
  separate checks. No four-study FASTQ run has been submitted by these tests.

## Native method references

- [rMATS paired statistics and input ordering](https://github.com/Xinglab/rmats-turbo)
- [SUPPA paired differential splicing](https://github.com/comprna/SUPPA)
- [DEXSeq alternative exon-usage models](https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html)

## Updating an HPC application image

When dependency locks match a verified release, `scripts/update_downstream_sif.py`
can install a new wheel into a separate writable Singularity sandbox and build
a new immutable SIF. The input archive is checksum-checked, and the old image
is unchanged. The build receipt identifies both the base and the new wheel.
Run `scripts/check_downstream_environment.py` in the new image before use;
`status=built` alone does not establish installation acceptance. A changed
dependency lock requires a full environment build instead.
