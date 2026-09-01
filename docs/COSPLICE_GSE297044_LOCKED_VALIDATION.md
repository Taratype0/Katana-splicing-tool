# GSE297044 locked external validation

## Design

GSE297044 is an independent human clear-cell renal-cell-carcinoma RNA-seq
study. The primary validation contrast is fixed before alignment as 786-O
SETD2 KO DMSO versus matched 786-O SETD2 WT DMSO, with three paired-end
replicates per group. A1331852 treatment and A498 SETD2 add-back samples are
secondary sensitivity analyses and cannot change primary thresholds.

The GEO supplementary featureCounts tables aggregate exons by gene name. They
are not splice-junction IJC/SJC evidence and are excluded from the CoSplice
network. The primary run therefore starts from ENA FASTQ, aligns against
GRCh38, and derives competing-junction usage plus read-count uncertainty from
either coordinate-sorted BAM or audited aligner-stream junction counts. The
Windows validation uses the official GRCh38 HISAT2 genome index, `--add-chrname`
for GENCODE-compatible coordinates, and never writes a large intermediate SAM.

## Locked endpoint

Frozen TCGA KIRC modules are projected into GSE297044 without selecting or
retuning modules on the validation phenotype. The main test asks whether genes
from each frozen KIRC program are enriched for GSE297044 junctions passing the
replicate-aware beta-binomial condition test and the predeclared absolute
delta-PSI cutoff. Event-model eligibility requires at least 10 competing-
junction reads in all three replicates of both conditions. It does not require
a minimum observed delta PSI, variance or dynamic range before model fitting.

TCGA SpliceSeq event IDs in the local panel lack genomic coordinates that can
be matched rigorously to STAR junctions. Results are therefore labeled
gene/program-level replication unless an exact coordinate map is available.
They must not be described as exact event preservation or event-direction
concordance.

## Reproduction

```powershell
python scripts/prepare_gse297044_external_validation.py `
  --source-root 'E:\Project\RNASEQ\Public_validation\GSE297044_SETD2_ccRCC' `
  --work-root 'E:\Project\RNASEQ\Public_validation\GSE297044_SETD2_ccRCC_locked_validation_v1' `
  --stage all --subset primary --threads 12 --download-workers 6

# Windows low-memory fallback after FASTQ/reference preparation:
python scripts/run_gse297044_hisat2_junctions.py `
  --work-root 'E:\Project\RNASEQ\Public_validation\GSE297044_SETD2_ccRCC_locked_validation_v1' `
  --workers 2 --threads-per-sample 6

# Rewrites the config to select the completed precomputed-junction backend.
python scripts/prepare_gse297044_external_validation.py `
  --source-root 'E:\Project\RNASEQ\Public_validation\GSE297044_SETD2_ccRCC' `
  --work-root 'E:\Project\RNASEQ\Public_validation\GSE297044_SETD2_ccRCC_locked_validation_v1' `
  --stage config

python -m katana_cosplice.cli stage --name matrix `
  --config 'E:\Project\RNASEQ\Public_validation\GSE297044_SETD2_ccRCC_locked_validation_v1\cosplice_config.yaml'

python -m katana_cosplice.cli stage --name stable_qc `
  --config 'E:\Project\RNASEQ\Public_validation\GSE297044_SETD2_ccRCC_locked_validation_v1\cosplice_config.yaml'

python scripts/validate_gse297044_setd2_programs.py `
  --locked-root 'E:\Project\RNASEQ\Public_validation\TCGA_production_release_v10_paired' `
  --cosplice-output 'E:\Project\RNASEQ\Public_validation\GSE297044_SETD2_ccRCC_locked_validation_v1\katana_cosplice_output' `
  --output 'E:\Project\RNASEQ\Public_validation\GSE297044_SETD2_ccRCC_locked_validation_v1\locked_program_validation'
```

The validator fits the preregistered beta-binomial event test directly from
the full PSI/IJC/SJC matrices after replicate/read measurement QC. The network
`stable_as` output remains a separate discovery diagnostic and is not used to
preselect differential events. The validator does not discover local modules
in the six validation samples.

`precomputed_junctions` means sample-level genomic introns and read counts from
the locked aligner/filter settings. It is not a precomputed PSI matrix: CoSplice
still builds competing-junction clusters, PSI, IJC/SJC, stable-event QC and
measurement uncertainty itself.

## Observed result

The independent event model fitted 45,710 replicate/read-eligible junctions;
44,340 models converged. Among 39,437 gene-annotated tested junctions, 331
passed the global beta-binomial FDR cutoff of 0.05, absolute delta PSI of 0.10,
and convergence requirement, representing 251 genes.

No frozen TCGA KIRC program was enriched at FDR 0.05. The best nominal result
was RI C01 for CoSplice (hypergeometric p = 0.01315, FDR = 0.138); the matched
event-level WGCNA module had the same result. This run validates an independent
SETD2-responsive junction signature, but it does not validate preservation of
the frozen KIRC programs and does not demonstrate a CoSplice advantage over
event-level WGCNA.
