# Nextflow `nf-core` DSL2 Integration Kit for KatanaAS

Contains 4 standalone `nf-core` DSL2 modules (`modules/local/katana_*`), 2 wired end-to-end subworkflows (`subworkflows/local/katana_bulk_as` and `katana_smartseq3_as`), and a standalone runnable `main.nf` + `nextflow.config`.

## 1. Run Directly with Nextflow (`-profile conda | docker | singularity`)

```bash
# Bulk RNA-seq end-to-end (SpliceNet -> SpliceTrace + Module Bridge + GTF/UniProt Proteome & 5 Plots)
nextflow run main.nf -profile conda \
  --mode bulk --species mouse \
  --usage psi_matrix.tsv --samples samples.tsv --support counts.tsv --annotation event2gene.tsv \
  --rmats ./rmats_out --fasta GRCm39.fa --gtf gencode.vM34.annotation.gtf \
  --outdir ./results_bulk

# Single-Cell Smart-seq3 end-to-end (SC-SpliceNet -> SC-SpliceTrace + Module Bridge + GTF/UniProt & 5 Plots)
nextflow run main.nf -profile conda \
  --mode smartseq3 --species human \
  --events EVENTS.tsv --psi PSI.tsv.gz --support SUPPORT.tsv.gz --metadata CELLS.tsv \
  --fasta GRCh38.fa --gtf gencode.v46.annotation.gtf.gz \
  --numerator Stim --denominator Ctrl --outdir ./results_sc
```

## 2. Embed Subworkflows into Any Existing `nf-core` Pipeline (e.g. `nf-core/rnasplice`)

Copy `modules/local/katana_*` and `subworkflows/local/katana_*` into your pipeline repository and include:
```groovy
include { KATANA_BULK_AS      } from './subworkflows/local/katana_bulk_as/main'
include { KATANA_SMARTSEQ3_AS } from './subworkflows/local/katana_smartseq3_as/main'
```
