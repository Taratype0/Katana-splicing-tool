# Katana CoSplice v1.2.3

## BAM-only external validation

- `rmats_root` is now optional when `matrix_input.input_mode` is
  `bam_junctions`.
- BAM junction discovery retains the same GTF annotation and competing-junction
  read-count semantics without inventing a placeholder rMATS directory.
- Optional rMATS outputs still add canonical event annotations when available.
- Matrix cache fingerprints and `parameters.json` record whether rMATS
  annotation was present.
- Added a fingerprinted `precomputed_junctions` backend for audited sample-level
  intron counts. It enters the same CoSplice clustering, PSI/IJC/SJC, stable-QC
  and uncertainty pipeline as BAM-derived counts.

## Disease-matched SETD2 validation

- Added a preregistered GSE297044 preparation workflow for 786-O SETD2 WT/KO
  DMSO paired-end RNA-seq with three replicates per group.
- Added resumable ENA FASTQ download, GENCODE v41 reference preparation, STAR
  alignment, BAM indexing, manifest/audit generation and input provenance.
- Added a Windows low-memory HISAT2 streaming path that writes junction counts
  directly and deletes transient uncompressed FASTQ; no giant SAM is retained.
- Added locked projection of TCGA KIRC splicing programs into GSE297044 using
  replicate/read-supported beta-binomial junction effects.
- Separated event-model measurement eligibility from network dynamic-range QC,
  preventing high-variance preselection before differential testing.
- Cross-platform outputs explicitly distinguish exact event matching from
  gene/program-level projection; unavailable coordinates remain unavailable.
- External event-model caches are parameter fingerprinted and invalidated when
  stable PSI/IJC/SJC or sample metadata changes.

## Validation

- Full repository regression suite passed before the external run.
- The GSE297044 manifest confirms balanced WT/KO batches, paired-end libraries
  and three primary replicates per condition.
- GEO supplementary featureCounts files are recorded as gene-assigned exon
  aggregates and excluded from junction uncertainty analysis.
- The independent GSE297044 run fitted 45,710 read-eligible beta-binomial
  models; 44,340 converged and 331 gene-annotated junctions passed global FDR
  0.05 plus absolute delta PSI 0.10.
- No frozen TCGA KIRC program passed external enrichment FDR 0.05. The best
  nominal CoSplice result had FDR 0.138 and was matched by event-level WGCNA.
  The release therefore makes no external-preservation or superiority claim.
