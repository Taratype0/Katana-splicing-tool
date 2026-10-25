# Upstream DEXSeq tools image

The [downloadable release](../../docs/RELEASE_INSTALL.md#linux-nextflow-with-docker)
includes the accepted `kanata-dexseq-tools:dev-20261024` Docker archive. Load it
and verify its image ID against the release manifest, or build locally below.

The raw FASTQ/BAM route needs DEXSeq annotation and exon counting in addition to
the Katana downstream image. The previous GHCR address was an unpublished
placeholder. Use the verified release archive or an explicitly built image;
downstream-only analysis of existing
results does not need this image.

From a checkout on a Linux x86-64 Docker host:

```sh
python scripts/build_dexseq_container.py \
  --context /tmp/kanata-dexseq-build \
  --tag kanata-dexseq-tools:0.1.1
```

Use a new build directory outside the checkout. The helper downloads the two
unchanged upstream scripts at a fixed commit, verifies their SHA256 values,
builds on an immutable HTSeq image and runs actual single-end and paired-end
GTF/BAM counting fixtures offline as an unprivileged user. A successful
`acceptance.receipt.json` records the image ID, native versions and exact counts.
There is no registry publication and no study data in the image.

Pass the image explicitly to the raw-input workflow:

```sh
nextflow run pipeline -profile docker \
  --mode complete --dexseq_tools_container kanata-dexseq-tools:0.1.1 \
  --katana_container "$KATANA_IMAGE" \
  --input /data/study/samples.csv --contrasts /data/study/contrasts.csv \
  --fasta /references/genome.fa --gtf /references/annotation.gtf \
  --project_template /data/study/project-template.json \
  --flow_config /data/study/flow.json --outdir /data/study/results
```

`--mode upstream-only` requires the same DEXSeq image option with container
profiles. For Apptainer/Singularity, supply a site-built and verified SIF at a
path visible to the compute nodes; Docker acceptance alone does not certify
that conversion or the site's Slurm deployment.

[Source and license notices](THIRD_PARTY.md) describe the separate GPL tools.
Their sources are downloaded only into the external build context. The
workflow never substitutes new scripts or statistical parameters for them.
