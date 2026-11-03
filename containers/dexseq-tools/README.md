# Katana DEXSeq tools image

The FASTQ/BAM workflow automatically obtains the published Katana DEXSeq tools
Docker image or SIF. It contains HTSeq 2.0.2 and the unchanged upstream DEXSeq
annotation/count scripts. Statistical fitting uses the existing Bioconductor
DEXSeq image. No separate download or local image build is required.

- Docker identity: `j387735517/katana-dexseq-tools@sha256:20efe62e4679d6437c21938c624f7fdb18cc3c7fa8289456e1a73d9086b29717`
- SIF: `kanata-dexseq-tools-20261024-linux-amd64.sif`, retained in the
  `v2.0.0-beta.5.dev20261028` release; SHA256
  `ccab8d61d15efb22da72de26da5dde9f61915399ddc1b1494a6d3fba96fd4564`.

The existing image is reused without rebuilding. The pipeline stages the same
unchanged scripts explicitly, and cloud checks exercise both the image-bundled
and staged copies with generated single-end and paired-end BAM counts.
[Source and license notices](THIRD_PARTY.md) preserve their GPL-3.0 attribution.

For developers rebuilding the unchanged tool environment, use
`python scripts/build_dexseq_container.py --context /tmp/katana-dexseq-build --tag katana-dexseq-tools:local`
in a new directory on a Linux Docker host. Normal pipeline users do not need
this command. `--dexseq_tools_container` remains an optional compatible override.
