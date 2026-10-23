FROM quay.io/biocontainers/htseq@sha256:ada151195080772235556fa8678c5e213e4e98b10b8c0e1c64ec9209db730222
USER root
COPY gpl-scripts/dexseq_prepare_annotation.py /usr/local/bin/dexseq_prepare_annotation.py
COPY gpl-scripts/dexseq_count.py /usr/local/bin/dexseq_count.py
COPY LICENSE.GPL-3.0 THIRD_PARTY.md sources.lock.json check_runtime.py /usr/share/doc/kanata-dexseq-tools/
RUN chmod 755 /usr/local/bin/dexseq_prepare_annotation.py /usr/local/bin/dexseq_count.py
LABEL org.opencontainers.image.title="kanata DEXSeq Python tools" \
      org.opencontainers.image.licenses="GPL-3.0-only" \
      org.opencontainers.image.version="0.1.1" \
      org.opencontainers.image.source="https://github.com/nf-core/rnasplice/tree/1d0494ae3402d1a46e0adadad24f81a0ff855c77" \
      org.opencontainers.image.description="Unmodified GPL DEXSeq Python scripts from nf-core/rnasplice 1.0.4 with HTSeq 2.0.2" \
      io.kanata.dexseq.upstream-revision="1d0494ae3402d1a46e0adadad24f81a0ff855c77"
