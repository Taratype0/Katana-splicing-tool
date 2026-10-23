# DEXSeq tool sources

The two Python tools are downloaded at build time from nf-core/rnasplice commit
`1d0494ae3402d1a46e0adadad24f81a0ff855c77`. Their complete source files and
original copyright/GPLv3 notices are included unchanged in the resulting image.
They are not included in the MIT workflow or Katana wheel. The download URLs
and exact SHA256 values are recorded in [sources.lock.json](sources.lock.json).

- `dexseq_count.py`: Simon Anders and Alejandro Reyes, EMBL, copyright 2010-2013.
- `dexseq_prepare_annotation.py`: Simon Anders, EMBL, copyright 2010.

The complete [GPL version 3 text](LICENSE.GPL-3.0) is supplied alongside them.
Its SHA256 is `3972dc9744f6499f0f9b2dbf76696f2ae7ad8af9b23dde66d6af86c9dfb36986`.
The base is the immutable BioContainers HTSeq 2.0.2 image named in the source
lock and Dockerfile. Its constituent software retains its own licenses.
Building this image does not publish a registry package or imply nf-core approval.
