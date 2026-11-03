# DEXSeq command-line scripts

These unmodified scripts come from nf-core/rnasplice 1.0.4, commit
`1d0494ae3402d1a46e0adadad24f81a0ff855c77`. They retain their original
copyright notices and **GPL-3.0** license; they are not relicensed under Katana's
MIT license. See `LICENSE.GPL-3.0` and `sources.lock.json` here.

Nextflow stages the scripts as explicit task inputs into the upstream HTSeq
2.0.2 BioContainer. No Katana-specific DEXSeq container build is needed. The
scripts and HTSeq version are identical to the former extra image. Statistical
DEXSeq fitting continues in its existing upstream Bioconductor container.
