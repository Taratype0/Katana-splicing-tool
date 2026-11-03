# Automatic analysis databases

For `--species Mus_musculus --reference_build GRCm39`, the default
`--resource_bundle auto` selects `mouse-grcm39-v1`. Only resources needed by
enabled modules download. The database cache occupies about 158 MB after
extraction; initial compressed transfers are about 19 MB. Container images,
FASTA/GTF, indexes and scientific results are separate from these sizes.

| Resource | Version and use | Source |
| --- | --- | --- |
| GO MGI mouse | Frozen 2026-09-09 derived GMT; annotations propagated through is_a/part_of | [Gene Ontology](https://www.geneontology.org/docs/download-go-annotations/) |
| KEGG mouse | KEGG 2019 Mouse via Enrichr; unique case-insensitive MGI symbol mapping. This is a historical library, not current KEGG. | [Enrichr library](https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse) |
| RBP identity | Frozen mouse RBP/splicing-factor membership, not regulator target genes | GO annotations and [ATtRACT](https://attract.cnic.es/) |
| Motif models | Pinned ATtRACT archive, filtered by species by the existing analysis | [ATtRACT download](https://attract.cnic.es/download) |
| NOVA1 / NOVA2 CLIP | GSE69710, CD1 E18.5 cortex; author peaks, mm9 converted to GRCm39 | [GEO GSE69710](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69710) |
| PTBP1 CLIP | GSM1828887 / GSM1828888, published mouse reference biosamples; converted to GRCm39 | [GSM1828887](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1828887), [GSM1828888](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1828888) |

These CLIP data are external reference binding evidence, not measurements of
your experimental animals. Separate target, biosample, assay, original build,
conversion and source identities are retained in the generated project. Supply
your own matched CLIP with `--clip_peaks` when available, or use
`--skip_motif_binding` to omit this analysis.

The machine-readable catalog is
[pipeline/assets/resources/mouse-grcm39-v1.json](../pipeline/assets/resources/mouse-grcm39-v1.json).
Versioned derived GO/RBP/CLIP files download from the matching Katana release;
ATtRACT and the KEGG source download directly from their providers. The KEGG
conversion uses a pinned symbol map and preserves the existing result bytes.
Each first download is checked against its pinned SHA256. Later uses share the
cached file. Interrupted downloads do not become accepted cache entries.
A changed upstream snapshot produces a clear download error instead of silently
changing the scientific reference. Use custom files or a newer reviewed catalog
when changing database versions.

## Launch-host downloads

`KATANA_FLOW_CONFIGURE` uses the local executor on the host running Nextflow.
This keeps database downloads off network-isolated compute nodes. It still
runs as a native Nextflow task in the selected Katana container. Scientific
tasks continue through the configured scheduler; no external coordinator or
additional user configuration file is needed. The launch host must support
the selected container runtime and network access, or have a populated cache.

If a source cannot be reached, the error identifies its filename and URL.
Completed cached resources are retained. Correct connectivity, then use
`-resume`; an aborted run with no DEXSeq counts does not imply duplicate samples.

## Overrides and offline use

Set `--resource_cache /project/cache/katana/resources` to share downloads between
studies. It is not removed by `--clean_work`. Pre-populate this directory on a
network-enabled host for offline clusters. Alternatively pass explicit files:

```bash
--gene_sets '/data/GO.gmt,/data/KEGG.gmt' \
--rbp_gene_sets /data/RBP.gmt \
--motif_database /data/ATtRACT.zip \
--clip_peaks '/data/NOVA1.bed,/data/PTBP1.bed'
```

A custom file replaces the catalog entries of that resource kind, not the other
kinds. Use `--resource_bundle none` for custom-only input. Other species/builds
currently require explicit compatible resources. FASTA/GTF are always supplied
by the user so the annotation and alignment reference are chosen together.

## Attribution

GO annotations are distributed with attribution under the
[GO citation and licensing policy](https://www.geneontology.org/docs/go-citation-policy/).
Retain the Gene Ontology Consortium and original annotation-provider citations.
For ATtRACT, cite Giudice et al., *Database* (2016), baw035. Cite Enrichr and KEGG
when using the KEGG-derived library, and each GEO study when using its CLIP data.
GEO's [data disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html)
and the original submitters' rights apply. Katana's software license does not
replace database licenses. All source accessions and processing metadata remain
in the catalog and exported analysis project.
