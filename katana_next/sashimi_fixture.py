"""Small real indexed BAMs for software acceptance, never biological validation."""
from pathlib import Path
import struct

import pandas as pd
from Bio.bgzf import BgzfWriter

from katana_next.domain.models import AnnotationReference, InputAsset, ProjectManifest, Sample
from katana_next.infrastructure.manifest_store import ManifestStore


def write_indexed_bam(path: Path, cell="c1"):
    text = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"
    header = b"BAM\x01" + struct.pack("<i", len(text)) + text + struct.pack("<ii", 1, 5) + b"chr1\0" + struct.pack("<i", 1000)
    with BgzfWriter(str(path), "wb") as writer:
        writer.write(header)
        begin = writer.tell()
        for i, (gap, umi) in enumerate([(100, "u1"), (100, "u1"), (100, "u2"), (200, "u3")]):
            name = f"r{i}".encode() + b"\0"
            cigar = struct.pack("<III", 20 << 4, gap << 4 | 3, 20 << 4)
            core = struct.pack("<iiIIiiii", 0, 100, (4681 << 16) | (60 << 8) | len(name), 3, 40, -1, -1, 0)
            tags = b"CBZ"+cell.encode()+b"\0UBZ"+umi.encode()+b"\0XSA+"
            block = core + name + cigar + bytes([0x11])*20 + bytes([30])*40 + tags
            writer.write(struct.pack("<i", len(block)) + block)
        end = writer.tell()
    index = b"BAI\x01" + struct.pack("<iiIiQQiQQ", 1, 1, 4681, 1, begin, end, 1, begin, 0)
    Path(str(path)+".bai").write_bytes(index)


def create_sashimi_fixture(root: Path):
    inputs = root / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    metadata, samples = [], []
    for condition in ("WT", "MUT"):
        for replicate in range(1, 4):
            sample = f"{condition}-{replicate}"
            cell = f"{sample}:cell"
            bam = inputs / f"{sample}.bam"
            write_indexed_bam(bam, cell)
            samples.append(Sample(sample_id=sample, biological_sample_id=sample, condition=condition, assets={"single_cell_bam": bam}))
            metadata.append(dict(cell_id=cell, sample_id=sample, donor_id=sample, condition=condition, cell_type="T", assay="10x", count_mode="umi"))
    meta = inputs / "cells.tsv"
    pd.DataFrame(metadata).to_csv(meta, sep="\t", index=False)
    gtf = inputs / "annotation.gtf"
    lines = []
    for tx, end in (("T1", 240), ("T2", 340)):
        for a, b in ((101, 120), (end-19, end)):
            lines.append(f'chr1\tfixture\texon\t{a}\t{b}\t.\t+\t.\tgene_id "G1"; transcript_id "{tx}"; gene_name "FixtureGene";')
    gtf.write_text("\n".join(lines)+"\n", encoding="ascii")
    manifest = ProjectManifest(name="Synthetic BAM bridge acceptance - NOT biological evidence", input_root=inputs,
                    output_root=root / "results", annotation=AnnotationReference(genome_build="synthetic", gtf=gtf),
                    assets=(InputAsset(asset_id="cells", kind="single_cell_metadata", path=meta),), samples=tuple(samples),
                    metadata={"purpose": "Software acceptance with generated indexed BAMs; not public-cohort or biological validation."})
    return ManifestStore().save(manifest, root / "katana.project.json")
