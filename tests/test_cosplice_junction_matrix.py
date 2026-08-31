from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from katana_cosplice.config import CoSpliceConfig, MatrixInputConfig, load_config
from katana_cosplice.junction_matrix import (
    _filter_and_cluster,
    _record_junctions,
    annotate_junctions_with_gtf,
    map_rmats_to_junctions,
)


def test_record_junctions_uses_zero_based_half_open_intron_coordinates() -> None:
    record = SimpleNamespace(
        is_unmapped=False,
        is_secondary=False,
        is_supplementary=False,
        is_qcfail=False,
        is_duplicate=False,
        mapping_quality=60,
        reference_name="chr1",
        reference_start=100,
        cigartuples=[(0, 10), (3, 90), (0, 12)],
    )
    settings = MatrixInputConfig(min_anchor_bp=8)

    assert _record_junctions(record, settings) == [("chr1", 110, 200)]


def test_competing_junction_clusters_drop_singletons_and_compute_membership() -> None:
    index = pd.MultiIndex.from_tuples(
        [
            ("chr1", 100, 200),
            ("chr1", 100, 250),
            ("chr1", 150, 250),
            ("chr2", 500, 600),
        ],
        names=["chrom", "intron_start_0based", "intron_end_0based"],
    )
    counts = pd.DataFrame(
        {
            "S1": [20, 10, 5, 100],
            "S2": [18, 12, 7, 100],
        },
        index=index,
    )
    settings = MatrixInputConfig(
        min_total_junction_reads=1,
        min_samples_with_reads=1,
        min_cluster_fraction=0.0,
    )

    filtered, clusters = _filter_and_cluster(counts, settings)

    assert len(filtered) == 3
    assert ("chr2", 500, 600) not in filtered.index
    assert len(set(clusters.values())) == 1


def test_gtf_annotation_recognizes_direct_and_novel_junctions(tmp_path: Path) -> None:
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text(
        "\n".join(
            [
                'chr1\ttest\texon\t1\t100\t.\t+\t.\tgene_id "ENSG1.1"; gene_name "Gene1"; transcript_id "TX1";',
                'chr1\ttest\texon\t201\t300\t.\t+\t.\tgene_id "ENSG1.1"; gene_name "Gene1"; transcript_id "TX1";',
                'chr1\ttest\texon\t251\t350\t.\t+\t.\tgene_id "ENSG1.1"; gene_name "Gene1"; transcript_id "TX2";',
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    table = annotate_junctions_with_gtf(
        {("chr1", 100, 200), ("chr1", 100, 250)},
        gtf,
    ).set_index(["chrom", "intron_start_0based", "intron_end_0based"])

    assert table.loc[("chr1", 100, 200), "gene_id"] == "ENSG1"
    assert bool(table.loc[("chr1", 100, 200), "annotated_direct_junction"])
    assert not bool(table.loc[("chr1", 100, 250), "annotated_direct_junction"])


def test_rmats_se_mapping_links_all_three_event_junctions(tmp_path: Path) -> None:
    post = tmp_path / "A-B" / "rmats_post"
    post.mkdir(parents=True)
    table = pd.DataFrame(
        [
            {
                "ID": 1,
                "GeneID": "ENSG1",
                "geneSymbol": "Gene1",
                "chr": "chr1",
                "strand": "+",
                "exonStart_0base": 200,
                "exonEnd": 250,
                "upstreamES": 0,
                "upstreamEE": 100,
                "downstreamES": 300,
                "downstreamEE": 400,
                "FDR": 0.01,
                "IncLevelDifference": 0.2,
            }
        ]
    )
    table.to_csv(post / "SE.MATS.JCEC.txt", sep="\t", index=False)

    mapping = map_rmats_to_junctions(
        {
            ("chr1", 100, 200),
            ("chr1", 250, 300),
            ("chr1", 100, 300),
        },
        tmp_path,
        "JCEC",
    )

    assert set(mapping["junction_role"]) == {
        "upstream_inclusion",
        "downstream_inclusion",
        "skipping",
    }
    assert set(mapping["comparison_id"]) == {"A-B"}
    assert set(mapping["rmats_gene_id"]) == {"ENSG1"}
    assert set(mapping["rmats_gene_symbol"]) == {"Gene1"}


def test_bam_matrix_source_round_trips_through_yaml(tmp_path: Path) -> None:
    rmats = tmp_path / "rmats"
    rmats.mkdir()
    metadata = tmp_path / "samples.csv"
    metadata.write_text("sample_id,condition\nS1,A\n", encoding="utf-8")
    config = CoSpliceConfig(
        project_name="junction-test",
        rmats_root=rmats,
        sample_metadata=metadata,
        output_dir=tmp_path / "output",
        matrix_input=MatrixInputConfig(
            input_mode="bam_junctions",
            backend="bamnostic",
            min_anchor_bp=12,
            max_parallel_samples=2,
        ),
    )

    path = config.write(tmp_path / "cosplice_config.yaml")
    loaded = load_config(path)

    assert loaded.matrix_input.input_mode == "bam_junctions"
    assert loaded.matrix_input.backend == "bamnostic"
    assert loaded.matrix_input.min_anchor_bp == 12
    assert loaded.matrix_input.max_parallel_samples == 2
