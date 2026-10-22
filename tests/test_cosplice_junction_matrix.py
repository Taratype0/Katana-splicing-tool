from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from katana_cosplice.audit import audit_inputs
from katana_cosplice.config import CoSpliceConfig, MatrixInputConfig, load_config
from katana_cosplice.junction_matrix import (
    _filter_and_cluster,
    _record_junctions,
    annotate_junctions_with_gtf,
    load_precomputed_junctions,
    map_rmats_to_junctions,
    resolve_precomputed_junction_sources,
)
from scripts.align_hisat2_stream_junctions import _junctions


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


def test_missing_rmats_annotation_returns_stable_empty_schema() -> None:
    mapping = map_rmats_to_junctions({("chr1", 100, 200)}, None, "JCEC")

    assert mapping.empty
    assert mapping.columns.tolist() == [
        "event_key",
        "comparison_id",
        "rmats_event_type",
        "rmats_event_id",
        "rmats_gene_id",
        "rmats_gene_symbol",
        "junction_role",
        "rMATS_FDR",
        "raw_dPSI",
    ]


def test_hisat2_sam_cigar_parser_matches_bam_coordinates() -> None:
    assert _junctions(
        "chr1",
        101,
        "20M80N30M",
        min_anchor=8,
        max_intron=500_000,
    ) == [("chr1", 120, 200)]


def test_precomputed_junction_source_round_trip(tmp_path: Path) -> None:
    source = tmp_path / "A_1.junction_counts.tsv.gz"
    pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "intron_start_0based": 120,
                "intron_end_0based": 200,
                "read_count": 7,
            }
        ]
    ).to_csv(source, sep="\t", index=False, compression="gzip")
    samples = pd.DataFrame([{"sample_id": "A_1", "condition": "A"}])

    resolved = resolve_precomputed_junction_sources(samples, tmp_path)
    counts, parameters = load_precomputed_junctions(source, "A_1")

    assert resolved.loc[0, "source_path"] == str(source)
    assert counts[("chr1", 120, 200)] == 7
    assert parameters["backend"] == "precomputed_junctions"


def test_precomputed_junction_source_sums_duplicate_coordinates(
    tmp_path: Path,
) -> None:
    source = tmp_path / "A_1.junction_counts.tsv.gz"
    pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "intron_start_0based": [120, 120],
            "intron_end_0based": [200, 200],
            "read_count": [3, 4],
        }
    ).to_csv(source, sep="\t", index=False, compression="gzip")

    counts, parameters = load_precomputed_junctions(source, "A_1")

    assert counts[("chr1", 120, 200)] == 7
    assert parameters["unique_junctions"] == 1
    assert parameters["junction_read_observations"] == 7


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


def test_bam_only_audit_does_not_require_rmats_layout(tmp_path: Path) -> None:
    metadata = tmp_path / "samples.csv"
    metadata.write_text(
        "sample_id,condition\nA_1,A\nA_2,A\nA_3,A\nB_1,B\nB_2,B\nB_3,B\n",
        encoding="utf-8",
    )
    bam_root = tmp_path / "bam"
    bam_root.mkdir()
    for sample in ("A_1", "A_2", "A_3", "B_1", "B_2", "B_3"):
        bam = bam_root / f"{sample}.bam"
        bam.write_bytes(b"bam")
        (bam_root / f"{sample}.bam.bai").write_bytes(b"bai")
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text("", encoding="utf-8")
    config = CoSpliceConfig(
        project_name="bam-only",
        rmats_root=None,
        sample_metadata=metadata,
        output_dir=tmp_path / "output",
        bam_root=bam_root,
        gtf=gtf,
        matrix_input=MatrixInputConfig(input_mode="bam_junctions"),
    )

    result = audit_inputs(config)
    config_path = config.write(tmp_path / "bam_only_config.yaml")
    loaded = load_config(config_path)

    assert result.passed
    assert loaded.rmats_root is None
    rmats_check = result.checks.set_index("check").loc["rmats_layout"]
    assert rmats_check["severity"] == "info"
    assert not bool(rmats_check["passed"])


def test_precomputed_junction_audit_accepts_all_samples(tmp_path: Path) -> None:
    metadata = tmp_path / "samples.tsv"
    metadata.write_text(
        "sample_id\tcondition\nA_1\tA\nA_2\tA\nA_3\tA\nB_1\tB\nB_2\tB\nB_3\tB\n",
        encoding="utf-8",
    )
    source_root = tmp_path / "junctions"
    source_root.mkdir()
    for sample in ("A_1", "A_2", "A_3", "B_1", "B_2", "B_3"):
        pd.DataFrame(
            [{"chrom": "chr1", "intron_start_0based": 10, "intron_end_0based": 20, "read_count": 2}]
        ).to_csv(
            source_root / f"{sample}.junction_counts.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text("", encoding="utf-8")
    config = CoSpliceConfig(
        project_name="precomputed",
        rmats_root=None,
        sample_metadata=metadata,
        output_dir=tmp_path / "output",
        bam_root=source_root,
        gtf=gtf,
        matrix_input=MatrixInputConfig(
            input_mode="bam_junctions",
            backend="precomputed_junctions",
        ),
    )

    result = audit_inputs(config)

    assert result.passed
    coverage = result.checks.set_index("check").loc["local_bam_coverage"]
    assert bool(coverage["passed"])
    assert "junction_tables=6" in str(coverage["detail"])
