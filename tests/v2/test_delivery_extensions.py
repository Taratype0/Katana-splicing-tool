from pathlib import Path
from types import SimpleNamespace
import pandas as pd
import pytest

from katana_next.science.isoform_consequences import TranscriptModel, transcript_features, load_models, load_sequences, consequences
from katana_next.science.single_cell_sashimi import regional_evidence, cigar_evidence
from katana_next.sashimi_fixture import write_indexed_bam


def model(strand="+"):
    tx = TranscriptModel("T1", "G1", "chr1", strand)
    tx.features = {"exon": [(100, 112, "."), (200, 300, ".")], "CDS": [(100, 106, "0")],
                   "start_codon": [(100, 103, ".")], "stop_codon": [(106, 109, ".")]}
    if strand == "-":
        tx.features = {"exon": [(100, 200, "."), (288, 300, ".")], "CDS": [(294, 300, "0")],
                       "start_codon": [(297, 300, ".")], "stop_codon": [(291, 294, ".")]}
    return tx


@pytest.mark.parametrize("strand", ["+", "-"])
def test_orf_mapping_respects_transcription_order(strand):
    tx = model(strand)
    result = transcript_features("T1", tx, "ATGAAATAA" + "A"*103)
    assert result["protein_sequence"] == "MK"
    assert result["utr3_length"] == 103
    assert result["stop_to_last_junction_nt"] == 3
    assert result["nmd_prediction"] == "rule_negative_not_proof_of_escape"
    tx.tags.add("cds_end_NF")
    assert transcript_features("T1", tx, "ATGAAATAA" + "A"*103)["protein_status"] == "partial_or_missing_codon_annotation"


def test_nmd_rule_is_prediction_and_strict_boundary():
    tx = model()
    tx.features["exon"] = [(100, 164, "."), (200, 248, ".")]
    assert transcript_features("T1", tx, "ATGAAATAA" + "A"*103)["nmd_prediction"].startswith("rule_negative")
    tx.features["exon"] = [(100, 165, "."), (200, 247, ".")]
    assert transcript_features("T1", tx, "ATGAAATAA" + "A"*103)["nmd_prediction"] == "rule_positive"


def test_missing_or_incompatible_sequence_never_means_normal():
    assert transcript_features("T", None, None)["protein_status"] == "not_evaluable"
    assert transcript_features("T", model(), "ATG")["sequence_status"] == "length_mismatch"
    assert transcript_features("T", model(), "ATGNNNTAA" + "A"*103)["protein_status"] == "invalid_or_ambiguous_CDS_sequence"


def test_version_collisions_rejected(tmp_path):
    fasta = tmp_path / "t.fa"
    fasta.write_text(">T1.1\nAAA\n>T1.2\nAAA\n")
    with pytest.raises(ValueError, match="Ambiguous"):
        load_sequences(fasta, {"T1"})


def test_fasta_annotation_version_mismatch_rejected(tmp_path):
    fasta = tmp_path / "t.fa"
    fasta.write_text(">T1.2\nAAA\n")
    tx = model()
    tx.transcript_id = "T1.1"
    with pytest.raises(ValueError, match="version mismatch"):
        load_sequences(fasta, {"T1"}, expected_models={"T1": tx})


def test_bam_bridge_workflow_runs_in_isolation(tmp_path):
    from scripts.verify_delivery_extensions import verify
    verify(tmp_path / "integration")


def test_actual_bam_fetch_umi_dedup_and_coverage(tmp_path):
    bam = tmp_path / "reads.bam"
    write_indexed_bam(bam)
    metadata = pd.DataFrame([dict(cell_id="c1", sample_id="s1", donor_id="d1", cell_type="T", condition="WT", assay="10x")])
    region, coverage, counts, audit = regional_evidence([SimpleNamespace(sample_id="s1", assets={"bam": bam})], metadata,
                                       ["chr1:120-220:+", "chr1:120-320:+"], flank=20)
    assert counts.junction_molecules.tolist() == [2, 1]
    assert counts.junction_usage.tolist() == [2/3, 1/3]
    assert coverage[("d1", "WT", "T")][0] == 4
    assert coverage[("d1", "WT", "T")][30] == 0
    assert audit.accepted_reads.iloc[0] == 4
    assert counts.cells_with_regional_reads.iloc[0] == 1


def test_deletions_do_not_become_coverage_or_junctions():
    blocks, junctions = cigar_evidence(SimpleNamespace(reference_start=100, cigartuples=[(0, 10), (2, 5), (0, 10), (3, 100), (0, 10)]), 0, 1000)
    assert blocks == [(100, 110), (115, 125), (225, 235)]
    assert junctions == [(125, 225)]


def test_desktop_runtime_stays_explicit(tmp_path):
    from katana_next.desktop import desktop_runtime_root
    assert desktop_runtime_root(tmp_path / "runtime") == (tmp_path / "runtime").resolve()
    assert (tmp_path / "runtime/cache").is_dir()
