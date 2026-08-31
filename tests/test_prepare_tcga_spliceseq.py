from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


def test_panel_downloader_imports_as_packaged_module() -> None:
    from scripts.download_tcga_spliceseq_panel import SUPPORTED_SPLICE_TYPES

    assert {"ES", "RI", "AA", "AD", "AP", "AT", "ME"} <= SUPPORTED_SPLICE_TYPES


def _module():
    script = Path(__file__).parents[1] / "scripts" / "prepare_tcga_spliceseq.py"
    spec = importlib.util.spec_from_file_location("prepare_tcga_spliceseq", script)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_prepare_tcga_spliceseq_removes_clinical_rows_and_normal_samples(
    tmp_path: Path,
) -> None:
    source = tmp_path / "download.tsv"
    source.write_text(
        "\t".join(
            [
                "symbol",
                "as_id",
                "splice_type",
                "exons",
                "from_exon",
                "to_exon",
                "novel_splice",
                "pct_with_values",
                "psi_range",
                "std_psi",
                "TCGA_A",
                "TCGA_A_Norm",
                "avg_exp_pct",
                "annotation",
                "",
            ]
        )
        + "\n"
        + "\t".join(["Gender"] + [""] * 9 + ["FEMALE", "FEMALE", "", "", ""])
        + "\n"
        + "\t".join(
            [
                "GENE1",
                "15",
                "ES",
                "9",
                "8",
                "10.1",
                "0",
                "1",
                "0.5",
                "0.2",
                "0.25",
                "0.75",
                "1",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    module = _module()
    result = module.prepare_tcga_spliceseq(source, tmp_path / "out", cohort="TEST")

    full = pd.read_csv(tmp_path / "out" / "TEST_ES_PSI_matrix.tsv", sep="\t", index_col=0)
    tumor = pd.read_csv(
        tmp_path / "out" / "TEST_ES_PSI_tumor_only.tsv",
        sep="\t",
        index_col=0,
    )
    sample_metadata = pd.read_csv(
        tmp_path / "out" / "TEST_sample_metadata.tsv",
        sep="\t",
    )
    assert full.index.tolist() == ["GENE1|ES|15.0"]
    assert full.columns.tolist() == ["TCGA_A", "TCGA_A_Norm"]
    assert tumor.columns.tolist() == ["TCGA_A"]
    assert result["n_events"] == 1
    assert result["n_normal_samples"] == 1
    assert sample_metadata["phenotype"].tolist() == ["tumor", "normal"]
    assert sample_metadata["Gender"].tolist() == ["FEMALE", "FEMALE"]


def test_prepare_tcga_spliceseq_writes_separate_event_type_files(
    tmp_path: Path,
) -> None:
    source = tmp_path / "download.tsv"
    header = [
        "symbol", "as_id", "splice_type", "exons", "from_exon", "to_exon",
        "novel_splice", "pct_with_values", "psi_range", "std_psi", "TCGA_A",
    ]
    source.write_text(
        "\t".join(header)
        + "\n"
        + "\t".join(["GENE1", "1", "RI", "1", "", "", "0", "1", "0.4", "0.2", "0.3"])
        + "\n"
        + "\t".join(["GENE2", "2", "ES", "2", "", "", "0", "1", "0.5", "0.2", "0.4"])
        + "\n",
        encoding="utf-8",
    )
    module = _module()
    result = module.prepare_tcga_spliceseq(
        source, tmp_path / "out", cohort="TEST", event_type="RI"
    )
    matrix = pd.read_csv(
        tmp_path / "out" / "TEST_RI_PSI_matrix.tsv", sep="\t", index_col=0
    )
    assert matrix.index.tolist() == ["GENE1|RI|1.0"]
    assert result["event_type"] == "RI"
    assert (tmp_path / "out" / "parameters_RI.json").exists()
