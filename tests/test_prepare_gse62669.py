from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


def _load_prepare_function():
    script = Path(__file__).parents[1] / "scripts" / "prepare_gse62669_exon_usage.py"
    spec = importlib.util.spec_from_file_location("prepare_gse62669_exon_usage", script)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.prepare_exon_usage


def test_gse62669_exon_usage_is_normalized_within_gene(tmp_path: Path) -> None:
    prepare_exon_usage = _load_prepare_function()

    source = tmp_path / "counts.tsv"
    pd.DataFrame(
        {
            "ENSEMBL Exon ID": ["GENE1_1", "GENE1_2", "GENE1_3", "GENE2_1"],
            "S1": [10, 20, 30, 50],
            "S2": [20, 20, 20, 50],
            "S3": [30, 20, 10, 50],
        }
    ).to_csv(source, sep="\t", index=False)
    matrix_path = tmp_path / "usage.tsv"
    prepare_exon_usage(
        source,
        matrix_path,
        tmp_path / "metadata.tsv",
        tmp_path / "audit.tsv",
        min_exons_per_gene=3,
        min_gene_total=1,
        min_feature_sd=0,
    )
    matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)
    assert len(matrix) == 3
    assert matrix.groupby(matrix.index.to_series().str.split("|", regex=False).str[0]).sum().iloc[
        0
    ].round(8).eq(1.0).all()
