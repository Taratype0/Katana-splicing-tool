from pathlib import Path

import pandas as pd

from src.services.jutils_benchmark_service import JutilsBenchmarkService


def test_jutils_benchmark_service_writes_auditable_outputs(tmp_path: Path) -> None:
    output_root = tmp_path / "katana_output"
    legacy_data = output_root / "08_jutils" / "data"
    legacy_plots = output_root / "08_jutils" / "plots"
    native_root = output_root / "06_native_splicing_patterns" / "pattern_explorer"
    native_plots = native_root / "plots"
    legacy_data.mkdir(parents=True)
    legacy_plots.mkdir(parents=True)
    native_plots.mkdir(parents=True)

    pd.DataFrame({"comparison_name": ["A", "B"]}).to_csv(
        legacy_data / "jutils_manifest.tsv",
        sep="\t",
        index=False,
    )
    (legacy_plots / "heatmap.pdf").write_bytes(b"%PDF")
    (native_plots / "heatmap.png").write_bytes(b"PNG")
    tables = {
        "dpsi_matrix": pd.DataFrame({"event_id": ["e1"], "A": [0.2], "B": [0.1]}),
        "pairwise_as_similarity": pd.DataFrame(
            {
                "comparison_A": ["A"],
                "comparison_B": ["B"],
                "shared_significant_events": [1],
                "shared_event_support_fraction": [0.5],
                "confidence_weighted_cosine_similarity": [0.75],
                "similarity_confidence": ["moderate"],
            }
        ),
    }

    generated = JutilsBenchmarkService().write(native_root, output_root, tables)

    assert generated == [
        native_root / "native_vs_legacy_jutils_benchmark.tsv",
        native_root / "native_vs_legacy_jutils_benchmark.md",
    ]
    benchmark = pd.read_csv(generated[0], sep="\t")
    input_row = benchmark.loc[benchmark["benchmark_item"] == "input_comparisons"].iloc[0]
    assert input_row["legacy_jutils"] == "A, B"
    similarity_row = benchmark.loc[
        benchmark["benchmark_item"] == "comparison_similarity_quantification"
    ].iloc[0]
    assert "confidence_weighted_cosine=0.750" in similarity_row["native_jutils"]
    assert "Native Jutils vs legacy Jutils benchmark" in generated[1].read_text(encoding="utf-8")
