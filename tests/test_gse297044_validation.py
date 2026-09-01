from __future__ import annotations

import pandas as pd

from scripts.validate_gse297044_setd2_programs import evaluate


def test_program_replication_requires_converged_read_supported_events() -> None:
    frozen = pd.DataFrame(
        {
            "run_id": ["run", "run"],
            "event_type": ["ES", "ES"],
            "method": ["CoSplice", "CoSplice"],
            "module": ["M1", "M1"],
            "event_key": ["A|ES|1", "B|ES|2"],
            "gene_symbol": ["A", "B"],
        }
    )
    models = pd.DataFrame(
        {
            "event_key": ["J1", "J2", "J3"],
            "mapped_genes": [{"A"}, {"B"}, {"C"}],
            "gene_symbol": ["A", "B", "C"],
            "event_coordinates": ["chr1:1-2", "chr1:3-4", "chr1:5-6"],
            "condition_effect_FDR": [0.01, 0.01, 0.01],
            "KO_minus_WT_delta_PSI": [0.20, 0.30, -0.20],
            "converged": ["True", "False", "True"],
            "median_total_junction_reads": [50, 60, 70],
        }
    )

    summary, membership, audit = evaluate(
        frozen,
        models,
        fdr_cutoff=0.05,
        dpsi_cutoff=0.10,
    )

    assert summary.loc[0, "responsive_genes"] == "A"
    assert membership["GSE297044_event_key"].tolist() == ["J1"]
    assert audit.set_index("metric").loc["GSE297044_responsive_junctions", "value"] == 2
