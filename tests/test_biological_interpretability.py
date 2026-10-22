from __future__ import annotations

import pandas as pd

from katana_cosplice.biological_interpretability import (
    evaluate_biological_interpretability,
    gene_from_event_key,
)


def test_gene_from_event_key_handles_spliceseq_and_katana_keys() -> None:
    assert gene_from_event_key("SRSF7|ES|123.0") == "SRSF7"
    assert gene_from_event_key("ENSG000001;SE|chr1|10|20") == "ENSG000001"


def test_interpretability_audit_detects_gene_dominated_module() -> None:
    assignments = {
        "cohort_a": {
            "CoSplice": pd.DataFrame(
                {
                    "event_key": [
                        "GENE1|ES|1",
                        "GENE1|ES|2",
                        "GENE1|ES|3",
                        "GENE2|ES|4",
                        "SRSF7|ES|5",
                        "GENE3|ES|6",
                    ],
                    "consensus_module": ["M1", "M1", "M1", "M1", "M2", "grey"],
                }
            )
        }
    }
    result = evaluate_biological_interpretability(
        assignments,
        run_functional_enrichment=False,
    )
    audit = result.module_gene_audit.set_index("module")
    assert audit.loc["M1", "n_events"] == 4
    assert audit.loc["M1", "n_unique_genes"] == 2
    assert audit.loc["M1", "top_gene_event_fraction"] == 0.75
    assert audit.loc["M2", "rbp_splicing_factor_gene_count"] == 1
    assert result.enrichment_sources.iloc[0]["status"] == "disabled"
    summary = result.method_summary.iloc[0]
    assert summary["n_modules_audited"] == 2
    assert summary["fraction_modules_not_gene_dominated"] == 0.0
