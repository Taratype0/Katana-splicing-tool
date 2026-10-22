"""Reference CLIP enrichment on independent, matched event-gene pairs."""
import numpy as np
import pandas as pd
from scipy.stats import binomtest

from katana_next.science.statistics import bh_adjust


def paired_binding_enrichment(evidence, pairs):
    matched = pairs.loc[pairs.status.eq("matched")]
    rows = []
    for resource, hits in evidence.groupby("resource_id", sort=True):
        if hits.event_id.duplicated().any():
            raise ValueError("Each event/resource must occur once")
        if hits.comparison_id.nunique() != 1:
            raise ValueError("CLIP enrichment must be computed separately per comparison")
        statuses = hits.set_index("event_id").status
        a = matched.foreground_event.map(statuses)
        b = matched.background_event.map(statuses)
        valid = a.isin(["overlap_in_reference_biosample", "no_overlap_in_resource"]) & b.isin(["overlap_in_reference_biosample", "no_overlap_in_resource"])
        a, b = a.loc[valid].eq("overlap_in_reference_biosample"), b.loc[valid].eq("overlap_in_reference_biosample")
        only_a, only_b = int((a & ~b).sum()), int((b & ~a).sum())
        discordant = only_a + only_b
        p = float(binomtest(only_a, discordant, .5, alternative="greater").pvalue) if discordant else 1.
        if len(a) < 2:
            p = np.nan
        row = hits.iloc[0]
        rows.append({"resource_id": resource, "comparison_id": row.comparison_id, "RBP": row.RBP,
                     "biosample": row.biosample, "accession": row.accession, "assay": row.assay,
                     "matched_gene_pairs": len(matched), "evaluable_gene_pairs": len(a),
                     "unevaluable_gene_pairs": len(matched) - len(a),
                     "foreground_hit_genes": int(a.sum()), "background_hit_genes": int(b.sum()),
                     "foreground_only_pairs": only_a, "background_only_pairs": only_b,
                     "hit_fraction_difference": float(a.mean() - b.mean()) if len(a) else np.nan,
                     "enrichment_pvalue": p,
                     "status": "tested" if len(a) >= 2 else "insufficient_evaluable_gene_pairs",
                     "test": "one-sided exact paired McNemar/binomial; one event per gene",
                     "fdr_family": "all declared experimental CLIP resource tests within this comparison"})
    result = pd.DataFrame(rows)
    if result.empty:
        raise ValueError("No experimental CLIP event/resource observations")
    result["enrichment_fdr"] = bh_adjust(result.enrichment_pvalue)
    result["row_id"] = result.resource_id
    return result
