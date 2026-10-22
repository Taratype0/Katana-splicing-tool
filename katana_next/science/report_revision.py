"""Contrast-specific report statistics; no pooled-condition inference."""
import numpy as np
import pandas as pd
from scipy.stats import t, ttest_ind
from katana_next.science.statistics import bh_adjust


def dtu_background_genes(features):
    """Include screened genes even when stageR did not confirm their transcripts."""
    confirmation=pd.to_numeric(features.dtu_fdr,errors="coerce").between(0,1)
    screening=pd.to_numeric(features.get("gene_screen_fdr",pd.Series(np.nan,index=features.index)),errors="coerce").between(0,1)
    return set(features.loc[features.source_kind.eq("dtu_result") & (confirmation|screening),"gene_id"])


def module_contrasts(scores, contrasts):
    rows = []
    for contrast in contrasts:
        for module in scores.columns:
            a = scores.reindex(contrast.numerator.sample_ids)[module].dropna().to_numpy(float)
            b = scores.reindex(contrast.denominator.sample_ids)[module].dropna().to_numpy(float)
            row = dict(comparison_id=contrast.contrast_id, module=module,
                       numerator=contrast.numerator.label, denominator=contrast.denominator.label,
                       n_numerator=len(a), n_denominator=len(b), effect=np.nan,
                       ci_low=np.nan, ci_high=np.nan, pvalue=np.nan, status="insufficient_samples")
            if min(len(a), len(b)) >= 2:
                va, vb = a.var(ddof=1)/len(a), b.var(ddof=1)/len(b)
                se = np.sqrt(va+vb)
                row["effect"] = a.mean()-b.mean()
                if se > 0:
                    df = (va+vb)**2/(va**2/(len(a)-1)+vb**2/(len(b)-1))
                    margin = t.ppf(.975, df)*se
                    row.update(ci_low=row["effect"]-margin, ci_high=row["effect"]+margin,
                               pvalue=ttest_ind(a,b,equal_var=False).pvalue, status="exploratory_same_data")
                else:
                    row["status"] = "zero_variance_not_tested"
            rows.append(row)
    result = pd.DataFrame(rows)
    result["FDR"] = bh_adjust(result.pvalue)
    result["FDR_family"] = "all declared pairwise module contrasts in this project"
    return result


def transcript_usage_reliability(transcripts, contrasts, dtu):
    """No gene-screen fallback and no invented stability pass/fail cutoff."""
    frame = transcripts.copy()
    totals = frame.groupby(["sample_id", "gene_id"]).tpm.transform("sum")
    frame["usage_fraction"] = frame.tpm.div(totals.where(totals.gt(0)))
    annotation = frame[["transcript_id", "gene_id", "gene_symbol"]].drop_duplicates()
    matrix = frame.pivot(index="transcript_id", columns="sample_id", values="usage_fraction")
    outputs = []
    for contrast in contrasts:
        a = matrix.reindex(columns=contrast.numerator.sample_ids)
        b = matrix.reindex(columns=contrast.denominator.sample_ids)
        effect = a.mean(axis=1)-b.mean(axis=1)
        loo = []
        for col in a:
            loo.append(a.drop(columns=col).mean(axis=1)-b.mean(axis=1))
        for col in b:
            loo.append(a.mean(axis=1)-b.drop(columns=col).mean(axis=1))
        sensitivity = pd.concat(loo,axis=1)
        part = pd.DataFrame(dict(transcript_id=matrix.index, delta_fraction=effect.values,
            numerator_mean_fraction=a.mean(axis=1).values, denominator_mean_fraction=b.mean(axis=1).values,
            numerator_valid_replicates=a.notna().sum(axis=1).values,
            denominator_valid_replicates=b.notna().sum(axis=1).values,
            loo_min_delta=sensitivity.min(axis=1).values, loo_max_delta=sensitivity.max(axis=1).values,
            loo_direction_fraction=sensitivity.apply(lambda x: np.sign(x)).eq(np.sign(effect),axis=0).where(sensitivity.notna()).mean(axis=1).values))
        evidence = dtu[dtu.comparison_id.eq(contrast.contrast_id) & dtu.source_kind.eq("dtu_result")]
        evidence = evidence.groupby("transcript_id",as_index=False).dtu_fdr.min()
        part = part.merge(annotation,on="transcript_id",validate="one_to_one").merge(evidence,on="transcript_id",how="left",validate="one_to_one")
        part["comparison_id"] = contrast.contrast_id
        part["DTU_supported"] = part.dtu_fdr.le(.05)
        part["quantification_uncertainty"] = "not assessed: bootstrap/Gibbs input not used"
        part["reliability_classification"] = "diagnostics only; no validated binary stability cutoff"
        outputs.append(part)
    return frame, pd.concat(outputs,ignore_index=True)
