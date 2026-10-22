"""Partition recovery with explicit grey-event and retained-universe denominators."""
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

UNASSIGNED = {"grey", "gray", "0", "unassigned", ""}


def partition_recovery(reference, candidate):
    a, b = np.asarray(reference, dtype=object), np.asarray(candidate, dtype=object)
    if a.ndim != 1 or b.shape != a.shape or not len(a):
        raise ValueError("Partition recovery requires equal nonempty event axes")
    a = pd.Series(a).fillna("grey").astype(str)
    b = pd.Series(b).fillna("grey").astype(str)
    a = a.mask(a.isin(UNASSIGNED), "grey")
    b = b.mask(b.isin(UNASSIGNED), "grey")
    ka, kb = a.ne("grey"), b.ne("grey")
    common = ka & kb
    n = len(a)
    def pairs(counts):
        return int(sum(int(v)*(int(v)-1)//2 for v in counts))
    reference_pairs, candidate_pairs = pairs(a[ka].value_counts()), pairs(b[kb].value_counts())
    recovered_pairs = pairs(pd.DataFrame({"a":a[common], "b":b[common]}).value_counts())
    eligible = bool(common.sum() >= 2 and a[common].nunique() >= 2 and b[common].nunique() >= 2)
    return {"ari_all_events": float(adjusted_rand_score(a,b)),
        "ari_jointly_assigned": float(adjusted_rand_score(a[common],b[common])) if eligible else np.nan,
        "assigned_ari_status": "available" if eligible else "fewer_than_two_modules_on_common_assigned_events",
        "event_universe":n, "reference_assigned_events":int(ka.sum()), "candidate_assigned_events":int(kb.sum()),
        "reference_coverage":float(ka.mean()), "candidate_coverage":float(kb.mean()),
        "jointly_assigned_events":int(common.sum()), "jointly_assigned_fraction":float(common.mean()),
        "assignment_jaccard":float(common.sum()/(ka|kb).sum()) if (ka|kb).any() else np.nan,
        "reference_modules":int(a[ka].nunique()), "candidate_modules":int(b[kb].nunique()),
        "reference_coassigned_pairs":reference_pairs, "candidate_coassigned_pairs":candidate_pairs,
        "recovered_coassigned_pairs":recovered_pairs,
        "coassignment_recall":recovered_pairs/reference_pairs if reference_pairs else np.nan,
        "coassignment_precision":recovered_pairs/candidate_pairs if candidate_pairs else np.nan,
        "scope":"Grey is not a biological module; common-assigned ARI is conditional on retention, read alongside coverage and pair recovery"}
