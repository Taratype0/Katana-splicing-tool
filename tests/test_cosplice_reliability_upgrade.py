import numpy as np
import pandas as pd
import pytest

from katana_cosplice.partition_metrics import partition_recovery
from katana_cosplice.resampling import bootstrap_plan
from katana_cosplice.reliability_audit import _partial_corr, event_junctions, fixed_module_audit, structural_components


def paired_samples():
    return pd.DataFrame({"sample_id": [f"d{d}_{c}" for d in range(4) for c in "AB"],
                         "condition": list("AB") * 4, "donor_id": np.repeat(["d0", "d1", "d2", "d3"], 2)})


def test_paired_units_stay_together_and_metadata_order_does_not_matter():
    samples = paired_samples()
    plan = bootstrap_plan(samples, samples.sample_id)
    other = bootstrap_plan(samples.sample(frac=1, random_state=10), samples.sample_id)
    rng, rng2 = np.random.default_rng(123), np.random.default_rng(123)
    for _ in range(30):
        indices = plan.draw(rng)
        np.testing.assert_array_equal(indices, other.draw(rng2))
        drawn = samples.iloc[indices]
        counts = drawn.groupby(["donor_id", "condition"]).size().unstack(fill_value=0)
        assert (counts.A == counts.B).all()
        assert drawn.condition.value_counts().to_dict() == {"A": 4, "B": 4}
    assert plan.audit()["resampling_units"] == 4
    assert len(plan.draw_hash(indices)) == 64


def test_bootstrap_does_not_guess_pairing_from_replicate_numbers():
    samples = paired_samples().drop(columns="donor_id").assign(replicate=np.repeat(range(4), 2))
    assert bootstrap_plan(samples, samples.sample_id).audit()["resampling_unit_field"] == "sample_id"
    samples["biological_sample_id"] = samples.sample_id.str[:2]
    assert bootstrap_plan(samples, samples.sample_id).audit()["resampling_units"] == 4
    samples.loc[0, "biological_sample_id"] = None
    with pytest.raises(ValueError, match="complete"):
        bootstrap_plan(samples, samples.sample_id)


def test_singleton_profiles_flag_no_between_unit_resampling():
    samples = pd.DataFrame({"sample_id": ["a", "b"], "condition": ["A", "B"], "donor_id": ["d1", "d2"]})
    assert bootstrap_plan(samples, samples.sample_id).audit()["resampling_status"] == "no_between_unit_resampling"
    with pytest.raises(ValueError, match="not estimable"):
        bootstrap_plan(samples, samples.sample_id).require_resampleable()


def test_discovery_and_validation_use_identical_block_draws():
    from katana_cosplice.consensus import _bootstrap_indices as consensus_draw
    from katana_cosplice.validation import _bootstrap_indices as validation_draw
    samples = paired_samples()
    expected = bootstrap_plan(samples, samples.sample_id).draw(np.random.default_rng(17))
    np.testing.assert_array_equal(consensus_draw(samples, samples.sample_id, np.random.default_rng(17)), expected)
    np.testing.assert_array_equal(validation_draw(samples, samples.sample_id, np.random.default_rng(17)), expected)


def test_incomplete_donor_profiles_are_explicit_and_preserve_sample_counts():
    samples = paired_samples().iloc[:-1].copy()
    plan = bootstrap_plan(samples, samples.sample_id)
    assert plan.audit()["resampling_status"] == "partially_resampleable_fixed_strata"
    assert plan.audit()["fixed_singleton_strata"] == 1
    selected = samples.iloc[plan.draw(np.random.default_rng(71))]
    assert selected.condition.value_counts().to_dict() == samples.condition.value_counts().to_dict()
    assert selected.donor_id.eq("d3").sum() == 1


def test_grey_is_not_evidence_of_stable_modules():
    result = partition_recovery(["grey"]*10, ["grey"]*10)
    assert result["ari_all_events"] == 1
    assert np.isnan(result["ari_jointly_assigned"])
    assert result["candidate_coverage"] == 0
    assert np.isnan(result["coassignment_recall"])


def test_perfect_conditional_ari_does_not_hide_lost_events():
    baseline = ["M1"]*10 + ["M2"]*10
    candidate = ["C1"]*2 + ["grey"]*8 + ["C2"]*2 + ["grey"]*8
    result = partition_recovery(baseline, candidate)
    assert result["ari_jointly_assigned"] == 1
    assert result["candidate_coverage"] == .2
    assert result["coassignment_recall"] == pytest.approx(2/90)
    assert result["coassignment_precision"] == 1


def test_partial_correlation_refits_shared_samples():
    rng = np.random.default_rng(23)
    design = np.column_stack([np.ones(20), np.repeat([0, 1], 10)])
    a, b = rng.normal(size=(2, 20)) + design[:, 1]*4
    a[[1, 3, 11]], b[[5, 14]] = np.nan, np.nan
    mask = np.isfinite(a) & np.isfinite(b)
    x = design[mask]
    y = np.column_stack([a[mask], b[mask]])
    residuals = y - x @ np.linalg.lstsq(x, y, rcond=None)[0]
    r, df, status = _partial_corr(a, b, design)
    assert r == pytest.approx(np.corrcoef(residuals.T)[0, 1])
    assert df == 13 and status == "available"
    b[10:] = np.nan
    assert _partial_corr(a, b, design)[2] == "incomplete_pair_design"


def event(key="e1", **changes):
    return {"event_key": key, "gene_symbol": "Gene", "event_type": "SE", "chr": "chr1", "strand": "+",
            "upstreamEE": 10, "exonStart_0base": 20, "exonEnd": 30, "downstreamES": 40, **changes}


def test_all_event_structure_audit_has_no_same_gene_shortcut():
    records = [event(), event("e2", exonStart_0base=21),
               event("e3", upstreamEE=100, exonStart_0base=120, exonEnd=130, downstreamES=140),
               event("unknown", strand=".")]
    structural = [event_junctions(r) for r in records]
    result = structural_components([r["event_key"] for r in records], structural, records).set_index("event_key")
    assert result.loc["e1", "structural_component"] == result.loc["e2", "structural_component"]
    assert result.loc["e3", "structural_status"] == "isolated_known"
    assert result.loc["unknown", "structural_status"] == "not_assessable"
    assert not event_junctions(event(exonStart_0base=20.5))
    assert not event_junctions(event(chr=np.nan))


def test_denominator_labels_do_not_link_different_chromosomes():
    records = [{"event_key": "a", "chr": "chr1", "strand": "+", "cluster_id": "1"},
               {"event_key": "b", "chr": "chr2", "strand": "+", "cluster_id": "1"},
               {"event_key": "c", "chr": "chr1", "strand": "+", "cluster_id": "1"}]
    result = structural_components([r["event_key"] for r in records], [set()]*3, records).set_index("event_key")
    assert result.loc["a", "structural_component"] == result.loc["c", "structural_component"]
    assert result.loc["b", "structural_component"] != result.loc["a", "structural_component"]


def test_components_are_invariant_to_event_order():
    records = [event("a"), event("b", exonStart_0base=22), event("c", strand="-")]
    def get(rows):
        return structural_components([r["event_key"] for r in rows], [event_junctions(r) for r in rows], rows).set_index("event_key").sort_index()
    pd.testing.assert_frame_equal(get(records), get(records[::-1]))


def test_fixed_audit_returns_complete_structure_and_pair_df():
    samples = paired_samples()
    psi = pd.DataFrame(np.random.default_rng(12).uniform(.1, .9, (3, 8)),
                       index=["e1", "e2", "e3"], columns=samples.sample_id)
    events = pd.DataFrame([event(k) for k in psi.index])
    assignments = pd.DataFrame({"event_key": psi.index, "module": "M1"})
    result = fixed_module_audit(psi, samples, events, assignments, max_pairs=1)
    assert len(result["events"]) == 3
    assert len(result["pairs"]) == 1
    assert result["modules"].iloc[0].all_event_structurally_linked_events == 3
    assert result["pairs"].iloc[0].pair_residual_df == 6
    assert result["coverage"].iloc[0].resampling_units == 4
    psi.iloc[0, 0] = np.inf
    with pytest.raises(ValueError, match="PSI"):
        fixed_module_audit(psi, samples, events, assignments)
