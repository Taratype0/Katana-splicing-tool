import numpy as np
import pandas as pd
import pytest

from katana_next.science.protein_domains import domain_consequences
from katana_next.science.replicate_prediction import heldout_arm, train_dispersion
from katana_next.science.single_cell import donor_celltype_pseudobulk
from katana_next.science.single_cell_landscape import build_single_cell_landscape
from katana_next.science.single_cell import _replicate_effect_diagnostics
from katana_next.science.experimental_design import bootstrap_blocks
from katana_next.infrastructure.model_checkpoint import ModelCheckpoint, model_groups
from katana_next.science.statistics import _bb_value_gradient, beta_binomial_group_test
from scipy.optimize._numdiff import approx_derivative


def test_binary_analytic_gradient_and_fit_match_original_likelihood():
    rng = np.random.default_rng(41)
    groups = np.repeat([0, 1], 20)
    total = np.full(40, 100.)
    success = rng.binomial(100, rng.beta(5+groups*3, 9-groups*3))
    design = np.column_stack([np.ones(40), groups])
    parameters = np.array([-.2, .6, 3.])
    value, gradient = _bb_value_gradient(parameters, design, success, total)
    numeric = approx_derivative(lambda p: _bb_value_gradient(p, design, success, total)[0], parameters).ravel()
    np.testing.assert_allclose(gradient, numeric, rtol=2e-5, atol=1e-7)
    original = beta_binomial_group_test(success, total, groups)
    analytic = beta_binomial_group_test(success, total, groups, analytic_gradient=True)
    assert original['status'] == analytic['status'] == 'fit'
    np.testing.assert_allclose([analytic['dpsi'], analytic['rho'], analytic['pvalue']],
                               [original['dpsi'], original['rho'], original['pvalue']], rtol=1e-3, atol=1e-5)


def test_model_checkpoint_resumes_only_integrity_checked_prefix(tmp_path):
    checkpoint = ModelCheckpoint(tmp_path/"models.json", "frozen-input-and-parameters")
    grouped = pd.DataFrame({"g": range(120)}).groupby("g")
    rows = []
    with pytest.raises(RuntimeError):
        for key, _ in model_groups(grouped, rows, checkpoint):
            if key == 105:
                raise RuntimeError("simulated interruption")
            rows.append({"g":int(key)})
    assert checkpoint.load()[0] == 100
    rows = []; visited = []
    for key, _ in model_groups(grouped, rows, checkpoint):
        visited.append(key); rows.append({"g":int(key)})
    assert visited == list(range(100, 120)) and len(rows) == 120
    assert checkpoint.load() == (120, rows)
    with pytest.raises(ValueError, match="fingerprint"):
        ModelCheckpoint(checkpoint.path, "changed-input").load()
    checkpoint.path.write_text(checkpoint.path.read_text().replace('"completed_groups": 120','"completed_groups": 119'))
    with pytest.raises(ValueError, match="integrity"):
        checkpoint.load()


@pytest.mark.parametrize("paired", [True, False])
def test_fast_donor_bootstrap_matches_original_resampling(paired):
    group = pd.DataFrame([
        dict(condition=condition, donor_id=f"d{j}" if paired else f"{condition}{j}",
             inclusion=5+j if condition == "A" else 20-j, total=40)
        for condition in ("A", "B") for j in range(3)
    ])
    result = _replicate_effect_diagnostics(group, numerator="A", denominator="B", iterations=50, seed=19)
    donor = group.groupby(["condition", "donor_id"], as_index=False).agg(inclusion=("inclusion", "sum"), total=("total", "sum"))
    rng = np.random.default_rng(19)
    effects = []
    for _ in range(50):
        frame = donor.iloc[bootstrap_blocks(donor, rng)]
        values = (frame.inclusion+.5)/(frame.total+1)
        effects.append(values[frame.condition.eq("A")].mean()-values[frame.condition.eq("B")].mean())
    assert np.allclose(result[1:], np.quantile(effects, [.025, .975]), atol=1e-14)


def test_paired_donors_remain_distinct_in_pseudobulk_and_landscape():
    metadata = pd.DataFrame([
        dict(cell_id=f"{donor}-{condition}", sample_id=f"{donor}-{condition}",
             donor_id=donor, condition=condition, cell_type="T", batch="", assay="smart_seq2")
        for donor in ("d1", "d2", "d3") for condition in ("A", "B")
    ])
    counts = pd.DataFrame([
        dict(cell_id=cell.cell_id, event_id=f"e{event}", cluster_id=f"k{event}",
             gene_id=f"g{event}", gene_symbol="", event_type="junction_choice", annotated=True,
             inclusion=10+event if cell.condition == "A" else 80-event,
             exclusion=90-event if cell.condition == "A" else 20+event, total=100)
        for cell in metadata.itertuples() for event in range(8)
    ])
    bulk, coverage = donor_celltype_pseudobulk(metadata, counts)
    assert len(bulk) == 48 and bulk.row_id.is_unique
    assert len(coverage) == 6 and coverage.row_id.is_unique
    result = build_single_cell_landscape(bulk, min_psi_variance=0, max_events=8)
    assert result.embedding.unit_id.nunique() == 6
    assert set(result.embedding.condition) == {"A", "B"}
    with pytest.raises(ValueError, match="Duplicate donor-condition"):
        build_single_cell_landscape(pd.concat([bulk, bulk.iloc[:1]]))


def prediction_fixture():
    rng=np.random.default_rng(19)
    q=rng.beta(3,3,(50,3)); i=rng.binomial(50,q)
    return pd.DataFrame({"event_id":[f"e{x}" for x in range(50)],"gene_id":[f"g{x}" for x in range(50)],
        "IJC_SAMPLE_1":[",".join(map(str,row)) for row in i],"SJC_SAMPLE_1":[",".join(map(str,row)) for row in 50-i],
        "IncFormLen":200.,"SkipFormLen":100.})


def test_heldout_counts_never_fit_their_own_predictions():
    frame=prediction_fixture()
    detail,fits=heldout_arm(frame,["r1","r2","r3"],["b1","b2","b3"],index=1,draws=100)
    altered=frame.copy()
    altered["IJC_SAMPLE_1"]=altered.IJC_SAMPLE_1.map(lambda s:"49,"+s.split(",",1)[1])
    altered["SJC_SAMPLE_1"]=altered.SJC_SAMPLE_1.map(lambda s:"1,"+s.split(",",1)[1])
    other,other_fits=heldout_arm(altered,["r1","r2","r3"],["b1","b2","b3"],index=1,draws=100)
    pd.testing.assert_frame_equal(fits[fits.held_sample_id.eq("r1")],other_fits[other_fits.held_sample_id.eq("r1")])
    columns=["event_id","method","predictive_lower","predictive_upper","predicted_PSI"]
    pd.testing.assert_frame_equal(detail.loc[detail.held_sample_id.eq("r1"),columns],other.loc[other.held_sample_id.eq("r1"),columns])
    assert set(detail.method)=={"pooled_count_only","training_gene_crossfit_dispersion"}
    assert detail.predictive_lower.between(0,1).all() and detail.predictive_upper.between(0,1).all()
    with pytest.raises(ValueError,match="technical"):
        heldout_arm(frame,["a","b","c"],["a","a","c"],index=1)


def test_missingness_and_gene_unknown_are_not_zero_or_negative_results():
    frame=prediction_fixture()
    frame.loc[0,"gene_id"]=""
    frame.loc[1,"IJC_SAMPLE_1"]="NA,NA,NA"
    result,_=heldout_arm(frame,["a","b","c"],["a","b","c"],index=1,draws=100)
    assert result[result.event_id.isin(["e0","e1"])].method.eq("not_evaluable").all()
    assert train_dispersion(np.ones((3,2)),np.ones((3,2))) is None


def test_domains_require_both_scanned_sequences_and_keep_unscanned_unknown():
    switches=pd.DataFrame([{"row_id":"s","comparison_id":"c","gene_id":"g","annotation_gene_match":True,
        "numerator_dominant_transcript":"a.1","denominator_dominant_transcript":"b.1"}])
    scan={"proteins":[{"transcript_id":"a"},{"transcript_id":"b"}],
          "domains":[{"transcript_id":"b","domain_accession":"PF00076.1"}]}
    result,details=domain_consequences(switches,scan)
    assert result.domain_consequence.iloc[0]=="predicted_domain_loss"
    assert details.delta_domain_count.iloc[0]==-1
    scan["proteins"]=[{"transcript_id":"b"}]
    result,details=domain_consequences(switches,scan)
    assert result.domain_consequence.iloc[0]=="not_evaluable_unscanned_or_gene_mismatch"
    assert details.empty
