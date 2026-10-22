import json

import numpy as np
import pandas as pd
import pytest

from katana_next.science import motiftrace_usage as usage
from katana_next.science.motiftrace_events import score_positions


def example_usage():
    return usage.usage_fractions(pd.DataFrame([
        {"sample_id": s, "gene_id": "g.1", "transcript_id": tx, "tpm": value}
        for s, values in (("N1", [2., 8.]), ("N2", [4., 16.]), ("D1", [8., 2.]), ("D2", [16., 4.]))
        for tx, value in zip(("t1.1", "t2.1"), values)
    ]), ["N1", "N2"], ["D1", "D2"])


def test_fractions_are_per_biological_sample_and_gene():
    frame = example_usage()
    assert frame.groupby("sample_id").transcript_fraction.sum().eq(1).all()
    assert frame.loc[frame.sample_id.eq("N1") & frame.transcript_id.eq("t1"), "transcript_fraction"].iloc[0] == .2


@pytest.mark.parametrize("mode", ["missing_sample", "missing_transcript", "duplicate", "negative", "mapping"])
def test_missing_and_ambiguous_inputs_are_not_zero_filled(mode):
    frame = example_usage()
    if mode == "missing_sample":
        frame = frame[frame.sample_id.ne("N1")]
    elif mode == "missing_transcript":
        frame = frame.drop(frame.index[0])
    elif mode == "duplicate":
        frame = pd.concat([frame, frame.iloc[:1]])
    elif mode == "negative":
        frame.loc[frame.index[0], "tpm"] = -1
    else:
        frame.loc[frame.index[0], "gene_id"] = "other"
    with pytest.raises(ValueError):
        usage.usage_fractions(frame, ["N1", "N2"], ["D1", "D2"])


def test_zero_gene_tpm_is_na_not_fraction_zero():
    frame = example_usage()
    frame.loc[frame.sample_id.eq("N1"), "tpm"] = 0
    frame = usage.usage_fractions(frame, ["N1", "N2"], ["D1", "D2"])
    assert frame.loc[frame.sample_id.eq("N1"), "transcript_fraction"].isna().all()


@pytest.mark.parametrize("kind", ["consensus", "PWM"])
@pytest.mark.parametrize("word", ["TTTTTTTACGATAGAGACCT", "ACGT", "ANNTTTTTTT", ""])
def test_fast_candidate_adapter_equals_existing_motiftrace(kind, word):
    model = {"kind": "consensus", "pattern": "TTT[TC]", "width": 4}
    if kind == "PWM":
        model = {"kind": "PWM", "matrix": [[.1, .1, .1, .7]]*5, "width": 5}
    scores, selected, _ = score_positions(word, model)
    starts, values, n = usage.candidate_positions(word, model)
    np.testing.assert_array_equal(starts, np.flatnonzero(selected))
    np.testing.assert_allclose(values, scores[starts], rtol=0, atol=0)
    assert n == np.isfinite(scores).sum()


def test_native_stager_gene_q_never_substitutes_transcript_confirmation():
    frame = usage.canonical_dtu(pd.DataFrame({"geneID": ["g.1"], "txID": ["t.1"], "gene": [.001], "transcript": [np.nan]}))
    assert frame.stageR_gene_screen_q.iloc[0] == .001
    assert np.isnan(frame.stageR_transcript_confirmation_q.iloc[0])


def test_dexseq_one_based_conversion_and_multigene_preservation():
    frame = pd.DataFrame({"groupID": ["g.2", "g.2+h.1"], "featureID": ["E1", "E2"], "pvalue": [.001, .2], "padj": [.01, .5],
                          "genomicData.seqnames": ["chr1", "chr1"], "genomicData.start": [101, 201], "genomicData.end": [110, 210], "genomicData.strand": ["-", "-"]})
    result = usage.canonical_dexseq(frame, "contrast")
    assert result.start0.tolist() == [100, 200]
    assert result.end0.tolist() == [110, 210]
    assert result.unambiguous_gene.tolist() == [True, False]
    assert result.DEXSeq_bin_FDR.tolist() == [.01, .5]


def test_exact_expression_usage_decomposition_with_replicate_covariance():
    frame = example_usage()
    frame.loc[frame.sample_id.eq("N2") & frame.transcript_id.eq("t1"), "tpm"] = 8
    frame = usage.usage_fractions(frame, ["N1", "N2"], ["D1", "D2"])
    counts = pd.DataFrame({"transcript_id": ["t1", "t2"], "model_id": ["m", "m"], "RBP": ["R", "R"], "candidate_instances": [2., 0.]})
    samples, summary = usage.weighted_model_opportunities(frame, counts)
    np.testing.assert_allclose(summary.decomposition_residual, 0, atol=1e-12)
    assert samples.motif_P.isna().all() and samples.motif_FDR.isna().all()
    assert summary.delta_mean_instances_per_RNA.iloc[0] < 0


def test_negative_strand_junction_locus_keeps_split_genomic_coordinates():
    rna = pd.DataFrame([{"transcript_id": "t1", "gene_id": "g", "chr": "chr1", "strand": "-", "blocks_json": json.dumps([[200, 203], [100, 103]]),
                         "junctions_json": "[3]", "sequence": "TTTCCC", "sequence_RNA": "UUUCCC"}])
    hits, counts = usage.scan_transcript_models(rna, {"m": {"kind": "consensus", "pattern": "TTCC", "width": 4, "RBPs": "R"}})
    assert len(hits) == 1 and hits.crosses_splice_join.iloc[0]
    assert json.loads(hits.genomic_blocks_json.iloc[0]) == [[200, 202], [101, 103]]
    assert hits.hit_RNA.iloc[0] == "UUCC"
    assert counts.candidate_instances.iloc[0] == 1


def test_locus_weighting_shared_sites_are_not_counted_twice():
    frame = example_usage()
    hits = pd.DataFrame([{**{k: "v" for k in ("model_id", "RBP", "model_kind", "chr", "strand", "genomic_blocks_json", "hit_RNA", "evidence")},
                           "site_key": "shared", "transcript_id": tx, "crosses_splice_join": False} for tx in ("t1", "t2")])
    sample, summary = usage.weighted_locus_opportunities(frame, hits, {"increased_transcript": "t2", "decreased_transcript": "t1"})
    assert sample.RNA_fraction_with_locus.eq(1).all()
    assert summary.pair_state.iloc[0] == "retained_in_pair"
    assert summary.delta_fraction_availability.iloc[0] == 0


def test_short_transcript_not_called_no_motif():
    rna = pd.DataFrame([{"transcript_id": "t", "gene_id": "g", "chr": "chr1", "strand": "+", "blocks_json": "[[0,2]]", "junctions_json": "[]", "sequence": "TT", "sequence_RNA": "UU"}])
    hits, counts = usage.scan_transcript_models(rna, {"m": {"kind": "consensus", "pattern": "TTT", "width": 3}})
    assert hits.empty and counts.candidate_instances.isna().all()


def test_selection_is_motif_blind_and_retains_unconfirmed_increase():
    frame = example_usage()
    dtu = usage.canonical_dtu(pd.DataFrame({"geneID": ["g", "g"], "txID": ["t1", "t2"], "gene": [.01, .01], "transcript": [.02, .8]}))
    bins = pd.DataFrame({"gene_id": ["g"], "unambiguous_gene": [True], "DEXSeq_bin_FDR": [.04]})
    case, all_rows = usage.select_case(frame, dtu, bins, eligible_genes={"g"})
    assert case["decreased_transcript"] == "t1"
    assert case["increased_transcript"] == "t2"
    assert all_rows.loc[all_rows.transcript_id.eq("t2"), "stageR_transcript_confirmation_q"].iloc[0] == .8


def test_gene_abundance_change_does_not_imply_composition_change():
    frame = example_usage()
    frame.loc[frame.arm.eq("numerator"), "tpm"] = [8, 2, 16, 4]
    frame = usage.usage_fractions(frame, ["N1", "N2"], ["D1", "D2"])
    frame.loc[frame.arm.eq("numerator"), "tpm"] *= 2
    frame = usage.usage_fractions(frame, ["N1", "N2"], ["D1", "D2"])
    counts = pd.DataFrame({"transcript_id": ["t1", "t2"], "model_id": ["m", "m"], "RBP": ["R", "R"], "candidate_instances": [2., 0.]})
    _, summary = usage.weighted_model_opportunities(frame, counts)
    assert summary.delta_mean_instances_per_RNA.iloc[0] == 0
    assert summary.composition_component_TPM_weighted.iloc[0] == 0
    assert summary.gene_abundance_component_TPM_weighted.iloc[0] > 0


def test_new_workflow_is_registered_without_replacing_old_entry():
    from katana_next.workflows.registry import WorkflowRegistry
    registry = WorkflowRegistry()
    new = registry.get("prioritize.motiftrace_usage")
    assert new.definition.supports_contrast_selection
    assert registry.get("prioritize.motiftrace_consequences") is not new
    assert "single-case motif P/FDR" in new.definition.method


def test_dexseq_sample_counts_must_be_present():
    bins = pd.DataFrame({"featureID": ["E1"]})
    with pytest.raises(ValueError, match="selected sample"):
        usage.dexseq_bin_support(bins, example_usage(), pd.DataFrame())


@pytest.mark.parametrize("change", ["none", "swapped_arms", "reference", "species"])
def test_prepared_case_cannot_be_silently_reused_for_another_context(tmp_path, change):
    from types import SimpleNamespace as S
    from katana_next.workflows.motiftrace_usage import validate_prepared_context
    paths = [tmp_path/n for n in ("a.gtf", "g.fa", "t.fa")]
    context = S(manifest=S(annotation=S(gtf=paths[0], genome_fasta=paths[1], transcript_fasta=paths[2], genome_build="build"), metadata={"species": "mouse"}))
    contrast = S(contrast_id="NvsD", numerator=S(sample_ids=["n"]), denominator=S(sample_ids=["d"]))
    plan = {"comparison_id": "NvsD", "numerator": {"sample_ids": ["n"]}, "denominator": {"sample_ids": ["d"]}, "genome_build": "build", "species": "mouse"}
    complete = {"samples": ["n", "d"], "reference_sources": [{"path": str(p)} for p in paths]}
    if change == "swapped_arms":
        contrast.numerator.sample_ids, contrast.denominator.sample_ids = ["d"], ["n"]
    elif change == "reference":
        context.manifest.annotation.gtf = tmp_path/"other.gtf"
    elif change == "species":
        context.manifest.metadata["species"] = "human"
    if change == "none":
        validate_prepared_context(context, contrast, plan, complete)
    else:
        with pytest.raises(ValueError):
            validate_prepared_context(context, contrast, plan, complete)
