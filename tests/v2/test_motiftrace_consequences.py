import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from katana_cosplice.motif_scanner import reverse_complement
from katana_next.science.bulk import READ_SUPPORT_POLICY_VERSION
from katana_next.science.motiftrace_events import (
    EVENT_TYPES, event_forms, extract_event_rna, match_event_controls, psi_values,
    score_positions, trace_instances,
)
from katana_next.science.motiftrace_clusters import (
    CLUSTER_Z, _profile, build_paired_datasets, cluster_runs, paired_cluster_family,
)


def reference(tmp_path):
    rng = np.random.default_rng(23)
    sequence = "".join(rng.choice(list("ACGT"), 1500))
    path = tmp_path/"reference.fa"
    path.write_text(">chr1\n"+sequence+"\n", encoding="ascii", newline="\n")
    Path(str(path)+".fai").write_text(f"chr1\t{len(sequence)}\t6\t{len(sequence)}\t{len(sequence)+1}\n")
    return path, sequence


def event(kind="SE", strand="+", **changes):
    row = dict(event_id="event", gene_id="g", gene_symbol="G", comparison_id="PBS", event_type=kind,
        strand=strand, chr="chr1", dpsi=.3, fdr=.001,
        significant=True, stable_significant=True, measurement_eligible=True,
        read_support_policy=READ_SUPPORT_POLICY_VERSION, criterion_min_reads=10,
        IJC_SAMPLE_1="60,70", SJC_SAMPLE_1="40,30",
        IJC_SAMPLE_2="30,40", SJC_SAMPLE_2="70,60",
        numerator_mean_reads=100, denominator_mean_reads=100, IncLevel1=".6,.7", IncLevel2=".3,.4",
        upstreamES=0, upstreamEE=40, downstreamES=300, downstreamEE=340,
        exonStart_0base=100, exonEnd=135, **{"1stExonStart_0base": 100, "1stExonEnd": 135,
        "2ndExonStart_0base": 180, "2ndExonEnd": 210}, riExonStart_0base=0, riExonEnd=340)
    if kind in {"A3SS", "A5SS"}:
        suffix = (kind == "A5SS") == (strand == "+")
        row.update(longExonStart_0base=100, longExonEnd=180,
                   shortES=100 if suffix else 130, shortEE=150 if suffix else 180,
                   flankingES=300 if suffix else 0, flankingEE=340 if suffix else 40)
    return {**row, **changes}


MODELS = {"test": {"kind": "consensus", "pattern": "AC", "width": 2, "RBPs": "test_only"}}


@pytest.mark.parametrize("kind,strand", itertools.product(EVENT_TYPES, ("+", "-")))
def test_five_types_both_strands_roundtrip_and_short_exons(tmp_path, kind, strand):
    fasta, genome = reference(tmp_path)
    row = event(kind, strand)
    data = extract_event_rna(pd.DataFrame([row]), fasta)
    assert {"inclusion", "exclusion", "affected_inclusion"} <= set(data.region)
    for r in data.itertuples():
        blocks = json.loads(r.blocks_json)
        expected = "".join(genome[a:b] if strand == "+" else reverse_complement(genome[a:b]) for a, b in blocks)
        assert expected == r.sequence
        assert r.sequence_RNA == expected.replace("T", "U")
        assert hashlib.sha256(expected.encode()).hexdigest() == r.sequence_sha256
    if kind == "SE":
        assert data.loc[data.region.eq("affected_inclusion"), "length"].iloc[0] == 35
    if kind == "RI":
        inc = data[data.region.eq("inclusion")].iloc[0]
        assert json.loads(inc.junctions_json) == []
        assert json.loads(inc.blocks_json) == [[0, 340]]


def test_negative_mxe_inclusion_is_second_genomic_exon_not_first():
    plus, minus = event_forms(event("MXE", "+")), event_forms(event("MXE", "-"))
    assert plus["affected_inclusion"] == [(100, 135)]
    assert minus["affected_inclusion"] == [(180, 210)]
    assert minus["inclusion"] == [(300, 340), (180, 210), (0, 40)]


@pytest.mark.parametrize("kind,strand", itertools.product(("A3SS", "A5SS"), ("+", "-")))
def test_wrong_alternative_site_geometry_fails(kind, strand):
    row = event(kind, strand)
    row["event_type"] = "A3SS" if kind == "A5SS" else "A5SS"
    with pytest.raises(ValueError, match="geometry"):
        event_forms(row)


def test_no_silent_coordinate_or_legacy_budget_truncation(tmp_path):
    fasta, _ = reference(tmp_path)
    expected = extract_event_rna(pd.DataFrame([event()]), fasta)
    pd.testing.assert_frame_equal(expected, extract_event_rna(pd.DataFrame([event()]), fasta, max_total_bases=10))
    with pytest.raises(ValueError):
        event_forms(event(exonEnd=100.5))
    with pytest.raises(ValueError, match="RI"):
        event_forms(event("RI", riExonEnd=341))
    with pytest.raises(ValueError, match="Flanking"):
        event_forms(event(downstreamES=30))


def test_exact_loci_not_net_counts_and_genuine_join_hits(tmp_path):
    fasta, genome = reference(tmp_path)
    rows = pd.DataFrame([event()])
    rna = extract_event_rna(rows, fasta)
    models = {**MODELS, "junction": {"kind": "consensus", "pattern": genome[38:40]+genome[300:302], "width": 4}}
    hits, changes, weighted = trace_instances(rna, rows, models)
    assert {"inclusion_only", "exclusion_only", "retained"} <= set(changes.instance_state)
    for h in hits.itertuples():
        form = rna[rna.region.eq(h.form)].iloc[0]
        assert form.sequence_RNA[h.RNA_start0:h.RNA_end0] == h.hit_RNA
        blocks = json.loads(h.genomic_blocks_json)
        assert sum(b-a for a, b in blocks) == h.RNA_end0-h.RNA_start0
        assert h.crosses_splice_join == (len(blocks) > 1)
    assert changes.loc[changes.instance_state.eq("retained"), "delta_sequence_opportunity"].eq(0).all()
    assert weighted.loc[weighted.instance_state.eq("retained"), "sequence_opportunity"].eq(1).all()
    assert changes[changes.instance_state.eq("inclusion_only")].delta_sequence_opportunity.eq(.3).all()


def test_missing_psi_stays_na_and_reverse_group_is_explicit():
    rows = psi_values(event(IncLevel1="NA,.8", IncLevel2=".3,.4", reverse_direction=True))
    assert rows[0]["group"] == "numerator" and rows[0]["PSI"] == .3
    assert np.isnan(rows[2]["PSI"])
    with pytest.raises(ValueError, match="PSI"):
        psi_values(event(IncLevel1="1.2"))


def test_short_and_ambiguous_opportunities_not_zero_and_methods_separate():
    model = {"kind": "consensus", "pattern": "AAA", "width": 3}
    assert np.isnan(_profile("AA", model, "finite_nonoverlap_counts", 0)).all()
    assert np.isnan(_profile("NNN", model, "overlapping_union_coverage", 0)).all()
    assert _profile("AAAAA", model, "finite_nonoverlap_counts", 0)[0] == 1
    assert _profile("AAAAA", model, "overlapping_union_coverage", 0)[0] == 1
    assert len(score_positions("AAAAA", model)[0]) == 3
    pwm = {"kind": "PWM", "matrix": [[1, 0, 0, 0]]*3, "width": 3}
    assert np.isfinite(_profile("AAA", pwm, "PWM_mean_odds", 0)).all()
    assert np.isnan(_profile("AA", pwm, "PWM_mean_odds", 0)).all()


def test_form_scores_keep_zero_hits_and_use_valid_start_denominators(tmp_path):
    from katana_next.science.motiftrace_events import paired_form_scores
    fasta, _ = reference(tmp_path)
    frame = pd.DataFrame([event()])
    rna = extract_event_rna(frame, fasta)
    models = {"none": {"kind": "consensus", "pattern": "A"*25, "width": 25}}
    result = paired_form_scores(rna, frame, models).iloc[0]
    assert result.inclusion_normalized_score == result.exclusion_normalized_score == 0
    assert result.inclusion_valid_starts == 115-25+1
    assert result.exclusion_valid_starts == 80-25+1
    assert result.delta_PSI_weighted_normalized_score == 0


def test_matching_is_type_specific_unique_genes_and_motif_blind(tmp_path):
    fasta, _ = reference(tmp_path)
    rows = [event(event_id="f1", gene_id="fg"), event(event_id="f2", gene_id="fg"),
            event(event_id="b", gene_id="bg", dpsi=0, fdr=.9, stable_significant=False),
            event("MXE", event_id="wrong_type", gene_id="other", dpsi=0, fdr=.9, stable_significant=False)]
    frame = pd.DataFrame(rows)
    rna = extract_event_rna(frame, fasta, feature_only=True)
    matches = match_event_controls(frame, {"f1", "f2"}, rna)
    assert matches.iloc[0].background_event == "b"
    assert matches.iloc[1].status == "same_gene_extra_event_trace_only"


def data(differences, **meta):
    return {"pair_ids": [str(i) for i in range(len(differences))], "positions": list(range(differences.shape[1])),
            "differences": differences, "meta": {"model_id": "m", "method": "test", "region": "affected_inclusion", "window_nt": 10, **meta}}


def test_cluster_does_not_bridge_na_or_nonadjacent_grids():
    assert [(a, b) for a, b, _ in cluster_runs(np.array([3, 4, np.nan, 3, 0]))] == [(0, 2), (3, 4)]
    bad = data(np.ones((20, 2)))
    bad["positions"] = [0, 2]
    with pytest.raises(ValueError, match="step-one"):
        paired_cluster_family([bad], bad["pair_ids"], permutations=99)


def test_whole_pair_swaps_cluster_p_and_new_window_bh_reproducible():
    diff = np.ones((24, 8))
    diff[:, 4] = np.nan
    first = data(diff)
    win, clusters, null = paired_cluster_family([first], first["pair_ids"], permutations=99)
    win2, clusters2, null2 = paired_cluster_family([first], first["pair_ids"], permutations=99)
    pd.testing.assert_frame_equal(win, win2)
    pd.testing.assert_frame_equal(clusters, clusters2)
    np.testing.assert_array_equal(null, null2)
    assert len(clusters) == 2
    assert win.loc[4, "status"] == "insufficient_independent_pairs"
    assert np.isnan(win.loc[4, "raw_P_new_paired_swaps"])
    assert clusters.cluster_P_FWER_new_family.ge(.01).all()
    assert clusters.cluster_P_FWER_new_family.lt(.05).all()
    assert win.raw_P_new_paired_swaps.dropna().ge(.01).all()
    assert win.raw_P_new_paired_swaps.isna().equals(win.q_new_window_family.isna())
    # Duplicating a perfectly correlated model does not create new independent evidence.
    _, other, null3 = paired_cluster_family([first, data(diff, model_id="duplicate")], first["pair_ids"], permutations=99)
    np.testing.assert_array_equal(null, null3)
    assert other.cluster_P_FWER_new_family.eq(clusters.cluster_P_FWER_new_family.iloc[0]).all()


def test_insufficient_pairs_never_invent_pvalue():
    item = data(np.ones((5, 3)))
    windows, clusters, _ = paired_cluster_family([item], item["pair_ids"], permutations=99)
    assert windows.raw_P_new_paired_swaps.isna().all() and clusters.empty


def test_cluster_null_reproduces_whole_vector_signs_independently():
    diff = np.arange(1, 25.)[:, None] * np.ones((24, 3))
    item = data(diff)
    _, _, actual = paired_cluster_family([item], item["pair_ids"], permutations=99, seed=7)
    signs = np.random.default_rng(7).choice(np.array([-1., 1.]), size=(99, 24))
    ranks = np.arange(1, 25.)
    z = signs @ ranks/np.sqrt(np.sum(ranks*ranks))
    np.testing.assert_allclose(actual, np.maximum(0, z-CLUSTER_Z)*3)


def test_full_pipeline_preserves_inputs_and_exports_all_types(tmp_path):
    from katana_next.science.motiftrace_consequences import analyze_consequences
    fasta, _ = reference(tmp_path)
    frame = pd.DataFrame([event(kind, event_id=kind, gene_id=kind) for kind in EVENT_TYPES])
    before, fa_before = frame.copy(deep=True), fasta.read_bytes()
    result = analyze_consequences(frame, set(frame.event_id), MODELS, fasta, tmp_path/"run", permutations=99)
    assert result["selected_events"] == 5
    assert result["inference_status"].startswith("not_estimable")
    assert result["foreground_events_including_unmatched"] == 5
    actual = pd.read_csv(tmp_path/"run/ACTUAL_RNA_FORMS.tsv.gz", sep="\t")
    assert set(actual.event_type) == set(EVENT_TYPES)
    pd.testing.assert_frame_equal(frame, before)
    assert fa_before == fasta.read_bytes()
    seal = json.loads((tmp_path/"run/COMPLETE.json").read_text())
    for name, record in seal["members"].items():
        assert hashlib.sha256((tmp_path/"run"/name).read_bytes()).hexdigest() == record["sha256"]
    with pytest.raises(FileExistsError):
        analyze_consequences(frame, set(frame.event_id), MODELS, fasta, tmp_path/"run", permutations=99)
    html = (tmp_path/"run/MOTIFTRACE_EVENT_CONSEQUENCES.html").read_text(encoding="utf-8")
    assert "__MOTIFTRACE_CONSEQUENCES__" not in html
    assert "cluster_P_FWER_new_family" in html and "same" not in result["inference_status"]


def test_product_registration_default_and_legacy_modes_retained():
    from katana_next.workflows.registry import WorkflowRegistry
    registry = WorkflowRegistry()
    item = registry.get("prioritize.motiftrace_consequences")
    assert item.normalize_parameters({})["cluster_permutations"] == 999
    primary = registry.get("prioritize.rbp_mechanism")
    assert primary.normalize_parameters({})["sequence_scope"] == "boundary_association"
    scope = next(p for p in primary.definition.parameters if p.key == "sequence_scope")
    assert {"completed_sequence_maps", "se_sequence_impact", "event_consequences"} <= set(scope.options)


def test_portable_report_preserves_full_precision_and_missing_statistics(tmp_path):
    from katana_next.science.motiftrace_consequences import render_consequences

    tiny_p, small_q = 1.234567890123456e-42, 1.234567890123456e-12
    tests = pd.DataFrame({"P": [tiny_p, np.nan], "q": [small_q, np.nan]})
    clusters = pd.DataFrame({"cluster_P_FWER_new_family": [tiny_p]})
    scores = pd.DataFrame({"delta": [1.234567890123456, np.nan]})
    original = tests.copy(deep=True)
    path = tmp_path / "precision.html"
    render_consequences(pd.DataFrame(), pd.DataFrame(), tests, clusters,
                        {"untrusted": "</script><script>alert(1)</script>"}, path,
                        paired_scores=scores)
    html = path.read_text(encoding="utf-8")
    data = json.JSONDecoder().raw_decode(html.split("const D=", 1)[1])[0]
    assert data["tests"][0] == {"P": tiny_p, "q": small_q}
    assert data["tests"][1] == {"P": None, "q": None}
    assert data["clusters"][0]["cluster_P_FWER_new_family"] == tiny_p
    assert data["paired_scores"][0]["delta"] == scores.iloc[0].delta
    assert data["paired_scores"][1]["delta"] is None
    assert "</script><script>alert(1)" not in html
    pd.testing.assert_frame_equal(tests, original)


def test_real_product_api_new_route_downloads_and_cache(tmp_path):
    from test_supervisor_feedback import _api_catalog_fixture
    from katana_next.domain.models import JobStatus, SelectionSet
    from katana_next.workflows.base import ArtifactWriter
    runtime, client, context, project, _, comparisons = _api_catalog_fixture(tmp_path, clip=True)
    try:
        # The shared Windows fixture writes CRLF but declares LF faidx offsets.
        # Normalize this synthetic FASTA, not any original scientific input.
        fasta = Path(context.manifest.annotation.genome_fasta)
        fasta.write_text(fasta.read_text(), encoding="ascii", newline="\n")
        rows = pd.DataFrame([event(kind, event_id=kind, gene_id=kind, comparison_id=comparisons[0],
            row_id=f"{comparisons[0]}::{kind}") for kind in EVENT_TYPES])
        writer = ArtifactWriter(context, "test.consequence_source", "existing-only", {}, "1.0")
        source = writer.table("Test five-type existing events", "events.parquet", rows, entity_type="splicing_event")
        writer.finish()
        saved = runtime.projects.save_selection(SelectionSet(project_id=project, name="Five test events", entity_type="splicing_event",
            source_artifact_id=source.artifact_id, entity_ids=tuple(rows.row_id)))
        motif = tmp_path/"custom.tsv"
        motif.write_text("motif_id\tRBP\tspecies\tevidence\tsource\tmotif\nM\tR\tMus_musculus\tuser_defined\ttest only\tAC\n")
        parameters = {"selection_id": saved.selection_id, "sequence_scope": "event_consequences",
            "custom_motif_path": str(motif), "cluster_permutations": 99}
        url = f"/api/projects/{project}/workflows/prioritize.rbp_mechanism/run"
        started = client.post(url, json={"parameters": parameters})
        assert started.status_code == 200, started.text
        job = runtime.jobs.wait(started.json()["job_id"], timeout=90)
        assert job.status == JobStatus.COMPLETE, job.error
        artifacts = [context.catalog.get_artifact(i) for i in job.artifact_ids]
        assert all(a.metadata.get("new_inference") is False for a in artifacts)
        report = next(a for a in artifacts if a.path.name == "MOTIFTRACE_EVENT_CONSEQUENCES_APP.html")
        html = client.get(f"/api/projects/{project}/artifacts/{report.artifact_id}/file")
        assert html.status_code == 200 and "ATTACHMENT_URLS=" in html.text
        encoded_urls = json.JSONDecoder().raw_decode(html.text.split("const ATTACHMENT_URLS=", 1)[1])[0]
        assert "ACTUAL_LOCAL_RNA.fa" in encoded_urls
        assert all(client.get(u).status_code == 200 for u in encoded_urls.values())
        assert client.post(url, json={"parameters": parameters}).json()["status"] == "cache_hit"
        pd.testing.assert_frame_equal(pd.read_parquet(source.path), rows)
    finally:
        runtime.jobs.shutdown()
