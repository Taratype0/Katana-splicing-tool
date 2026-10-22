import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from katana_next.science.se_rbp import (
    REGIONS, _blocks, extract_se_sequences, mapped_blocks, region_enrichment, tx_slice,
)


def genome(tmp_path):
    sequence = ("ACGTTGCA" * 30)
    path = tmp_path/"genome.fa"
    path.write_text(">chr1\n"+sequence+"\n", encoding="ascii", newline="\n")
    path.with_suffix(".fa.fai").write_text(f"chr1\t{len(sequence)}\t6\t{len(sequence)}\t{len(sequence)+1}\n")
    return path, sequence


def event(strand):
    return {"event_id": "e", "gene_id": "g", "gene_symbol": "G", "chr": "chr1", "strand": strand,
            "upstreamES": 0, "upstreamEE": 24, "exonStart_0base": 64, "exonEnd": 88,
            "downstreamES": 136, "downstreamEE": 160, "dpsi": -.2, "fdr": .01}


@pytest.mark.parametrize("strand", ["+", "-"])
def test_real_sequence_and_junction_block_round_trip(tmp_path, strand):
    from katana_cosplice.motif_scanner import reverse_complement
    fasta, reference = genome(tmp_path)
    windows = extract_se_sequences(pd.DataFrame([event(strand)]), fasta, 7)
    assert set(windows.region) == set(REGIONS)
    for row in windows.itertuples():
        blocks = json.loads(row.blocks_json)
        expected = "".join(reference[a:b] if strand == "+" else reverse_complement(reference[a:b]) for a, b in blocks)
        assert expected == row.sequence
        for a, b in [(0, 1), (0, len(expected)), (2, min(10, len(expected)))]:
            projected = mapped_blocks(blocks, strand, a, b)
            actual = "".join(reference[x:y] if strand == "+" else reverse_complement(reference[x:y]) for x, y in projected)
            assert actual == row.sequence[a:b]
    junction = windows[windows.region.eq("exclusion_junction")].iloc[0]
    assert len(mapped_blocks(json.loads(junction.blocks_json), strand, 4, 8)) == 2
    assert junction.sequence_rna == junction.sequence.replace("T", "U")


def test_negative_transcript_order_not_rmats_genomic_upstream():
    blocks = _blocks(SimpleNamespace(**event("-")))
    assert blocks == [(136, 160), (64, 88), (0, 24)]
    assert tx_slice((100, 120), "-", 0, 6) == (114, 120)
    with pytest.raises(ValueError):
        mapped_blocks([(100, 105)], "-", 0, 6)


def test_invalid_exons_are_not_silently_scanned():
    data = event("+")
    data["exonStart_0base"] = 20
    with pytest.raises(ValueError, match="Overlapping"):
        _blocks(SimpleNamespace(**data))


def test_all_regions_directions_and_evidence_share_declared_bh_family():
    catalog = pd.DataFrame({"motif_token": ["a", "b"], "motif_id": ["M1", "M2"],
                            "RBP": ["R", "R"], "evidence_class": ["direct_experiment", "homology_inferred"]})
    pairs = pd.DataFrame({"status": ["matched"]*4, "direction": ["inclusion_up"]*2+["inclusion_down"]*2,
                          "foreground_event": ["f1", "f2", "f3", "f4"], "background_event": ["b1", "b2", "b3", "b4"]})
    windows = pd.DataFrame([{"event_id": e, "region": r} for e in ["f1", "f2", "f3", "f4", "b1", "b2", "b3", "b4"] for r in REGIONS])
    sites = pd.DataFrame({"motif_token": ["a", "a"], "event_id": ["f1", "f2"], "region": ["cassette_exon"]*2})
    result = region_enrichment(sites, catalog, pairs, windows)
    assert len(result) == 2*2*len(REGIONS)
    row = result[(result.evidence_class == "direct_experiment") & (result.region == "cassette_exon") & (result.direction == "inclusion_up")].iloc[0]
    assert row.enrichment_pvalue == .25
    assert row.enrichment_fdr == 1
    assert row.foreground_hit_genes == 2
    assert np.isclose(row.hit_fraction_difference, 1)


def test_background_excludes_opposite_significant_genes_and_reuses_no_gene(tmp_path):
    from katana_next.science.se_rbp import prepare_se_background
    fasta, _ = genome(tmp_path)
    rows = []
    for identifier, gene, change, fdr, sig in [("a", "FG", .2, .01, True), ("b", "FG", -.2, .01, True),
                                               ("c", "OPPOSITE", -.3, .02, True), ("d", "BG", .01, .9, False)]:
        row = event("+")
        row.update(event_id=identifier, gene_id=gene, dpsi=change, fdr=fdr, stable_significant=sig,
                   event_type="SE", comparison_id="A_vs_B", measurement_eligible=True,
                   numerator_mean_reads=50, denominator_mean_reads=50, IncLevel1=".4,.4", IncLevel2=".4,.4")
        rows.append(row)
    pairs, _, _ = prepare_se_background(pd.DataFrame(rows), {"a", "b"}, fasta)
    assert pairs.iloc[0].background_gene == "BG"
    assert pairs.iloc[1].status == "same_gene_additional_event_sites_only"
    assert "OPPOSITE" not in set(pairs.background_gene)


def test_crossing_filter_sequence_and_indistinguishable_rbp_assignments(tmp_path, monkeypatch):
    from katana_cosplice.motif_scanner import Motif
    import katana_next.science.se_rbp as module
    fasta, _ = genome(tmp_path)
    windows = extract_se_sequences(pd.DataFrame([event("-")]), fasta, 7)
    windows["role"] = "foreground"
    windows.loc[windows.region.eq("upstream_intron"), "role"] = "background"
    def fake_scan(motifs, frame, directory, runtime, threshold):
        assert len(motifs) == 1
        junction = frame[frame.region.eq("exclusion_junction")].iloc[0].to_dict()
        hits = []
        for a, b in [(1, 4), (5, 8)]:
            hits.append({**junction, "motif_token": "m000000", "start_1based": a, "stop_1based": b,
                         "score": 10, "site_pvalue": 1e-5})
        return pd.DataFrame(hits), None, {"tested_positions_all_motifs": 20}
    monkeypatch.setattr(module, "calibrated_scan", fake_scan)
    pwm = np.full((4, 4), .25)
    motifs = [Motif(rbp, "M", "source", pwm=pwm, evidence_class="direct_experiment", target_species="Mus_musculus") for rbp in ("R1", "R2")]
    sites, catalog, provenance = module.scan_se_sequences(motifs, windows, tmp_path, "unused")
    assert len(sites) == 2 and len(catalog) == 2
    assert set(sites.start_1based) == {5}
    assert set(sites.indistinguishable_PWM_RBPs) == {"R1;R2"}
    assert all(len(json.loads(x)) == 2 for x in sites.hit_genomic_blocks_json)
    assert provenance["unique_pwms_scanned"] == 1


def test_registered_sequence_route_keeps_legacy_available(monkeypatch):
    from katana_next.workflows.registry import WorkflowRegistry
    from katana_next.workflows.se_rbp_impact import SeRbpImpactWorkflow
    registry = WorkflowRegistry()
    main = registry.get("prioritize.rbp_mechanism")
    assert main.normalize_parameters({})["sequence_scope"] == "boundary_association"
    assert registry.get("prioritize.se_rbp_impact").version == "1.0.1"
    def fake(self, context, job, fingerprint, parameters):
        assert fingerprint == "test-se"
        assert parameters["selection_id"] == "chosen_SE"
        return ["new_sequence_evidence"]
    monkeypatch.setattr(SeRbpImpactWorkflow, "run", fake)
    assert main.run(None, None, "test", main.normalize_parameters({"sequence_scope": "se_sequence_impact", "selection_id": "chosen_SE"})) == ["new_sequence_evidence"]
