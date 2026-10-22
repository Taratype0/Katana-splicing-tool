from types import SimpleNamespace

from katana_next.infrastructure.fingerprints import fingerprint_file, fingerprint_value
from katana_next.workflows.programs import CoSpliceWorkflow, configure_production_membership


def test_product_default_annotates_instead_of_rejecting():
    workflow = CoSpliceWorkflow()
    assert workflow.normalize_parameters({})["membership_policy"] == "backbone_confidence_only"
    assert workflow.normalize_parameters({"membership_policy": "backbone_refined_consensus"})["membership_policy"] == "backbone_refined_consensus"


def test_old_size_cutoff_cannot_silently_prune_annotation_only_modules():
    config = SimpleNamespace(consensus=SimpleNamespace(min_module_size=100),
                             network=SimpleNamespace(min_module_size=20))
    configure_production_membership(config, {})
    assert config.consensus.clustering_method == "backbone_confidence_only"
    assert config.consensus.min_module_size == 20


def test_project_reference_overrides_are_fingerprinted(tmp_path, monkeypatch):
    config = tmp_path / "config.yaml"
    config.write_text("config", encoding="utf-8")
    references = [tmp_path / name for name in ("genes.gtf", "genome.fa", "transcripts.fa")]
    for path in references:
        path.write_text("original", encoding="utf-8")
    monkeypatch.setattr("katana_next.workflows.programs.load_config", lambda _: object())
    monkeypatch.setattr("katana_next.workflows.programs.configuration_input_files", lambda _: [])
    context = SimpleNamespace(manifest=SimpleNamespace(
        assets=[SimpleNamespace(kind="cosplice_config", path=config)], samples=[],
        annotation=SimpleNamespace(gtf=references[0], genome_fasta=references[1],
                                   transcript_fasta=references[2])))
    workflow = CoSpliceWorkflow()
    paths = workflow.input_paths(context, {})
    assert set(paths) == {config, *references}
    before = fingerprint_value([fingerprint_file(p) for p in paths])
    references[0].write_text("updated annotation", encoding="utf-8")
    after = fingerprint_value([fingerprint_file(p) for p in workflow.input_paths(context, {})])
    assert before != after
    references[1].unlink()
    assert references[1] in workflow.input_paths(context, {})
