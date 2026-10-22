import pandas as pd
import pytest

from katana_next.api.schemas import ConfigureGeneFilterRequest
from katana_next.domain.models import AnnotationReference, InputAsset, ProjectManifest
from katana_next.infrastructure.manifest_store import ManifestStore
from katana_next.science.gene_filter import apply_gene_filter, filter_config
from katana_next.services.project_service import ProjectService
from src.domain.models.gene_filter import GeneFilterConfig


def _inputs(tmp_path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text('chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "G1"; gene_biotype "protein_coding";\n'
                   'MT\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "G2"; gene_biotype "protein_coding";\n'
                   'chr1\ttest\tgene\t101\t200\t.\t+\t.\tgene_id "G3"; gene_biotype "processed_pseudogene";\n', encoding="utf-8")
    result = tmp_path / "native.tsv"
    pd.DataFrame({"gene_id": ["G1", "G2", "G3"], "log2FoldChange": [1.23, -.75, 2.5],
                  "padj": [.04, .009, .95]}).to_csv(result, sep="\t", index=False)
    manifest = ProjectManifest(name="new scope", input_root=tmp_path,
                               output_root=tmp_path / "results", annotation=AnnotationReference(genome_build="test", gtf=gtf),
                               assets=(InputAsset(asset_id="native", kind="deg_result", path=result),))
    return manifest, result


def test_new_request_defaults_on_but_legacy_decoder_does_not_change():
    requested = ConfigureGeneFilterRequest()
    assert requested.exclude_mitochondrial and requested.exclude_pseudogenes
    assert GeneFilterConfig.for_new_analysis().exclude_pseudogenes
    assert not GeneFilterConfig.from_mapping(None).active
    assert not ConfigureGeneFilterRequest(exclude_mitochondrial=False, exclude_pseudogenes=False).exclude_pseudogenes


def test_new_bulk_creation_excludes_before_workflows_and_preserves_imported_stats(tmp_path):
    manifest, native = _inputs(tmp_path)
    original = native.read_bytes()
    projects = ProjectService(application_root=tmp_path / "source")
    target = tmp_path / "new" / "katana.project.json"
    context = projects.create_bulk_project(manifest, target)
    config = filter_config(context.manifest)
    assert config.exclude_pseudogenes and config.exclude_mitochondrial
    filtered = apply_gene_filter(context, pd.read_csv(native, sep="\t"))
    assert filtered.gene_id.tolist() == ["G1"]
    assert filtered.padj.tolist() == [.04]
    assert filtered.log2FoldChange.tolist() == [1.23]
    assert native.read_bytes() == original
    policy = context.manifest.metadata["gene_filter"]
    assert policy["application_stage"] == "before_package_workflows"
    assert policy["upstream_model_exclusions"] == "not_verified"
    assert policy["upstream_refit_performed"] is False
    with pytest.raises(FileExistsError, match="cannot overwrite"):
        projects.create_bulk_project(manifest, target)


def test_opening_completed_legacy_manifest_does_not_inject_default_policy(tmp_path):
    manifest, native = _inputs(tmp_path)
    path = ManifestStore().save(manifest, tmp_path / "legacy.json")
    before = path.read_bytes()
    context = ProjectService(application_root=tmp_path / "source").open_project(path)
    assert "gene_filter" not in context.manifest.metadata
    assert not filter_config(context.manifest).active
    assert len(apply_gene_filter(context, pd.read_csv(native, sep="\t"))) == 3
    assert path.read_bytes() == before


def test_new_bulk_creation_validates_annotation_before_writing(tmp_path):
    manifest, _ = _inputs(tmp_path)
    no_annotation = manifest.model_copy(update={"annotation": AnnotationReference(genome_build="test")})
    target = tmp_path / "unwritten" / "katana.project.json"
    with pytest.raises(ValueError, match="Gene-type filtering requires"):
        ProjectService(application_root=tmp_path / "source").create_bulk_project(no_annotation, target)
    assert not target.exists()


def test_explicit_saved_false_policy_is_preserved_for_new_copy(tmp_path):
    manifest, _ = _inputs(tmp_path)
    manifest = manifest.model_copy(update={"metadata": {"gene_filter": {"config": {
        "exclude_pseudogenes": False, "exclude_mitochondrial": False}}}})
    context = ProjectService(application_root=tmp_path / "source").create_bulk_project(manifest, tmp_path / "explicit.json")
    assert not filter_config(context.manifest).active
