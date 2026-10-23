"""Reusable actual export fixture for installed-wheel/browser acceptance.

Generated indexed BAMs are software test data, never research evidence. The
helper creates a real ProjectService catalog and complete native-style tables,
then exports through the public exporter (no hand-written request bypass).
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


def create_target_runner_fixture(root, *, sashimi=True):
    from katana_next.domain.models import AnnotationReference, Contrast, ContrastArm, ProjectManifest, Sample
    from katana_next.infrastructure.manifest_store import ManifestStore
    from katana_next.science.bulk import READ_SUPPORT_POLICY_VERSION
    from katana_next.services.project_service import ProjectService
    from katana_next.sashimi_fixture import write_indexed_bam
    from katana_next.targets.bundle import export_targets
    from katana_next.workflows.base import ArtifactWriter
    from katana_next.workflows.registry import WorkflowRegistry

    root = Path(root).resolve()
    source = root / "source-project"
    source.mkdir(parents=True, exist_ok=True)
    resource_root = root / "bound-resources"
    resource_root.mkdir()
    gtf = resource_root / "annotation.gtf"
    gtf.write_text('chr1\tfixture\texon\t101\t120\t.\t+\t.\tgene_id "G2"; transcript_id "T2";\n'
                   'chr1\tfixture\texon\t221\t240\t.\t+\t.\tgene_id "G2"; transcript_id "T2";\n'
                   'chr1\tfixture\texon\t321\t340\t.\t+\t.\tgene_id "G2"; transcript_id "T2";\n', encoding="utf-8")
    samples = []
    for group in ("A", "B"):
        for rep in (1, 2):
            sid = f"{group}{rep}"
            directory = resource_root / sid
            directory.mkdir()
            bam = directory / "alignment.bam"
            write_indexed_bam(bam)
            samples.append(Sample(sample_id=sid, biological_sample_id=f"mouse-{sid}", condition=group,
                                  attributes={"donor_id": f"donor-{sid}"},
                                  assets={"bam": bam, "bai": Path(str(bam) + ".bai")}))
    manifest = ProjectManifest(project_id="portable-target-fixture", name="Generated target workflow acceptance",
                               input_root=source, output_root=source / "results",
                               annotation=AnnotationReference(genome_build="synthetic-v1", gtf=gtf),
                               samples=tuple(samples), contrasts=(Contrast(contrast_id="B_vs_A", label="B versus A",
                                   numerator=ContrastArm(label="B", sample_ids=("B2", "B1")),
                                   denominator=ContrastArm(label="A", sample_ids=("A1", "A2"))),))
    manifest_path = ManifestStore().save(manifest, source / "katana.project.json")
    context = ProjectService().open_project(manifest_path)
    rows = []
    for i in range(1, 6):
        rows.append(dict(row_id=f"event-row-{i}", event_id=f"SE:e{i}", event_type="SE",
                         gene_id="G1" if i < 5 else "G2", gene_symbol="ExampleOne" if i < 5 else "ExampleTwo",
                         comparison_id="B_vs_A", chr="chr1", strand="+", upstreamES=100, upstreamEE=120,
                         exonStart_0base=220, exonEnd=240, downstreamES=320, downstreamEE=340,
                         fdr=1e-200 if i == 1 else .9, pvalue=3e-205 if i == 1 else .4, dpsi=.75,
                         IncLevelDifference=.75, IncLevel1="0.9,0.8", IncLevel2="0.1,0.1",
                         IJC_SAMPLE_1="10,11", SJC_SAMPLE_1="1,1", IJC_SAMPLE_2="1,2", SJC_SAMPLE_2="10,10"))
    frame = pd.DataFrame(rows)
    frame["exonEnd"] = frame["exonEnd"].astype("int32")
    writer = ArtifactWriter(context, "differential.as", "target-fixture-events", {"count_mode": "JC", "fdr": .05}, "3.0.0",
                            common_metadata={"read_support_policy": READ_SUPPORT_POLICY_VERSION})
    artifact = writer.table("All measured splice events", "events.parquet", frame, entity_type="splice_event")
    writer.finish()
    # Source coefficients are native saved results, not newly inferred values.
    extra = [
        ("differential.deg", "Gene expression", pd.DataFrame([
            dict(row_id="deg1", gene_id="G1", gene_symbol="ExampleOne", comparison_id="B_vs_A", log2FoldChange=1.5, padj=.013, pvalue=.002),
            dict(row_id="deg2", gene_id="G2", gene_symbol="ExampleTwo", comparison_id="B_vs_A", log2FoldChange=-.2, padj=.81, pvalue=.6)])),
        ("differential.dtu", "Transcript usage", pd.DataFrame([
            dict(row_id="dtu1", gene_id="G1", transcript_id="T1", comparison_id="B_vs_A", usage_difference=.3, FDR=.021),
            dict(row_id="dtu2", gene_id="G2", transcript_id="T2", comparison_id="B_vs_A", usage_difference=-.1, FDR=.71)])),
        ("differential.dexseq", "Relative exon usage", pd.DataFrame([
            dict(row_id="exon1", gene_id="G1", exon_id="E1", comparison_id="B_vs_A", log2fold=.8, padj=.032),
            dict(row_id="exon2", gene_id="G2", exon_id="E2", comparison_id="B_vs_A", log2fold=-.4, padj=.42)])),
        ("programs.cosplice", "Module membership and hubs", pd.DataFrame([
            dict(row_id=f"member-{i}", event_id=f"SE:e{i}", gene_id="G1" if i < 5 else "G2", module_id="M1" if i < 5 else "M2",
                 membership=.91 - i * .02, hub_rank=i) for i in range(1, 6)])),
        ("motiftrace.consequences", "Saved motif sites and whole-family statistics", pd.DataFrame([
            dict(row_id=f"motif-{i}", event_id=f"SE:e{i}", gene_id="G1" if i < 5 else "G2", comparison_id="B_vs_A",
                 model_id="MODEL1", motif_sequence="UGCAUG", motif_start=223, motif_end=229, consequence="lost",
                 pvalue=.000123456789012345, FDR=.03123456789012345, statistical_family="whole-event-model-window-family") for i in range(1, 6)])),
        ("landscape.bulk", "Complete sample PSI matrix", pd.DataFrame([
            {"row_id": f"psi-{i}", "event_id": f"SE:e{i}", "A1": .1, "A2": .1, "B1": .8, "B2": .9} for i in range(1, 6)])),
    ]
    extra_artifacts = []
    registry = WorkflowRegistry()
    for i, (analysis, name, table) in enumerate(extra):
        try:
            version = registry.get(analysis).version
        except KeyError:
            version = "fixture-1.0"
        writer = ArtifactWriter(context, analysis, f"target-fixture-extra-{i}", {}, version,
                                common_metadata={"read_support_policy": READ_SUPPORT_POLICY_VERSION,
                                                 "purpose": "Generated software acceptance fixture"})
        extra_artifacts.append(writer.table(name, "table.parquet", table,
                                           entity_type="splice_event" if "event_id" in table.columns else "gene"))
        writer.finish()
    selections = [{"source_artifact_id": artifact.artifact_id, "row_ids": ["event-row-1"], "target_kind": "gene",
                   "filters": {"fdr": .01, "min_abs_dpsi": .1, "event_set": "stable_significant"}},
                  {"source_artifact_id": artifact.artifact_id, "row_ids": ["event-row-5"], "target_kind": "event",
                   "filters": {"fdr": .05, "min_abs_dpsi": .1, "event_set": "all_events"}}]
    exported = export_targets(context, selections, root / "export",
                              modules=["gene_card", "event_card"] + (["sashimi"] if sashimi else []))
    bindings_document = json.loads(Path(exported["local_bindings_template"]).read_text(encoding="utf-8"))
    evidence = json.loads(Path(exported["request_path"]).read_text(encoding="utf-8"))["evidence"]
    assert {e["name"] for e in evidence} == {"All measured splice events", *(name for _, name, _ in extra)}
    return {"context": context, "manifest_path": manifest_path, "source_dir": source,
            "resource_root": resource_root, "selections": selections, "source_artifact": artifact,
            "evidence_artifacts": extra_artifacts, "events": frame, "export": exported,
            "bindings": bindings_document, "targets": exported["targets"]}
