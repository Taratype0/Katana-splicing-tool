"""Product entry for verified completed short-exon and sliding-motif evidence."""
from pathlib import Path

import pandas as pd

from katana_next.science.rbp_sequence_results import ASSET_KIND, CompletedRbpBundle, render_browser
from katana_next.science.motiftrace_custom import NAME, load_custom_library, add_custom_results
from katana_next.science.motiftrace_windows import ASSET_KIND as WINDOW_KIND, MatchedWindowBundle, render_windows
from katana_next.workflows.base import ArtifactWriter, ParameterDefinition, WorkflowBase, WorkflowDefinition
from katana_next.workflows.sources import selected_source, selection_rows


class CompletedRbpSequenceWorkflow(WorkflowBase):
    version = "1.2.0"
    definition = WorkflowDefinition(
        workflow_id="prioritize.rbp_sequence_results", workspace_id="prioritize",
        label="Katana MotifTrace: SE RNA, motifs and enrichment",
        question="Where do known motifs occur in the actual SE-affected reference RNA, and what did the completed enrichment test show?",
        method="Reuse SHA-bound public results with exact DAS events, coordinates, effects and background. Preserve original models, short-exon NA, all windows and the public BH family. Optional user models default to position preview; explicit cached-RNA enrichment has a separate, labelled BH family.",
        output_summary="Searchable RNA/motif position browser, original numeric or schematic logos, source/evidence mappings, unchanged complete P/FDR exports, actual RNA FASTA and split-coordinate hit attachments.",
        analysis_unit="frozen splice event, local RNA form and native motif model; not a new biological replicate",
        required_asset_kinds=(ASSET_KIND,),
        optional_asset_kinds=(WINDOW_KIND,),
        parameters=(
            ParameterDefinition("selection_id", "Full frozen DAS SE selection", "selection", "", "Use the original complete foreground; search/filter events inside the viewer without refitting enrichment."),
            ParameterDefinition("bundle_asset_id", "Completed result asset ID", "text", "", "Leave empty only when exactly one registered, SHA-bound completed RBP result bundle exists."),
            ParameterDefinition("window_bundle_asset_id", "Matched count / coverage result asset ID (optional)", "text", "", "Reuse a registered original rMAPS3-count/coverage result bundle and separate external-reference annotations. Empty uses the unique registered bundle; no new statistical test."),
            ParameterDefinition("view_scope", "Display comparisons", "select", "selected_comparison", "Display-only: the original multi-comparison FDR family remains unchanged.", ("selected_comparison", "all_frozen_comparisons")),
            ParameterDefinition("custom_motif_path", "User motif library path (optional)", "text", "", "Consensus TSV or Minimal MEME .meme with .meme.metadata.tsv; exact species/source/evidence required. Original public-library results stay unchanged."),
            ParameterDefinition("custom_motif_mode", "User motif analysis", "select", "preview_only", "Preview adds RNA positions/logos only. Explicit enrichment uses cached RNA and a separate, labelled custom-library BH family; no quantification or original-test rerun.", ("preview_only", "fixed_enrichment")),
        ),
        presentation_role="advanced",
        limitations="Reference sequence matches are not measured binding, methylation, affinity or causality. User-library follow-up is not independent validation or improved sensitivity. This entry requires verified completed SE evidence, not arbitrary cohort inputs. Missing CLIP remains not computed.",
    )

    def _bundle(self, context, parameters):
        assets = [asset for asset in context.manifest.assets if asset.kind == ASSET_KIND]
        if parameters.get("bundle_asset_id"):
            assets = [asset for asset in assets if asset.asset_id == parameters["bundle_asset_id"]]
        if len(assets) != 1:
            raise ValueError("Select exactly one registered completed RBP result asset")
        return CompletedRbpBundle(assets[0], species=context.manifest.metadata.get("species", ""),
                                  genome_build=context.manifest.annotation.genome_build)

    def parent_artifacts(self, context, parameters):
        source = selected_source(context, parameters)
        if source is None:
            raise ValueError("Save the full frozen DAS SE selection from its complete source catalog")
        return [source]

    def _windows(self, context, parameters, bundle):
        assets = [a for a in context.manifest.assets if a.kind == WINDOW_KIND]
        if parameters.get("window_bundle_asset_id"):
            assets = [a for a in assets if a.asset_id == parameters["window_bundle_asset_id"]]
            if not assets:
                raise ValueError("Requested matched count/coverage asset is not registered")
        if not assets:
            return None
        if len(assets) != 1:
            raise ValueError("Select exactly one matched count/coverage result asset")
        return MatchedWindowBundle(assets[0], species=context.manifest.metadata.get("species", ""),
            genome_build=context.manifest.annotation.genome_build, source_manifest_sha256=bundle.asset.metadata["sha256"])

    def input_paths(self, context, parameters):
        if parameters.get("custom_motif_mode") == "fixed_enrichment" and not parameters.get("custom_motif_path"):
            raise ValueError("An explicit user motif file is required for custom enrichment")
        bundle = self._bundle(context, parameters)
        source = self.parent_artifacts(context, parameters)[0]
        complete = pd.read_parquet(source.path)
        bundle.bind(complete, selection_rows(context, parameters, source, complete), view_scope=parameters["view_scope"])
        science = Path(__file__).resolve().parents[1] / "science"
        custom = load_custom_library(parameters["custom_motif_path"], species=context.manifest.metadata.get("species", "")) if parameters.get("custom_motif_path") else None
        paths = [source.path, *bundle.input_paths(), Path(__file__), science / "rbp_sequence_results.py",
                 science / "rbp_sequence_browser.html", science / "motiftrace_custom.py", science / "rbp_short_extension.py",
                 science / "motiftrace_windows.py", science / "motiftrace_windows.html"]
        windows = self._windows(context, parameters, bundle)
        if windows:
            paths.extend(windows.input_paths())
        if custom:
            paths.extend(Path(row["path"]) for row in custom["inputs"])
        return paths

    def run(self, context, job, run_fingerprint, parameters):
        bundle = self._bundle(context, parameters)
        source = self.parent_artifacts(context, parameters)[0]
        complete = pd.read_parquet(source.path)
        selected = selection_rows(context, parameters, source, complete)
        cases = bundle.bind(complete, selected, view_scope=parameters["view_scope"])
        job.check_cancelled()
        job.report(.1, "Verifying native matrices, all local RNA forms and original scalar P/FDR; no rescoring")
        payload = bundle.payload(cases)
        windows = self._windows(context, parameters, bundle)
        has_custom = bool(parameters.get("custom_motif_path"))
        new_custom_inference = has_custom and parameters["custom_motif_mode"] == "fixed_enrichment"
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={"reuse_completed_results": True, "displayed_comparisons": cases,
                             "frozen_FDR_family_comparisons": list(bundle.cases), "new_inference": new_custom_inference,
                             "public_inference_unchanged": True, "tool_name": NAME,
                             "quantification_rerun": False, "CLIP_validation": "reused_external_reference_only" if windows else "not_computed",
                             "genome_build": context.manifest.annotation.genome_build,
                             "limitations": self.definition.limitations})
        if has_custom:
            custom = load_custom_library(parameters["custom_motif_path"], species=context.manifest.metadata["species"])
            payload = add_custom_results(payload, custom, mode=parameters["custom_motif_mode"],
                                         output=writer.directory / "custom_library", check_cancelled=job.check_cancelled)
        parents = (source.artifact_id,)
        summary = pd.DataFrame([row for row in payload["summary"] if row["contrast"] in cases])
        writer.table(NAME + ": original enrichment scopes including negatives and NA", "summary.parquet", summary,
                     entity_type="rbp_enrichment_scope", parents=parents)
        events = pd.DataFrame([{key: value for key, value in event.items() if key != "forms"} | {"comparison_id": case, "row_id": case + "::" + event["event_id"]}
                               for case, rows in payload["events"].items() for event in rows])
        writer.table(NAME + ": complete SE event index (short events retained)", "events.parquet", events,
                     entity_type="rbp_reference_event", parents=parents)
        rna = pd.DataFrame([{key: value for key, value in event.items() if key != "forms"} |
                            {"comparison_id": case, "region": form["region"], "reference_RNA": form["RNA"],
                             "genomic_blocks_RNA_order": str(form["blocks"]), "junction_offsets0": str(form["junctions"])}
                            for case, rows in payload["events"].items() for event in rows for form in event["forms"]])
        writer.table(NAME + ": actual reference RNA and split genomic coordinates", "rna_forms.parquet", rna,
                     entity_type="rbp_reference_rna", parents=parents)
        attachments = []
        if windows:
            job.report(.35, "Verifying saved count/coverage P/FDR and external-reference identities; no new inference")
            window_data = windows.payload(payload, check_cancelled=job.check_cancelled)
            window_attachments = []
            for name in ("tests", "plan", "external_windows", "external_actual", "external_hits", "quality", "external_delivery_acceptance"):
                artifact = writer.existing_file(NAME + ": matched windows / " + name, windows.path(name),
                    entity_type="motiftrace_window_evidence", parents=parents,
                    metadata={"separate_method_FDR": True, "experimental_support_not_significance": True})
                window_attachments.append({"label": name, "url": f"/api/projects/{context.manifest.project_id}/artifacts/{artifact.artifact_id}/file?download=true"})
            window_report = render_windows(window_data, writer.directory / "MOTIFTRACE_MATCHED_SLIDING_WINDOWS.html", window_attachments)
            artifact = writer.existing_file(NAME + ": original rMAPS3 counts and coverage, separate P/FDR", window_report,
                entity_type="motiftrace_matched_windows", metadata={"mime_type": "text/html", "new_inference": False}, parents=parents)
            payload["matched_windows"] = {"url": f"/api/projects/{context.manifest.project_id}/artifacts/{artifact.artifact_id}/file",
                "offline_file": window_report.name, "summary": window_data["summary"], "families": window_data["families"],
                "experimental_support_not_significance": True}
            context.results.write_plot_spec({"verified_members": windows.verified, "registration": str(windows.manifest),
                "original_method_families": window_data["families"], "new_inference": False,
                "evidence_status_counts": pd.Series([s for e in window_data["evidence"].values() for s in e["status"]]).value_counts().to_dict()},
                writer.directory / "MATCHED_WINDOW_REUSE_RECEIPT.json")
            writer.existing_file(NAME + ": matched-window identity and reuse receipt", writer.directory / "MATCHED_WINDOW_REUSE_RECEIPT.json",
                entity_type="rbp_provenance", parents=parents)
        if has_custom:
            for path in sorted((writer.directory / "custom_library").glob("*")):
                artifact = writer.existing_file(NAME + ": user-library " + path.name, path,
                    entity_type="motiftrace_custom_evidence", parents=parents,
                    metadata={"separate_custom_family": True, "not_independent_validation": True})
                attachments.append({"label": "自定义库独立结果 / " + path.name,
                    "url": f"/api/projects/{context.manifest.project_id}/artifacts/{artifact.artifact_id}/file?download=true"})
        for relative, label in (
            ("ALL_NATIVE_ASSIGNMENTS_AND_DEDUPLICATION.tsv", "全部motif/RBP/物种/证据来源映射"),
            ("DELIVERY/ALL_COMPLETE_SEQUENCE_AND_SHORT_STRATUM_TESTS.tsv.gz", "全部完整序列、短序列与接头P/FDR（含NA）"),
            ("ALL_TESTS_WITH_NA_DENOMINATORS.tsv.gz", "全部逐窗及完整序列原P/FDR与分母"),
            ("ACTUAL_AFFECTED_RNA_MOTIF_HITS.tsv.gz", "全部原始序列匹配、候选及真实分段坐标"),
            ("NATIVE_NUMERIC_MATRICES.npz", "原生数值矩阵"),
            ("NATIVE_CONSENSUS_PATTERNS.json", "原生有限模式"),
            ("ANALYSIS_PLAN_AND_FREEZE.json", "原统计设计与冻结家族"),
        ):
            job.check_cancelled()
            path = bundle.path(relative)
            artifact = writer.existing_file(label, path, entity_type="rbp_frozen_attachment", parents=parents,
                                            metadata={"frozen_scope": list(bundle.cases), "not_display_filtered": True})
            attachments.append({"label": label, "url": f"/api/projects/{context.manifest.project_id}/artifacts/{artifact.artifact_id}/file?download=true"})
        for name in sorted(bundle.members):
            if name.startswith("DELIVERY/figures/") and name.endswith(".png") and ("__logo" in name or any(case in name for case in cases)):
                writer.existing_file(NAME + ": original frozen figure: " + Path(name).stem, bundle.path(name), entity_type="rbp_frozen_figure",
                                     metadata={"mime_type": "image/png", "display_only": True}, parents=parents)
        job.report(.85, "Rendering portable event/motif browser; original enrichment family unchanged")
        report = render_browser(payload, writer.directory / "SE_RNA_MOTIF_BROWSER.html", attachments=attachments)
        writer.existing_file(NAME + ": RNA, motif positions, logos and original P/FDR", report,
                             entity_type="rbp_sequence_browser", metadata={"mime_type": "text/html"}, parents=parents)
        context.results.write_plot_spec({"verified_members": bundle.verified, "new_inference": new_custom_inference,
                                         "tool_name": NAME, "public_inference_unchanged": True,
                                         "displayed_comparisons": cases, "frozen_family": list(bundle.cases)},
                                        writer.directory / "COMPLETED_RESULT_REUSE_RECEIPT.json")
        writer.existing_file(NAME + ": completed-result reuse verification", writer.directory / "COMPLETED_RESULT_REUSE_RECEIPT.json",
                             entity_type="rbp_provenance", parents=parents)
        job.report(.98, "MotifTrace evidence is browsable/exportable; public results unchanged, no quantification rerun")
        return writer.finish()
