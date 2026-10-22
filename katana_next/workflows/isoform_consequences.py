"""Explicit downstream annotation/sequence analysis of cached transcript switches."""
import hashlib

import pandas as pd
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

from katana_next.science.isoform_consequences import load_models, load_sequences, consequences, normalized_id
from katana_next.workflows.base import WorkflowBase, WorkflowDefinition, ParameterDefinition, ArtifactWriter
from katana_next.workflows.sources import resolve_table
from katana_next.workflows.selected_gene_views import save_figure
from katana_next.workflows.figure_exports import _RENDER_LOCK


class IsoformConsequencesWorkflow(WorkflowBase):
    version = "1.2.1"
    definition = WorkflowDefinition(
        workflow_id="prioritize.isoform_consequences", workspace_id="prioritize",
        label="Isoform sequence and functional-consequence support",
        question="Do observed dominant transcript changes alter the annotated coding sequence, UTRs, exon structure or NMD-rule status?",
        method="Match cached transcript IDs to GTF and transcript FASTA; validate annotated ORFs, compare protein sequences and predict the final-junction NMD rule (>55 nt). When versioned Pfam and PyHMMER assets are configured, scan every validated protein at family-specific gathering thresholds and compare detected domain copy counts.",
        output_summary="All switch consequences, transcript annotation audit, validated protein FASTA, optional Pfam hits and gain/loss predictions, summary and paired transcript schematics; source DTU evidence stays separate.",
        analysis_unit="annotated transcript pair within a gene and comparison",
        supports_contrast_selection=True,
        parameters=(ParameterDefinition("max_schematics", "Transcript schematics", "number", 12,
                     "Plot the highest-abundance switches; all switches remain in the table. This is a display limit only.", minimum=0, maximum=100, step=1),),
        limitations="Annotation-based prioritization, not measured protein function, confirmed NMD or a new switch-significance test. Missing/partial CDS stays not evaluable. A missing Pfam hit is not experimental domain loss; identical domain counts do not prove identical function. Signal peptides and localization are not inferred from length.",
    )

    def parent_artifacts(self, context, parameters):
        names = ["Descriptive dominant isoform switches", "Isoform delta TPM by contrast"]
        return [a for name in names if (a := resolve_table(context, "prioritize.isoform_support", name, parameters))]

    def input_paths(self, context, parameters):
        annotation = context.manifest.annotation
        from pathlib import Path
        from katana_next.science.versioned_evidence import verified_resource
        resources=[a for a in context.manifest.assets if a.kind=="uniprot_protein_json"]
        if resources and not context.manifest.metadata.get("species"):
            raise ValueError("Declare project species before using versioned protein annotations")
        for asset in resources:
            verified_resource(asset,species=context.manifest.metadata["species"])
        return [a.path for a in self.parent_artifacts(context, parameters)] + [p for p in (annotation.gtf, annotation.transcript_fasta) if p] + [a.path for a in context.manifest.assets if a.kind in {"uniprot_protein_json","pfam_hmm","hmmer_runtime"}] + [Path(__file__).resolve().parents[2]/"scripts"/"pfam_worker.py"]

    def run(self, context, job, run_fingerprint, parameters):
        source = resolve_table(context, "prioritize.isoform_support", "Descriptive dominant isoform switches", parameters)
        if source is None:
            raise ValueError("Run Isoform support explicitly first. This workflow does not rebuild upstream results.")
        annotation = context.manifest.annotation
        if not annotation.gtf or not annotation.transcript_fasta or not annotation.gtf.is_file() or not annotation.transcript_fasta.is_file():
            raise ValueError("Configure a matching GTF and transcript FASTA in the project annotation.")
        switches = pd.read_parquet(source.path)
        ids = parameters.get("contrast_ids") or [c.contrast_id for c in context.manifest.contrasts]
        switches = switches[switches.comparison_id.isin(ids)].copy()
        if switches.empty:
            raise ValueError("No cached dominant transcript switches in these comparisons. No consequence is inferred.")
        requested = set(switches.numerator_dominant_transcript.map(normalized_id)) | set(switches.denominator_dominant_transcript.map(normalized_id))
        job.report(.08, f"Reading GTF once for {len(requested)} unique transcripts")
        models = load_models(annotation.gtf, requested, job.check_cancelled)
        job.report(.24, "Reading transcript FASTA and checking annotation compatibility")
        sequences = load_sequences(annotation.transcript_fasta, requested, job.check_cancelled, expected_models=models)
        result, features = consequences(switches, models, sequences)
        features["row_id"] = features.transcript_id
        support = resolve_table(context, "prioritize.isoform_support", "Isoform delta TPM by contrast", parameters)
        if support:
            usage = pd.read_parquet(support.path)
            columns = [c for c in ("best_dtu_fdr", "best_dexseq_fdr") if c in usage]
            if columns:
                evidence = usage.groupby(["comparison_id", "gene_id"], as_index=False)[columns].min()
                result = result.merge(evidence, on=["comparison_id", "gene_id"], how="left", validate="many_to_one")
        labels = {c.contrast_id: c.label for c in context.manifest.contrasts}
        result["comparison_label"] = result.comparison_id.map(labels)
        result["direction"] = "denominator transcript -> numerator transcript"
        parents = tuple(a.artifact_id for a in self.parent_artifacts(context, parameters))
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
                                common_metadata={"method": self.definition.method, "limitations": self.definition.limitations,
                                    "NMD_rule": "distance from stop-codon end to final exon junction >55 nt; prediction with exceptions",
                                    "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2697156/",
                                    "DTU_evidence": "best_dtu_fdr is the inherited minimum gene-associated DTU q-value, not this switch-pair FDR"})
        from katana_next.science.protein_domains import resource_paths, scan_domains, domain_consequences
        domain_resources=resource_paths(context)
        if domain_resources:
            job.report(.4,"Scanning all validated isoform proteins against the declared Pfam release")
            scan,key=scan_domains(features,*domain_resources,writer.directory/"pfam",job,
                cache_root=context.manifest.output_root/".katana"/"pfam",
                previous_directories=sorted(writer.directory.parent.glob("*/pfam")))
            result,domain_changes=domain_consequences(result,scan)
            writer.table("Isoform domain gain/loss predictions","domain_changes.parquet",domain_changes,entity_type="domain_consequence",parents=parents)
            writer.table("Per-isoform Pfam domain hits","pfam_hits.parquet",pd.DataFrame(scan["domains"],columns=["transcript_id","domain_accession","domain_name","start_1based","end_inclusive","hmm_start","hmm_end","hmm_length","domain_score","domain_i_evalue"]),entity_type="protein_domain",parents=parents)
            writer.table("Pfam complete protein scan audit","pfam_proteins.parquet",pd.DataFrame(scan["proteins"],columns=["transcript_id","protein_sha256"]),entity_type="protein_scan",parents=parents,
                metadata={"profiles_scanned":scan["profiles_scanned"],"scan_fingerprint":key,"database":domain_resources[0].metadata,"limitations":"GA-threshold prediction; same domain counts do not prove equal function; absent hit is not experimental domain loss"})
            if (writer.directory/"pfam"/"domains.json").exists():
                writer.existing_file("Pfam scan provenance",writer.directory/"pfam"/"domains.json",entity_type="resource_audit",parents=parents)
            if (writer.directory/"pfam"/"reuse.json").exists():
                writer.existing_file("Verified Pfam cache reuse",writer.directory/"pfam"/"reuse.json",entity_type="resource_audit",parents=parents)
        writer.table("Isoform consequence support", "isoform_consequences.parquet", result, entity_type="transcript_switch", parents=parents)
        writer.table("Transcript sequence and annotation audit", "transcript_annotation_audit.parquet", features.drop(columns="protein_sequence"), entity_type="transcript", parents=parents)
        resources = [a for a in context.manifest.assets if a.kind == "uniprot_protein_json"]
        if resources:
            from katana_next.science.versioned_evidence import protein_feature_support
            annotations, provenance = protein_feature_support(features, resources, species=context.manifest.metadata.get("species", ""))
            writer.table("Exact-sequence protein functional annotations", "protein_features.parquet", annotations,
                entity_type="protein_feature", parents=parents,
                metadata={"scope": "Exact complete protein sequence matching only; an unmatched isoform cannot be called domain-lost. UniProt evidence codes retain experimental/predicted distinction."})
            writer.table("Protein annotation provenance", "protein_resources.parquet", provenance, entity_type="resource_audit")
        fasta = writer.directory / "validated_annotated_proteins.fa"
        with fasta.open("w", encoding="ascii") as handle:
            for row in features[features.protein_status.eq("validated_annotated_ORF")].itertuples():
                handle.write(f">{row.transcript_id} annotated_ORF_only\n{row.protein_sequence}\n")
        writer.existing_file("Validated annotated protein sequences", fasta, entity_type="protein_sequence", parents=parents)
        with _RENDER_LOCK:
            if resources:
                statuses=annotations.groupby("transcript_id").mapping_status.agg(set)
                valid_ids=set(features.loc[features.protein_status.eq("validated_annotated_ORF"),"transcript_id"])
                coverage={"Exact sequence + annotated features":0,"Exact sequence, no retained feature":0,
                    "Validated protein, no exact reference match":0,"No validated complete ORF":0}
                for transcript, states in statuses.items():
                    category=("Exact sequence + annotated features" if "exact_protein_sequence" in states else
                        "Exact sequence, no retained feature" if "exact_sequence_no_supported_feature" in states else
                        "Validated protein, no exact reference match" if transcript in valid_ids else "No validated complete ORF")
                    coverage[category]+=1
                figure=Figure(figsize=(11,4.7),layout="constrained")
                axis=figure.subplots()
                bars=axis.barh(list(coverage),list(coverage.values()),color=["#0072B2","#56B4E9","#E69F00","#999999"])
                axis.bar_label(bars,padding=4)
                axis.set_xlim(0,max(coverage.values(),default=1)*1.15)
                axis.set_xlabel("Unique input transcripts (each counted once)")
                axis.set_title("Protein annotation coverage | exact complete-sequence matches only")
                axis.spines[["top","right"]].set_visible(False)
                save_figure(writer,figure,"protein_annotation_coverage","Protein reference coverage",parents)
            if domain_resources:
                figure=Figure(figsize=(11,5),layout="constrained")
                axis=figure.subplots()
                groups=result.groupby(["comparison_label","domain_consequence"]).size().unstack(fill_value=0)
                bottom=[0]*len(groups)
                for i,column in enumerate(groups):
                    axis.barh(groups.index,groups[column],left=bottom,label=column.replace("_"," "),color=["#0072B2","#D55E00","#009E73","#999999","#E69F00"][i%5])
                    bottom=groups[column].to_numpy()+bottom
                axis.set_xlabel("Descriptive transcript pairs; all retained, including not evaluable")
                axis.set_title("Sequence-scanned domain evidence, not measured protein function")
                axis.legend(loc="upper center",bbox_to_anchor=(.5,-.16),ncol=2,frameon=False,fontsize=8)
                axis.spines[["top","right"]].set_visible(False)
                save_figure(writer,figure,"domain_consequence_summary","Isoform domain evidence overview",parents)
            figure = Figure(figsize=(10, 5), layout="constrained")
            axis = figure.subplots()
            summary = result.groupby(["comparison_label", "protein_change"]).size().unstack(fill_value=0)
            bottom = [0] * len(summary)
            for i, column in enumerate(summary):
                values = summary[column].to_numpy()
                axis.barh(summary.index, values, left=bottom, label=column.replace("_", " "), color=["#0072B2", "#8b98a0", "#009E73", "#D55E00"][i % 4])
                bottom = values + bottom
            axis.set_xlabel("Descriptive dominant-transcript switches (all retained)")
            axis.set_title("Annotated protein consequences | support, not functional proof")
            axis.legend(loc="upper center", bbox_to_anchor=(.5, -.16), ncol=2, frameon=False)
            axis.spines[["top", "right"]].set_visible(False)
            save_figure(writer, figure, "isoform_consequence_summary", "Isoform consequence overview", parents)
            shown = result.assign(display_abundance=result.numerator_dominant_tpm + result.denominator_dominant_tpm).sort_values("display_abundance", ascending=False).head(int(parameters["max_schematics"]))
            for i, row in enumerate(shown.itertuples()):
                job.check_cancelled()
                job.report(.55 + .4*i/max(1, len(shown)), f"Transcript consequence schematic {i+1}/{len(shown)}")
                figure = Figure(figsize=(11, 4.3), layout="constrained")
                axis = figure.subplots()
                for y, tx, side in ((1, row.denominator_dominant_transcript, "Reference"), (0, row.numerator_dominant_transcript, "Experiment")):
                    model = models.get(normalized_id(tx))
                    if not model:
                        continue
                    offset = 0
                    for a, b in model.exons:
                        axis.add_patch(Rectangle((offset, y-.1), b-a, .2, facecolor="#d5dee4", edgecolor="#334155", linewidth=.7))
                        offset += b-a
                    for a, b, _ in model.mapped("CDS"):
                        axis.add_patch(Rectangle((a, y-.1), b-a, .2, facecolor="#0072B2", edgecolor="none"))
                    axis.plot([0, offset], [y, y], color="#334155", lw=.5, zorder=0)
                axis.autoscale_view()
                axis.set_ylim(-.5, 1.5)
                axis.set_yticks([0, 1], [f"Experiment\n{row.numerator_dominant_transcript}", f"Reference\n{row.denominator_dominant_transcript}"], fontsize=8)
                axis.set_xlabel("Spliced transcript position (nt), 5' to 3' | blue: annotated CDS; gray: UTR/noncoding")
                axis.set_title(f"{row.gene_symbol} | {row.comparison_label}\n{row.protein_change.replace('_', ' ')}; not a functional assay", fontsize=11)
                axis.spines[["top", "right", "left"]].set_visible(False)
                stem = hashlib.sha256(str(row.row_id).encode()).hexdigest()[:14]
                save_figure(writer, figure, f"isoform_structure_{stem}", f"Transcript structure: {row.gene_symbol} | {row.comparison_label}", parents)
        return writer.finish()
