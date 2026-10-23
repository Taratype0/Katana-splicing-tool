from pathlib import Path
import subprocess
import time
import pandas as pd

from katana_next.workflows.base import WorkflowBase, WorkflowDefinition, ParameterDefinition, ArtifactWriter
from katana_next.science.dexseq_interaction import read_design


class DexseqInteractionWorkflow(WorkflowBase):
    version = "3.0.1"
    definition = WorkflowDefinition(
        workflow_id="differential.replicate_interaction", workspace_id="differential",
        label="Design-driven DEXSeq interaction",
        question="Does the effect of factor A on exon usage depend on factor B?",
        method="Official DEXSeq negative-binomial full-vs-reduced likelihood-ratio test; one declared 2x2 interaction; BH across tested exon bins.",
        output_summary="Exon-bin interaction FDR, descriptive bin fractions, four-group plots, filtering audit and fitted R object.",
        analysis_unit="DEXSeq exon counting bin, independent biological sample",
        parameters=(
            ParameterDefinition("design_path", "Explicit four-group design TSV", "text", "", "Columns: sample_id, factor_a, factor_b, count_file; optional biological_sample_id and batch. No filename inference."),
            ParameterDefinition("reference_a", "Factor A reference", "text", "", "Explicit baseline, e.g. WT."),
            ParameterDefinition("reference_b", "Factor B reference", "text", "", "Explicit baseline, e.g. BS."),
            ParameterDefinition("r_library", "Optional R package library", "text", "", "Directory containing the installed official DEXSeq package."),
            ParameterDefinition("min_reads", "Minimum bin counts", "number", 10, "Eligibility only; no significance-based selection.", minimum=1),
            ParameterDefinition("min_replicates", "Samples supporting bin eligibility", "number", 2, "Across the declared samples; actual per-cell replication is checked separately.", minimum=2),
        ),
        limitations="Exon-bin usage, not rMATS event PSI. Published modeling strategy does not validate this implementation automatically. Independent samples only; no donor GLMM. Raw bin fractions are descriptive, not model coefficients. Requires full design and count files; never subtracts existing result tables.",
    )

    def _inputs(self, context, parameters):
        path = parameters.get("design_path")
        if not path:
            raise ValueError("Declare a four-group design TSV; all-pair inference is no longer the default")
        frame = read_design(path, parameters.get("reference_a"), parameters.get("reference_b"),
                            [s.sample_id for s in context.manifest.samples])
        return Path(path), frame

    def input_paths(self, context, parameters):
        path, frame = self._inputs(context, parameters)
        return [path, *map(Path, frame.count_file), Path(__file__).parents[1]/"science/dexseq_interaction.R"]

    def run(self, context, job, run_fingerprint, parameters):
        if parameters.get("selection_id") or parameters.get("contrast_ids"):
            raise ValueError("Use the explicit four-group design, not event selection or arbitrary comparison pairs")
        _, frame = self._inputs(context, parameters)
        from katana_next.flow.stages import _executable, _r_environment
        rscript = _executable(context.manifest.metadata.get("rscript_path"), "Rscript")
        r_environment = _r_environment(rscript, parameters)
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={"method":self.definition.method,"limitations":self.definition.limitations,
                "analysis_unit":"exon_bin","reference":"https://doi.org/10.1038/ncomms14060",
                "implementation_status":"requires real-data and simulation acceptance; not automatically validated by citation"})
        directory = writer.directory
        directory.mkdir(parents=True, exist_ok=True)
        frame.to_csv(directory/"experimental_design.tsv",sep="\t",index=False)
        command = [str(rscript),str(Path(__file__).parents[1]/"science/dexseq_interaction.R"),
            str(directory/"experimental_design.tsv"),str(directory),parameters["reference_a"],parameters["reference_b"],
            str(parameters.get("r_library", "")),str(parameters["min_reads"]),str(parameters["min_replicates"])]
        job.report(.05,"DEXSeq: validating design and fitting count model")
        with (directory/"DEXSeq.stdout.log").open("w") as stdout, (directory/"DEXSeq.stderr.log").open("w") as stderr:
            process = subprocess.Popen(command,stdout=stdout,stderr=stderr,env=r_environment)
            try:
                while process.poll() is None:
                    job.check_cancelled()
                    time.sleep(1)
            finally:
                if process.poll() is None:
                    process.terminate()
                    process.wait()
        if process.returncode:
            raise RuntimeError(f"DEXSeq failed; see {directory/'DEXSeq.stderr.log'}")
        result = pd.read_csv(directory/"dexseq_interaction.tsv",sep="\t")
        if context.manifest.metadata.get("gene_filter"):
            from katana_next.science.gene_filter import filter_exon_bins
            (directory/"dexseq_interaction.tsv").replace(directory/"unfiltered_model_dexseq_interaction.tsv")
            result = filter_exon_bins(context, result)
            result.to_csv(directory/"dexseq_interaction.tsv",sep="\t",index=False)
        writer.table("DEXSeq factorial exon-usage interactions","dexseq_interaction.parquet",result,entity_type="exon_bin")
        fractions = pd.read_csv(directory/"exon_bin_fractions.tsv",sep="\t")
        if context.manifest.metadata.get("gene_filter"):
            (directory/"exon_bin_fractions.tsv").replace(directory/"unfiltered_model_exon_bin_fractions.tsv")
            fractions = filter_exon_bins(context, fractions)
            fractions.to_csv(directory/"exon_bin_fractions.tsv",sep="\t",index=False)
        fractions["row_id"] = fractions.bin_id
        writer.table("Descriptive exon-bin fractions, not PSI","exon_bin_fractions.parquet",fractions,entity_type="exon_bin")
        self._plot(result,fractions,frame,directory)
        for name in ["experimental_design.tsv","filter_audit.tsv","DEXSeq_interaction.rds","R_sessionInfo.txt","model.txt","four_group_exon_usage.png","four_group_exon_usage.svg","DEXSeq.stdout.log","DEXSeq.stderr.log"]:
            writer.existing_file(name,directory/name,entity_type="dexseq_interaction_provenance")
        return writer.finish()

    @staticmethod
    def _plot(results, fractions, design, directory):
        import matplotlib.pyplot as plt
        import numpy as np
        chosen = results.dropna(subset=["interaction_FDR"]).nsmallest(12,"interaction_FDR")
        labels=[]
        if chosen.empty:
            fig, ax = plt.subplots(figsize=(8,3))
            ax.text(.5,.5,"No testable exon-bin interaction results",ha="center",va="center"); ax.axis("off")
        else:
            fig, axes = plt.subplots(int(np.ceil(len(chosen)/3)),3,figsize=(14,3.5*int(np.ceil(len(chosen)/3))),squeeze=False)
            groups = list(design.groupby(["factor_a","factor_b"],sort=True))
            indexed = fractions.set_index("bin_id")
            for ax, (_, row) in zip(axes.flat,chosen.iterrows()):
                for x,(group,meta) in enumerate(groups):
                    values = pd.to_numeric(indexed.loc[row.bin_id,meta.sample_id],errors="coerce").to_numpy(dtype=float)
                    ax.scatter(x+np.linspace(-.09,.09,len(values)),values,s=20)
                    ax.plot([x-.15,x+.15],[np.nanmean(values)]*2,color="black")
                ax.set_xticks(range(4),[" / ".join(g) for g,_ in groups],rotation=35,ha="right")
                identifier=str(row.bin_id)
                group,feature=identifier.rsplit(':',1) if ':' in identifier else (identifier,'')
                genes=group.split('+')
                short=genes[0]+(f" + {len(genes)-1} genes" if len(genes)>1 else '')+':'+feature
                labels.append({'display_label':short,'bin_id':identifier})
                ax.set_title(f"{short}\nInteraction FDR={row.interaction_FDR:.3g}",fontsize=9)
                ax.set_ylabel("Exon-bin / gene-bin count fraction")
            for ax in list(axes.flat)[len(chosen):]: ax.axis("off")
        fig.suptitle("Four-group exon usage: descriptive fractions, not rMATS PSI")
        fig.tight_layout()
        for extension in ["png","svg"]: fig.savefig(directory/f"four_group_exon_usage.{extension}",dpi=160,bbox_inches="tight")
        plt.close(fig)
        pd.DataFrame(labels,columns=['display_label','bin_id']).to_csv(directory/'figure_bin_labels.tsv',sep='\t',index=False)
