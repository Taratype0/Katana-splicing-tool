"""Function annotation of one frozen event-module lineage."""
from pathlib import Path
import json
import pandas as pd
from scipy.stats import hypergeom
from katana_cosplice.utils import bh_adjust
from katana_next.workflows.base import WorkflowBase, WorkflowDefinition, ParameterDefinition, ArtifactWriter
from katana_next.workflows.sources import validate_source


def enrich_modules(assignments, annotation, libraries):
    label = next(c for c in ("analysis_module", "consensus_module", "module") if c in assignments)
    if not assignments.event_key.is_unique or not annotation.event_key.is_unique:
        raise ValueError("Unique event identifiers are required")
    symbol = next(c for c in ("gene_symbol", "gene_name") if c in annotation)
    members = assignments[["event_key", label]].merge(annotation[["event_key", symbol]], on="event_key", how="left")
    members[symbol] = members[symbol].fillna("").astype(str)
    background = set(members.loc[members[symbol].ne(""), symbol])
    if not background:
        raise ValueError("No gene names map to the measured network events")
    rows = []
    for library, path in libraries.items():
        terms = []
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            fields = line.split("\t")
            if len(fields) >= 3:
                genes = set(fields[2:]) & background
                if genes:
                    terms.append((fields[0], genes))
        for module, group in members.loc[members[label].ne("grey")].groupby(label):
            genes = set(group.loc[group[symbol].ne(""), symbol])
            for term, annotated in terms:
                hits = genes & annotated
                rows.append(dict(module=module, library=library, term=term, input_genes=len(genes),
                    input_gene_list=";".join(sorted(genes)), background_genes=len(background),
                    annotated_background_genes=len(annotated), overlap_count=len(hits),
                    overlap_genes=";".join(sorted(hits)), gene_ratio=len(hits)/len(genes),
                    pvalue=float(hypergeom.sf(len(hits)-1, len(background), len(annotated), len(genes)))))
    result = pd.DataFrame(rows, columns=["module", "library", "term", "input_genes", "input_gene_list", "background_genes", "annotated_background_genes", "overlap_count", "overlap_genes", "gene_ratio", "pvalue"])
    result["FDR"] = result.groupby("library").pvalue.transform(bh_adjust) if len(result) else pd.Series(dtype=float)
    return result, pd.DataFrame({"gene_symbol": sorted(background)})


class ModuleFunctionsWorkflow(WorkflowBase):
    version = "1.0.0"
    definition = WorkflowDefinition("programs.module_functions", "programs", "Module GO / KEGG functions",
        "Which functions are enriched among the genes containing a module's events?",
        "Hypergeometric test; all measured network genes, including unassigned events, form the background. BH across all modules and terms within each library.",
        "All tested terms, module input genes, contributing genes and measured background", "module gene set",
        parameters=(ParameterDefinition("source_artifact_id", "Completed event-module assignments", "artifact", "", "Select the immutable assignment table from one completed discovery run.", options=("programs.cosplice|programs.psi_network", "consensus_modules__assignments")),
                    ParameterDefinition("libraries", "GMT resources (JSON mapping)", "text", "{}", 'Map each library to its species-matched GMT path, e.g. {"GO_BP":"resources/GO_BP.gmt"}.')),
        limitations="Whole-module gene annotation, not condition-specific DEG enrichment or experimental validation of an event's function.")

    def parent_artifacts(self, context, parameters):
        source = context.catalog.get_artifact(parameters.get("source_artifact_id", ""))
        validate_source(context, source)
        if source.analysis_id not in {"programs.cosplice", "programs.psi_network"} or Path(source.path).stem != "consensus_modules__assignments":
            raise ValueError("Select a consensus event-module assignment table")
        from katana_next.workflows.sources import source_complete
        annotation = next((a for a in context.catalog.list_artifacts(context.manifest.project_id)
            if a.run_fingerprint == source.run_fingerprint and a.analysis_id == source.analysis_id
            and Path(a.path).stem == "matrix__events" and source_complete(context, a)), None)
        validate_source(context, annotation)
        return [source, annotation]

    def input_paths(self, context, parameters):
        libraries = parameters.get("libraries")
        if isinstance(libraries, str):
            libraries = json.loads(libraries)
        if not isinstance(libraries, dict) or not libraries:
            raise ValueError("Supply explicit species-matched GMT libraries")
        paths = [Path(p).resolve() for p in libraries.values()]
        if any(not p.is_file() for p in paths):
            raise ValueError("A configured GMT resource is missing")
        return [a.path for a in self.parent_artifacts(context, parameters)] + paths

    def run(self, context, job, run_fingerprint, parameters):
        source, annotation = self.parent_artifacts(context, parameters)
        job.report(.2, "Mapping measured event genes and explicit GMT libraries")
        libraries = parameters["libraries"]
        if isinstance(libraries, str):
            libraries = json.loads(libraries)
        result, background = enrich_modules(pd.read_parquet(source.path), pd.read_parquet(annotation.path), libraries)
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={"source_discovery_run": source.run_fingerprint, "background_scope": "All measured network event genes, including grey", "method": self.definition.method})
        for name, frame in (("module_functions__enrichment", result), ("module_functions__background", background)):
            frame.insert(0, "row_id", [f"{name}:{i}" for i in range(len(frame))])
            writer.table(name, name+".parquet", frame, entity_type="module_function", parents=(source.artifact_id, annotation.artifact_id))
        return writer.finish()
