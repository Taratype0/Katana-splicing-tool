from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from katana_next.domain.models import (
    AnalysisArtifact,
    ProjectManifest,
    ProjectSnapshot,
    SelectionSet,
)
from katana_next.infrastructure.catalog import ArtifactCatalog
from katana_next.infrastructure.fingerprints import fingerprint_file, fingerprint_value
from katana_next.infrastructure.manifest_store import ManifestStore
from katana_next.infrastructure.result_store import ResultStore


@dataclass(frozen=True)
class ProjectContext:
    manifest_path: Path
    manifest: ProjectManifest
    snapshot: ProjectSnapshot
    catalog: ArtifactCatalog
    results: ResultStore


class ProjectService:
    def __init__(self, application_root: str | Path | None = None) -> None:
        self.application_root = (
            Path(application_root).resolve()
            if application_root is not None
            else Path(__file__).resolve().parents[2]
        )
        self.manifests = ManifestStore()
        self._projects: dict[str, ProjectContext] = {}

    def open_project(self, manifest_path: str | Path) -> ProjectContext:
        path, manifest = self.manifests.load(manifest_path)
        if manifest.metadata.get("gene_filter_identity"):
            from katana_next.science.gene_filter import gene_filter_details
            if gene_filter_details(manifest).get("identity") != manifest.metadata["gene_filter_identity"]:
                raise ValueError("Gene-filter annotation or expression source changed. Reconfigure exclusions from the original project.")
        self._validate_output_root(manifest.output_root)
        required_missing = [
            asset.path
            for asset in manifest.assets
            if asset.required and not asset.path.is_file()
        ]
        for sample in manifest.samples:
            required_missing.extend(
                path for path in sample.assets.values() if not path.is_file()
            )
        if required_missing:
            formatted = ", ".join(str(item) for item in required_missing[:5])
            suffix = (
                ""
                if len(required_missing) <= 5
                else f" (+{len(required_missing) - 5} more)"
            )
            raise FileNotFoundError(
                f"Required project input is missing: {formatted}{suffix}"
            )

        input_paths = self._input_paths(manifest)
        fingerprints = tuple(fingerprint_file(input_path) for input_path in input_paths)
        warnings = tuple(
            f"Optional input is missing: {item.path}"
            for item in manifest.assets
            if not item.required and not item.path.is_file()
        )
        available_kinds = {item.kind for item in manifest.assets if item.path.is_file()}
        warnings += tuple(
            f"{label} is not configured. Attach a species-matched resource bundle explicitly; no database is downloaded when opening a project."
            for kind, label in (("gene_sets_gmt", "GO / KEGG enrichment"),
                                ("rbp_motif_database", "RBP sequence motifs"))
            if kind not in available_kinds
        )
        if "rbp_motif_database" in available_kinds and "fimo_runtime" not in available_kinds:
            warnings += ("RBP motif enrichment needs a configured FIMO runtime; relative motif scores do not replace calibrated enrichment.",)
        manifest_file_fingerprint = fingerprint_file(path, force_full=True)
        manifest_fingerprint = fingerprint_value(
            {
                "manifest": manifest.model_dump(mode="json"),
                "manifest_file": manifest_file_fingerprint.model_dump(mode="json"),
                "inputs": [item.model_dump(mode="json") for item in fingerprints],
            }
        )

        manifest.output_root.mkdir(parents=True, exist_ok=True)
        snapshot = ProjectSnapshot(
            project_id=manifest.project_id,
            manifest_path=path,
            name=manifest.name,
            input_root=manifest.input_root,
            output_root=manifest.output_root,
            genome_build=manifest.annotation.genome_build,
            manifest_fingerprint=manifest_fingerprint,
            input_fingerprints=fingerprints,
            sample_count=len(manifest.samples),
            contrast_count=len(manifest.contrasts),
            warnings=warnings,
        )
        catalog = ArtifactCatalog(manifest.output_root)
        catalog.upsert_project(snapshot)
        context = ProjectContext(
            manifest_path=path,
            manifest=manifest,
            snapshot=snapshot,
            catalog=catalog,
            results=ResultStore(manifest.output_root),
        )
        self._projects[manifest.project_id] = context
        return context

    def get(self, project_id: str) -> ProjectContext:
        try:
            return self._projects[project_id]
        except KeyError as exc:
            raise KeyError(
                f"Project is not open in this session: {project_id}"
            ) from exc

    def list_artifacts(
        self,
        project_id: str,
        *,
        analysis_id: str | None = None,
        kind: str | None = None,
    ) -> list[AnalysisArtifact]:
        context = self.get(project_id)
        return context.catalog.list_artifacts(
            project_id, analysis_id=analysis_id, kind=kind
        )

    def save_selection(self, selection: SelectionSet) -> SelectionSet:
        context = self.get(selection.project_id)
        artifact = context.catalog.get_artifact(selection.source_artifact_id)
        if artifact is None or artifact.project_id != selection.project_id:
            raise ValueError(
                "Selection source artifact does not belong to the project."
            )
        if (
            artifact.metadata.get("project_snapshot_fingerprint")
            != context.snapshot.manifest_fingerprint
        ):
            raise ValueError(
                "Selection source artifact is stale for the current project inputs."
            )
        from katana_next.workflows.sources import validate_source
        import pyarrow.dataset as ds
        validate_source(context, artifact)
        if artifact.kind != "table" or selection.entity_type != artifact.entity_type:
            raise ValueError("Selection must retain its source table entity type")
        if not selection.entity_ids or len(set(selection.entity_ids)) != len(selection.entity_ids):
            raise ValueError("Selection requires unique nonempty row identities")
        dataset = ds.dataset(artifact.path, format="parquet")
        if "row_id" not in dataset.schema.names:
            raise ValueError("Source table has no stable row identity")
        selected = dataset.to_table(columns=["row_id"], filter=ds.field("row_id").cast("string").isin(list(selection.entity_ids)))
        found = [str(x) for x in selected.column("row_id").to_pylist()]
        if len(found) != len(selection.entity_ids) or set(found) != set(selection.entity_ids):
            raise ValueError("Selected row IDs are missing or duplicated in their original source")
        return context.catalog.save_selection(selection)

    def _validate_output_root(self, output_root: Path) -> None:
        resolved = output_root.resolve()
        try:
            resolved.relative_to(self.application_root)
        except ValueError:
            return
        raise ValueError(
            f"Output root cannot be inside the Katana source tree: {resolved}. "
            "Choose a project-specific results directory."
        )

    @staticmethod
    def _input_paths(manifest: ProjectManifest) -> list[Path]:
        paths: list[Path] = [asset.path for asset in manifest.assets]
        for sample in manifest.samples:
            paths.extend(sample.assets.values())
        annotation = manifest.annotation
        paths.extend(
            item
            for item in (
                annotation.gtf,
                annotation.genome_fasta,
                annotation.transcript_fasta,
            )
            if item is not None
        )
        return list(dict.fromkeys(path.resolve() for path in paths))
