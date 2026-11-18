"""Read-only exploration of imported complete report tables."""
import json
import duckdb
from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse
from pathlib import Path
from pydantic import BaseModel, Field
from typing import Literal
from katana_next.services.report_query import result_table, browse as browse_table, plot, select_rows, query_table, clauses_for

class ImportReport(BaseModel):
    report_path: str
    output_root: str = ''

class BrowseReport(BaseModel):
    comparison: str = ''
    identifiers: list[str] = Field(default_factory=list)
    search: str = ''
    fdr_column: str = ''
    fdr_max: float | None = Field(default=None, ge=0, le=1)
    effect_column: str = ''
    effect_min: float | None = Field(default=None, ge=0)
    fdr_min: float | None = Field(default=None, ge=0, le=1)
    effect_low: float | None = None
    effect_high: float | None = None
    direction: Literal['both', 'positive', 'negative'] = 'both'
    scopes: dict[str, str] = Field(default_factory=dict)
    row_ids: list[str] | None = None
    sort_column: str = ''
    sort_desc: bool = False
    offset: int = Field(default=0, ge=0)
    limit: int = Field(default=100, ge=1, le=1000)

class EnrichReport(BaseModel):
    filters: BrowseReport
    library_path: str
    library_species: Literal['human', 'mouse', 'other']
    gene_column: str = 'gene_id'
    test_column: str
    minimum_term_size: int = Field(default=1, ge=1)
    maximum_term_size: int = Field(default=2000, ge=1)

class ContrastRescuePayload(BaseModel):
    contrast_a: str
    contrast_b: str
    min_abs_dpsi: float = Field(default=0.10, ge=0.0, le=1.0)
    max_fdr: float = Field(default=0.05, ge=0.0, le=1.0)
    relation: Literal['reversed', 'concordant', 'specific_to_a'] = 'reversed'
    rs_only: bool = True

class GeneCardPreviewPayload(BaseModel):
    gene_query: str = ''
    event_id: str = ''
    comparison: str = ''
    custom_sequence: str = ''
    custom_motif: str = ''

class PsiLandscapePayload(BaseModel):
    comparison: str = ''
    mode: Literal['unsupervised', 'supervised'] = 'unsupervised'
    top_n: int = Field(default=150, ge=15, le=500)
    manifold_top_n: int = Field(default=2000, ge=50, le=10000)
    prior_kappa: float = Field(default=10.0, ge=0.5, le=100.0)
    z_clip: float = Field(default=2.5, ge=1.0, le=5.0)
    metric: Literal['euclidean', 'cosine', 'correlation'] = 'euclidean'
    event_types: list[str] = Field(default_factory=list)
    umap_neighbors: int = Field(default=12, ge=3, le=50)
    umap_min_dist: float = Field(default=0.18, ge=0.01, le=0.95)

class LiveOraPayload(BaseModel):
    comparison: str = ''
    preset: str = 'ps_das'
    custom_genes: list[str] = Field(default_factory=list)
    fdr_cutoff: float = Field(default=0.05, ge=0.0, le=1.0)
    effect_cutoff: float = Field(default=0.10, ge=0.0, le=5.0)
    min_term_size: int = Field(default=5, ge=2, le=500)
    max_term_size: int = Field(default=1500, ge=10, le=5000)

class LiveMotifScanPayload(BaseModel):
    comparison: str = ''
    query_motif_or_rbp: str = 'UGGACU'
    custom_sequence: str = ''

class LiveTargetRankerPayload(BaseModel):
    comparison: str = ''
    w_dpsi: float = Field(default=1.0, ge=0.0, le=5.0)
    w_isoform: float = Field(default=1.0, ge=0.0, le=5.0)
    w_splicenet: float = Field(default=1.0, ge=0.0, le=5.0)
    w_motif: float = Field(default=1.0, ge=0.0, le=5.0)
    w_proteome: float = Field(default=1.5, ge=0.0, le=5.0)
    min_abs_dpsi: float = Field(default=0.10, ge=0.0, le=1.0)
    max_fdr: float = Field(default=0.05, ge=0.0, le=1.0)
    deg_independent_only: bool = False
    require_domain_or_nmd: bool = False


def create_report_router(runtime):
    router = APIRouter()

    @router.get('/api/report-workspace/version')
    def version():
        return {'release':'katana-v2.7.0', 'application_root':str(Path(__file__).resolve().parents[2]),
                'base_release':'local-app-rebuild-20261113', 'base_native_nextflow_execution_verified':True,
                'base_verified_executor':'Linux local', 'slurm_execution_verified':False,
                'acceptance_run':'https://github.com/Taratype0/Katana-splicing-tool/actions/runs/42395577378'}

    @router.post('/api/projects/import-report')
    def load(payload: ImportReport):
        from katana_next.services.report_import import import_report, member, sha
        from katana_next.api.app import _project_payload
        from katana_next.api.local_files import workspace_home
        from uuid import uuid4
        try:
            from katana_next.services.report_archive import report_source
            source = report_source(payload.report_path, workspace_home())
            recorded = json.loads((source / 'report-manifest.json').read_text('utf-8-sig'))
            try:
                existing = runtime.projects.get(recorded['project_id'])
            except KeyError:
                existing = None
            if existing and not payload.output_root.strip():
                metadata = existing.manifest.metadata
                if not metadata.get('report_workspace'):
                    raise ValueError('This project is already open from a different source. Open its existing workspace instead.')
                if metadata['source_report_sha256'] != sha(source / 'report-manifest.json'):
                    raise ValueError('The source report has changed since import')
                receipt = existing.manifest_path.parent / 'report-import.json'
                if not receipt.is_file():
                    raise ValueError('Previous report import did not complete. Import into a new workspace.')
                for artifact in existing.catalog.list_artifacts(existing.manifest.project_id):
                    if artifact.analysis_id == 'imported.report':
                        src_member = member(source, artifact.metadata['source_relative_path'])
                        if not artifact.path.is_file() or not src_member.is_file():
                            raise ValueError('Report evidence has changed since import')
                        expected_src = artifact.metadata.get('source_table_sha256')
                        expected_imp = artifact.metadata.get('imported_table_sha256')
                        if (expected_src and sha(src_member) != expected_src) or (expected_imp and sha(artifact.path) != expected_imp):
                            raise ValueError('Report evidence has changed since import')
                return _project_payload(existing)
            destination = payload.output_root.strip() or str(workspace_home() / f'report-{uuid4().hex[:12]}')
            context = import_report(source, destination, runtime.projects)
            runtime.jobs.bind_project(context.manifest.project_id, context.manifest.output_root)
            return _project_payload(context)
        except (ValueError, OSError, KeyError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post('/api/projects/{project_id}/report-tables/{artifact_id}')
    def browse(project_id: str, artifact_id: str, payload: BrowseReport):
        try:
            context = runtime.projects.get(project_id)
            return browse_table(result_table(context, artifact_id), payload)
        except (ValueError, OSError, KeyError, duckdb.Error) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.get('/api/projects/{project_id}/report-views')
    def report_views(project_id: str):
        from katana_next.services.report_import import member
        try:
            context = runtime.projects.get(project_id)
            root = Path(context.manifest.metadata['source_report_directory']).resolve()
            manifest = json.loads(member(root, 'report-manifest.json').read_text('utf-8-sig'))
            views = []
            for item in manifest.get('network_browsers', []):
                relative = item.get('source', '')
                member(root, relative)
                views.append({'title': item.get('title', 'Module network'), 'source': relative, 'workspace': 'programs'})
            return views
        except (ValueError, OSError, KeyError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.get('/api/projects/{project_id}/report-presentation')
    def presentation(project_id: str):
        from katana_next.services.report_presentation import presentation
        try: return presentation(runtime.projects.get(project_id))
        except (ValueError, OSError, KeyError) as exc: raise HTTPException(422,str(exc)) from exc

    @router.get('/api/projects/{project_id}/report-figures/{panel_id}')
    def figure(project_id: str,panel_id: str):
        from katana_next.services.report_presentation import figure
        try: return figure(runtime.projects.get(project_id),panel_id)
        except (ValueError, OSError, KeyError) as exc: raise HTTPException(422,str(exc)) from exc

    @router.get('/api/projects/{project_id}/report-network')
    def network(project_id: str):
        from katana_next.services.mckinsey_interactive import unified_network
        try: return unified_network(runtime.projects.get(project_id))
        except (ValueError, OSError, KeyError) as exc: raise HTTPException(422,str(exc)) from exc

    @router.get('/api/projects/{project_id}/splicetrace-workbench')
    def splicetrace_view(project_id: str, comparison: str = ''):
        from katana_next.services.mckinsey_interactive import splicetrace_workbench
        try: return splicetrace_workbench(runtime.projects.get(project_id), comparison=comparison)
        except (ValueError, OSError, KeyError) as exc: raise HTTPException(422,str(exc)) from exc

    @router.get('/api/projects/{project_id}/quadrant-explorer')
    def quadrant_view(project_id: str, comparison: str = ''):
        from katana_next.services.mckinsey_interactive import quadrant_data
        try: return quadrant_data(runtime.projects.get(project_id), comparison=comparison)
        except (ValueError, OSError, KeyError) as exc: raise HTTPException(422,str(exc)) from exc

    @router.post('/api/projects/{project_id}/contrast-rescue')
    def rescue_view(project_id: str, payload: ContrastRescuePayload):
        from katana_next.services.mckinsey_interactive import contrast_rescue_filter
        try: return contrast_rescue_filter(runtime.projects.get(project_id), **payload.model_dump())
        except (ValueError, OSError, KeyError) as exc: raise HTTPException(422,str(exc)) from exc

    @router.get('/api/projects/{project_id}/single-cell-manifold')
    def sc_manifold_view(project_id: str, feature: str = ''):
        from katana_next.services.mckinsey_interactive import single_cell_manifold
        try: return single_cell_manifold(runtime.projects.get(project_id), feature=feature)
        except (ValueError, OSError, KeyError) as exc: raise HTTPException(422,str(exc)) from exc

    @router.post('/api/projects/{project_id}/genecard-preview')
    def genecard_view(project_id: str, payload: GeneCardPreviewPayload):
        from katana_next.services.mckinsey_interactive import genecard_preview
        try: return genecard_preview(runtime.projects.get(project_id), **payload.model_dump())
        except (ValueError, OSError, KeyError) as exc: raise HTTPException(422,str(exc)) from exc

    @router.post('/api/projects/{project_id}/compute/psi-landscape')
    def live_psi_landscape_view(project_id: str, payload: PsiLandscapePayload):
        from katana_next.services.mckinsey_compute import compute_live_psi_landscape
        try: return compute_live_psi_landscape(runtime.projects.get(project_id), **payload.model_dump())
        except (ValueError, OSError, KeyError) as exc: raise HTTPException(422, str(exc)) from exc

    @router.post('/api/projects/{project_id}/compute/live-ora')
    def live_ora_view(project_id: str, payload: LiveOraPayload):
        from katana_next.services.mckinsey_compute import compute_live_go_kegg_ora
        try: return compute_live_go_kegg_ora(runtime.projects.get(project_id), **payload.model_dump())
        except (ValueError, OSError, KeyError) as exc: raise HTTPException(422, str(exc)) from exc

    @router.post('/api/projects/{project_id}/compute/motif-scan')
    def live_motif_scan_view(project_id: str, payload: LiveMotifScanPayload):
        from katana_next.services.mckinsey_compute import compute_live_motif_scan
        try: return compute_live_motif_scan(runtime.projects.get(project_id), **payload.model_dump())
        except (ValueError, OSError, KeyError) as exc: raise HTTPException(422, str(exc)) from exc

    @router.post('/api/projects/{project_id}/compute/target-ranker')
    def live_target_ranker_view(project_id: str, payload: LiveTargetRankerPayload):
        from katana_next.services.mckinsey_compute import compute_live_target_ranker
        try: return compute_live_target_ranker(runtime.projects.get(project_id), **payload.model_dump())
        except (ValueError, OSError, KeyError) as exc: raise HTTPException(422, str(exc)) from exc

    @router.get('/api/projects/{project_id}/live-compute/{run_id}/{filename}')
    def live_compute_file(project_id: str, run_id: str, filename: str):
        import re
        allowed = {
            'psi_clustermap.png', 'psi_clustermap.svg', 'psi_pca.png', 'psi_pca.svg',
            'dual_psi_umap.png', 'dual_psi_umap.svg', 'landscape_coordinates.tsv',
            'live_go_kegg_ora.png', 'live_go_kegg_ora.svg', 'live_go_kegg_complete.tsv',
            'foreground_genes.tsv', 'background_genes.tsv',
        }
        if not re.fullmatch('[a-f0-9]{32}', run_id) or filename not in allowed:
            raise HTTPException(404, 'Unknown live compute artifact')
        fpath = runtime.projects.get(project_id).manifest.output_root / 'live-compute' / run_id / filename
        if not fpath.is_file():
            raise HTTPException(404, 'Live compute artifact not found')
        return FileResponse(fpath, filename=filename)
    @router.post('/api/projects/{project_id}/report-tables/{artifact_id}/plot')
    def plot_results(project_id: str, artifact_id: str, payload: BrowseReport):
        try:
            return plot(result_table(runtime.projects.get(project_id), artifact_id), payload)
        except (ValueError, OSError, KeyError, duckdb.Error) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post('/api/projects/{project_id}/report-tables/{artifact_id}/selection')
    def selection(project_id: str, artifact_id: str, payload: BrowseReport):
        try:
            return select_rows(result_table(runtime.projects.get(project_id), artifact_id), payload)
        except (ValueError, OSError, KeyError, duckdb.Error) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post('/api/projects/{project_id}/report-tables/{artifact_id}/export')
    def export_rows(project_id: str, artifact_id: str, payload: BrowseReport):
        from uuid import uuid4
        from starlette.background import BackgroundTask
        try:
            context = runtime.projects.get(project_id)
            artifact = result_table(context, artifact_id)
            destination = context.manifest.output_root / 'browse-exports'
            destination.mkdir(parents=True, exist_ok=True)
            path = destination / (uuid4().hex + '.tsv')
            with query_table(artifact) as (db, columns):
                where, values = clauses_for(payload, columns)
                # DuckDB streams to disk; the full result never becomes a browser preview.
                db.execute('CREATE TEMP TABLE export_rows AS SELECT * FROM evidence' + where, values)
                db.sql('SELECT * FROM export_rows').write_csv(str(path), sep='\t', header=True)
            return FileResponse(path, filename='filtered-results.tsv', background=BackgroundTask(path.unlink, missing_ok=True))
        except (ValueError, OSError, KeyError, duckdb.Error) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post('/api/projects/{project_id}/report-tables/{artifact_id}/enrichment')
    def enrichment(project_id: str, artifact_id: str, payload: EnrichReport):
        from katana_next.services.report_enrichment import enrich
        try:
            context = runtime.projects.get(project_id)
            result = enrich(context, result_table(context, artifact_id), **payload.model_dump(exclude={'filters'}), filters=payload.filters)
            result['files'] = {name: f'/api/projects/{project_id}/local-enrichment/{result["audit"]["run_id"]}/{name}'
                               for name in ('enrichment.tsv', 'background.tsv', 'foreground.tsv', 'audit.json')}
            return result
        except (ValueError, OSError, KeyError, duckdb.Error) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.get('/api/projects/{project_id}/local-enrichment/{run_id}/{filename}')
    def enrichment_file(project_id: str, run_id: str, filename: str):
        import re
        if not re.fullmatch('[a-f0-9]{32}', run_id) or filename not in {'enrichment.tsv', 'background.tsv', 'foreground.tsv', 'audit.json'}:
            raise HTTPException(404, 'Unknown enrichment output')
        context = runtime.projects.get(project_id)
        return FileResponse(context.manifest.output_root / 'local-enrichment' / run_id / filename, filename=filename)

    @router.get('/api/projects/{project_id}/report-source/{asset_path:path}')
    def report_source(project_id: str, asset_path: str):
        from katana_next.services.report_import import member
        try:
            context = runtime.projects.get(project_id)
            if not context.manifest.metadata.get('report_workspace'):
                raise ValueError('No imported report in this workspace')
            root = Path(context.manifest.metadata['source_report_directory']).resolve()
            return FileResponse(member(root, asset_path))
        except (KeyError, ValueError, OSError) as exc:
            raise HTTPException(404, str(exc)) from exc

    @router.post('/api/projects/{project_id}/report-compact')
    def compact_report(project_id: str):
        from uuid import uuid4
        from katana_next.services.report_archive import compact_copy
        try:
            context = runtime.projects.get(project_id)
            source = context.manifest.metadata.get('source_report_directory')
            if not source:
                raise ValueError('Import a completed report before creating a compact copy')
            run_id = uuid4().hex
            destination = context.manifest.output_root / 'compact-reports' / run_id
            compact_copy(source, destination)
            return {'directory': str(destination), 'archive_url': f'/api/projects/{project_id}/compact-reports/{run_id}'}
        except (ValueError, OSError, KeyError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.get('/api/projects/{project_id}/compact-reports/{run_id}')
    def compact_archive(project_id: str, run_id: str):
        import re
        if not re.fullmatch('[a-f0-9]{32}', run_id):
            raise HTTPException(404, 'Unknown report')
        path = runtime.projects.get(project_id).manifest.output_root / 'compact-reports' / run_id / 'REPORT.zip'
        if not path.is_file():
            raise HTTPException(404, 'Report archive is not ready')
        return FileResponse(path, filename='REPORT.zip')
    return router
