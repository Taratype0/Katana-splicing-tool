"""Read-only exploration of imported complete report tables."""
import json
import duckdb
from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse
from pathlib import Path
from pydantic import BaseModel, Field

class ImportReport(BaseModel):
    report_path: str
    output_root: str

class BrowseReport(BaseModel):
    comparison: str = ''
    identifiers: list[str] = Field(default_factory=list)
    search: str = ''
    fdr_column: str = ''
    fdr_max: float | None = Field(default=None, ge=0, le=1)
    effect_column: str = ''
    effect_min: float | None = Field(default=None, ge=0)
    offset: int = Field(default=0, ge=0)
    limit: int = Field(default=100, ge=1, le=1000)

def create_report_router(runtime):
    router = APIRouter()

    @router.get('/api/report-workspace/version')
    def version():
        return {'release':'report-workspace-20261108','native_nextflow_execution_verified':False}

    @router.post('/api/projects/import-report')
    def load(payload: ImportReport):
        from katana_next.services.report_import import import_report
        from katana_next.api.app import _project_payload
        try:
            context = import_report(payload.report_path, payload.output_root, runtime.projects)
            runtime.jobs.bind_project(context.manifest.project_id, context.manifest.output_root)
            return _project_payload(context)
        except (ValueError, OSError, KeyError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.post('/api/projects/{project_id}/report-tables/{artifact_id}')
    def browse(project_id: str, artifact_id: str, payload: BrowseReport):
        try:
            context = runtime.projects.get(project_id)
            artifact = context.catalog.get_artifact(artifact_id)
            if not artifact or artifact.project_id != project_id or artifact.analysis_id != 'imported.report':
                raise ValueError('Imported table not found')
            def q(name):
                if name not in columns:
                    raise ValueError(f'Unknown column: {name}')
                return '"' + name.replace('"', '""') + '"'
            with duckdb.connect() as db:
                db.read_parquet(str(artifact.path)).create_view('evidence')
                columns = [r[0] for r in db.execute('DESCRIBE evidence').fetchall()]
                clauses, values = [], []
                if payload.comparison:
                    clauses.append(q('comparison_id') + ' = ?'); values.append(payload.comparison)
                if payload.identifiers:
                    names = [c for c in ('gene_id','gene_symbol','event_id') if c in columns]
                    if not names:
                        raise ValueError('This table has no gene/event identifiers')
                    clauses.append('(' + ' OR '.join(f'CAST({q(c)} AS VARCHAR) IN (SELECT unnest(?))' for c in names) + ')')
                    values.extend([payload.identifiers] * len(names))
                if payload.search:
                    clauses.append('(' + ' OR '.join(f'contains(lower(CAST({q(c)} AS VARCHAR)), lower(?))' for c in columns) + ')')
                    values.extend([payload.search] * len(columns))
                for column, bound, operator in [(payload.fdr_column,payload.fdr_max,'<='),(payload.effect_column,payload.effect_min,'>=')]:
                    if bound is not None:
                        numeric = f'TRY_CAST({q(column)} AS DOUBLE)'
                        if operator == '>=':
                            numeric = f'abs({numeric})'
                        clauses.append(f'{numeric} {operator} ?'); values.append(bound)
                where = ' WHERE ' + ' AND '.join(clauses) if clauses else ''
                total = db.execute('SELECT count(*) FROM evidence' + where, values).fetchone()[0]
                frame = db.execute('SELECT * FROM evidence' + where + ' ORDER BY row_id LIMIT ? OFFSET ?', values + [payload.limit,payload.offset]).fetchdf()
                return {'columns': columns, 'rows': json.loads(frame.to_json(orient='records', double_precision=15)),
                        'total': total, 'offset': payload.offset, 'limit': payload.limit}
        except (ValueError, OSError, KeyError, duckdb.Error) as exc:
            raise HTTPException(422, str(exc)) from exc
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
    return router
