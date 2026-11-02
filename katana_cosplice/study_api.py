"""Whole-study CoSplice through the same workflow used by Nextflow."""
from os import PathLike


class CoSpliceStudy:
    workflow_id = 'programs.cosplice'

    def __init__(self, project):
        from katana_next.client import KatanaProject
        self._owns_project = isinstance(project, (str, PathLike))
        self.project = KatanaProject(project) if self._owns_project else project

    def plan(self, parameters=None):
        from katana_cosplice.study_design import project_study_design
        from katana_next.workflows.cosplice_study import CoSpliceWorkflow
        parameters = CoSpliceWorkflow().normalize_parameters(parameters)
        design, contrasts = project_study_design(self.project.context.manifest, paired=parameters['paired'], association_contrast_ids=parameters.get('association_contrast_ids'))
        return {**self.project.plan(self.workflow_id, parameters), 'study_design': design, 'contrasts': contrasts}

    def run(self, parameters=None, *, force=False, timeout=None):
        result = self.plan(parameters)
        if not result['ready']:
            raise ValueError('; '.join(str(e) for e in result['errors']))
        from .study_runtime import check_runtime
        check_runtime()
        return self.project.run(self.workflow_id, parameters or {}, force=force, timeout=timeout)

    def artifacts(self):
        return self.project.artifacts(self.workflow_id)

    def close(self):
        if self._owns_project:
            self.project.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()
