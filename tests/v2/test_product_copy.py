from katana_next.workflows.registry import WorkflowRegistry
from katana_next.workflows.catalog import WORKSPACES
from katana_next.workflows.presentation import WORKFLOW_COPY


def test_all_modules_have_concrete_purpose_and_output():
    registry = WorkflowRegistry()
    assert set(WORKFLOW_COPY) == {w.definition.workflow_id for w in registry.all()}
    for workflow in registry.all():
        data = workflow.definition.to_dict()
        assert "?" not in data["question"]
        assert data["label"] and data["output_summary"]
        assert data["workflow_id"] == workflow.definition.workflow_id
        assert data["analysis_unit"] == workflow.definition.analysis_unit
        assert [p["default"] for p in data["parameters"]] == [p.default for p in workflow.definition.parameters]
        assert [p["key"] for p in data["parameters"]] == [p.key for p in workflow.definition.parameters]


def test_navigation_is_analysis_specific_without_changing_workspace_ids():
    labels = {w.workspace_id: w.label for w in WORKSPACES}
    assert labels["prioritize"] == "Functional Analysis"
    assert labels["programs"] == "Splicing Networks"
    assert labels["differential"] == "Differential Analysis"


def test_primary_rbp_and_single_cell_addon_remain_available():
    registry = WorkflowRegistry()
    assert registry.get("prioritize.rbp_mechanism").definition.to_dict()["presentation_role"] == "primary"
    assert "Add-on" in registry.get("single_cell.cell_embedding").definition.to_dict()["label"]
