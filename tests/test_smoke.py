from src.adapters.native_sashimi_adapter import NativeSashimiAdapter
from src.services.project_service import ProjectService


def test_project_service_initializes_without_project() -> None:
    service = ProjectService()

    assert service.current_project is None
    assert service.module_state("candidate_ranking")["status"] == "not_run"


def test_native_sashimi_adapter_initializes() -> None:
    adapter = NativeSashimiAdapter()

    assert adapter is not None
