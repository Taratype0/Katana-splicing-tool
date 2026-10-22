from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CrossComparisonPageContract:
    subbranch: str
    title: str
    run_label: str
    result_level: str
    needs_pair_tabs: bool
    action: str


_CONTRACTS = {
    "matrix": CrossComparisonPageContract(
        "matrix",
        "5.4.1 Pairwise Candidate Matrix",
        "Build Pairwise Candidate Matrix",
        "gene",
        True,
        "candidate_matrix",
    ),
    "shared": CrossComparisonPageContract(
        "shared",
        "5.4.2 Gene-level Candidate Pattern",
        "Build Pairwise Gene Patterns",
        "gene",
        True,
        "candidate_matrix",
    ),
    "interaction": CrossComparisonPageContract(
        "interaction",
        "5.4.3 Replicate-aware AS Interaction",
        "Apply / Recalculate 5.4.3",
        "event",
        True,
        "interaction",
    ),
    "venn": CrossComparisonPageContract(
        "venn",
        "5.4.4 Strict-event Venn / Overlap",
        "Run 5.4.4 Event Venn",
        "event",
        True,
        "venn",
    ),
    "venn_multi": CrossComparisonPageContract(
        "venn_multi",
        "5.4.5 Selected-comparison Event Venn",
        "Run 5.4.5 Selected Venn",
        "event",
        False,
        "venn_multi",
    ),
    "venn_multi_similarity": CrossComparisonPageContract(
        "venn_multi_similarity",
        "5.4.6 Selected-comparison Event Similarity",
        "Run 5.4.6 Similarity Heatmap",
        "event",
        False,
        "venn_multi",
    ),
    "pattern": CrossComparisonPageContract(
        "pattern",
        "Legacy Pairwise dPSI Pattern Table",
        "Apply / Recalculate legacy dPSI table",
        "event",
        True,
        "pattern",
    ),
    "reversal": CrossComparisonPageContract(
        "reversal",
        "Legacy Gene-level Direction Change",
        "Build Gene-level Direction Change View",
        "gene",
        True,
        "candidate_matrix",
    ),
    "comparison": CrossComparisonPageContract(
        "comparison",
        "Legacy Group Comparison",
        "Run Legacy Group Comparison",
        "gene",
        True,
        "program_comparison",
    ),
    "heatmap": CrossComparisonPageContract(
        "heatmap",
        "Legacy Program Heatmap",
        "Run Legacy Program Heatmap",
        "gene",
        True,
        "program_comparison",
    ),
    "direction": CrossComparisonPageContract(
        "direction",
        "Legacy Direction vs Strength",
        "Run Legacy Direction vs Strength",
        "gene",
        True,
        "program_comparison",
    ),
}


class CrossComparisonController:
    """Central 5.4 UI contract and independent action routing."""

    @staticmethod
    def contract(subbranch: str | None) -> CrossComparisonPageContract:
        return _CONTRACTS.get(str(subbranch or ""), _CONTRACTS["matrix"])

    @classmethod
    def run_label(cls, subbranch: str | None) -> str:
        return cls.contract(subbranch).run_label

    @classmethod
    def needs_pair_tabs(cls, subbranch: str | None, *, interaction_view: str = "") -> bool:
        contract = cls.contract(subbranch)
        if contract.subbranch == "interaction" and interaction_view == "batch_summary":
            return False
        return contract.needs_pair_tabs

    @classmethod
    def action(cls, subbranch: str | None) -> str:
        return cls.contract(subbranch).action
