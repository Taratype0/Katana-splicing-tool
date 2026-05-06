from __future__ import annotations

from pathlib import Path

from src.domain.models.project import ProjectConfig
from src.io.loaders.cards_loader import CardsLoader
from src.io.loaders.counts_loader import CountsLoader
from src.io.loaders.deg_loader import DEGLoader
from src.io.loaders.dexseq_loader import DEXSeqLoader
from src.io.loaders.dtu_loader import DTULoader
from src.io.loaders.quant_loader import QuantLoader
from src.io.loaders.rmats_loader import RMATSLoader
from src.io.loaders.sashimi_loader import SashimiLoader


class ProjectValidator:
    def __init__(self) -> None:
        self.rmats_loader = RMATSLoader()
        self.deg_loader = DEGLoader()
        self.dexseq_loader = DEXSeqLoader()
        self.dtu_loader = DTULoader()
        self.quant_loader = QuantLoader()
        self.counts_loader = CountsLoader()
        self.cards_loader = CardsLoader()
        self.sashimi_loader = SashimiLoader()

    def validate_project(self, project: ProjectConfig) -> list[str]:
        messages = list(project.scan_messages)

        for comparison in project.available_comparisons:
            if comparison.rmats_path is not None:
                rmats_post = comparison.rmats_path / "rmats_post"
                for mode in comparison.detected_rmats_modes:
                    sample = rmats_post / f"SE.MATS.{mode}.txt"
                    messages.extend(self.rmats_loader.validate(sample))
            if comparison.deg_path is not None:
                messages.extend(self.deg_loader.validate(comparison.deg_path))
            if comparison.quant_dir is not None:
                messages.extend(self.quant_loader.validate(comparison.quant_dir / "quant.sf"))

        if project.counts_path is not None:
            messages.extend(self.counts_loader.validate(project.counts_path))

        messages.extend(self._validate_first_match(project.dexseq_root, "perGeneQValue.*.csv", self.dexseq_loader.validate))
        messages.extend(self._validate_first_match(project.dtu_root, "DEXSeqResults.*.tsv", self.dtu_loader.validate_result))
        messages.extend(self._validate_first_match(project.dtu_root, "getAdjustedPValues.*.tsv", self.dtu_loader.validate_stager))
        sashimi_root = project.output_root / "06_sashimi" / "data" if project.output_root is not None else None
        shortlist_root = project.output_root / "04b_shortlist" / "data" if project.output_root is not None else None
        cards_root = project.output_root / "05_cards" / "data" if project.output_root is not None else None
        messages.extend(self._validate_first_match(sashimi_root, "rmats2sashimi_manifest.tsv", self.sashimi_loader.validate_manifest))
        messages.extend(self._validate_first_match(shortlist_root, "*_shortlist*.tsv", self.cards_loader.validate_shortlist))
        messages.extend(self._validate_first_match(cards_root, "expression_support_group_means.tsv", self.cards_loader.validate_expression_support))
        return messages

    def _validate_first_match(self, root: Path | None, pattern: str, validator) -> list[str]:
        if root is None:
            return []
        matches = list(root.rglob(pattern))
        if not matches:
            return []
        return validator(matches[0])
