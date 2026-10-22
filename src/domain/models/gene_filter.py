from __future__ import annotations

from dataclasses import asdict, dataclass, field
import math
import re


@dataclass
class GeneFilterConfig:
    """Downstream scope only: never refits upstream differential models."""

    exclude_mitochondrial: bool = False
    exclude_pseudogenes: bool = False
    exclude_rrna: bool = False
    protein_coding_only: bool = False
    excluded_genes: list[str] = field(default_factory=list)
    min_expression: float = 0.0
    min_samples: int = 1

    @property
    def active(self) -> bool:
        return any((self.exclude_mitochondrial, self.exclude_pseudogenes,
                    self.exclude_rrna, self.protein_coding_only,
                    self.excluded_genes, self.min_expression > 0))

    @classmethod
    def from_mapping(cls, values: dict | None) -> GeneFilterConfig:
        values = values or {}
        boolean_names = ("exclude_mitochondrial", "exclude_pseudogenes",
                         "exclude_rrna", "protein_coding_only")
        flags = {}
        for name in boolean_names:
            value = values.get(name, False)
            if not isinstance(value, bool):
                raise ValueError(f"{name} must be true or false.")
            flags[name] = value
        genes = values.get("excluded_genes", [])
        if isinstance(genes, str):
            genes = re.split(r"[\s,;]+", genes)
        if not isinstance(genes, (list, tuple)):
            raise ValueError("Excluded genes must be a list of gene IDs or symbols.")
        genes = sorted({str(gene).strip() for gene in genes if str(gene).strip()})
        minimum = float(values.get("min_expression", 0))
        samples = float(values.get("min_samples", 1))
        if not math.isfinite(minimum) or minimum < 0:
            raise ValueError("Minimum expression must be finite and non-negative.")
        if not math.isfinite(samples) or not samples.is_integer() or samples < 1:
            raise ValueError("Minimum samples must be a positive integer.")
        return cls(**flags, excluded_genes=genes, min_expression=minimum, min_samples=int(samples))

    def to_dict(self) -> dict:
        return asdict(self)
