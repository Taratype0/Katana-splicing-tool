from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.domain.models.project import IsoformSampleDefinition


class IsoformSwitchAdapter:
    def build_import_manifest(self, quant_files: list[Path], output_dir: Path) -> Path:
        output_dir.mkdir(parents=True, exist_ok=True)
        out = output_dir / "isoform_quant_manifest.tsv"
        self.build_import_manifest_frame(quant_files).to_csv(out, sep="\t", index=False)
        return out

    def build_import_manifest_from_samples(self, samples: list[IsoformSampleDefinition], output_dir: Path) -> Path:
        output_dir.mkdir(parents=True, exist_ok=True)
        out = output_dir / "isoform_quant_manifest.tsv"
        self.build_import_manifest_frame(samples).to_csv(out, sep="\t", index=False)
        return out

    def build_import_manifest_frame(self, quant_inputs) -> pd.DataFrame:
        rows = []
        for item in quant_inputs:
            if isinstance(item, IsoformSampleDefinition):
                rows.append(
                    {
                        "sample_id": item.sample_id,
                        "quant_sf": str(item.quant_sf),
                        "comparison_id": item.comparison_id or "",
                        "condition": item.condition or "",
                        "batch": item.batch or "",
                    }
                )
            else:
                quant_file = Path(item)
                rows.append(
                    {
                        "sample_id": quant_file.parent.name,
                        "quant_sf": str(quant_file),
                    }
                )
        return pd.DataFrame(rows)
