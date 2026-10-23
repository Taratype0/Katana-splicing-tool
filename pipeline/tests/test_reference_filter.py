from __future__ import annotations
import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "bin/kanata_filter_analysis_gtf.py"
SPEC = importlib.util.spec_from_file_location("analysis_gtf_filter", SCRIPT)
FILTER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FILTER)


def fixture_gtf(path: Path) -> None:
    genes = [("nuclear", "chr1", "protein_coding", "Normal"),
             ("nuclear_mt_name", "chr2", "protein_coding", "MT-Nuclear"),
             ("pseudo", "chr3", "processed_pseudogene", "Pseudogene"),
             ("mt", "MT", "protein_coding", "Mito"),
             ("mt_alias", "chrMT", "Mt_rRNA", "MitoAlias"),
             ("unknown", "chr4", "", "Unknown")]
    lines = ["# fixture with complete retained and excluded feature blocks\n"]
    for gene, chrom, kind, name in genes:
        for feature in ("gene", "transcript", "exon"):
            attrs = f'gene_id "{gene}"; gene_name "{name}";'
            if kind:
                attrs += f' gene_biotype "{kind}";'
            if feature != "gene":
                attrs += f' transcript_id "{gene}.t1";'
            lines.append(f"{chrom}\tfixture\t{feature}\t1\t10\t.\t+\t.\t{attrs}\n")
    path.write_text("".join(lines), encoding="utf-8")


class ReferenceFilterTests(unittest.TestCase):
    def test_all_features_of_excluded_genes_removed_and_source_unchanged(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "raw.gtf"
            fixture_gtf(source)
            before = source.read_bytes()
            receipt = FILTER.filter_gtf(source, root / "out")
            self.assertEqual(source.read_bytes(), before)
            self.assertEqual(receipt["input_genes"], 6)
            self.assertEqual(receipt["excluded_genes"], 3)
            self.assertEqual(receipt["retained_features"], 9)
            self.assertEqual(receipt["unknown_biotype_genes"], 1)
            output = (root / "out/annotation.analysis.gtf").read_text()
            self.assertIn('gene_id "nuclear_mt_name"', output)
            self.assertIn('gene_id "unknown"', output)
            for excluded in ("pseudo", "mt", "mt_alias"):
                self.assertNotIn(f'gene_id "{excluded}";', output)
            with (root / "out/annotation_decisions.tsv").open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 6)
            self.assertEqual(json.loads((root / "out/filtering.receipt.json").read_text()), receipt)

    def test_missing_biotypes_rejected_instead_of_claiming_exclusion(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "missing.gtf"
            source.write_text('chr1\tx\tgene\t1\t10\t.\t+\t.\tgene_id "a";\n')
            with self.assertRaisesRegex(ValueError, "requires actual"):
                FILTER.filter_gtf(source, root / "out")
            self.assertFalse((root / "out").exists())

    def test_existing_output_and_conflicting_annotations_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "raw.gtf"
            fixture_gtf(source)
            output = root / "out"
            output.mkdir()
            with self.assertRaisesRegex(ValueError, "Output exists"):
                FILTER.filter_gtf(source, output)
            source.write_text(source.read_text() + 'chr9\tx\tgene\t1\t10\t.\t+\t.\tgene_id "nuclear"; gene_type "protein_coding";\n')
            with self.assertRaisesRegex(ValueError, "multiple chromosomes"):
                FILTER.filter_gtf(source, root / "other")


if __name__ == "__main__":
    unittest.main()
