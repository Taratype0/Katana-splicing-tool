"""Read-only command line entry for exchange-contract validation."""
import argparse
import json

from .validation import validate_file
from .sample_identity import validate_samplesheets


def main(argv=None):
    parser = argparse.ArgumentParser(description="Validate Katana exchange contracts without creating/opening a project or starting analysis.")
    commands = parser.add_subparsers(dest="command", required=True)
    validate = commands.add_parser("validate")
    validate.add_argument("manifest")
    validate.add_argument("--check-files", action="store_true", help="Check only local assets declared present; no remote fetch.")
    validate.add_argument("--verify-hashes", action="store_true", help="Also verify local size/SHA256; implies --check-files.")
    sheets = commands.add_parser("samplesheets", help="Check upstream sample IDs against explicit experiment metadata; nf-core retains processing ownership.")
    sheets.add_argument("--input", required=True)
    sheets.add_argument("--metadata", required=True)
    sheets.add_argument("--format", choices=["fastq", "bam", "sample_ids"], default="fastq")
    sheets.add_argument("--assay", choices=["bulk_rna", "single_cell_rna", "single_nucleus_rna"], default="bulk_rna")
    sheets.add_argument("--project", help="Optionally compare identity with an existing T09 project; never rewrite it.")
    args = parser.parse_args(argv)
    if args.command == "samplesheets":
        report = validate_samplesheets(args.input, args.metadata, input_format=args.format, assay_mode=args.assay, project_path=args.project)
    else:
        report = validate_file(args.manifest, check_files=args.check_files, verify_hashes=args.verify_hashes)
    print(json.dumps(report, ensure_ascii=False, allow_nan=False))
    return 0 if report["valid"] else 2
