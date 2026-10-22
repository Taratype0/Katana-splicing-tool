"""Read-only command line entry for exchange-contract validation."""
import argparse
import json

from .validation import validate_file


def main(argv=None):
    parser = argparse.ArgumentParser(description="Validate Katana exchange contracts without creating/opening a project or starting analysis.")
    commands = parser.add_subparsers(dest="command", required=True)
    validate = commands.add_parser("validate")
    validate.add_argument("manifest")
    validate.add_argument("--check-files", action="store_true", help="Check only local assets declared present; no remote fetch.")
    validate.add_argument("--verify-hashes", action="store_true", help="Also verify local size/SHA256; implies --check-files.")
    args = parser.parse_args(argv)
    report = validate_file(args.manifest, check_files=args.check_files, verify_hashes=args.verify_hashes)
    print(json.dumps(report, ensure_ascii=False, allow_nan=False))
    return 0 if report["valid"] else 2
