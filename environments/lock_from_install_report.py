"""Lock the exact additional pip distributions from an actual clean install.

The Conda explicit file locks the native base. This companion records pip's
download hashes, excluding the locally built Katana wheel installed separately.
"""
import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    report = json.loads(args.report.read_text(encoding="utf-8"))
    rows = []
    for item in report["install"]:
        name = item["metadata"]["name"].lower().replace("_", "-")
        if name == "katana-splicing-tool":
            continue
        digest = item.get("download_info", {}).get("archive_info", {}).get("hashes", {}).get("sha256")
        if not digest or len(digest) != 64:
            raise ValueError(f"Missing actual download SHA256: {name}")
        rows.append(f"{name}=={item['metadata']['version']} --hash=sha256:{digest}")
    if not rows or len(rows) != len(set(rows)):
        raise ValueError("Empty or duplicate package installation evidence")
    args.destination.parent.mkdir(parents=True, exist_ok=True)
    args.destination.write_text(
        "# Generated from a successful clean pip installation report.\n"
        "# Apply after the matching Conda lock, with --no-deps --require-hashes --no-build-isolation.\n"
        "# Source archives retain their upstream hash; native builds need the locked compiler base.\n"
        + "\n".join(sorted(rows)) + "\n", encoding="utf-8")
    print(f"Locked {len(rows)} additional pip distributions")


if __name__ == "__main__":
    main()
