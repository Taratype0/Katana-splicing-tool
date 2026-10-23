"""Build and test the unchanged upstream DEXSeq tools without publishing them.

GPL program sources are downloaded into a separate build context, not copied
into the MIT workflow or Python wheel. Every download is bound by SHA256.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]


def prepare_context(context, download=None):
    context = Path(context).resolve()
    if context.exists() or context.is_relative_to(ROOT):
        raise ValueError("Use a new build context outside the checkout")
    source = ROOT / "containers/dexseq-tools"
    lock = json.loads((source / "sources.lock.json").read_text(encoding="utf-8"))
    if (lock.get("schema") != "katana.dexseq_sources/v1"
            or (source / "Dockerfile").read_text(encoding="utf-8").splitlines()[0] != "FROM " + lock["base_image"]
            or "@sha256:" not in lock["base_image"]):
        raise ValueError("Dockerfile must match the immutable source lock base")
    context.mkdir(parents=True)
    (context / "gpl-scripts").mkdir()
    for name in ("Dockerfile", "LICENSE.GPL-3.0", "THIRD_PARTY.md", "sources.lock.json", "check_runtime.py"):
        shutil.copyfile(source / name, context / name)
    for item in lock["files"]:
        if Path(item["name"]).name != item["name"] or not item["url"].startswith(lock["source_url_prefix"]):
            raise ValueError("Invalid pinned upstream source identity")
        if download is None:
            with urlopen(item["url"], timeout=60) as response:
                data = response.read()
        else:
            data = download(item["url"])
        if hashlib.sha256(data).hexdigest() != item["sha256"]:
            raise ValueError(f"Upstream source checksum mismatch: {item['name']}")
        (context / "gpl-scripts" / item["name"]).write_bytes(data)
    return lock


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--context", type=Path, required=True)
    parser.add_argument("--tag", default="kanata-dexseq-tools:0.1.1")
    parser.add_argument("--docker", default="docker")
    args = parser.parse_args()
    context = args.context.resolve()
    lock = prepare_context(context)
    with (context / "build.log").open("w", encoding="utf-8") as output:
        subprocess.run([args.docker, "build", "--platform", "linux/amd64", "--tag", args.tag, str(context)],
                       stdout=output, stderr=subprocess.STDOUT, check=True)
    image = json.loads(subprocess.check_output([args.docker, "image", "inspect", args.tag], text=True))[0]
    command = [args.docker, "run", "--rm", "--network", "none", "--user", "1000:1000", args.tag,
               "python", "/usr/share/doc/kanata-dexseq-tools/check_runtime.py"]
    result = subprocess.run(command, capture_output=True, text=True)
    (context / "runtime.stdout.log").write_text(result.stdout, encoding="utf-8")
    (context / "runtime.stderr.log").write_text(result.stderr, encoding="utf-8")
    receipt = {"schema": "katana.dexseq_image_acceptance/v1", "status": "failed", "tag": args.tag,
               "image_id": image["Id"], "platform": "linux/amd64", "source_lock": lock,
               "native_exit_code": result.returncode, "published": False, "study_data_used": False}
    if result.returncode == 0:
        observed = json.loads(result.stdout)
        if observed.get("status") == "passed_real_dexseq_packaged_software_fixture":
            receipt.update(status="passed", runtime=observed)
    (context / "acceptance.receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    if receipt["status"] != "passed":
        raise RuntimeError(f"DEXSeq runtime check failed; see {context}")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
