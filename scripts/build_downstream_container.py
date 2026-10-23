"""Build the locked Linux runtime from a wheel in an isolated Docker context.

The context contains software only. An optional Conda archive cache speeds the
install; environments, project files and results are never copied into it.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
from urllib.parse import urlsplit
from urllib.request import urlopen
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def download_system_packages(lock_path, destination, base):
    """Download only official, individually hashed system packages for this base."""
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    if (lock.get("schema") != "katana.container_system_packages/v1"
            or lock.get("platform") != "linux/amd64"
            or lock.get("base_image") != base):
        raise ValueError("System package lock must match the immutable Linux base image")
    packages = lock.get("packages", [])
    names = [package["name"] for package in packages]
    if not packages or len(names) != len(set(names)):
        raise ValueError("System package lock is empty or contains duplicate package names")
    destination.mkdir(parents=True, exist_ok=False)
    received = []
    for package in packages:
        parsed = urlsplit(package["url"])
        filename = Path(parsed.path).name
        if (parsed.scheme != "https" or parsed.netloc != "deb.debian.org"
                or not parsed.path.startswith("/debian/pool/")
                or ".." in parsed.path.split("/") or not filename.endswith("_amd64.deb")
                or package.get("architecture") != "amd64"
                or not re.fullmatch(r"[0-9a-f]{64}", package.get("sha256", ""))
                or not isinstance(package.get("size_bytes"), int) or package["size_bytes"] <= 0):
            raise ValueError("Invalid official system package URL, architecture or checksum")
        target = destination / filename
        temporary = target.with_suffix(target.suffix + ".part")
        with urlopen(package["url"], timeout=60) as response, temporary.open("xb") as output:
            shutil.copyfileobj(response, output)
        if temporary.stat().st_size != package["size_bytes"] or digest(temporary) != package["sha256"]:
            raise ValueError(f"System package checksum/size mismatch: {package['name']}")
        temporary.rename(target)
        received.append({key: package[key] for key in ("name", "version", "architecture", "sha256")})
    return received


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--context", type=Path, required=True, help="New directory outside the checkout")
    parser.add_argument("--tag", required=True, help="Local image name; this command does not publish")
    parser.add_argument("--base", default="mambaorg/micromamba:2.9.0@sha256:4f29a2b62a6e65ffa3dfbccb62d8a90514e95954aae52a90afaa8e7dae347bbf")
    parser.add_argument("--docker", default="docker")
    parser.add_argument("--network", default="default", choices=("default", "host", "none"))
    parser.add_argument("--conda-cache", type=Path)
    args = parser.parse_args()
    context = args.context.resolve()
    if context.exists() or context.is_relative_to(ROOT.resolve()):
        parser.error("Use a new context directory outside the source checkout")
    wheel = args.wheel.resolve(strict=True)
    with zipfile.ZipFile(wheel) as archive:
        if "katana_next/flow/runtime_checks.py" not in archive.namelist():
            parser.error("Wheel does not contain the downstream runtime implementation")
    locks = [ROOT / "environments/locks" / name for name in
             ("linux-64.explicit.txt", "pip-linux-64.lock")]
    for path in locks:
        if not path.is_file():
            parser.error(f"Accepted installation lock missing: {path.name}")
    subprocess.run([args.docker, "pull", "--platform", "linux/amd64", args.base], check=True)
    inspected = json.loads(subprocess.check_output([args.docker, "image", "inspect", args.base], text=True))[0]
    candidates = inspected.get("RepoDigests", [])
    base = args.base if "@sha256:" in args.base else next(iter(candidates), None)
    if not base or "@sha256:" not in base:
        raise RuntimeError("Base image has no immutable registry digest")
    context.mkdir(parents=True)
    (context / "locks").mkdir(); (context / "wheel").mkdir(); (context / "conda-cache").mkdir()
    for path in locks:
        shutil.copy2(path, context / "locks" / path.name)
    shutil.copy2(wheel, context / "wheel" / wheel.name)
    shutil.copy2(ROOT / "containers/downstream/Dockerfile", context / "Dockerfile")
    shutil.copy2(ROOT / "scripts/check_downstream_environment.py", context / "check_downstream_environment.py")
    system_lock = ROOT / "containers/downstream/system-packages-linux-amd64.json"
    system_packages = download_system_packages(system_lock, context / "system-packages", base)
    shutil.copy2(system_lock, context / "locks" / system_lock.name)
    if args.conda_cache:
        # Only archives named in this explicit lock may enter the build context.
        urls = {line.split("#")[0] for line in locks[0].read_text().splitlines()
                if line.startswith("https://")}
        for url in sorted(urls):
            parsed = urlsplit(url)
            relative = Path(parsed.scheme) / parsed.netloc / parsed.path.lstrip("/")
            if ".." in relative.parts:
                raise ValueError("Unsafe archive path in explicit lock")
            path = args.conda_cache / relative
            if not path.is_file():
                path = args.conda_cache / relative.name
            if path.is_file() and not path.is_symlink():
                destination = context / "conda-cache" / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, destination)
    inputs = {"schema": "katana.container_build/v1", "base": base,
              "wheel": {"name": wheel.name, "sha256": digest(wheel)},
              "locks": {p.name: digest(p) for p in locks}, "platform": "linux/amd64",
              "system_packages_lock_sha256": digest(system_lock), "system_packages": system_packages,
              "study_data": False, "published": False}
    (context / "build-inputs.json").write_text(json.dumps(inputs, indent=2) + "\n")
    argv = [args.docker, "build", "--platform", "linux/amd64", "--network", args.network,
            "--build-arg", f"BASE_IMAGE={base}", "--tag", args.tag, str(context)]
    with (context / "build.log").open("w", encoding="utf-8") as log:
        result = subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"Image build failed; inspect {context / 'build.log'}")
    image = json.loads(subprocess.check_output([args.docker, "image", "inspect", args.tag], text=True))[0]
    receipt = {**inputs, "status": "built", "image_id": image["Id"], "tag": args.tag,
               "runtime_acceptance": "not_run"}
    (context / "build.receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
