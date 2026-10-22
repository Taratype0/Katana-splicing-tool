"""Attach an explicitly installed FIMO runtime to a new project manifest."""
import argparse
import json
from pathlib import Path
from katana_next.domain.models import InputAsset
from katana_next.infrastructure.manifest_store import ManifestStore


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--runtime-json", type=Path, required=True)
    parser.add_argument("--python", required=True)
    parser.add_argument("--wsl-distribution", default="")
    args = parser.parse_args()
    if args.manifest.resolve() == args.destination.resolve():
        raise ValueError("Use a new destination manifest to preserve the original")
    args.runtime_json.parent.mkdir(parents=True, exist_ok=True)
    args.runtime_json.write_text(json.dumps({"python": args.python, "wsl_distribution": args.wsl_distribution,
                                             "engine": "pymemesuite.FIMO", "version": "0.1.0a4"}, indent=2), encoding="utf-8")
    if args.destination.exists():
        raise FileExistsError(args.destination)
    store = ManifestStore()
    _, manifest = store.load(args.manifest)
    asset = InputAsset(asset_id="motif::fimo_runtime", kind="fimo_runtime", path=args.runtime_json.resolve(), metadata={"role": "calibrated FIMO site significance"})
    updated = manifest.model_copy(update={"assets": tuple(a for a in manifest.assets if a.kind != "fimo_runtime") + (asset,)})
    args.destination.parent.mkdir(parents=True, exist_ok=True)
    store.save(updated, args.destination)
    print(args.destination)


if __name__ == "__main__":
    main()
