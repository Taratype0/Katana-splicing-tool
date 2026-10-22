"""Prepare and explicitly attach versioned mouse GO/RBP resources on the data drive."""
import argparse
from pathlib import Path
from katana_next.resources import prepare_mouse_resources, attach_resources

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("resource_root", type=Path)
    parser.add_argument("--attract", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--destination", type=Path)
    parser.add_argument("--output-root", type=Path)
    args = parser.parse_args()
    bundle = prepare_mouse_resources(args.resource_root, args.attract)
    print(bundle, flush=True)
    if args.manifest:
        if not args.destination:
            parser.error("--destination is required to preserve the original manifest")
        print(attach_resources(args.manifest, bundle, args.destination, args.output_root), flush=True)
