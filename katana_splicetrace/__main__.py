"""Allow `python -m katana_splicetrace` execution."""
from katana_splicetrace.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
