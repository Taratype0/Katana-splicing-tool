from __future__ import annotations

import argparse

import uvicorn

from katana_next.api.app import create_app


app = create_app()


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the Katana Next local API.")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", default=8765, type=int)
    parser.add_argument("--reload", action="store_true")
    args = parser.parse_args()
    uvicorn.run(
        "katana_next.api.main:app",
        host=args.host,
        port=args.port,
        reload=args.reload,
    )


if __name__ == "__main__":
    main()
