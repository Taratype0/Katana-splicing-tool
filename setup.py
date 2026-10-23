"""Include pre-built browser assets in wheels; consumers do not need Node.js."""
from pathlib import Path
import os
import shutil

from setuptools import setup
from setuptools.command.build_py import build_py


class BuildWithBrowser(build_py):
    def run(self):
        super().run()
        pipeline = Path(__file__).parent / "pipeline"
        if (pipeline / "main.nf").is_file():
            shutil.copytree(pipeline, Path(self.build_lib) / "katana_next" / "flow" / "_pipeline",
                            dirs_exist_ok=True, ignore=shutil.ignore_patterns("__pycache__", "*.pyc", ".nextflow*", "work"))
        source = Path(__file__).parent / "apps" / "web" / "dist"
        if not (source / "index.html").is_file():
            if os.environ.get("KATANA_REQUIRE_WEB") == "1":
                raise RuntimeError("Build apps/web with pnpm build before building a browser wheel.")
            self.announce("Browser not built; CLI-only wheel. Build apps/web first for a browser wheel.", level=2)
            return
        destination = Path(self.build_lib) / "katana_next" / "_web"
        # Only clean this build-owned destination, never source data or results.
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(source, destination)


setup(cmdclass={"build_py": BuildWithBrowser})
