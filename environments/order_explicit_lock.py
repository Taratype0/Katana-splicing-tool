"""Preserve a successful solver transaction's linking order in explicit locks.

Micromamba 2.9 links explicit URLs in file order. Alphabetical exports can run
Bioconductor post-link scripts before their Python/R/yq/curl prerequisites exist.
Only a complete actual successful transaction log is accepted here. The two
downloaded annotation packages are linked last, after their installer tools.
"""
import argparse
import re
from pathlib import Path
from urllib.parse import urlsplit


POST_LINK_DATA = {"bioconductor-genomeinfodbdata", "bioconductor-go.db"}
POST_LINK_TOOLS = {"python", "r-base", "yq", "curl", "bioconductor-data-packages"}


def package_identity(url):
    filename = urlsplit(url).path.rsplit("/", 1)[-1]
    if filename.endswith(".conda"):
        return filename[:-6]
    if filename.endswith(".tar.bz2"):
        return filename[:-8]
    raise ValueError(f"Unsupported package archive: {filename}")


def reorder(lock_text, transaction_text):
    if "Transaction finished" not in transaction_text:
        raise ValueError("Successful transaction completion is missing")
    urls = [line.strip() for line in lock_text.splitlines() if line.startswith("https://")]
    by_identity = {package_identity(url): url for url in urls}
    if not urls or len(by_identity) != len(urls):
        raise ValueError("Empty or duplicate explicit package identities")
    linked = re.findall(r"^Linking (\S+)\s*$", transaction_text, re.MULTILINE)
    if len(linked) != len(set(linked)) or set(linked) != set(by_identity):
        raise ValueError("Transaction must cover every locked package exactly once")
    names = {identity.rsplit("-", 2)[0] for identity in linked}
    if POST_LINK_DATA & names and not POST_LINK_TOOLS <= names:
        raise ValueError(f"Missing Bioconductor post-link tools: {sorted(POST_LINK_TOOLS - names)}")
    normal = [identity for identity in linked if identity.rsplit("-", 2)[0] not in POST_LINK_DATA]
    data = [identity for identity in linked if identity.rsplit("-", 2)[0] in POST_LINK_DATA]
    return ("# platform: linux-64\n"
            "# Actual successful solver link order; annotation post-link packages last.\n"
            "@EXPLICIT\n" + "\n".join(by_identity[identity] for identity in normal + data) + "\n")


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("transaction_log", type=Path)
    parser.add_argument("locks", nargs="+", type=Path)
    args = parser.parse_args()
    log = args.transaction_log.read_text(encoding="utf-8")
    outputs = [(path, reorder(path.read_text(encoding="utf-8"), log)) for path in args.locks]
    for path, text in outputs:
        path.write_text(text, encoding="utf-8")
        print(f"Ordered {len(text.splitlines()) - 3} unchanged package URLs: {path.name}")


if __name__ == "__main__":
    main()
