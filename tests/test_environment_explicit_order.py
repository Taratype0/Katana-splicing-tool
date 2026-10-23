import importlib.util
from pathlib import Path

import pytest


spec = importlib.util.spec_from_file_location(
    "order_explicit_lock", Path(__file__).parents[1] / "environments" / "order_explicit_lock.py")
ordering = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ordering)


def fixture():
    names = ["python", "r-base", "yq", "curl", "bioconductor-data-packages",
             "bioconductor-genomeinfodbdata", "bioconductor-go.db", "r-wgcna"]
    identities = [f"{name}-1.0-0" for name in names]
    urls = [f"https://example.invalid/linux-64/{identity}.conda#{'a' * 32}" for identity in identities]
    return urls, "\n".join(f"Linking {identity}" for identity in identities) + "\nTransaction finished\n"


def test_preserves_every_locked_identity_and_hash_but_orders_postlink_last():
    urls, log = fixture()
    output = ordering.reorder("\n".join(reversed(urls)), log)
    ordered_urls = [line for line in output.splitlines() if line.startswith("https:")]
    assert set(ordered_urls) == set(urls)
    assert ordered_urls[-2:] == urls[5:7]
    assert ordered_urls.index(urls[2]) < ordered_urls.index(urls[5])
    assert "@EXPLICIT" in output


def test_refuses_partial_or_failed_transaction():
    urls, log = fixture()
    with pytest.raises(ValueError, match="completion"):
        ordering.reorder("\n".join(urls), log.replace("Transaction finished", ""))
    with pytest.raises(ValueError, match="every locked package"):
        ordering.reorder("\n".join(urls), log.replace("Linking curl-1.0-0\n", ""))


def test_refuses_postlink_without_declared_installation_tools():
    urls, log = fixture()
    urls = [url for url in urls if "/yq-" not in url]
    log = log.replace("Linking yq-1.0-0\n", "")
    with pytest.raises(ValueError, match="post-link tools"):
        ordering.reorder("\n".join(urls), log)
