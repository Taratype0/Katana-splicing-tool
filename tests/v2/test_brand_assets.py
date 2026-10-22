from pathlib import Path

from katana_next.desktop import desktop_icon_path


def test_original_icon_is_available_to_desktop():
    from PIL import Image

    path = desktop_icon_path()
    assert path == Path(__file__).resolve().parents[2] / "assets/katana_icon.ico"
    with Image.open(path) as image:
        assert image.format == "ICO"
        assert {(16, 16), (32, 32), (48, 48), (128, 128)} <= image.ico.sizes()


def test_frontend_and_packager_bind_original_artwork():
    repo = Path(__file__).resolve().parents[2]
    assert 'href="/katana_icon.ico"' in (repo / "apps/web/index.html").read_text()
    assert 'publicDir: "../../assets"' in (repo / "apps/web/vite.config.ts").read_text()
    packager = (repo / "scripts/build_next_desktop.py").read_text()
    assert '"--icon", str(repo / "assets/katana_icon.ico")' in packager
