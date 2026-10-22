from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from threading import RLock
from typing import Any

_RENDER_LOCK = RLock()


def _static_layout(spec: dict[str, Any]) -> None:
    encoding = spec.get("encoding", {})
    column = encoding.get("column", {})
    values = spec.get("data", {}).get("values", [])
    field = column.get("field")
    columns = max(1, len({str(row.get(field)) for row in values})) if field else 1
    if spec.get("width") == "container":
        # Vega-Lite width is per facet, not the combined width of the figure.
        spec["width"] = max(260, (1000 - 70 * (columns - 1)) // columns)
    if isinstance(spec.get("spec"), dict):
        _static_layout(spec["spec"])
    for key in ("hconcat", "vconcat", "concat"):
        for child in spec.get(key, []):
            _static_layout(child)


def export_vega(specification: dict[str, Any], directory: Path, stem: str) -> dict[str, str]:
    import vl_convert as vlc

    spec = deepcopy(specification)
    def static_fonts(value):
        if isinstance(value, dict):
            for key, item in value.items():
                if key in {"font", "labelFont", "titleFont"}:
                    value[key] = "Arial"
                else:
                    static_fonts(item)
        elif isinstance(value, list):
            for item in value:
                static_fonts(item)
    # Web fonts are bundled with the UI, not installed system-wide for rendering.
    static_fonts(spec)
    spec["background"] = "white"
    _static_layout(spec)
    spec["autosize"] = {"type": "pad", "contains": "padding"}
    directory.mkdir(parents=True, exist_ok=True)
    outputs = {extension: directory / f"{stem}.{extension}" for extension in ("png", "svg")}
    with _RENDER_LOCK:
        outputs["png"].write_bytes(vlc.vegalite_to_png(spec, scale=1.5))
        outputs["svg"].write_text(vlc.vegalite_to_svg(spec), encoding="utf-8")
    return {extension: str(path) for extension, path in outputs.items()}
