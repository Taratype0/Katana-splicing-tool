"""Public module name, separate from durable engine and result identifiers."""

import re


MODULE_NAME = "Katana-Splice-Network"
LEGACY_MODULE_NAME = "CoSplice"
_LEGACY_DISPLAY_NAME = re.compile(
    r"(?<![A-Za-z0-9_./\\-])(?:Katana-Splice(?:-Network)?|katana-splice-network|(?:Katana-)?CoSplice|cosplice)(?![A-Za-z0-9_./\\-])"
)


def display_module_name(value: str) -> str:
    """Rename presentation text without renaming paths, keys or peer tools."""
    return _LEGACY_DISPLAY_NAME.sub(MODULE_NAME, value)
