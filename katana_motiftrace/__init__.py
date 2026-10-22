"""Public MotifTrace integration over the established Katana algorithms."""
from .api import MotifTrace
from .catalog import PACKAGE_VERSION, SCIENTIFIC_RELEASE, mode_catalog

__version__ = PACKAGE_VERSION
scientific_release = SCIENTIFIC_RELEASE
__all__ = ["MotifTrace", "mode_catalog", "SCIENTIFIC_RELEASE", "scientific_release", "__version__"]
