"""River Avulsion Module."""

from ._version import get_versions
from .riverbmi import BmiRiverModule
from .rivermodule import RiverModule

__all__ = ["BmiRiverModule", "RiverModule"]

__version__ = get_versions()["version"]
del get_versions
