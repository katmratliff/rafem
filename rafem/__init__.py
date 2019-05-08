"""River Avulsion Module."""

from .riverbmi import BmiRiverModule
from .rivermodule import RiverModule


__all__ = ['BmiRiverModule', 'RiverModule']

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
