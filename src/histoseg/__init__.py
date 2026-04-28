"""HistoSeg: HE image analysis and spatial contour analysis tools."""

from . import contour, he

try:
    from ._version import __version__
except Exception:
    __version__ = "0+unknown"

__all__ = ["__version__", "contour", "he"]
