"""Package initialization for :mod:`link_budget`."""

# Re-export the ``setup_gui`` function for external use.
from .gui import setup_gui
__all__ = ["setup_gui"]
