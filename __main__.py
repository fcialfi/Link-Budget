"""Entry point for the Link Budget GUI.

This module supports execution both as ``python -m link_budget`` and by
invoking ``python __main__.py`` directly (useful when building installers
that execute the bundled script without package context).
"""

import os
import sys

if __package__ is None or __package__ == "":
    # Running as a script (e.g., PyInstaller entry-point): ensure local imports work
    sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
    from gui import setup_gui  # type: ignore  # pylint: disable=import-error
else:
    # Running as a package (``python -m link_budget``)
    from .gui import setup_gui

if __name__ == "__main__":
    setup_gui()
