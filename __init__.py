"""Package initialization for :mod:`link_budget`."""

if __package__:
    # Re-export the ``setup_gui`` function for external use.
    from .gui import setup_gui

    __all__ = ["setup_gui"]
else:
    # Imported without a parent package -- e.g. by pytest's package-aware
    # test collection, which imports every __init__.py it finds above the
    # tests/ directory purely to build its collection tree. Relative imports
    # aren't available in that context, and gui.py pulls in Tkinter, which
    # need not be installed in a headless test environment. Skip the
    # re-export rather than failing; nothing in this situation needs it.
    __all__ = []
