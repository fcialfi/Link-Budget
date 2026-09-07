import os
import sys

# calculations.py lives at the repository root and is imported as a plain
# top-level module (see gui.py), not as part of the ``link_budget`` package.
# Make sure that directory is importable regardless of the current working
# directory pytest is invoked from.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
