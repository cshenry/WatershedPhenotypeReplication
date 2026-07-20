"""pytest bootstrap for WatershedPhenotypeReplication.

The project's ``util.py`` lives under a PRJ-prefixed directory (not a valid
Python package name, per the work-notebook convention -- see
notebooks/PRJ-watershed_phenotype_replication/util.py's docstring), so it is
imported here via an explicit sys.path insert rather than a package import.
This mirrors the notebook's ``%run util.py`` usage without requiring a live
KBase session: util.py's module-level code only builds path constants and a
local NotebookSession cache -- it does not touch the network at import time.
"""
from __future__ import annotations

import sys
from pathlib import Path

_PRJ_DIR = Path(__file__).resolve().parent.parent / "notebooks" / "PRJ-watershed_phenotype_replication"

if str(_PRJ_DIR) not in sys.path:
    sys.path.insert(0, str(_PRJ_DIR))
