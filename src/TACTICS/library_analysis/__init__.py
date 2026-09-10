"""Analysis and visualisation utilities for TACTICS benchmark output.

Plotting requires the optional ``viz`` extra (``pip install chem-tactics[viz]``).
``TS_Benchmarks`` is resolved lazily so importing this package -- or
``TACTICS.library_analysis.diagnostic_plots`` -- does not load altair.
"""

import sys
from typing import TYPE_CHECKING

from .._lazy import install as _install

if TYPE_CHECKING:
    from .visualization import TS_Benchmarks

_install(sys.modules[__name__], {"TS_Benchmarks": ".visualization"})
