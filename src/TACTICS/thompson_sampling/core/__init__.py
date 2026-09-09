"""
Core Thompson Sampling functionality.

This package contains the main Thompson Sampling implementation and evaluators.
Names are resolved lazily on first access (see ``TACTICS._lazy``) so that
importing configuration modules does not load RDKit, sqlitedict or OpenEye.
"""

import sys
from typing import TYPE_CHECKING

from ..._lazy import install as _install

if TYPE_CHECKING:  # static analyzers see the real symbols
    from .sampler import ThompsonSampler
    from .reagent import Reagent
    from .evaluators import (
        ROCSEvaluator,
        LookupEvaluator,
        DBEvaluator,
        FredEvaluator,
        FPEvaluator,
        MWEvaluator,
        MLClassifierEvaluator,
        CustomEvaluator,
    )

_install(sys.modules[__name__], {
    "ThompsonSampler": ".sampler",
    "Reagent": ".reagent",
    "ROCSEvaluator": ".evaluators",
    "LookupEvaluator": ".evaluators",
    "DBEvaluator": ".evaluators",
    "FredEvaluator": ".evaluators",
    "FPEvaluator": ".evaluators",
    "MWEvaluator": ".evaluators",
    "MLClassifierEvaluator": ".evaluators",
    "CustomEvaluator": ".evaluators",
})
