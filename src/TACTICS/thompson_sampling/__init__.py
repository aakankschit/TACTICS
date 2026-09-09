"""
Thompson Sampling package for TACTICS.

This package provides Thompson Sampling implementations for combinatorial
library screening. Names are resolved lazily on first access (see
``TACTICS._lazy``), so ``from TACTICS.thompson_sampling.config import ...``
stays cheap.
"""

import sys
from typing import TYPE_CHECKING

from .._lazy import install as _install

if TYPE_CHECKING:
    from .config import ThompsonSamplingConfig
    from .presets import ConfigPresets, get_preset
    from .core import (
        ThompsonSampler,
        Reagent,
        ROCSEvaluator,
        LookupEvaluator,
        DBEvaluator,
        FredEvaluator,
        FPEvaluator,
        MWEvaluator,
        MLClassifierEvaluator,
    )
    from .strategies import (
        SelectionStrategy,
        GreedySelection,
        RouletteWheelSelection,
        UCBSelection,
        EpsilonGreedySelection,
        BayesUCBSelection,
        TopTwoSelection,
    )
    from .warmup import WarmupStrategy, BalancedWarmup, EnhancedWarmup
    from .utils import get_logger, read_reagents, create_reagents

_install(sys.modules[__name__], {
    # Configuration / presets
    "ThompsonSamplingConfig": ".config",
    "ConfigPresets": ".presets",
    "get_preset": ".presets",
    # Core
    "ThompsonSampler": ".core.sampler",
    "Reagent": ".core.reagent",
    # Selection strategies
    "SelectionStrategy": ".strategies.base_strategy",
    "GreedySelection": ".strategies.greedy_selection",
    "RouletteWheelSelection": ".strategies.roulette_wheel",
    "UCBSelection": ".strategies.ucb_selection",
    "EpsilonGreedySelection": ".strategies.epsilon_greedy",
    "BayesUCBSelection": ".strategies.bayes_ucb_selection",
    "TopTwoSelection": ".strategies.top_two_selection",
    # Warmup strategies
    "WarmupStrategy": ".warmup.base",
    "BalancedWarmup": ".warmup.balanced",
    "EnhancedWarmup": ".warmup.enhanced",
    # Evaluators
    "ROCSEvaluator": ".core.evaluators",
    "LookupEvaluator": ".core.evaluators",
    "DBEvaluator": ".core.evaluators",
    "FredEvaluator": ".core.evaluators",
    "FPEvaluator": ".core.evaluators",
    "MWEvaluator": ".core.evaluators",
    "MLClassifierEvaluator": ".core.evaluators",
    # Utilities
    "get_logger": ".utils.ts_logger",
    "read_reagents": ".utils.ts_utils",
    "create_reagents": ".utils.ts_utils",
})
