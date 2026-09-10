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
    from .core.evaluators import Evaluator, CustomEvaluator
    from .strategies.config import (
        GreedyConfig,
        RouletteWheelConfig,
        UCBConfig,
        EpsilonGreedyConfig,
        BayesUCBConfig,
        TopTwoConfig,
    )
    from .warmup.config import (
        EnhancedWarmupConfig,
        BalancedWarmupConfig,
    )
    from .core.evaluator_config import (
        LookupEvaluatorConfig,
        DBEvaluatorConfig,
        FPEvaluatorConfig,
        MWEvaluatorConfig,
        ROCSEvaluatorConfig,
        FredEvaluatorConfig,
        MLClassifierEvaluatorConfig,
        CustomEvaluatorConfig,
    )

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
    # Evaluator classes (ABC + user-function adapter)
    "Evaluator": ".core.evaluators",
    "CustomEvaluator": ".core.evaluators",
    # Pydantic component configs
    "GreedyConfig": ".strategies.config",
    "RouletteWheelConfig": ".strategies.config",
    "UCBConfig": ".strategies.config",
    "EpsilonGreedyConfig": ".strategies.config",
    "BayesUCBConfig": ".strategies.config",
    "TopTwoConfig": ".strategies.config",
    "EnhancedWarmupConfig": ".warmup.config",
    "BalancedWarmupConfig": ".warmup.config",
    "LookupEvaluatorConfig": ".core.evaluator_config",
    "DBEvaluatorConfig": ".core.evaluator_config",
    "FPEvaluatorConfig": ".core.evaluator_config",
    "MWEvaluatorConfig": ".core.evaluator_config",
    "ROCSEvaluatorConfig": ".core.evaluator_config",
    "FredEvaluatorConfig": ".core.evaluator_config",
    "MLClassifierEvaluatorConfig": ".core.evaluator_config",
    "CustomEvaluatorConfig": ".core.evaluator_config",
    # Utilities
    "get_logger": ".utils.ts_logger",
    "read_reagents": ".utils.ts_utils",
    "create_reagents": ".utils.ts_utils",
})
