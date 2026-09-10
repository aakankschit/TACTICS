"""
TACTICS: Thompson Sampling-Assisted Chemical Targeting and Iterative Compound Selection.

A Python library for Thompson Sampling-based optimization of chemical combinatorial
libraries for drug discovery.

Recommended usage::

    from TACTICS import ThompsonSampler, get_preset
    from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef
    from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

    pipeline = SynthesisPipeline(ReactionConfig(...))
    config = get_preset(synthesis_pipeline=pipeline,
                        evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv"))
    sampler = ThompsonSampler.from_config(config)
    results = sampler.search(num_cycles=config.num_ts_iterations)
    sampler.close()

Selection strategies (recommended):
    - :class:`TopTwoSelection`: Best overall performance (Top-Two Thompson Sampling)
    - :class:`RouletteWheelSelection`: CATS with thermal cycling and GMIC rotation

Selection strategies (baselines):
    - :class:`GreedySelection`: Pure argmax Thompson Sampling
    - :class:`UCBSelection`: Upper Confidence Bound
    - :class:`EpsilonGreedySelection`: Epsilon-greedy with decay
    - :class:`BayesUCBSelection`: Bayesian UCB
"""

import sys
from typing import TYPE_CHECKING

from ._lazy import install as _install

try:
    from importlib.metadata import version as _version
    __version__ = _version("chem-tactics")
except Exception:  # not installed (e.g. running from a source checkout)
    __version__ = "0.0.0"

if TYPE_CHECKING:
    from .thompson_sampling.config import ThompsonSamplingConfig
    from .thompson_sampling.presets import ConfigPresets, get_preset
    from .thompson_sampling.core import ThompsonSampler, Reagent
    from .thompson_sampling.strategies import (
        SelectionStrategy,
        TopTwoSelection,
        RouletteWheelSelection,
        GreedySelection,
        UCBSelection,
        EpsilonGreedySelection,
        BayesUCBSelection,
    )
    from .thompson_sampling.warmup import WarmupStrategy, EnhancedWarmup, BalancedWarmup
    from .thompson_sampling.core import (
        ROCSEvaluator,
        LookupEvaluator,
        DBEvaluator,
        FredEvaluator,
        FPEvaluator,
        MWEvaluator,
        MLClassifierEvaluator,
    )
    from .thompson_sampling.utils import get_logger, read_reagents, create_reagents
    from .thompson_sampling.core.evaluators import Evaluator, CustomEvaluator
    from .thompson_sampling.strategies.config import (
        GreedyConfig,
        RouletteWheelConfig,
        UCBConfig,
        EpsilonGreedyConfig,
        BayesUCBConfig,
        TopTwoConfig,
    )
    from .thompson_sampling.warmup.config import (
        EnhancedWarmupConfig,
        BalancedWarmupConfig,
    )
    from .thompson_sampling.core.evaluator_config import (
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
    "ThompsonSamplingConfig": ".thompson_sampling.config",
    "ConfigPresets": ".thompson_sampling.presets",
    "get_preset": ".thompson_sampling.presets",
    # Core
    "ThompsonSampler": ".thompson_sampling.core.sampler",
    "Reagent": ".thompson_sampling.core.reagent",
    # Selection strategies -- recommended
    "SelectionStrategy": ".thompson_sampling.strategies.base_strategy",
    "TopTwoSelection": ".thompson_sampling.strategies.top_two_selection",
    "RouletteWheelSelection": ".thompson_sampling.strategies.roulette_wheel",
    # Selection strategies -- baselines
    "GreedySelection": ".thompson_sampling.strategies.greedy_selection",
    "UCBSelection": ".thompson_sampling.strategies.ucb_selection",
    "EpsilonGreedySelection": ".thompson_sampling.strategies.epsilon_greedy",
    "BayesUCBSelection": ".thompson_sampling.strategies.bayes_ucb_selection",
    # Warmup strategies
    "WarmupStrategy": ".thompson_sampling.warmup.base",
    "EnhancedWarmup": ".thompson_sampling.warmup.enhanced",
    "BalancedWarmup": ".thompson_sampling.warmup.balanced",
    # Evaluators
    "ROCSEvaluator": ".thompson_sampling.core.evaluators",
    "LookupEvaluator": ".thompson_sampling.core.evaluators",
    "DBEvaluator": ".thompson_sampling.core.evaluators",
    "FredEvaluator": ".thompson_sampling.core.evaluators",
    "FPEvaluator": ".thompson_sampling.core.evaluators",
    "MWEvaluator": ".thompson_sampling.core.evaluators",
    "MLClassifierEvaluator": ".thompson_sampling.core.evaluators",
    # Evaluator classes (ABC + user-function adapter)
    "Evaluator": ".thompson_sampling.core.evaluators",
    "CustomEvaluator": ".thompson_sampling.core.evaluators",
    # Pydantic component configs
    "GreedyConfig": ".thompson_sampling.strategies.config",
    "RouletteWheelConfig": ".thompson_sampling.strategies.config",
    "UCBConfig": ".thompson_sampling.strategies.config",
    "EpsilonGreedyConfig": ".thompson_sampling.strategies.config",
    "BayesUCBConfig": ".thompson_sampling.strategies.config",
    "TopTwoConfig": ".thompson_sampling.strategies.config",
    "EnhancedWarmupConfig": ".thompson_sampling.warmup.config",
    "BalancedWarmupConfig": ".thompson_sampling.warmup.config",
    "LookupEvaluatorConfig": ".thompson_sampling.core.evaluator_config",
    "DBEvaluatorConfig": ".thompson_sampling.core.evaluator_config",
    "FPEvaluatorConfig": ".thompson_sampling.core.evaluator_config",
    "MWEvaluatorConfig": ".thompson_sampling.core.evaluator_config",
    "ROCSEvaluatorConfig": ".thompson_sampling.core.evaluator_config",
    "FredEvaluatorConfig": ".thompson_sampling.core.evaluator_config",
    "MLClassifierEvaluatorConfig": ".thompson_sampling.core.evaluator_config",
    "CustomEvaluatorConfig": ".thompson_sampling.core.evaluator_config",
    # Utilities
    "get_logger": ".thompson_sampling.utils.ts_logger",
    "read_reagents": ".thompson_sampling.utils.ts_utils",
    "create_reagents": ".thompson_sampling.utils.ts_utils",
})
