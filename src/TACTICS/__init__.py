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

# Library enumeration
from .library_enumeration import LibraryEnumerator, initializer

# Configuration
from .thompson_sampling.config import ThompsonSamplingConfig, RandomBaselineConfig

# Presets
from .thompson_sampling.presets import ConfigPresets, get_preset

# Core
from .thompson_sampling.core import ThompsonSampler, Reagent

# Selection strategies — recommended
from .thompson_sampling.strategies import (
    SelectionStrategy,
    TopTwoSelection,
    RouletteWheelSelection,
)

# Selection strategies — baselines
from .thompson_sampling.strategies import (
    GreedySelection,
    UCBSelection,
    EpsilonGreedySelection,
    BayesUCBSelection,
)

# Warmup strategies
from .thompson_sampling.warmup import (
    WarmupStrategy,
    EnhancedWarmup,
    BalancedWarmup,
    StandardWarmup,
)

# Evaluators
from .thompson_sampling.core import (
    ROCSEvaluator,
    LookupEvaluator,
    DBEvaluator,
    FredEvaluator,
    FPEvaluator,
    MWEvaluator,
    MLClassifierEvaluator,
)

# Baseline functions
from .thompson_sampling.baseline import run_random_baseline, run_exhaustive_baseline

# Utilities
from .thompson_sampling.utils import get_logger, read_reagents, create_reagents

__all__ = [
    # Configuration
    "ThompsonSamplingConfig",
    "RandomBaselineConfig",

    # Presets
    "ConfigPresets",
    "get_preset",

    # Core classes
    "ThompsonSampler",
    "Reagent",

    # Selection strategies — recommended
    "SelectionStrategy",
    "TopTwoSelection",
    "RouletteWheelSelection",

    # Selection strategies — baselines
    "GreedySelection",
    "UCBSelection",
    "EpsilonGreedySelection",
    "BayesUCBSelection",

    # Warmup strategies
    "WarmupStrategy",
    "EnhancedWarmup",
    "BalancedWarmup",
    "StandardWarmup",

    # Evaluators
    "ROCSEvaluator",
    "LookupEvaluator",
    "DBEvaluator",
    "FredEvaluator",
    "FPEvaluator",
    "MWEvaluator",
    "MLClassifierEvaluator",

    # Baseline functions
    "run_random_baseline",
    "run_exhaustive_baseline",

    # Library enumeration
    "LibraryEnumerator",
    "initializer",

    # Utilities
    "get_logger",
    "read_reagents",
    "create_reagents",
]
