"""
Thompson Sampling package for TACTICS.

This package provides Thompson Sampling implementations for combinatorial library screening.
"""

# Import configuration
from .config import ThompsonSamplingConfig

# Import core functionality
from .core import (
    ThompsonSampler,
    Reagent,
    ROCSEvaluator,
    LookupEvaluator,
    DBEvaluator,
    FredEvaluator,
    FPEvaluator,
    MWEvaluator,
    MLClassifierEvaluator
)

# Import strategies
from .strategies import (
    SelectionStrategy,
    GreedySelection,
    RouletteWheelSelection,
    UCBSelection,
    EpsilonGreedySelection,
    BayesUCBSelection,
    TopTwoSelection
)

# Import warmup strategies
from .warmup import (
    WarmupStrategy,
    BalancedWarmup,
    StandardWarmup,
    EnhancedWarmup
)

# Import utilities
from .utils import get_logger, read_reagents, create_reagents


# Import presets for easy access
from .presets import ConfigPresets, get_preset

__all__ = [
    # Configuration
    'ThompsonSamplingConfig',

    # Presets
    'ConfigPresets',
    'get_preset',


    # Core classes
    'ThompsonSampler',
    'Reagent',

    # Selection strategies
    'SelectionStrategy',
    'GreedySelection',
    'RouletteWheelSelection',
    'UCBSelection',
    'EpsilonGreedySelection',
    'BayesUCBSelection',
    'TopTwoSelection',

    # Warmup strategies
    'WarmupStrategy',
    'BalancedWarmup',
    'StandardWarmup',
    'EnhancedWarmup',

    # Evaluators
    'ROCSEvaluator',
    'LookupEvaluator',
    'DBEvaluator',
    'FredEvaluator',
    'FPEvaluator',
    'MWEvaluator',
    'MLClassifierEvaluator',

    # Utilities
    'get_logger',
    'read_reagents',
    'create_reagents',
]