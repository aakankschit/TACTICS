"""
Thompson Sampling package for TACTICS.

This package provides Thompson Sampling implementations for combinatorial library screening.
"""

# Import configuration
from .config import ThompsonSamplingConfig, RandomBaselineConfig

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
    EpsilonGreedySelection
)

# Import warmup strategies
from .warmup import (
    WarmupStrategy,
    StandardWarmup,
    StratifiedWarmup,
    EnhancedWarmup,
    LatinHypercubeWarmup
)

# Import utilities
from .utils import get_logger, read_reagents, create_reagents

# Import main function
from .main import run_ts

# Import baseline functionality
from .baseline import run_random_baseline, run_exhaustive_baseline

# Import presets for easy access
from .presets import ConfigPresets, get_preset

__all__ = [
    # Configuration
    'ThompsonSamplingConfig',
    'RandomBaselineConfig',

    # Main function
    'run_ts',

    # Presets
    'ConfigPresets',
    'get_preset',

    # Baseline functions
    'run_random_baseline',
    'run_exhaustive_baseline',

    # Core classes
    'ThompsonSampler',
    'Reagent',

    # Selection strategies
    'SelectionStrategy',
    'GreedySelection',
    'RouletteWheelSelection',
    'UCBSelection',
    'EpsilonGreedySelection',

    # Warmup strategies
    'WarmupStrategy',
    'StandardWarmup',
    'StratifiedWarmup',
    'EnhancedWarmup',
    'LatinHypercubeWarmup',

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