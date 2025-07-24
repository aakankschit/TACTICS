"""
Thompson Sampling package for PRISMS.

This package provides Thompson Sampling implementations for combinatorial library screening.
"""

# Import configuration
from .config import StandardSamplerConfig, EnhancedSamplerConfig, RandomBaselineConfig, TSConfig

# Import core functionality
from .core import (
    StandardThompsonSampler,
    EnhancedThompsonSampler,
    ROCSEvaluator,
    LookupEvaluator,
    DBEvaluator,
    FredEvaluator,
    FPEvaluator,
    MWEvaluator,
    MLClassifierEvaluator
)

# Import utilities
from .utils import get_logger, read_reagents, create_reagents

# Import main function
from .main import run_ts

# Import baseline functionality
from .baseline import run_random_baseline, run_exhaustive_baseline

__all__ = [
    # Configuration
    'StandardSamplerConfig',
    'EnhancedSamplerConfig',
    'RandomBaselineConfig',
    'TSConfig',
    
    # Main function
    'run_ts',
    
    # Baseline functions
    'run_random_baseline',
    'run_exhaustive_baseline',
    
    # Core samplers
    'StandardThompsonSampler',
    'EnhancedThompsonSampler',
    
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