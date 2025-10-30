"""
Warmup strategies for Thompson Sampling.

This module provides different strategies for the warmup phase of Thompson Sampling,
which initializes reagent posteriors before the main search begins.

Available strategies:
- StandardWarmup: Random partner selection with replacement (baseline)
- StratifiedWarmup: Random selection without replacement (recommended)
- EnhancedWarmup: Parallel pairing with shuffling (legacy)
- LatinHypercubeWarmup: Space-filling stratified sampling (advanced)
"""

from .base import WarmupStrategy
from .standard import StandardWarmup
from .stratified import StratifiedWarmup
from .enhanced import EnhancedWarmup
from .latin_hypercube import LatinHypercubeWarmup

__all__ = [
    'WarmupStrategy',
    'StandardWarmup',
    'StratifiedWarmup',
    'EnhancedWarmup',
    'LatinHypercubeWarmup',
]
