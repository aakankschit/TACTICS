"""
Warmup strategies for Thompson Sampling.

This module provides different strategies for the warmup phase of Thompson Sampling,
which initializes reagent posteriors before the main search begins.

Available strategies:
- BalancedWarmup: Exactly K observations per reagent with stratified partners (recommended)
- StandardWarmup: Random partner selection with replacement (baseline)
- EnhancedWarmup: Parallel pairing with shuffling (legacy)
"""

from .base import WarmupStrategy
from .standard import StandardWarmup
from .enhanced import EnhancedWarmup
from .balanced import BalancedWarmup

__all__ = [
    'WarmupStrategy',
    'BalancedWarmup',
    'StandardWarmup',
    'EnhancedWarmup',
]
