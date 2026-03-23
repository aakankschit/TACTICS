"""
Warmup strategies for Thompson Sampling.

This module provides different strategies for the warmup phase of Thompson Sampling,
which initializes reagent posteriors before the main search begins.

Recommended:
    EnhancedWarmup: Stochastic parallel pairing. Universally optimal across both
    balanced and imbalanced libraries (114,450+ trials, 21 libraries). On imbalanced
    libraries, its natural over-sampling of the small component pre-solves that
    component's ranking, which GMIC-weighted rotation then exploits.

Alternatives:
    BalancedWarmup: Exactly K observations per reagent with stratified partners.
    Good for isolating framework gains (e.g., Balanced-Greedy vs Legacy-Greedy),
    but consistently slightly below EnhancedWarmup when paired with GMIC rotation.

Baseline:
    StandardWarmup: Random partner selection with replacement. For comparison only.
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
