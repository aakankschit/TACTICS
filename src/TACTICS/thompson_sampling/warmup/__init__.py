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
    Good for isolating framework gains (Balanced + Greedy vs random warmup + Greedy),
    but consistently slightly below EnhancedWarmup when paired with GMIC rotation.
"""

from .base import WarmupStrategy
from .enhanced import EnhancedWarmup
from .balanced import BalancedWarmup

__all__ = [
    'WarmupStrategy',
    'BalancedWarmup',
    'EnhancedWarmup',
]
