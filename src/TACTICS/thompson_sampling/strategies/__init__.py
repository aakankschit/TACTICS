"""
Selection strategies for Thompson Sampling.

This package contains different strategies for selecting reagents during sampling.
"""

from .base_strategy import SelectionStrategy
from .greedy_selection import GreedySelection
from .roulette_wheel import RouletteWheelSelection
from .ucb_selection import UCBSelection
from .epsilon_greedy import EpsilonGreedySelection

__all__ = [
    'SelectionStrategy',
    'GreedySelection',
    'RouletteWheelSelection',
    'UCBSelection',
    'EpsilonGreedySelection',
]
