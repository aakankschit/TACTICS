"""
Selection strategies for Thompson Sampling.

This package provides pluggable strategies for selecting reagents during
Thompson Sampling search. All strategies inherit from :class:`SelectionStrategy`.

Recommended strategies:
    TopTwoSelection
        Best overall performance (86.1% mean recovery, 114,450 trials, 21 libraries).
        Uses Top-Two Thompson Sampling with adaptive per-component thermal cycling
        and GMIC-weighted component rotation. Targets best-arm identification.
    RouletteWheelSelection
        Close second (85.5% mean recovery). Uses CATS (Component-Aware Thompson
        Sampling) with Boltzmann thermal cycling and GMIC criticality-weighted
        rotation. The original TACTICS method.

Baseline strategies (for comparison and ablation):
    GreedySelection
        Pure argmax Thompson Sampling. Simplest baseline — sample posteriors,
        pick the best reagent at each component.
    UCBSelection
        Upper Confidence Bound. Consistently outperformed by TopTwoSelection
        and RouletteWheelSelection in benchmarks.
    EpsilonGreedySelection
        Epsilon-greedy with decay. Consistently outperformed by TopTwoSelection
        and RouletteWheelSelection in benchmarks.
    BayesUCBSelection
        Bayesian UCB variant. Consistently outperformed by TopTwoSelection
        and RouletteWheelSelection in benchmarks.
"""

from .base_strategy import SelectionStrategy
from .greedy_selection import GreedySelection
from .roulette_wheel import RouletteWheelSelection
from .ucb_selection import UCBSelection
from .epsilon_greedy import EpsilonGreedySelection
from .bayes_ucb_selection import BayesUCBSelection
from .top_two_selection import TopTwoSelection

__all__ = [
    "SelectionStrategy",
    "TopTwoSelection",
    "RouletteWheelSelection",
    "GreedySelection",
    "UCBSelection",
    "EpsilonGreedySelection",
    "BayesUCBSelection",
]
