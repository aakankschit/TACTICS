"""Greedy selection strategy — pure argmax/argmin Thompson Sampling."""

import numpy as np
from .base_strategy import SelectionStrategy


class GreedySelection(SelectionStrategy):
    """Standard greedy selection: sample posteriors, pick argmax (or argmin)."""

    def __init__(self, mode="maximize"):
        super().__init__(mode)

    def select_reagent(self, reagent_list, disallow_mask=None, **kwargs):
        """Select the reagent with the best sampled posterior value."""
        rng = kwargs.get('rng', np.random.default_rng())

        stds = np.array([r.std for r in reagent_list])
        mu = np.array([r.mean for r in reagent_list])
        scores = rng.normal(size=len(reagent_list)) * stds + mu

        if disallow_mask:
            scores[np.array(list(disallow_mask))] = np.nan

        if self.mode in ["maximize", "maximize_boltzmann"]:
            return np.nanargmax(scores)
        else:
            return np.nanargmin(scores)
