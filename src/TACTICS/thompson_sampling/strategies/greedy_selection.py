import numpy as np
from .base_strategy import SelectionStrategy


class GreedySelection(SelectionStrategy):
    """Standard greedy selection (argmax/argmin)"""

    def select_reagent(self, reagent_list, disallow_mask=None, **kwargs):
        rng = kwargs.get('rng', np.random.default_rng())
        scores = self.prepare_scores(reagent_list, rng)

        if disallow_mask:
            scores[np.array(list(disallow_mask))] = np.nan

        if self.mode in ["maximize", "maximize_boltzmann"]:
            return np.nanargmax(scores)
        else:
            return np.nanargmin(scores)
