import numpy as np
from .base_strategy import SelectionStrategy


class UCBSelection(SelectionStrategy):
    """Upper Confidence Bound selection"""

    def __init__(self, mode="maximize", c=2.0):
        super().__init__(mode)
        self.c = c  # Exploration parameter

    def select_reagent(self, reagent_list, disallow_mask=None, **kwargs):
        total_samples = kwargs.get('iteration', 1)

        ucb_scores = np.zeros(len(reagent_list))
        for i, reagent in enumerate(reagent_list):
            if reagent.n_samples == 0:
                ucb_scores[i] = np.inf  # Force exploration of unsampled reagents
            else:
                exploration = self.c * np.sqrt(np.log(total_samples) / reagent.n_samples)
                if self.mode == "maximize":
                    ucb_scores[i] = reagent.mean + exploration
                else:
                    ucb_scores[i] = reagent.mean - exploration

        if disallow_mask:
            ucb_scores[np.array(list(disallow_mask))] = -np.inf if self.mode == "maximize" else np.inf

        return np.argmax(ucb_scores) if self.mode == "maximize" else np.argmin(ucb_scores)
