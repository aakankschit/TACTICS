import numpy as np
from .base_strategy import SelectionStrategy


class UCBSelection(SelectionStrategy):
    """Upper Confidence Bound selection"""

    def __init__(self, mode="maximize", c=2.0):
        super().__init__(mode)
        self.c = c  # Exploration parameter

    def select_reagent(self, reagent_list, disallow_mask=None, **kwargs):
        total_samples = max(kwargs.get('iteration', 1), 1)  # Ensure at least 1
        exploration_factor = kwargs.get('exploration_factor', 1.0)

        ucb_scores = np.zeros(len(reagent_list))
        for i, reagent in enumerate(reagent_list):
            if reagent.n_samples == 0:
                ucb_scores[i] = np.inf  # Force exploration of unsampled reagents
            else:
                # Add numerical stability: ensure log argument is > 1
                # and avoid division issues when n_samples is large
                log_term = max(np.log(max(total_samples, 2)), 0.01)
                ratio = log_term / max(reagent.n_samples, 1)
                base_exploration = self.c * np.sqrt(ratio)
                # Apply exploration_factor with minimum floor to maintain diversity
                # Even critical components (low factor) get 40% of base exploration
                exploration = base_exploration * (0.4 + 0.6 * exploration_factor)

                if self.mode == "maximize":
                    ucb_scores[i] = reagent.mean + exploration
                else:
                    ucb_scores[i] = reagent.mean - exploration

        if disallow_mask:
            ucb_scores[np.array(list(disallow_mask))] = -np.inf if self.mode == "maximize" else np.inf

        return np.argmax(ucb_scores) if self.mode == "maximize" else np.argmin(ucb_scores)
