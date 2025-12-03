import numpy as np
from .base_strategy import SelectionStrategy


class EpsilonGreedySelection(SelectionStrategy):
    """Epsilon-greedy selection with decaying epsilon"""

    def __init__(self, mode="maximize", epsilon=0.1, decay=0.995):
        super().__init__(mode)
        self.epsilon = epsilon
        self.initial_epsilon = epsilon
        self.decay = decay

    def select_reagent(self, reagent_list, disallow_mask=None, **kwargs):
        rng = kwargs.get('rng', np.random.default_rng())
        iteration = kwargs.get('iteration', 0)

        # Decay epsilon over time
        current_epsilon = self.initial_epsilon * (self.decay ** iteration)
        epsilon_scaling = kwargs.get('epsilon_scaling', 1.0) 
        current_epsilon = current_epsilon * epsilon_scaling  

        # Explore with probability epsilon
        if rng.random() < current_epsilon:
            # Random selection (exploration)
            valid_indices = [i for i in range(len(reagent_list))
                           if disallow_mask is None or i not in disallow_mask]
            return rng.choice(valid_indices)
        else:
            # Greedy selection (exploitation)
            scores = self.prepare_scores(reagent_list, rng)
            if disallow_mask:
                scores[np.array(list(disallow_mask))] = np.nan
            return np.nanargmax(scores) if self.mode == "maximize" else np.nanargmin(scores)
