import numpy as np
from .base_strategy import SelectionStrategy


class GreedySelection(SelectionStrategy):
    """Standard greedy selection (argmax/argmin) with optional CATS variance scaling"""

    def select_reagent(self, reagent_list, disallow_mask=None, **kwargs):
        """
        Select a reagent using greedy Thompson Sampling.
        
        Args:
            reagent_list: List of Reagent objects
            disallow_mask: Set of indices to exclude
            **kwargs: Additional parameters including:
                - rng: Random number generator
                - variance_scaling: Optional float from CATS (component-specific scaling)
        """
        rng = kwargs.get('rng', np.random.default_rng())
        
        # CATS MODIFICATION: Get variance scaling for this component
        variance_scaling = kwargs.get('variance_scaling', None)
        
        # Sample from posteriors
        stds = np.array([r.std for r in reagent_list])
        mu = np.array([r.mean for r in reagent_list])
        
        # CATS: Apply variance scaling if provided
        if variance_scaling is not None:
            # variance_scaling is a single float for this component
            # Lower values = less exploration (critical components)
            # Higher values = more exploration (flexible components)
            stds = stds * np.sqrt(variance_scaling)
        
        scores = rng.normal(size=len(reagent_list)) * stds + mu

        if disallow_mask:
            scores[np.array(list(disallow_mask))] = np.nan

        if self.mode in ["maximize", "maximize_boltzmann"]:
            return np.nanargmax(scores)
        else:
            return np.nanargmin(scores)