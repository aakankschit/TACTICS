import numpy as np
from .base_strategy import SelectionStrategy


class RouletteWheelSelection(SelectionStrategy):
    """
    Roulette wheel selection with adaptive thermal cycling.

    Implements the enhanced Thompson sampling approach from:
    Zhao, H., Nittinger, E. & Tyrchan, C. Enhanced Thompson Sampling by Roulette
    Wheel Selection for Screening Ultra-Large Combinatorial Libraries.
    bioRxiv 2024.05.16.594622 (2024)
    """

    def __init__(self, mode="maximize", alpha=0.1, beta=None, scaling=1.0,
                 alpha_increment=0.01, beta_increment=0.001, efficiency_threshold=0.1):
        super().__init__(mode)
        self.initial_alpha = alpha
        self.initial_beta = beta if beta is not None else alpha
        self.alpha = alpha
        self.beta = beta if beta is not None else alpha
        self.scaling = scaling
        self.alpha_increment = alpha_increment
        self.beta_increment = beta_increment
        self.efficiency_threshold = efficiency_threshold
        self.current_component_idx = 0

    def select_reagent(self, reagent_list, disallow_mask=None, **kwargs):
        """
        Select a single reagent using roulette wheel selection.

        For batch selection, use select_batch() instead for better efficiency.
        """
        rng = kwargs.get('rng', np.random.default_rng())
        component_idx = kwargs.get('component_idx', 0)

        # Sample base scores
        stds = np.array([r.std for r in reagent_list])
        mu = np.array([r.mean for r in reagent_list])
        scores = rng.normal(size=len(reagent_list)) * stds + mu

        # CRITICAL: Invert scores for minimize mode before thermal cycling
        # In minimize mode, better scores are more negative, so we negate them
        # to make better scores positive before the exp transform
        if self.mode not in ["maximize", "maximize_boltzmann"]:
            scores = -scores

        # Apply thermal cycling
        if component_idx == self.current_component_idx:
            # Heat up - increase exploration
            scores = np.exp((scores - np.mean(scores)) / np.std(scores) / self.alpha)
        else:
            # Cool down - focus on exploitation
            scores = np.exp((scores - np.mean(scores)) / np.std(scores) / self.beta)

        # Normalize to probabilities
        probs = scores / np.sum(scores)

        # Apply disallow mask
        if disallow_mask:
            probs[np.array(list(disallow_mask))] = 0
            if np.sum(probs) > 0:
                probs = probs / np.sum(probs)  # Renormalize
            else:
                # All reagents disallowed - uniform fallback
                probs = np.ones(len(reagent_list)) / len(reagent_list)

        return np.random.choice(len(reagent_list), p=probs)

    def select_batch(self, reagent_list, batch_size, disallow_mask=None, **kwargs):
        """
        Select multiple reagents using roulette wheel selection (batch mode).

        This is more efficient than calling select_reagent multiple times
        as it computes probabilities once and samples multiple times.
        """
        rng = kwargs.get('rng', np.random.default_rng())
        component_idx = kwargs.get('component_idx', 0)

        # Sample base scores
        stds = np.array([r.std for r in reagent_list])
        mu = np.array([r.mean for r in reagent_list])
        scores = rng.normal(size=len(reagent_list)) * stds + mu

        # CRITICAL: Invert scores for minimize mode before thermal cycling
        # In minimize mode, better scores are more negative, so we negate them
        # to make better scores positive before the exp transform
        if self.mode not in ["maximize", "maximize_boltzmann"]:
            scores = -scores

        # Apply thermal cycling
        if component_idx == self.current_component_idx:
            # Heat up - increase exploration
            scores = np.exp((scores - np.mean(scores)) / np.std(scores) / self.alpha)
        else:
            # Cool down - focus on exploitation
            scores = np.exp((scores - np.mean(scores)) / np.std(scores) / self.beta)

        # Normalize to probabilities
        probs = scores / np.sum(scores)

        # Apply disallow mask
        if disallow_mask:
            probs[np.array(list(disallow_mask))] = 0
            if np.sum(probs) > 0:
                probs = probs / np.sum(probs)  # Renormalize
            else:
                # All reagents disallowed - uniform fallback
                probs = np.ones(len(reagent_list)) / len(reagent_list)

        # Sample batch_size reagents with replacement
        return np.random.choice(len(reagent_list), size=batch_size, p=probs)

    def update_temperature(self, n_unique: int, batch_size: int):
        """
        Adaptively adjust temperature based on sampling efficiency.

        Parameters:
        -----------
        n_unique : int
            Number of unique compounds generated in this cycle
        batch_size : int
            Total number of compounds sampled
        """
        efficiency = n_unique / batch_size if batch_size > 0 else 0

        # If efficiency falls below threshold, increase exploration
        if efficiency < self.efficiency_threshold:
            self.alpha += self.alpha_increment

            # If no unique compounds at all, also adjust beta
            if n_unique == 0:
                self.beta += self.beta_increment

    def rotate_component(self, n_components: int):
        """
        Rotate to the next component for thermal cycling.

        Parameters:
        -----------
        n_components : int
            Total number of reagent components
        """
        self.current_component_idx = (self.current_component_idx + 1) % n_components

    def reset_temperature(self):
        """Reset temperature parameters to initial values."""
        self.alpha = self.initial_alpha
        self.beta = self.initial_beta
