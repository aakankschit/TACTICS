"""
Interaction detector for combinatorial libraries.

Detects whether interaction effects dominate the score landscape by fitting
an additive model and computing R². When the additive model explains poorly
(low R²), per-component criticality estimates are unreliable and CATS
influence should be dampened.
"""

import logging
from typing import List, Tuple

import numpy as np

logger = logging.getLogger(__name__)


class InteractionDetector:
    """
    Detect interaction effects in combinatorial score landscapes.

    Periodically fits an additive model to observed scores and reports
    the fraction of variance explained. High R² means the additive
    assumption holds (CATS is reliable); low R² means interaction
    effects dominate (CATS should be dampened).

    Parameters:
        n_components: Number of reagent components
        component_sizes: List of reagent counts per component
        activation_fraction: Don't start fitting until this fraction of
            total_cycles has passed (default: 0.3)
        refit_interval: Re-fit the model every N observations (default: 200)
        interaction_threshold: R² below this triggers dampening (default: 0.5)
    """

    def __init__(
        self,
        n_components: int,
        component_sizes: List[int],
        activation_fraction: float = 0.3,
        refit_interval: int = 200,
        interaction_threshold: float = 0.5,
    ):
        self.n_components = n_components
        self.component_sizes = component_sizes
        self.activation_fraction = activation_fraction
        self.refit_interval = refit_interval
        self.interaction_threshold = interaction_threshold

        # Observation storage
        self._combinations: List[Tuple[int, ...]] = []
        self._scores: List[float] = []

        # State
        self.interaction_strength: float = 0.0
        self.additive_r2: float = 1.0
        self._last_fit_count: int = 0
        self._min_observations = 50

    def add_observation(self, combination_indices: List[int], score: float) -> None:
        """
        Record an observation.

        Parameters:
            combination_indices: Reagent indices for each component
            score: Observed score
        """
        self._combinations.append(tuple(combination_indices))
        self._scores.append(score)

    def should_refit(self, current_cycle: int, total_cycles: int) -> bool:
        """
        Check whether the model should be re-fit.

        Returns True when:
        1. We've passed the activation fraction of total cycles
        2. We have at least min_observations
        3. At least refit_interval new observations since last fit

        Parameters:
            current_cycle: Current search cycle
            total_cycles: Total number of cycles

        Returns:
            True if we should re-fit
        """
        # Check activation fraction
        if total_cycles > 0 and current_cycle < self.activation_fraction * total_cycles:
            return False

        n_obs = len(self._scores)

        # Minimum observations
        if n_obs < self._min_observations:
            return False

        # Refit interval
        if n_obs - self._last_fit_count < self.refit_interval:
            return False

        return True

    def fit(self) -> float:
        """
        Fit additive model and compute R².

        For K components, the additive model predicts:
            score_hat = sum_k(mean_k[i_k]) - (K-1) * global_mean

        where mean_k[i_k] is the average score for reagent i_k in component k.

        Returns:
            R² of the additive model (0 to 1)
        """
        n_obs = len(self._scores)
        if n_obs < self._min_observations:
            return 1.0  # Insufficient data, assume additive

        scores = np.array(self._scores)
        global_mean = scores.mean()

        # Compute per-reagent means for each component
        component_means = []
        for k in range(self.n_components):
            reagent_sums = {}
            reagent_counts = {}
            for idx, (comb, score) in enumerate(zip(self._combinations, self._scores)):
                r_idx = comb[k]
                reagent_sums[r_idx] = reagent_sums.get(r_idx, 0.0) + score
                reagent_counts[r_idx] = reagent_counts.get(r_idx, 0) + 1

            means = {}
            for r_idx in reagent_sums:
                means[r_idx] = reagent_sums[r_idx] / reagent_counts[r_idx]
            component_means.append(means)

        # Predict using additive model
        predictions = np.zeros(n_obs)
        for i, comb in enumerate(self._combinations):
            pred = 0.0
            for k in range(self.n_components):
                r_idx = comb[k]
                pred += component_means[k].get(r_idx, global_mean)
            pred -= (self.n_components - 1) * global_mean
            predictions[i] = pred

        # Compute R²
        ss_res = np.sum((scores - predictions) ** 2)
        ss_tot = np.sum((scores - global_mean) ** 2)

        if ss_tot < 1e-10:
            r2 = 1.0  # All scores identical
        else:
            r2 = 1.0 - ss_res / ss_tot

        self.additive_r2 = max(r2, 0.0)
        self.interaction_strength = 1.0 - self.additive_r2
        self._last_fit_count = n_obs

        logger.info(
            f"InteractionDetector: R²={self.additive_r2:.3f}, "
            f"interaction_strength={self.interaction_strength:.3f} "
            f"(n_obs={n_obs})"
        )

        return self.additive_r2
