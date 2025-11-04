"""
Adaptive Bayes-UCB Selection Strategy with Thermal Cycling.

This module implements an adaptive Bayes Upper Confidence Bound (UCB) selection
strategy that integrates thermal cycling for combinatorial library screening.

The key innovation is using adaptive percentile parameters as an analog to
temperature in thermal cycling:
- Higher percentile (e.g., 0.95) → wider confidence bounds → more exploration
- Lower percentile (e.g., 0.60) → tighter bounds → more exploitation

References:
    Kaufmann, E., Cappé, O., & Garivier, A. (2012). On Bayesian upper confidence
    bounds for bandit problems. In AISTATS.

    Zhao, H., Nittinger, E. & Tyrchan, C. Enhanced Thompson Sampling by Roulette
    Wheel Selection for Screening Ultra-Large Combinatorial Libraries.
    bioRxiv 2024.05.16.594622 (2024)
"""

import numpy as np
from scipy import stats
from .base_strategy import SelectionStrategy


class BayesUCBSelection(SelectionStrategy):
    """
    Adaptive Bayes-UCB selection with thermal cycling.

    This strategy computes Upper Confidence Bound indices using Bayesian posteriors
    (Student-t quantiles) and adapts exploration/exploitation dynamically based on
    sampling efficiency.

    The algorithm maintains two percentile levels:
    - p_high: For "heated" component (more exploration)
    - p_low: For "cooled" components (more exploitation)

    One component is heated per cycle (thermal cycling), and percentiles adapt
    based on whether the search is stuck (low sampling efficiency) or making
    progress (high sampling efficiency).

    Parameters
    ----------
    mode : str, default="maximize"
        Either "maximize" or "minimize" for optimization objective
    initial_p_high : float, default=0.90
        Initial percentile for heated component
    initial_p_low : float, default=0.90
        Initial percentile for cooled components
    efficiency_threshold : float, default=0.10
        Sampling efficiency threshold (0.10 = 10% unique compounds)
    p_high_bounds : tuple, default=(0.85, 0.995)
        Lower and upper bounds for p_high
    p_low_bounds : tuple, default=(0.50, 0.90)
        Lower and upper bounds for p_low
    delta_high : float, default=0.01
        Step size for increasing p_high when stuck
    delta_low : float, default=0.005
        Step size for decreasing p_low when making progress

    Attributes
    ----------
    p_high : float
        Current high percentile (for heated component)
    p_low : float
        Current low percentile (for cooled components)
    current_component_idx : int
        Index of currently heated component

    Examples
    --------
    >>> strategy = BayesUCBSelection(
    ...     mode="maximize",
    ...     initial_p_high=0.90,
    ...     initial_p_low=0.90,
    ...     efficiency_threshold=0.10
    ... )
    >>> sampler = ThompsonSampler(selection_strategy=strategy)
    >>> # ... configure and run sampling ...
    """

    def __init__(self,
                 mode="maximize",
                 initial_p_high=0.90,
                 initial_p_low=0.90,
                 efficiency_threshold=0.10,
                 p_high_bounds=(0.85, 0.995),
                 p_low_bounds=(0.50, 0.90),
                 delta_high=0.01,
                 delta_low=0.005):
        super().__init__(mode)

        # Store initial values for reset
        self.initial_p_high = initial_p_high
        self.initial_p_low = initial_p_low

        # Current percentile values
        self.p_high = initial_p_high
        self.p_low = initial_p_low

        # Adaptation parameters
        self.efficiency_threshold = efficiency_threshold
        self.p_high_bounds = p_high_bounds
        self.p_low_bounds = p_low_bounds
        self.delta_high = delta_high
        self.delta_low = delta_low

        # Thermal cycling state
        self.current_component_idx = 0

    def select_reagent(self, reagent_list, disallow_mask=None, **kwargs):
        """
        Select a single reagent using Bayes-UCB indices.

        Computes UCB index for each reagent based on posterior distribution
        (mean, std, n_samples) and current percentile, then selects via argmax.

        Parameters
        ----------
        reagent_list : List[Reagent]
            List of reagent objects with posterior distributions
        disallow_mask : set, optional
            Indices to exclude from selection
        **kwargs : dict
            Additional parameters:
            - component_idx : int, which component is being selected
            - rng : np.random.Generator (not used but kept for API compatibility)
            - iteration : int (not used but kept for API compatibility)

        Returns
        -------
        int
            Index of selected reagent
        """
        component_idx = kwargs.get('component_idx', 0)

        # Determine percentile based on thermal cycling
        if component_idx == self.current_component_idx:
            percentile = self.p_high  # Heated component - more exploration
        else:
            percentile = self.p_low   # Cooled component - more exploitation

        # Compute UCB indices for all reagents
        ucb_indices = self._compute_ucb_indices(reagent_list, percentile)

        # Apply disallow mask
        if disallow_mask:
            ucb_indices = ucb_indices.copy()
            if self.mode == "maximize":
                ucb_indices[np.array(list(disallow_mask))] = -np.inf
            else:
                ucb_indices[np.array(list(disallow_mask))] = np.inf

        # Select reagent with best UCB index
        if self.mode == "maximize":
            return np.argmax(ucb_indices)
        else:
            return np.argmin(ucb_indices)

    def select_batch(self, reagent_list, batch_size, disallow_mask=None, **kwargs):
        """
        Select multiple reagents using Bayes-UCB indices (batch mode).

        Note: This implementation samples with replacement by calling select_reagent
        multiple times. The percentile and thermal state remain constant within
        a batch.

        Parameters
        ----------
        reagent_list : List[Reagent]
            List of reagent objects with posterior distributions
        batch_size : int
            Number of reagents to select
        disallow_mask : set, optional
            Indices to exclude from selection
        **kwargs : dict
            Additional parameters passed to select_reagent

        Returns
        -------
        np.ndarray
            Array of selected reagent indices
        """
        return np.array([
            self.select_reagent(reagent_list, disallow_mask, **kwargs)
            for _ in range(batch_size)
        ])

    def _compute_ucb_indices(self, reagent_list, percentile):
        """
        Compute Bayes-UCB indices for all reagents.

        Uses Student-t quantiles for proper Bayesian treatment:
        UCB_i = μ_i + σ_i * t_{df}(percentile) / sqrt(n_i)

        where:
        - μ_i: posterior mean
        - σ_i: posterior standard deviation
        - df = n_i - 1: degrees of freedom
        - t_{df}(percentile): Student-t quantile at given percentile
        - n_i: number of observations

        For reagents with n < 2, uses a conservative large bonus.

        Parameters
        ----------
        reagent_list : List[Reagent]
            List of reagent objects
        percentile : float
            Confidence level (e.g., 0.95 for 95th percentile)

        Returns
        -------
        np.ndarray
            UCB indices for all reagents
        """
        # Vectorized computation for speed
        n_reagents = len(reagent_list)
        ucb_indices = np.zeros(n_reagents)

        # Extract arrays
        means = np.array([r.mean for r in reagent_list])
        stds = np.array([r.std for r in reagent_list])
        n_samples = np.array([r.n_samples for r in reagent_list])

        # Handle under-explored reagents (n < 2)
        unexplored_mask = n_samples < 2
        if self.mode == "maximize":
            ucb_indices[unexplored_mask] = means[unexplored_mask] + 3.0 * stds[unexplored_mask]
        else:
            ucb_indices[unexplored_mask] = means[unexplored_mask] - 3.0 * stds[unexplored_mask]

        # Handle explored reagents (n >= 2)
        explored_mask = ~unexplored_mask
        if np.any(explored_mask):
            explored_n = n_samples[explored_mask]
            explored_means = means[explored_mask]
            explored_stds = stds[explored_mask]

            # Compute t-quantiles (still in loop but only for unique df values)
            # Group by degrees of freedom to minimize ppf calls
            unique_dfs = np.unique(explored_n - 1)
            t_quantiles = np.zeros(len(explored_n))

            for df in unique_dfs:
                mask = (explored_n - 1) == df
                t_quantiles[mask] = stats.t.ppf(percentile, df)

            # Compute UCB indices
            if self.mode == "maximize":
                ucb_indices[explored_mask] = explored_means + explored_stds * t_quantiles / np.sqrt(explored_n)
            else:
                ucb_indices[explored_mask] = explored_means - explored_stds * t_quantiles / np.sqrt(explored_n)

        return ucb_indices

    def update_percentiles(self, n_unique, batch_size):
        """
        Adaptively adjust percentiles based on sampling efficiency.

        When sampling efficiency drops below threshold (search is stuck):
        - Increase p_high (more exploration for heated component)
        - Slightly decrease p_low if no unique compounds at all

        When sampling efficiency is above threshold (making progress):
        - Slightly decrease p_high (can be more focused)
        - Decrease p_low (more exploitation)

        This adaptive mechanism helps escape local optima while efficiently
        exploiting promising regions.

        Parameters
        ----------
        n_unique : int
            Number of unique compounds generated in this cycle
        batch_size : int
            Total number of compounds sampled
        """
        efficiency = n_unique / batch_size if batch_size > 0 else 0

        if efficiency < self.efficiency_threshold:
            # Stuck - increase exploration
            self.p_high += self.delta_high

            # If completely stuck (no unique compounds), also adjust p_low
            if n_unique == 0:
                self.p_low = max(self.p_low - self.delta_low * 0.5, self.p_low_bounds[0])
        else:
            # Making progress - can be more greedy
            # Decrease exploration slowly
            self.p_high = max(self.p_high - self.delta_high * 0.2, self.p_high_bounds[0])
            self.p_low = max(self.p_low - self.delta_low, self.p_low_bounds[0])

        # Enforce bounds
        self.p_high = np.clip(self.p_high, self.p_high_bounds[0], self.p_high_bounds[1])
        self.p_low = np.clip(self.p_low, self.p_low_bounds[0], self.p_low_bounds[1])

    def rotate_component(self, n_components):
        """
        Rotate to the next component for thermal cycling.

        This cycles through components, heating one component at a time
        while keeping others cooled.

        Parameters
        ----------
        n_components : int
            Total number of reagent components
        """
        self.current_component_idx = (self.current_component_idx + 1) % n_components

    def reset_percentiles(self):
        """
        Reset percentile parameters to initial values.

        Useful for multi-run experiments or when restarting search.
        """
        self.p_high = self.initial_p_high
        self.p_low = self.initial_p_low
