"""
Bayesian Upper Confidence Bound Selection with Component-Aware Thompson Sampling (CATS).

This module implements Bayes-UCB with CATS for combinatorial library screening.
Uses Student-t quantiles for proper Bayesian treatment with percentile-based
thermal cycling and component criticality analysis.

References:
    Kaufmann, E., Cappé, O., & Garivier, A. (2012). On Bayesian upper confidence
    bounds for bandit problems. In AISTATS.

    Zhao, H., Nittinger, E. & Tyrchan, C. Enhanced Thompson Sampling by Roulette
    Wheel Selection for Screening Ultra-Large Combinatorial Libraries.
    bioRxiv 2024.05.16.594622 (2024)
"""

import warnings
import numpy as np
from scipy import stats
from .base_strategy import SelectionStrategy


class BayesUCBSelection(SelectionStrategy):
    """
    Bayesian Upper Confidence Bound selection with Component-Aware Thompson Sampling (CATS).

    Combines percentile-based thermal cycling with component criticality analysis
    for efficient exploration of ultra-large combinatorial libraries.

    Uses Student-t quantiles for proper Bayesian treatment of uncertainty.
    Percentile levels serve as analog to temperature in RWS:
    - Higher percentile → wider confidence bounds → more exploration
    - Lower percentile → tighter bounds → more exploitation
    """

    def __init__(
        self,
        mode="maximize",
        initial_p_high=0.90,
        initial_p_low=0.60,
        exploration_phase_end=0.20,
        transition_phase_end=0.60,
        min_observations=5,
        cats_exploration_fraction=0.3,
        **kwargs
    ):
        """
        Initialize Bayes-UCB Selection with CATS.

        Args:
            mode: "maximize" or "minimize" optimization mode
            initial_p_high: Base percentile for heated component (default: 0.90)
            initial_p_low: Base percentile for cooled component (default: 0.60)
            exploration_phase_end: Fraction of iterations before CATS starts (default: 0.20)
            transition_phase_end: Fraction of iterations when CATS is fully applied (default: 0.60)
            min_observations: Minimum observations per reagent before trusting criticality (default: 5)
            cats_exploration_fraction: Fraction of total cycles during which CATS explores
                at full strength. After this point, CATS influence decays linearly if
                criticality remains low. Set to None to disable decay (default: 0.5).
            **kwargs: Catches deprecated parameters with warnings
        """
        super().__init__(mode)

        # Core parameters
        self.initial_p_high = initial_p_high
        self.initial_p_low = initial_p_low
        self.p_high = initial_p_high
        self.p_low = initial_p_low

        # CATS parameters
        self.exploration_phase_end = exploration_phase_end
        self.transition_phase_end = transition_phase_end
        self.min_observations = min_observations
        self.cats_exploration_fraction = cats_exploration_fraction

        # Thermal cycling state
        self.current_component_idx = 0

        # Derive CATS range from p_high/p_low ratio
        self.ratio = self.p_high / self.p_low if self.p_low > 0 else 1.0
        self.cats_max_mult = self.ratio
        self.cats_min_mult = 1.0 / self.ratio if self.ratio > 0 else 1.0

        # Validation warnings
        if abs(self.p_high - self.p_low) < 1e-10:
            warnings.warn(
                f"initial_p_high ({self.p_high}) equals initial_p_low ({self.p_low}). "
                "Thermal cycling will have no effect. "
                "For thermal cycling to work, set initial_p_low < initial_p_high "
                "(e.g., initial_p_low=0.60, initial_p_high=0.90).",
                UserWarning,
                stacklevel=2
            )

        # Deprecated parameter warnings
        deprecated = {'efficiency_threshold', 'p_high_bounds', 'p_low_bounds', 'delta_high', 'delta_low'}
        found_deprecated = deprecated & set(kwargs.keys())
        if found_deprecated:
            warnings.warn(
                f"Parameters {found_deprecated} are deprecated and will be ignored. "
                f"CATS now derives all parameters from initial_p_high/initial_p_low ratio. "
                f"Remove these parameters from your configuration.",
                DeprecationWarning,
                stacklevel=2
            )

    def _calculate_criticality(self, reagent_list):
        """
        Calculate component criticality using Shannon entropy.

        Shannon entropy measures posterior distribution concentration:
        - High entropy (uniform) → Flexible component (many good options)
        - Low entropy (peaked) → Critical component (few good options)

        Criticality = 1 - (entropy / max_entropy) ∈ [0, 1]
        - criticality ≈ 0: Flexible component → should EXPLORE
        - criticality ≈ 1: Critical component → should EXPLOIT

        Args:
            reagent_list: List of Reagent objects with posterior distributions

        Returns:
            Criticality score in [0, 1]. Returns 0.5 (neutral) if insufficient data.
        """
        # Check if we have sufficient data
        observations = [r.n_samples for r in reagent_list]
        if min(observations) < self.min_observations:
            return 0.5  # Neutral criticality if insufficient data

        # Extract posterior means
        means = np.array([r.mean for r in reagent_list])

        # Handle edge case: all means identical
        if np.std(means) < 1e-10:
            return 0.0  # Perfectly flexible (all equally good)

        # Convert means to probability distribution using softmax
        # For minimization, negate means (want higher probability for lower scores)
        if self.mode == "minimize":
            means = -means

        # Numerical stability: subtract max before exp
        exp_means = np.exp(means - means.max())
        probabilities = exp_means / exp_means.sum()

        # Shannon entropy: H = -Σ p_i log(p_i)
        entropy = -np.sum(probabilities * np.log(probabilities + 1e-10))

        # Maximum entropy (uniform distribution)
        max_entropy = np.log(len(reagent_list))

        # Criticality: 1 - normalized_entropy
        if max_entropy < 1e-10:
            return 0.5  # Edge case: single reagent

        normalized_entropy = entropy / max_entropy
        criticality = 1.0 - normalized_entropy

        return np.clip(criticality, 0.0, 1.0)

    def _get_criticality_weight(self, reagent_list, current_cycle=None, total_cycles=None):
        """
        Calculate observation-gated criticality weight.

        Instead of a fixed three-phase schedule, weight is determined by how
        many observations we have — i.e. how much we can trust the criticality
        estimate. This is data-driven rather than schedule-driven.

        Weight ramps from 0 to 1 as min observations go from 0 to
        2 * min_observations, then stays at 1.0.

        Args:
            reagent_list: List of Reagent objects (used to check observation counts)
            current_cycle: Unused, kept for API compatibility
            total_cycles: Unused, kept for API compatibility

        Returns:
            Criticality weight in [0, 1]
        """
        min_obs = min(r.n_samples for r in reagent_list)
        weight = min_obs / (2.0 * self.min_observations) if self.min_observations > 0 else 1.0
        return float(np.clip(weight, 0.0, 1.0))

    def _get_cats_multiplier(self, criticality):
        """
        Map component criticality to percentile multiplier.

        CATS range is derived from p_high/p_low ratio:
        - cats_max = p_high / p_low (for flexible components, criticality ≈ 0)
        - cats_min = p_low / p_high (for critical components, criticality ≈ 1)

        Args:
            criticality: Component criticality in [0, 1]

        Returns:
            CATS multiplier in [cats_min, cats_max]
        """
        # Map criticality [0,1] to multiplier [cats_max, cats_min]
        # Higher criticality → lower multiplier (more exploitation)
        multiplier = self.cats_min_mult + (self.cats_max_mult - self.cats_min_mult) * (1 - criticality)
        return multiplier

    def _get_component_percentile(self, component_idx, reagent_list, current_cycle, total_cycles):
        """
        Get CATS-adjusted percentile for a component.

        Combines thermal cycling with component criticality:
        1. Thermal cycling sets base percentile (p_high for heated, p_low for cooled)
        2. CATS adjusts based on criticality
        3. Progressive weighting controls CATS influence
        4. Exploration decay reduces CATS when criticality stays low

        Args:
            component_idx: Which reaction component
            reagent_list: List of Reagent objects for this component
            current_cycle: Current search cycle
            total_cycles: Total number of cycles

        Returns:
            Final percentile value
        """
        # Step 1: Base percentile from thermal cycling
        base_percentile = self.p_high if component_idx == self.current_component_idx else self.p_low

        # Step 2: Calculate criticality
        criticality = self._calculate_criticality(reagent_list)

        # Step 3: Get observation-gated weight (ramps up with data)
        weight = self._get_criticality_weight(reagent_list)

        # Step 4: Exploration decay — reduce CATS when it hasn't found structure
        if (
            self.cats_exploration_fraction is not None
            and total_cycles > 0
            and current_cycle is not None
        ):
            exploration_end = self.cats_exploration_fraction * total_cycles
            if current_cycle > exploration_end:
                remaining = total_cycles - exploration_end
                progress = (current_cycle - exploration_end) / max(remaining, 1)
                decay = criticality + (1.0 - criticality) * (1.0 - progress)
                weight *= decay

        # Step 5: Get CATS multiplier
        cats_mult = self._get_cats_multiplier(criticality)

        # Step 6: Blend: weight=0 → no adjustment, weight=1 → full CATS
        effective_mult = (1.0 - weight) * 1.0 + weight * cats_mult

        # Step 7: Apply to base
        final_percentile = base_percentile * effective_mult

        # Clamp to valid percentile range; 0.999 avoids t.ppf(1.0)=inf
        return np.clip(final_percentile, 0.5, 0.999)

    def select_reagent(self, reagent_list, disallow_mask=None, **kwargs):
        """
        Select a single reagent using Bayes-UCB indices with CATS.

        Computes UCB index for each reagent based on posterior distribution
        (mean, std, n_samples) and CATS-adjusted percentile, then selects via argmax.

        Args:
            reagent_list: List of Reagent objects with posterior distributions
            disallow_mask: Optional set of indices to exclude from selection
            **kwargs: Additional context:
                - component_idx: Which reaction component
                - current_cycle: Current search cycle (for CATS)
                - total_cycles: Total number of cycles (for CATS)
                - rng: Random number generator (not used but kept for API compatibility)

        Returns:
            Index of selected reagent
        """
        component_idx = kwargs.get('component_idx', 0)
        current_cycle = kwargs.get('current_cycle', 0)
        total_cycles = kwargs.get('total_cycles', 1)

        # Get CATS-adjusted percentile
        percentile = self._get_component_percentile(
            component_idx, reagent_list, current_cycle, total_cycles
        )

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
        Select multiple reagents using Bayes-UCB indices with CATS (batch mode).

        Note: This implementation samples with replacement by calling select_reagent
        multiple times. The CATS-adjusted percentile and thermal state remain constant
        within a batch.

        Args:
            reagent_list: List of Reagent objects with posterior distributions
            batch_size: Number of reagents to select
            disallow_mask: Optional set of indices to exclude from selection
            **kwargs: Additional context passed to select_reagent:
                - component_idx: Which reaction component
                - current_cycle: Current search cycle (for CATS)
                - total_cycles: Total number of cycles (for CATS)

        Returns:
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
            ucb_indices[unexplored_mask] = means[unexplored_mask] + 3.0 * np.maximum(stds[unexplored_mask], 1e-6)
        else:
            ucb_indices[unexplored_mask] = means[unexplored_mask] - 3.0 * np.maximum(stds[unexplored_mask], 1e-6)

        # Handle explored reagents (n >= 2)
        explored_mask = ~unexplored_mask
        if np.any(explored_mask):
            explored_n = n_samples[explored_mask]
            explored_means = means[explored_mask]
            explored_stds = np.maximum(stds[explored_mask], 1e-6)  # Avoid zero std

            # Compute t-quantiles (still in loop but only for unique df values)
            # Group by degrees of freedom to minimize ppf calls
            unique_dfs = np.unique(explored_n - 1)
            t_quantiles = np.zeros(len(explored_n))

            for df in unique_dfs:
                mask = (explored_n - 1) == df
                # Clamp df to avoid numerical issues with very small degrees of freedom
                safe_df = max(df, 1)
                t_quantiles[mask] = stats.t.ppf(percentile, safe_df)

            # Compute UCB indices with numerical stability
            sqrt_n = np.sqrt(np.maximum(explored_n, 1))
            if self.mode == "maximize":
                ucb_indices[explored_mask] = explored_means + explored_stds * t_quantiles / sqrt_n
            else:
                ucb_indices[explored_mask] = explored_means - explored_stds * t_quantiles / sqrt_n

        return ucb_indices

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
