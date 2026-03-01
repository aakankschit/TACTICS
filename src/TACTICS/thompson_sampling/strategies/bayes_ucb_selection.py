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
from typing import Any, Dict, List, Optional
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
        criticality_metric="ipr",
        n_adaptive_sharpening=True,
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
            criticality_metric: "ipr" (Inverse Participation Ratio) or "shannon"
                (legacy Shannon entropy). (default: "ipr")
            n_adaptive_sharpening: If True and criticality_metric="ipr", apply
                sqrt(log(N)) sharpening to z-scores. (default: True)
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
        self.criticality_metric = criticality_metric
        self.n_adaptive_sharpening = n_adaptive_sharpening

        # Thermal cycling state
        self.current_component_idx = 0
        # Mean criticality across components (updated by rotate_component_weighted)
        self._mean_criticality: float = 0.5

        # Derive CATS range from p_high/p_low ratio
        self.ratio = self.p_high / self.p_low if self.p_low > 0 else 1.0
        self.cats_max_mult = self.ratio
        self.cats_min_mult = 1.0 / self.ratio if self.ratio > 0 else 1.0

        # Cache for last component state (avoids double computation)
        self._last_component_states: Dict[int, Dict[str, Any]] = {}

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
        Calculate component criticality using z-score softmax with SNR dampening.

        Aligned with RouletteWheelSelection: uses z-scores (scale-invariant),
        signal-to-noise dampening, N-adaptive sharpening, and supports both
        IPR and Shannon entropy metrics.

        Criticality in [0, 1]:
        - criticality ~ 0: Flexible component -> should EXPLORE
        - criticality ~ 1: Critical component -> should EXPLOIT

        Args:
            reagent_list: List of Reagent objects with posterior distributions

        Returns:
            Criticality score in [0, 1]. Returns 0.5 (neutral) if insufficient data.
        """
        # Filter out retired reagents
        active_reagents = [r for r in reagent_list if r.n_samples > 0]
        if len(active_reagents) < 2:
            return 0.5

        observations = [r.n_samples for r in active_reagents]
        if min(observations) < self.min_observations:
            return 0.5

        means = np.array([r.mean for r in active_reagents])

        if np.std(means) < 1e-10:
            return 0.0

        if self.mode == "minimize":
            means = -means

        mean_std = np.std(means)
        if mean_std < 1e-10:
            return 0.0

        # Signal-to-noise dampening (aligned with RWS)
        se_squared = np.array([
            r.std ** 2 / max(r.n_samples, 1) for r in active_reagents
        ])
        noise_std = np.sqrt(np.mean(se_squared))
        snr = mean_std / max(noise_std, 1e-10)
        imbalance_strength = 1.0 - np.exp(-max(snr - 1.0, 0.0))

        z_scores = (means - np.mean(means)) / mean_std
        z_scores *= imbalance_strength

        # N-adaptive sharpening
        N = len(active_reagents)
        if self.criticality_metric == "ipr" and self.n_adaptive_sharpening and N > 2:
            sharpening = max(1.0, np.sqrt(np.log(N)))
            z_scores *= sharpening

        exp_means = np.exp(z_scores - z_scores.max())  # numerical stability
        probabilities = exp_means / exp_means.sum()

        if self.criticality_metric == "ipr":
            ipr = np.sum(probabilities ** 2)
            effective_N = 1.0 / ipr
            criticality = 1.0 - (effective_N / N)
        else:  # shannon (legacy)
            entropy = -np.sum(probabilities * np.log(probabilities + 1e-10))
            max_entropy = np.log(N)
            if max_entropy < 1e-10:
                return 0.5
            criticality = 1.0 - (entropy / max_entropy)

        return float(np.clip(criticality, 0.0, 1.0))

    def _calculate_criticality_details(self, reagent_list):
        """Calculate criticality and return all intermediate values.

        Returns:
            Tuple of (criticality, details_dict) matching RouletteWheelSelection schema.
        """
        details = {
            "snr": float("nan"),
            "imbalance_strength": float("nan"),
            "normalized_entropy": float("nan"),
            "n_active_reagents": 0,
            "participation_ratio": float("nan"),
            "effective_n": float("nan"),
            "sharpening_factor": float("nan"),
        }

        active_reagents = [r for r in reagent_list if r.n_samples > 0]
        details["n_active_reagents"] = len(active_reagents)

        if len(active_reagents) < 2:
            return 0.5, details

        observations = [r.n_samples for r in active_reagents]
        if min(observations) < self.min_observations:
            return 0.5, details

        means = np.array([r.mean for r in active_reagents])

        if np.std(means) < 1e-10:
            details["snr"] = 0.0
            details["imbalance_strength"] = 0.0
            details["normalized_entropy"] = 1.0
            details["participation_ratio"] = 1.0 / len(active_reagents)
            details["effective_n"] = float(len(active_reagents))
            details["sharpening_factor"] = 1.0
            return 0.0, details

        if self.mode == "minimize":
            means = -means

        mean_std = np.std(means)
        if mean_std < 1e-10:
            details["snr"] = 0.0
            details["imbalance_strength"] = 0.0
            details["normalized_entropy"] = 1.0
            details["participation_ratio"] = 1.0 / len(active_reagents)
            details["effective_n"] = float(len(active_reagents))
            details["sharpening_factor"] = 1.0
            return 0.0, details

        se_squared = np.array([
            r.std ** 2 / max(r.n_samples, 1) for r in active_reagents
        ])
        noise_std = np.sqrt(np.mean(se_squared))
        snr = mean_std / max(noise_std, 1e-10)
        imbalance_strength = 1.0 - np.exp(-max(snr - 1.0, 0.0))

        details["snr"] = float(snr)
        details["imbalance_strength"] = float(imbalance_strength)

        z_scores = (means - np.mean(means)) / mean_std
        z_scores *= imbalance_strength

        N = len(active_reagents)
        sharpening = 1.0
        if self.criticality_metric == "ipr" and self.n_adaptive_sharpening and N > 2:
            sharpening = max(1.0, np.sqrt(np.log(N)))
            z_scores *= sharpening
        details["sharpening_factor"] = float(sharpening)

        exp_means = np.exp(z_scores - z_scores.max())
        probabilities = exp_means / exp_means.sum()

        ipr = float(np.sum(probabilities ** 2))
        effective_n = 1.0 / ipr
        details["participation_ratio"] = ipr
        details["effective_n"] = float(effective_n)

        if self.criticality_metric == "ipr":
            criticality = 1.0 - (effective_n / N)
            entropy = -np.sum(probabilities * np.log(probabilities + 1e-10))
            max_entropy = np.log(N)
            details["normalized_entropy"] = float(entropy / max_entropy) if max_entropy > 1e-10 else float("nan")
        else:
            entropy = -np.sum(probabilities * np.log(probabilities + 1e-10))
            max_entropy = np.log(N)
            if max_entropy < 1e-10:
                details["normalized_entropy"] = float("nan")
                return 0.5, details
            normalized_entropy = entropy / max_entropy
            details["normalized_entropy"] = float(normalized_entropy)
            criticality = 1.0 - normalized_entropy

        return float(np.clip(criticality, 0.0, 1.0)), details

    def get_component_criticality(self, reagent_list) -> float:
        """Return CATS criticality score for a component."""
        return self._calculate_criticality(reagent_list)

    def get_component_state(
        self,
        reagent_list: List,
        component_idx: int,
        current_cycle: int,
        total_cycles: int,
    ) -> Optional[Dict[str, Any]]:
        """Return full intermediate state for a component.

        Same schema as RouletteWheelSelection for DataFrame consistency.
        Now includes SNR dampening (aligned with RWS) and IPR details.
        ``final_temperature`` maps to the CATS-adjusted percentile value.

        Returns:
            Dict with all intermediate values.
        """
        # Step 1: Criticality with details (aligned with RWS)
        criticality, crit_details = self._calculate_criticality_details(reagent_list)

        # Step 2: Base percentile from thermal cycling
        is_heated = component_idx == self.current_component_idx
        base_temp = self.p_high if is_heated else self.p_low

        # Step 3: Observation-gated weight
        weight = self._get_criticality_weight(reagent_list)

        # Step 4: Exploration decay
        decay = 1.0
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

        # Step 5: CATS multiplier
        cats_mult = self._get_cats_multiplier(criticality)

        # Step 6: Blend
        effective_mult = (1.0 - weight) * 1.0 + weight * cats_mult

        # Step 7: Final percentile (mapped to final_temperature for schema)
        final_percentile = base_temp * effective_mult
        final_percentile = float(np.clip(final_percentile, 0.5, 0.999))

        state = {
            # Metadata
            "component_idx": component_idx,
            "current_cycle": current_cycle,
            "total_cycles": total_cycles,
            # Criticality details
            "criticality": criticality,
            "snr": crit_details["snr"],
            "imbalance_strength": crit_details["imbalance_strength"],
            "normalized_entropy": crit_details["normalized_entropy"],
            "n_active_reagents": crit_details["n_active_reagents"],
            # IPR details
            "participation_ratio": crit_details["participation_ratio"],
            "effective_n": crit_details["effective_n"],
            "sharpening_factor": crit_details["sharpening_factor"],
            # Temperature/percentile pipeline
            "base_temp": base_temp,
            "is_heated": is_heated,
            "criticality_weight": weight,
            "decay": decay,
            "cats_multiplier": cats_mult,
            "effective_multiplier": effective_mult,
            "final_temperature": final_percentile,
        }

        self._last_component_states[component_idx] = state
        return state

    def _get_criticality_weight(self, reagent_list, current_cycle=None, total_cycles=None):
        """
        Calculate observation-gated criticality weight.

        Instead of a fixed three-phase schedule, weight is determined by how
        many observations we have — i.e. how much we can trust the criticality
        estimate. This is data-driven rather than schedule-driven.

        Uses the 25th percentile of observation counts rather than the minimum.
        With large components (e.g., 3844 dipeptides), the least-sampled reagent
        may never be re-sampled after warmup, permanently pinning min_obs at the
        warmup count. The 25th percentile reflects the bulk of the distribution,
        allowing CATS to ramp up as the majority of reagents accumulate data.

        Weight ramps from 0 to 1 as p25 observations go from 0 to
        2 * min_observations, then stays at 1.0.

        Args:
            reagent_list: List of Reagent objects (used to check observation counts)
            current_cycle: Unused, kept for API compatibility
            total_cycles: Unused, kept for API compatibility

        Returns:
            Criticality weight in [0, 1]
        """
        # Filter out retired reagents (aligned with RWS)
        active_obs = [r.n_samples for r in reagent_list if r.n_samples > 0]
        if not active_obs:
            return 0.0
        p25_obs = float(np.percentile(active_obs, 25))
        weight = p25_obs / (2.0 * self.min_observations) if self.min_observations > 0 else 1.0
        return float(np.clip(weight, 0.0, 1.0))

    def _get_cats_multiplier(self, criticality):
        """
        Map component criticality to a bidirectional percentile multiplier
        using a relative neutral point.

        The neutral point (multiplier = 1.0) is at the mean criticality across
        all components, updated each cycle by ``rotate_component_weighted``.
        Components below the mean get an exploration boost (multiplier > 1);
        components above the mean get an exploitation reduction (multiplier < 1).

        Mapping (with mean_crit as the neutral point):
            criticality = 0.0        →  mult = cats_max_mult  (max exploration)
            criticality = mean_crit  →  mult = 1.0            (neutral)
            criticality = 1.0        →  mult = cats_min_mult  (max exploitation)

        Args:
            criticality: Component criticality in [0, 1]

        Returns:
            CATS multiplier in [cats_min_mult, cats_max_mult]
        """
        neutral = self._mean_criticality

        if criticality <= neutral:
            # Below mean → explore: interpolate [cats_max_mult, 1.0]
            t = criticality / max(neutral, 1e-10)
            multiplier = self.cats_max_mult + t * (1.0 - self.cats_max_mult)
        else:
            # Above mean → exploit: interpolate [1.0, cats_min_mult]
            t = (criticality - neutral) / max(1.0 - neutral, 1e-10)
            multiplier = 1.0 + t * (self.cats_min_mult - 1.0)

        return multiplier

    def _get_component_percentile(self, component_idx, reagent_list, current_cycle, total_cycles):
        """
        Get CATS-adjusted percentile for a component.

        Delegates to ``get_component_state()`` which computes and caches all
        intermediate values, then returns just the final percentile.

        Args:
            component_idx: Which reaction component
            reagent_list: List of Reagent objects for this component
            current_cycle: Current search cycle
            total_cycles: Total number of cycles

        Returns:
            Final percentile value
        """
        state = self.get_component_state(
            reagent_list, component_idx, current_cycle, total_cycles
        )
        return state["final_temperature"]

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

    def rotate_component_weighted(self, n_components, criticalities, rng=None):
        """
        Rotate to next component using criticality-weighted probabilities.

        Flexible components (low criticality) get heated more often.
        Also caches the mean criticality so that ``_get_cats_multiplier`` can
        use a relative neutral point instead of a fixed 0.5.

        Parameters
        ----------
        n_components : int
            Total number of reagent components
        criticalities : list of float
            Per-component criticality values in [0, 1]
        rng : numpy.random.Generator, optional
            Random number generator for reproducibility
        """
        if rng is None:
            rng = np.random.default_rng()
        crits = np.array(criticalities, dtype=float)
        self._mean_criticality = float(crits.mean())
        flexibility = np.maximum(1.0 - crits, 0.1)
        heat_probs = flexibility / flexibility.sum()
        self.current_component_idx = int(rng.choice(n_components, p=heat_probs))

    def reset_percentiles(self):
        """
        Reset percentile parameters to initial values.

        Useful for multi-run experiments or when restarting search.
        """
        self.p_high = self.initial_p_high
        self.p_low = self.initial_p_low
