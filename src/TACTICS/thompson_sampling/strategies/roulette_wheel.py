import warnings
from typing import Any, Dict, List, Optional
import numpy as np
from .base_strategy import SelectionStrategy


class RouletteWheelSelection(SelectionStrategy):
    """
    Roulette wheel selection with Component-Aware Thompson Sampling (CATS).

    Combines thermal cycling with component criticality analysis for efficient
    exploration of ultra-large combinatorial libraries.

    References:
        Zhao, H., Nittinger, E. & Tyrchan, C. Enhanced Thompson Sampling by Roulette
        Wheel Selection for Screening Ultra-Large Combinatorial Libraries.
        bioRxiv 2024.05.16.594622 (2024)
    """

    def __init__(
        self,
        mode="maximize",
        alpha=0.1,
        beta=0.1,
        exploration_phase_end=0.20,
        transition_phase_end=0.60,
        min_observations=5,
        adaptive_temperature=False,
        alpha_increment=0.01,
        beta_increment=0.001,
        efficiency_threshold=0.10,
        alpha_max=2.0,
        cats_exploration_fraction=0.3,
        cats_range=None,
        criticality_metric="ipr",
        n_adaptive_sharpening=True,
        **kwargs,
    ):
        """
        Initialize Roulette Wheel Selection with CATS.

        Args:
            mode: "maximize" or "minimize" optimization mode
            alpha: Base temperature for heated component (default: 0.1)
            beta: Base temperature for cooled component (default: 0.1)
            exploration_phase_end: Fraction of iterations before CATS starts (default: 0.20)
            transition_phase_end: Fraction of iterations when CATS is fully applied (default: 0.60)
            min_observations: Minimum observations per reagent before trusting criticality (default: 5)
            adaptive_temperature: Enable legacy-inspired adaptive temperature control (default: False)
            alpha_increment: Amount to increase alpha when efficiency drops (default: 0.01)
            beta_increment: Amount to increase beta when zero unique found (default: 0.001)
            efficiency_threshold: Efficiency below which alpha is incremented (default: 0.10)
            alpha_max: Maximum alpha value (default: 2.0)
            cats_exploration_fraction: Fraction of total cycles during which CATS explores
                at full strength. After this point, CATS influence decays linearly if
                criticality remains low. Set to None to disable decay (default: 0.5).
            cats_range: Override alpha/beta-derived CATS multiplier range. If set,
                cats_max = cats_range, cats_min = 1/cats_range. (default: None)
            criticality_metric: "ipr" (Inverse Participation Ratio) or "shannon"
                (legacy Shannon entropy). IPR is more sensitive to probability
                concentration, especially at large N. (default: "ipr")
            n_adaptive_sharpening: If True and criticality_metric="ipr", apply
                sqrt(log(N)) sharpening to z-scores before softmax to counteract
                flattening at large N. (default: True)
            **kwargs: Catches deprecated parameters with warnings
        """
        super().__init__(mode)

        # Core parameters
        self.initial_alpha = alpha
        self.initial_beta = beta
        self.alpha = alpha
        self.beta = beta

        # Adaptive temperature parameters
        self.adaptive_temperature = adaptive_temperature
        self.alpha_increment = alpha_increment
        self.beta_increment = beta_increment
        self.efficiency_threshold = efficiency_threshold
        self.alpha_max = alpha_max

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

        # Derive CATS range from alpha/beta ratio or explicit override
        if cats_range is not None:
            self.ratio = cats_range
            self.cats_max_mult = cats_range
            self.cats_min_mult = 1.0 / cats_range
        else:
            self.ratio = self.alpha / self.beta if self.beta > 0 else 1.0
            self.cats_max_mult = self.ratio
            self.cats_min_mult = 1.0 / self.ratio if self.ratio > 0 else 1.0

        # Cache for last component state (avoids double computation)
        self._last_component_states: Dict[int, Dict[str, Any]] = {}

        # Validation warnings
        if abs(self.alpha - self.beta) < 1e-10:
            warnings.warn(
                f"alpha ({self.alpha}) equals beta ({self.beta}). "
                "Thermal cycling will have no effect. "
                "For thermal cycling to work, set beta < alpha (e.g., beta=0.05, alpha=0.1).",
                UserWarning,
                stacklevel=2,
            )

    def _calculate_criticality(self, reagent_list, rng=None):
        """
        Calculate component criticality using z-score softmax with signal-to-noise
        dampening. Supports IPR (default) and Shannon entropy metrics.

        Converts posterior means to z-scores before softmax, making criticality
        scale-invariant. A signal-to-noise ratio (SNR) check dampens the z-scores
        when between-reagent spread is comparable to sampling noise, preventing
        false criticality signals on balanced libraries.

        SNR = std(posterior_means) / expected_noise_std
        - SNR <= 1: noise dominates -> z-scores suppressed -> criticality -> 0
        - SNR > 1: real signal -> z-scores preserved -> criticality reflects structure

        Criticality in [0, 1]:
        - criticality ~ 0: Flexible component -> should EXPLORE
        - criticality ~ 1: Critical component -> should EXPLOIT

        Args:
            reagent_list: List of Reagent objects with posterior distributions
            rng: Unused, kept for API compatibility

        Returns:
            Criticality score in [0, 1]. Returns 0.5 (neutral) if insufficient data.
        """
        # Filter out retired reagents (those that never got warmup observations)
        active_reagents = [r for r in reagent_list if r.n_samples > 0]
        if len(active_reagents) < 2:
            return 0.5  # Not enough active reagents for entropy

        # Check if we have sufficient data
        observations = [r.n_samples for r in active_reagents]
        if min(observations) < self.min_observations:
            return 0.5  # Neutral criticality if insufficient data

        # Extract posterior means
        means = np.array([r.mean for r in active_reagents])

        # Handle edge case: all means identical
        if np.std(means) < 1e-10:
            return 0.0  # Perfectly flexible (all equally good)

        # For minimization, negate means (want higher probability for lower scores)
        if self.mode == "minimize":
            means = -means

        # Normalize to z-scores so criticality is scale-invariant.
        mean_std = np.std(means)
        if mean_std < 1e-10:
            return 0.0  # All means identical after sign flip

        # Signal-to-noise dampening
        se_squared = np.array([
            r.std ** 2 / max(r.n_samples, 1) for r in active_reagents
        ])
        noise_std = np.sqrt(np.mean(se_squared))
        snr = mean_std / max(noise_std, 1e-10)

        imbalance_strength = 1.0 - np.exp(-max(snr - 1.0, 0.0))

        z_scores = (means - np.mean(means)) / mean_std
        z_scores *= imbalance_strength  # Dampen toward 0 when SNR is low

        # N-adaptive sharpening (counteracts softmax flattening for large N)
        N = len(active_reagents)
        if self.criticality_metric == "ipr" and self.n_adaptive_sharpening and N > 2:
            sharpening = max(1.0, np.sqrt(np.log(N)))
            z_scores *= sharpening

        exp_means = np.exp(z_scores - z_scores.max())  # numerical stability
        probabilities = exp_means / exp_means.sum()

        if self.criticality_metric == "ipr":
            # Inverse Participation Ratio: sensitive to probability concentration
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

    def _calculate_criticality_details(self, reagent_list, rng=None):
        """Calculate criticality and return all intermediate values.

        Returns:
            Tuple of (criticality, details_dict) where details_dict contains
            snr, imbalance_strength, normalized_entropy, n_active_reagents,
            participation_ratio, effective_n, sharpening_factor.
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

        # N-adaptive sharpening
        N = len(active_reagents)
        sharpening = 1.0
        if self.criticality_metric == "ipr" and self.n_adaptive_sharpening and N > 2:
            sharpening = max(1.0, np.sqrt(np.log(N)))
            z_scores *= sharpening
        details["sharpening_factor"] = float(sharpening)

        exp_means = np.exp(z_scores - z_scores.max())  # numerical stability
        probabilities = exp_means / exp_means.sum()

        # Always compute IPR-based metrics for details
        ipr = float(np.sum(probabilities ** 2))
        effective_n = 1.0 / ipr
        details["participation_ratio"] = ipr
        details["effective_n"] = float(effective_n)

        if self.criticality_metric == "ipr":
            criticality = 1.0 - (effective_n / N)
            # Compute normalized entropy for the details dict
            entropy = -np.sum(probabilities * np.log(probabilities + 1e-10))
            max_entropy = np.log(N)
            details["normalized_entropy"] = float(entropy / max_entropy) if max_entropy > 1e-10 else float("nan")
        else:  # shannon (legacy)
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

        Computes all values from the criticality + temperature pipeline and
        caches the result on ``self._last_component_states[component_idx]``.

        Returns:
            Dict with all intermediate values.
        """
        # Step 1: Criticality with details
        criticality, crit_details = self._calculate_criticality_details(reagent_list)

        # Step 2: Base temperature from thermal cycling
        is_heated = component_idx == self.current_component_idx
        base_temp = self.alpha if is_heated else self.beta

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

        # Step 7: Final temperature
        final_temp = base_temp * effective_mult

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
            # Temperature pipeline
            "base_temp": base_temp,
            "is_heated": is_heated,
            "criticality_weight": weight,
            "decay": decay,
            "cats_multiplier": cats_mult,
            "effective_multiplier": effective_mult,
            "final_temperature": final_temp,
        }

        # Cache for sampler to read without recomputing
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
        # Filter out retired reagents (those that never got warmup observations)
        active_obs = [r.n_samples for r in reagent_list if r.n_samples > 0]
        if not active_obs:
            return 0.0

        # Use 25th percentile instead of min to avoid a single never-resampled
        # reagent throttling the entire component's CATS weight.
        p25_obs = float(np.percentile(active_obs, 25))
        # Ramp from 0 → 1 as p25 observations go from 0 → 2*min_observations
        # Below min_observations: low confidence, partial weight
        # Above 2*min_observations: full confidence, weight = 1.0
        weight = p25_obs / (2.0 * self.min_observations) if self.min_observations > 0 else 1.0
        return float(np.clip(weight, 0.0, 1.0))

    def _get_cats_multiplier(self, criticality):
        """
        Map component criticality to a bidirectional temperature multiplier
        using a relative neutral point.

        The neutral point (multiplier = 1.0) is at the mean criticality across
        all components, updated each cycle by ``rotate_component_weighted``.
        Components below the mean get an exploration boost (multiplier > 1);
        components above the mean get an exploitation reduction (multiplier < 1).

        This ensures CATS always differentiates between components regardless
        of the absolute criticality scale. With IPR + N-adaptive sharpening,
        online criticalities cluster in [0.6, 1.0], so a fixed neutral point
        at 0.5 would put all components in the exploit zone, nullifying the
        exploration boost entirely.

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

    def _get_component_temperature(
        self, component_idx, reagent_list, current_cycle, total_cycles, rng=None
    ):
        """
        Get CATS-adjusted temperature for a component.

        Delegates to ``get_component_state()`` which computes and caches all
        intermediate values, then returns just the final temperature.

        Args:
            component_idx: Which reaction component
            reagent_list: List of Reagent objects for this component
            current_cycle: Current search cycle
            total_cycles: Total number of cycles
            rng: Optional numpy random generator (unused, kept for API compat)

        Returns:
            Final temperature value
        """
        state = self.get_component_state(
            reagent_list, component_idx, current_cycle, total_cycles
        )
        return state["final_temperature"]

    def select_reagent(self, reagent_list, disallow_mask=None, **kwargs):
        """
        Select a reagent using roulette wheel selection with CATS.

        Args:
            reagent_list: List of Reagent objects with posterior distributions
            disallow_mask: Optional set of indices to exclude from selection
            **kwargs: Additional context:
                - rng: Random number generator
                - component_idx: Which reaction component
                - current_cycle: Current search cycle (for CATS)
                - total_cycles: Total number of cycles (for CATS)

        Returns:
            Index of selected reagent
        """
        rng = kwargs.get("rng", np.random.default_rng())
        component_idx = kwargs.get("component_idx", 0)
        current_cycle = kwargs.get("current_cycle", 0)
        total_cycles = kwargs.get("total_cycles", 1)

        # Sample base scores
        stds = np.array([r.std for r in reagent_list])
        mu = np.array([r.mean for r in reagent_list])
        scores = rng.normal(size=len(reagent_list)) * stds + mu

        # Invert for minimize mode
        if self.mode not in ["maximize", "maximize_boltzmann"]:
            scores = -scores

        # Get CATS-adjusted temperature
        effective_temp = self._get_component_temperature(
            component_idx, reagent_list, current_cycle, total_cycles, rng=rng
        )

        # Apply temperature via Boltzmann distribution
        # Handle case where all scores are identical (std=0)
        score_std = np.std(scores)
        if score_std < 1e-10:
            # All scores identical, use uniform distribution
            probs = np.ones(len(reagent_list)) / len(reagent_list)
        else:
            scores = np.exp((scores - np.mean(scores)) / score_std / effective_temp)
            # Normalize to probabilities
            probs = scores / np.sum(scores)

        # Apply disallow mask
        if disallow_mask:
            probs[np.array(list(disallow_mask))] = 0
            if np.sum(probs) > 0:
                probs = probs / np.sum(probs)
            else:
                probs = np.ones(len(reagent_list)) / len(reagent_list)

        return rng.choice(len(reagent_list), p=probs)

    def select_batch(self, reagent_list, batch_size, disallow_mask=None, **kwargs):
        """
        Select multiple reagents using roulette wheel selection with CATS (batch mode).

        This is more efficient than calling select_reagent multiple times
        as it computes probabilities once and samples multiple times.

        Args:
            reagent_list: List of Reagent objects with posterior distributions
            batch_size: Number of reagents to select
            disallow_mask: Optional set of indices to exclude from selection
            **kwargs: Additional context:
                - rng: Random number generator
                - component_idx: Which reaction component
                - current_cycle: Current search cycle (for CATS)
                - total_cycles: Total number of cycles (for CATS)

        Returns:
            Array of selected reagent indices
        """
        rng = kwargs.get("rng", np.random.default_rng())
        component_idx = kwargs.get("component_idx", 0)
        current_cycle = kwargs.get("current_cycle", 0)
        total_cycles = kwargs.get("total_cycles", 1)

        # Sample base scores
        stds = np.array([r.std for r in reagent_list])
        mu = np.array([r.mean for r in reagent_list])
        scores = rng.normal(size=len(reagent_list)) * stds + mu

        # Invert scores for minimize mode
        if self.mode not in ["maximize", "maximize_boltzmann"]:
            scores = -scores

        # Get CATS-adjusted temperature
        effective_temp = self._get_component_temperature(
            component_idx, reagent_list, current_cycle, total_cycles, rng=rng
        )

        # Apply temperature via Boltzmann distribution
        # Handle case where all scores are identical (std=0)
        score_std = np.std(scores)
        if score_std < 1e-10:
            # All scores identical, use uniform distribution
            probs = np.ones(len(reagent_list)) / len(reagent_list)
        else:
            scores = np.exp((scores - np.mean(scores)) / score_std / effective_temp)
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
        return rng.choice(len(reagent_list), size=batch_size, p=probs)

    def adapt_temperatures(self, n_unique, n_attempted):
        """
        Adapt temperatures based on sampling efficiency (legacy RWS-inspired).

        When posteriors tighten, selection concentrates on a few reagents, leading
        to more duplicate combinations. Increasing temperatures counteracts this
        by broadening the selection distribution.

        Mirrors the adaptive mechanism from the ETS paper's RWSSampler:
        - alpha += alpha_increment when efficiency < threshold
        - beta += beta_increment when zero unique compounds found

        Args:
            n_unique: Number of unique compounds generated in this batch
            n_attempted: Number of compounds attempted in this batch

        Returns:
            True if temperatures were adjusted, False otherwise
        """
        if not self.adaptive_temperature:
            return False

        adjusted = False
        efficiency = n_unique / max(n_attempted, 1)

        if efficiency < self.efficiency_threshold and self.alpha < self.alpha_max:
            self.alpha += self.alpha_increment
            adjusted = True

        if n_unique == 0:
            self.beta += self.beta_increment
            adjusted = True

        return adjusted

    def rotate_component(self, n_components: int):
        """
        Rotate to the next component for thermal cycling.

        Parameters:
        -----------
        n_components : int
            Total number of reagent components
        """
        self.current_component_idx = (self.current_component_idx + 1) % n_components

    def rotate_component_weighted(self, n_components, criticalities, rng=None):
        """
        Rotate to next component using criticality-weighted probabilities.

        Flexible components (low criticality) get heated more often because
        they benefit more from exploration. Critical components (high criticality)
        already have strong signal and need less heating.

        Also caches the mean criticality so that ``_get_cats_multiplier`` can
        use a relative neutral point instead of a fixed 0.5.

        Parameters:
        -----------
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

    def reset_temperature(self):
        """Reset temperature parameters to initial values."""
        self.alpha = self.initial_alpha
        self.beta = self.initial_beta
