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
        divergence_threshold=0.1,
        cats_ema_decay=None,
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
            divergence_threshold: KL divergence threshold for switching from diversity
                to GMIC criticality mode (default: 0.1)
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

        # GMIC + divergence state
        self._divergence_threshold = divergence_threshold
        # Per-component posterior snapshots: {comp_idx: (means_array, stds_array)}
        self._posterior_snapshots: Dict[int, tuple] = {}
        # Per-component divergence values: {comp_idx: float}
        self._divergence_values: Dict[int, float] = {}
        # Cached GMIC per component (updated by rotate_component_weighted)
        self._cached_gmics: Dict[int, float] = {}
        # CATS EMA smoothing for relative GMIC
        self.cats_ema_decay = cats_ema_decay
        self._ema_relative_gmic: Dict[int, float] = {}

        # Validation warnings
        if abs(self.alpha - self.beta) < 1e-10:
            warnings.warn(
                f"alpha ({self.alpha}) equals beta ({self.beta}). "
                "Thermal cycling will have no effect. "
                "For thermal cycling to work, set beta < alpha (e.g., beta=0.05, alpha=0.1).",
                UserWarning,
                stacklevel=2,
            )

    def _calculate_gmic(self, reagent_list) -> float:
        """
        Calculate Gaussian Mutual Information Criticality (GMIC) for a component.

        GMIC measures how much information the posterior means carry about which
        reagent is best, relative to the noise level. High GMIC = critical component
        (clear winners), low GMIC = flexible component (all reagents similar).

        Returns:
            GMIC value >= 0. Returns 0.0 if insufficient data.
        """
        active = [r for r in reagent_list if r.n_samples > 0]
        if len(active) < 2:
            return 0.0
        means = np.array([r.mean for r in active])
        signal = np.var(means)
        noise = np.mean([r.std ** 2 for r in active])
        return 0.5 * np.log1p(signal / max(noise, 1e-10))

    def _calculate_gmic_details(self, reagent_list) -> tuple:
        """
        Calculate GMIC with full diagnostic details.

        Returns:
            Tuple of (gmic, details_dict) where details_dict contains
            signal_var, mean_noise_var, n_active_reagents.
        """
        active = [r for r in reagent_list if r.n_samples > 0]
        details = {
            "signal_var": float("nan"),
            "mean_noise_var": float("nan"),
            "n_active_reagents": len(active),
        }
        if len(active) < 2:
            return 0.0, details
        means = np.array([r.mean for r in active])
        signal = float(np.var(means))
        noise = float(np.mean([r.std ** 2 for r in active]))
        details["signal_var"] = signal
        details["mean_noise_var"] = noise
        gmic = 0.5 * np.log1p(signal / max(noise, 1e-10))
        return float(gmic), details

    def _update_divergence(self, component_idx: int, reagent_list) -> float:
        """
        Update posterior divergence tracker for a component.

        Computes the average KL divergence between current and previous posterior
        snapshots across all active reagents. Uses Gaussian KL:
            D_KL(N(mu1,s1^2) || N(mu2,s2^2)) = log(s2/s1) + (s1^2 + (mu1-mu2)^2)/(2*s2^2) - 0.5

        Args:
            component_idx: Index of the component
            reagent_list: List of Reagent objects

        Returns:
            Average KL divergence (inf if no prior snapshot).
        """
        active = [r for r in reagent_list if r.n_samples > 0]
        if len(active) < 2:
            return float("inf")

        current_means = np.array([r.mean for r in active])
        current_stds = np.array([np.maximum(r.std, 1e-10) for r in active])

        if component_idx not in self._posterior_snapshots:
            # First call: store snapshot, return inf (no comparison possible)
            self._posterior_snapshots[component_idx] = (
                current_means.copy(),
                current_stds.copy(),
            )
            self._divergence_values[component_idx] = float("inf")
            return float("inf")

        prev_means, prev_stds = self._posterior_snapshots[component_idx]

        # Handle size changes (e.g., new active reagents)
        n = min(len(current_means), len(prev_means))
        if n < 2:
            self._posterior_snapshots[component_idx] = (
                current_means.copy(),
                current_stds.copy(),
            )
            self._divergence_values[component_idx] = float("inf")
            return float("inf")

        mu1 = current_means[:n]
        s1 = current_stds[:n]
        mu2 = prev_means[:n]
        s2 = prev_stds[:n]

        # Gaussian KL: D_KL(current || previous)
        kl_per_reagent = (
            np.log(s2 / s1)
            + (s1 ** 2 + (mu1 - mu2) ** 2) / (2.0 * s2 ** 2)
            - 0.5
        )
        avg_kl = float(np.mean(kl_per_reagent))

        # Update snapshot
        self._posterior_snapshots[component_idx] = (
            current_means.copy(),
            current_stds.copy(),
        )
        self._divergence_values[component_idx] = avg_kl
        return avg_kl

    def _get_component_state_gmic(
        self,
        reagent_list: List,
        component_idx: int,
        current_cycle: int,
        total_cycles: int,
    ) -> Dict[str, Any]:
        """
        Compute full intermediate state for a component using GMIC + Divergence gate.

        1. Compute GMIC + details
        2. Update divergence tracker
        3. Check stability: is_stable = divergence < threshold
        4. Base temp from thermal cycling (alpha/beta)
        5. If stable: apply GMIC-based directional multiplier
        6. If not stable: use base temp (diversity fallback)

        Returns:
            Dict with all intermediate values.
        """
        # Step 1: GMIC with details
        gmic, gmic_details = self._calculate_gmic_details(reagent_list)

        # Step 2: Update divergence
        divergence = self._update_divergence(component_idx, reagent_list)

        # Step 3: Stability check
        is_stable = divergence < self._divergence_threshold

        # Step 4: Base temperature from thermal cycling
        is_heated = component_idx == self.current_component_idx
        base_temp = self.alpha if is_heated else self.beta

        # Step 5/6: GMIC-based multiplier or diversity fallback
        if is_stable and gmic > 0:
            # Use cached mean GMIC for relative comparison
            mean_gmic = np.mean(list(self._cached_gmics.values())) if self._cached_gmics else gmic
            mean_gmic = max(mean_gmic, 1e-10)

            # Relative GMIC: >1 = more critical than average, <1 = more flexible
            relative_gmic = gmic / mean_gmic

            # EMA smoothing: accumulate directional signal across cycles
            if self.cats_ema_decay is not None:
                if component_idx in self._ema_relative_gmic:
                    self._ema_relative_gmic[component_idx] = (
                        self.cats_ema_decay * relative_gmic
                        + (1.0 - self.cats_ema_decay) * self._ema_relative_gmic[component_idx]
                    )
                else:
                    self._ema_relative_gmic[component_idx] = relative_gmic
                relative_gmic = self._ema_relative_gmic[component_idx]

            if is_heated:
                # Flexible components (low relative GMIC) get amplified heating
                if relative_gmic < 1.0:
                    t = 1.0 - relative_gmic  # 0 to ~1
                    cats_mult = 1.0 + t * (self.cats_max_mult - 1.0)
                else:
                    cats_mult = 1.0
            else:
                # Critical components (high relative GMIC) get amplified cooling
                if relative_gmic > 1.0:
                    t = min(relative_gmic - 1.0, 1.0)  # 0 to 1, capped
                    cats_mult = 1.0 + t * (self.cats_min_mult - 1.0)
                else:
                    cats_mult = 1.0

            cats_mode = "criticality"
            effective_mult = cats_mult
        else:
            # Diversity fallback: base temperature only
            cats_mult = 1.0
            effective_mult = 1.0
            cats_mode = "diversity"

        # Step 7: Final temperature
        final_temp = base_temp * effective_mult

        state = {
            # Metadata
            "component_idx": component_idx,
            "current_cycle": current_cycle,
            "total_cycles": total_cycles,
            # GMIC details
            "gmic": gmic,
            "criticality": gmic,
            "signal_var": gmic_details["signal_var"],
            "mean_noise_var": gmic_details["mean_noise_var"],
            "n_active_reagents": gmic_details["n_active_reagents"],
            # Divergence gate
            "divergence": divergence,
            "divergence_threshold": self._divergence_threshold,
            "is_stable": is_stable,
            "cats_mode": cats_mode,
            # Temperature pipeline
            "base_temp": base_temp,
            "is_heated": is_heated,
            "cats_multiplier": cats_mult,
            "effective_multiplier": effective_mult,
            "final_temperature": final_temp,
            # EMA smoothing state
            "ema_relative_gmic": self._ema_relative_gmic.get(component_idx, float("nan")),
        }

        self._last_component_states[component_idx] = state
        return state

    def get_component_criticality(self, reagent_list) -> float:
        """Return GMIC criticality score for a component."""
        return self._calculate_gmic(reagent_list)

    def get_component_state(
        self,
        reagent_list: List,
        component_idx: int,
        current_cycle: int,
        total_cycles: int,
    ) -> Optional[Dict[str, Any]]:
        """Return full intermediate state for a component.

        Computes GMIC criticality with divergence-gated temperature adjustment.
        Caches the result on ``self._last_component_states[component_idx]``.

        Returns:
            Dict with all intermediate values.
        """
        return self._get_component_state_gmic(
            reagent_list, component_idx, current_cycle, total_cycles
        )

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

    def rotate_component_weighted(self, n_components, reagent_lists, rng=None):
        """
        Rotate to next component using GMIC-weighted probabilities.

        Flexible components (low GMIC) get heated more often. Weight = 1 / (1 + gmic).
        Also caches per-component GMIC values for use by get_component_state().

        Parameters:
        -----------
        n_components : int
            Total number of reagent components
        reagent_lists : list of list
            List of reagent lists, one per component
        rng : numpy.random.Generator, optional
            Random number generator for reproducibility
        """
        if rng is None:
            rng = np.random.default_rng()

        gmics = []
        for i in range(n_components):
            g = self._calculate_gmic(reagent_lists[i])
            self._cached_gmics[i] = g
            gmics.append(g)

        gmics_arr = np.array(gmics)
        # Weight by flexibility: low GMIC → more heating
        flexibility = 1.0 / (1.0 + gmics_arr)
        heat_probs = flexibility / flexibility.sum()

        self.current_component_idx = int(rng.choice(n_components, p=heat_probs))

        # Cache mean GMIC for relative comparison
        self._mean_criticality = float(gmics_arr.mean())

    def reset_temperature(self):
        """Reset temperature parameters to initial values."""
        self.alpha = self.initial_alpha
        self.beta = self.initial_beta
