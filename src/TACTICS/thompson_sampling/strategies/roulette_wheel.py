import warnings
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

        # Derive CATS range from alpha/beta ratio
        self.ratio = self.alpha / self.beta if self.beta > 0 else 1.0
        self.cats_max_mult = self.ratio
        self.cats_min_mult = 1.0 / self.ratio if self.ratio > 0 else 1.0

        # Validation warnings
        if abs(self.alpha - self.beta) < 1e-10:
            warnings.warn(
                f"alpha ({self.alpha}) equals beta ({self.beta}). "
                "Thermal cycling will have no effect. "
                "For thermal cycling to work, set beta < alpha (e.g., beta=0.05, alpha=0.1).",
                UserWarning,
                stacklevel=2,
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
        max_entropy = np.log(len(active_reagents))

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
        # Filter out retired reagents (those that never got warmup observations)
        active_obs = [r.n_samples for r in reagent_list if r.n_samples > 0]
        if not active_obs:
            return 0.0

        min_obs = min(active_obs)
        # Ramp from 0 → 1 as observations go from 0 → 2*min_observations
        # Below min_observations: low confidence, partial weight
        # Above 2*min_observations: full confidence, weight = 1.0
        weight = min_obs / (2.0 * self.min_observations) if self.min_observations > 0 else 1.0
        return float(np.clip(weight, 0.0, 1.0))

    def _get_cats_multiplier(self, criticality):
        """
        Map component criticality to temperature multiplier.

        CATS range is derived from alpha/beta ratio:
        - cats_max = alpha / beta (for flexible components, criticality ≈ 0)
        - cats_min = beta / alpha (for critical components, criticality ≈ 1)

        Args:
            criticality: Component criticality in [0, 1]

        Returns:
            CATS multiplier in [cats_min, cats_max]
        """
        # Map criticality [0,1] to multiplier [cats_max, cats_min]
        # Higher criticality → lower multiplier (more exploitation)
        multiplier = self.cats_min_mult + (self.cats_max_mult - self.cats_min_mult) * (
            1 - criticality
        )
        return multiplier

    def _get_component_temperature(
        self, component_idx, reagent_list, current_cycle, total_cycles
    ):
        """
        Get CATS-adjusted temperature for a component.

        Combines thermal cycling with component criticality:
        1. Thermal cycling sets base temperature (alpha for heated, beta for cooled)
        2. CATS adjusts based on criticality
        3. Progressive weighting controls CATS influence
        4. Exploration decay reduces CATS when criticality stays low

        Args:
            component_idx: Which reaction component
            reagent_list: List of Reagent objects for this component
            current_cycle: Current search cycle
            total_cycles: Total number of cycles

        Returns:
            Final temperature value
        """
        # Step 1: Base temperature from thermal cycling
        base_temp = (
            self.alpha if component_idx == self.current_component_idx else self.beta
        )

        # Step 2: Calculate criticality
        criticality = self._calculate_criticality(reagent_list)

        # Step 3: Get observation-gated weight (ramps up with data)
        weight = self._get_criticality_weight(reagent_list)

        # Step 4: Exploration decay — reduce CATS when it hasn't found structure
        # After cats_exploration_fraction of cycles, if criticality is still low,
        # linearly decay the CATS weight so posteriors drive selection directly.
        # High criticality → keep CATS active; low criticality → decay to neutral.
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
        final_temp = base_temp * effective_mult

        return final_temp

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
            component_idx, reagent_list, current_cycle, total_cycles
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
            component_idx, reagent_list, current_cycle, total_cycles
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

    def reset_temperature(self):
        """Reset temperature parameters to initial values."""
        self.alpha = self.initial_alpha
        self.beta = self.initial_beta
