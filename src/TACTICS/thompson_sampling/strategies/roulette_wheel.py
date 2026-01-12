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
            **kwargs: Catches deprecated parameters with warnings
        """
        super().__init__(mode)

        # Core parameters
        self.initial_alpha = alpha
        self.initial_beta = beta
        self.alpha = alpha
        self.beta = beta

        # CATS parameters
        self.exploration_phase_end = exploration_phase_end
        self.transition_phase_end = transition_phase_end
        self.min_observations = min_observations

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

        # Deprecated parameter warnings
        deprecated = {
            "alpha_max",
            "alpha_increment",
            "beta_increment",
            "efficiency_threshold",
            "scaling",
        }
        found_deprecated = deprecated & set(kwargs.keys())
        if found_deprecated:
            warnings.warn(
                f"Parameters {found_deprecated} are deprecated and will be ignored. "
                f"CATS now derives all parameters from alpha/beta ratio. "
                f"Remove these parameters from your configuration.",
                DeprecationWarning,
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

    def _get_criticality_weight(self, current_cycle, total_cycles):
        """
        Calculate progressive criticality weight for current iteration.

        Three-phase progression prevents early-iteration bias:
        - Phase 1 (0-20%): Pure exploration, weight = 0.0
        - Phase 2 (20-60%): Gradual introduction, weight increases linearly 0→1
        - Phase 3 (60-100%): Full CATS, weight = 1.0

        Args:
            current_cycle: Current search cycle (0-indexed)
            total_cycles: Total number of search cycles

        Returns:
            Criticality weight in [0, 1]
        """
        if total_cycles <= 0:
            return 1.0  # Edge case

        # Calculate progress as fraction of total
        progress = current_cycle / total_cycles

        # Phase 1: Pure exploration (0-20%)
        if progress < self.exploration_phase_end:
            return 0.0

        # Phase 2: Gradual introduction (20-60%)
        elif progress < self.transition_phase_end:
            phase_progress = (progress - self.exploration_phase_end) / (
                self.transition_phase_end - self.exploration_phase_end
            )
            return phase_progress

        # Phase 3: Full CATS (60-100%)
        else:
            return 1.0

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

        # Step 3: Get progressive weight
        weight = self._get_criticality_weight(current_cycle, total_cycles)

        # Step 4: Get CATS multiplier
        cats_mult = self._get_cats_multiplier(criticality)

        # Step 5: Progressive blending
        effective_mult = (1.0 - weight) * 1.0 + weight * cats_mult

        # Step 6: Apply to base
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

        return np.random.choice(len(reagent_list), p=probs)

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
        return np.random.choice(len(reagent_list), size=batch_size, p=probs)

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
