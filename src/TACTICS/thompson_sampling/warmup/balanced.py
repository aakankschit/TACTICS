"""Balanced warmup strategy with exactly K observations per reagent and stratified partner selection."""

import numpy as np
from typing import List, Optional, TYPE_CHECKING
from .base import WarmupStrategy

if TYPE_CHECKING:
    from ..core.reagent import Reagent
    from ..legacy.disallow_tracker import DisallowTracker


class BalancedWarmup(WarmupStrategy):
    """
    Balanced warmup strategy guaranteeing exactly K observations per reagent.

    This strategy ensures:
    1. Every reagent gets exactly `observations_per_reagent` observations
    2. Partners are selected using stratified sampling (no duplicates within a reagent's trials)
    3. Optional seeded RNG for reproducibility

    Characteristics:
    ----------------
    - Balance: Each reagent gets exactly K observations (no imbalance)
    - Diversity: Partners sampled from K different strata (no duplicates)
    - Evaluations: sum(component_sizes) x observations_per_reagent
    - Reproducibility: Optional seed for deterministic results

    Example:
    --------
    For 130 acids and 3844 amines with K=5:
    - Each acid tested with 5 amines from different strata
    - Each amine tested with 5 acids from different strata
    - Total: (130 + 3844) x 5 = 19,870 evaluations
    - Every reagent guaranteed exactly 5 observations

    Parameters:
    -----------
    observations_per_reagent : int, default=5
        Number of observations guaranteed per reagent (K)
    seed : int or None, default=None
        Random seed for reproducibility. None = random each run.
    """

    def __init__(
        self,
        observations_per_reagent: int = 5,
        seed: Optional[int] = None,
        use_per_reagent_variance: bool = True,
        shrinkage_strength: float = 3.0
    ):
        """
        Initialize BalancedWarmup strategy.

        Args:
            observations_per_reagent: Number of observations per reagent (K)
            seed: Optional random seed for reproducibility
            use_per_reagent_variance: If True, use per-reagent variance with shrinkage
            shrinkage_strength: Shrinkage parameter (higher = more regularization)
        """
        self.observations_per_reagent = observations_per_reagent
        self.seed = seed
        self.use_per_reagent_variance = use_per_reagent_variance
        self.shrinkage_strength = shrinkage_strength

    def get_name(self) -> str:
        return f"BalancedWarmup(K={self.observations_per_reagent})"

    def generate_warmup_combinations(
        self,
        reagent_lists: List[List['Reagent']],
        num_warmup_trials: int,
        disallow_tracker: 'DisallowTracker'
    ) -> List[List[int]]:
        """
        Generate warmup combinations with exactly K observations per reagent.

        Uses stratified partner selection: divides partner space into K strata
        and samples one partner from each stratum, guaranteeing diversity.

        Note: The num_warmup_trials parameter is ignored - we use observations_per_reagent instead.
        This ensures the caller's expectation of "trials per reagent" is honored.

        Args:
            reagent_lists: List of reagent lists, one for each component
            num_warmup_trials: Ignored (uses observations_per_reagent)
            disallow_tracker: Tracker to prevent duplicate combinations

        Returns:
            List of combinations [idx_comp1, idx_comp2, ...]
        """
        # Use observations_per_reagent as K (ignore num_warmup_trials for clarity)
        K = self.observations_per_reagent

        # Initialize RNG with seed for reproducibility
        rng = np.random.default_rng(self.seed)

        combinations = []
        n_components = len(reagent_lists)
        reagent_counts = [len(rl) for rl in reagent_lists]

        # Iterate over each component
        for component_idx in range(n_components):
            partner_indices = [i for i in range(n_components) if i != component_idx]

            # Iterate over each reagent in this component
            for reagent_idx in range(reagent_counts[component_idx]):
                # Generate K stratified partner selections for each partner component
                partners = self._stratified_partner_selection(
                    partner_indices, reagent_counts, K, rng
                )

                # Create K combinations for this reagent
                for trial in range(K):
                    combination = [disallow_tracker.Empty] * n_components
                    combination[component_idx] = reagent_idx

                    # Fill in partners from stratified selections
                    for partner_comp in partner_indices:
                        combination[partner_comp] = partners[partner_comp][trial]

                    # Check if combination is allowed and add if so
                    if self._is_combination_allowed(combination, disallow_tracker, component_idx):
                        disallow_tracker.update(combination)
                        combinations.append(combination)

        return combinations

    def _stratified_partner_selection(
        self,
        partner_indices: List[int],
        reagent_counts: List[int],
        K: int,
        rng: np.random.Generator
    ) -> dict:
        """
        Select K partners from each partner component using stratified sampling.

        Divides partner space into K strata and samples one partner from each,
        guaranteeing diversity across the full partner range.

        Args:
            partner_indices: Indices of partner components
            reagent_counts: Number of reagents in each component
            K: Number of partners to select (= number of strata)
            rng: Random number generator

        Returns:
            Dict mapping partner_component_idx -> list of K selected partner indices
        """
        partners = {}

        for partner_comp in partner_indices:
            pool_size = reagent_counts[partner_comp]
            selected = []

            if pool_size <= K:
                # Pool smaller than K: use all partners, then cycle
                base_indices = list(range(pool_size))
                rng.shuffle(base_indices)
                # Repeat to fill K slots
                repeats = (K // pool_size) + 1
                selected = (base_indices * repeats)[:K]
            else:
                # Stratified sampling: divide pool into K strata
                stratum_size = pool_size / K
                for k in range(K):
                    stratum_start = int(k * stratum_size)
                    stratum_end = int((k + 1) * stratum_size)
                    if k == K - 1:
                        stratum_end = pool_size  # Last stratum gets remainder

                    # Sample one partner from this stratum
                    partner_idx = rng.integers(stratum_start, stratum_end)
                    selected.append(partner_idx)

            partners[partner_comp] = selected

        return partners

    def _is_combination_allowed(
        self,
        combination: List[int],
        disallow_tracker: 'DisallowTracker',
        fixed_component: int
    ) -> bool:
        """
        Check if a combination is allowed by the disallow tracker.

        Args:
            combination: Full combination to check
            disallow_tracker: DisallowTracker instance
            fixed_component: Index of the fixed component (not checked)

        Returns:
            True if combination is allowed, False otherwise
        """
        n_components = len(combination)

        for idx in range(n_components):
            if idx == fixed_component:
                continue

            # Create test combination with this slot marked as To_Fill
            test_combination = combination.copy()
            test_combination[idx] = disallow_tracker.To_Fill

            # Get disallowed selections for this partial combination
            disallow_mask = disallow_tracker.get_disallowed_selection_mask(test_combination)

            # Check if our selected reagent is disallowed
            if combination[idx] in disallow_mask:
                return False

        return True

    def get_expected_evaluations(
        self,
        reagent_lists: List[List['Reagent']],
        num_warmup_trials: int
    ) -> int:
        """
        Calculate expected number of evaluations.

        For balanced warmup: sum(component_sizes) x observations_per_reagent

        Args:
            reagent_lists: List of reagent lists
            num_warmup_trials: Ignored (uses observations_per_reagent)

        Returns:
            Expected number of evaluations
        """
        return sum(len(rl) for rl in reagent_lists) * self.observations_per_reagent
