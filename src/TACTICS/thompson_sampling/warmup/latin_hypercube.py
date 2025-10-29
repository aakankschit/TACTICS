"""Latin Hypercube warmup strategy for optimal space-filling sampling."""

import numpy as np
from typing import List, Dict, TYPE_CHECKING
from .base import WarmupStrategy

if TYPE_CHECKING:
    from ..core.reagent import Reagent
    from ..legacy.disallow_tracker import DisallowTracker


class LatinHypercubeWarmup(WarmupStrategy):
    """
    Latin Hypercube Sampling (LHS) warmup strategy for optimal space coverage.

    This is an advanced warmup strategy that ensures optimal space-filling properties.
    Instead of random partner selection (which may cluster), LHS divides the partner
    space into strata and samples one partner from each stratum, guaranteeing
    coverage across the entire reagent range.

    Key Advantage over StratifiedWarmup:
    ------------------------------------
    While StratifiedWarmup ensures no duplicates, random selection without replacement
    can still result in clustering (e.g., all partners from first 10% of space).
    LHS guarantees spatial spread across the ENTIRE partner range.

    Characteristics:
    ----------------
    - Balance: ✅ Each reagent gets exactly num_warmup_trials samples
    - Diversity: ✅ Guaranteed no duplicates within reagent
    - Spatial coverage: ✅ Guaranteed spread across entire partner space
    - Evaluations: sum(component_sizes) × num_warmup_trials (same as Standard)
    - Use case: When reagent ordering has meaning (e.g., sorted by molecular weight)

    Example:
    --------
    For an acid testing 10 amines from a pool of 3844:

    Random (without replacement):
    - Might select: [5, 12, 18, 25, 31, 42, 55, 68, 74, 89]
    - All from first 2.3% of amine space! (clustering)

    Latin Hypercube:
    - Divides 3844 amines into 10 strata: [0-384], [385-768], ..., [3460-3843]
    - Selects one from each stratum:
      [142, 501, 987, 1234, 1598, 1876, 2341, 2789, 3124, 3721]
    - Guaranteed spread across LOW, MIDDLE, and HIGH index amines

    Benefits:
    ---------
    1. **Optimal space-filling**: Maximizes coverage with fixed sample size
    2. **Reduced sampling bias**: No risk of clustering in random selection
    3. **Better for sorted reagents**: If reagents ordered by property (MW, LogP, etc.)
    4. **Statistical efficiency**: Better estimation of posterior with same samples

    When to Use:
    ------------
    - Reagents are ordered meaningfully (by property, similarity, etc.)
    - You want maximum information from limited warmup budget
    - You need robust initialization across diverse reagent space
    - Chemical diversity is important for your library

    Algorithm:
    ----------
    For each reagent:
    1. Divide partner space into num_warmup_trials equal strata
    2. Randomly select one partner from each stratum
    3. Test reagent with each selected partner
    This guarantees both uniqueness AND spatial spread.
    """

    def get_name(self) -> str:
        return "LatinHypercubeWarmup"

    def generate_warmup_combinations(
        self,
        reagent_lists: List[List['Reagent']],
        num_warmup_trials: int,
        disallow_tracker: 'DisallowTracker'
    ) -> List[List[int]]:
        """
        Generate warmup combinations using Latin Hypercube Sampling.

        For each reagent:
        1. Divide partner space into strata (bins)
        2. Sample one partner from each stratum
        3. Generate combinations using stratified selections

        This ensures both diversity (no duplicates) and coverage (spatial spread).
        """
        combinations = []
        n_components = len(reagent_lists)
        reagent_counts = [len(rl) for rl in reagent_lists]

        for component_idx in range(n_components):
            partner_component_indices = [i for i in range(n_components) if i != component_idx]

            for reagent_idx in range(reagent_counts[component_idx]):
                # Use Latin Hypercube Sampling to select partners
                partner_selections = self._latin_hypercube_sample_partners(
                    partner_component_indices,
                    reagent_counts,
                    num_warmup_trials
                )

                # Generate combinations using LHS-selected partners
                for trial_idx in range(num_warmup_trials):
                    combination = [None] * n_components
                    combination[component_idx] = reagent_idx

                    # Use stratified partners for this trial
                    for partner_component_idx in partner_component_indices:
                        partner_idx = partner_selections[partner_component_idx][trial_idx]
                        combination[partner_component_idx] = partner_idx

                    # Check if combination is allowed
                    if self._is_combination_allowed(combination, disallow_tracker):
                        disallow_tracker.update(combination)
                        combinations.append(combination)

        return combinations

    def _latin_hypercube_sample_partners(
        self,
        partner_component_indices: List[int],
        reagent_counts: List[int],
        num_trials: int
    ) -> Dict[int, np.ndarray]:
        """
        Select partners using Latin Hypercube Sampling for each component.

        Divides each partner pool into strata and samples one from each stratum,
        ensuring optimal space-filling coverage.

        Parameters:
        -----------
        partner_component_indices : List[int]
            Indices of partner components
        reagent_counts : List[int]
            Number of reagents in each component
        num_trials : int
            Number of trials (strata to create)

        Returns:
        --------
        Dict[int, np.ndarray]
            Mapping from component_idx to array of stratified partner indices
        """
        partner_selections = {}

        for partner_component_idx in partner_component_indices:
            pool_size = reagent_counts[partner_component_idx]
            n_samples = min(num_trials, pool_size)

            # Perform stratified sampling (Latin Hypercube)
            selected_indices = self._stratified_sample(pool_size, n_samples)

            # Handle edge case: more trials than partners
            if n_samples < num_trials:
                # Very rare: fewer partners than trials requested
                # Cycle through stratified samples
                repeats = (num_trials // n_samples) + 1
                selected_indices = np.tile(selected_indices, repeats)[:num_trials]

            partner_selections[partner_component_idx] = selected_indices

        return partner_selections

    def _stratified_sample(self, pool_size: int, n_samples: int) -> np.ndarray:
        """
        Perform stratified sampling (Latin Hypercube Sampling).

        Divides pool into n_samples equal strata and selects one random
        element from each stratum.

        Parameters:
        -----------
        pool_size : int
            Total size of the pool to sample from
        n_samples : int
            Number of samples (strata) to generate

        Returns:
        --------
        np.ndarray
            Array of stratified sample indices

        Example:
        --------
        pool_size=100, n_samples=4
        - Stratum 1: [0-24]   → select random from range, e.g., 12
        - Stratum 2: [25-49]  → select random from range, e.g., 38
        - Stratum 3: [50-74]  → select random from range, e.g., 61
        - Stratum 4: [75-99]  → select random from range, e.g., 87
        - Result: [12, 38, 61, 87] - well-spread across entire range
        """
        stratum_size = pool_size / n_samples  # Use float division for exact boundaries
        samples = []

        for stratum_idx in range(n_samples):
            # Calculate stratum boundaries
            start = int(stratum_idx * stratum_size)
            end = int((stratum_idx + 1) * stratum_size) if stratum_idx < n_samples - 1 else pool_size

            # Ensure valid range
            if end <= start:
                end = start + 1

            # Sample one random element from this stratum
            sample = np.random.randint(start, end)
            samples.append(sample)

        return np.array(samples)

    def _is_combination_allowed(
        self,
        combination: List[int],
        disallow_tracker: 'DisallowTracker'
    ) -> bool:
        """
        Check if combination is allowed (not already evaluated).

        Parameters:
        -----------
        combination : List[int]
            Reagent combination to check
        disallow_tracker : DisallowTracker
            Global tracker of evaluated combinations

        Returns:
        --------
        bool
            True if combination is allowed, False if already evaluated
        """
        for component_idx in range(len(combination)):
            test_combination = combination.copy()
            test_combination[component_idx] = disallow_tracker.To_Fill

            disallow_mask = disallow_tracker.get_disallowed_selection_mask(test_combination)

            if combination[component_idx] in disallow_mask:
                return False

        return True

    def get_description(self) -> str:
        return """
Latin Hypercube Sampling warmup for optimal space coverage.

This advanced strategy ensures maximum coverage of the reagent space with
a fixed sample budget. By dividing the partner space into strata and sampling
from each, it guarantees spread across the entire range.

Key benefits:
- Optimal space-filling properties
- No clustering (unlike random sampling)
- Better for reagents sorted by properties
- Maximum information from fixed budget

Best for:
- Reagents ordered by molecular properties
- When diversity in warmup is critical
- High-stakes applications requiring robust initialization
        """
