"""Standard warmup strategy with random partner selection (with replacement)."""

import numpy as np
from typing import List, TYPE_CHECKING
from .base import WarmupStrategy

if TYPE_CHECKING:
    from ..core.reagent import Reagent
    from ..legacy.disallow_tracker import DisallowTracker


class StandardWarmup(WarmupStrategy):
    """
    Standard warmup strategy using random partner selection with replacement.

    This is the baseline/legacy warmup approach where each reagent is tested
    num_warmup_trials times with randomly selected partners. Partners are chosen
    independently for each trial, so the same partner may be selected multiple times.

    Characteristics:
    ----------------
    - Balance: ✅ Each reagent gets exactly num_warmup_trials samples
    - Diversity: ❌ ~30% chance of duplicate partners for small partner pools
    - Evaluations: sum(component_sizes) × num_warmup_trials
    - Use case: Backward compatibility, simple baseline

    Example:
    --------
    For 130 acids and 3844 amines with 10 trials:
    - Each acid tested with 10 random amines (may repeat)
    - Each amine tested with 10 random acids (may repeat)
    - Total: 39,740 evaluations
    """

    def get_name(self) -> str:
        return "StandardWarmup"

    def generate_warmup_combinations(
        self,
        reagent_lists: List[List['Reagent']],
        num_warmup_trials: int,
        disallow_tracker: 'DisallowTracker'
    ) -> List[List[int]]:
        """
        Generate warmup combinations using random partner selection.

        For each reagent in each component:
        1. Fix the reagent
        2. For each trial, randomly select partners from other components
        3. Add combination to list if not already evaluated

        Note: Partner selection is independent across trials, so duplicates may occur.
        """
        combinations = []
        n_components = len(reagent_lists)
        reagent_counts = [len(rl) for rl in reagent_lists]

        for component_idx in range(n_components):
            partner_indices = [i for i in range(n_components) if i != component_idx]

            for reagent_idx in range(reagent_counts[component_idx]):
                # Test this reagent num_warmup_trials times
                for trial in range(num_warmup_trials):
                    # Initialize combination: fixed reagent + Empty for partners
                    combination = [disallow_tracker.Empty] * n_components
                    combination[component_idx] = reagent_idx

                    # Select random partners for this trial
                    # Fill partners in randomized order to satisfy DisallowTracker constraint
                    partner_order = np.random.permutation(partner_indices).tolist()

                    for partner_component_idx in partner_order:
                        partner_pool_size = reagent_counts[partner_component_idx]

                        # Mark this position as "to fill" while keeping others as Empty
                        # This ensures exactly ONE To_Fill slot for DisallowTracker
                        combination[partner_component_idx] = disallow_tracker.To_Fill
                        disallow_mask = disallow_tracker.get_disallowed_selection_mask(combination)

                        # Generate random scores and select highest (not in disallow mask)
                        selection_scores = np.random.uniform(size=partner_pool_size)
                        selection_scores[list(disallow_mask)] = np.nan

                        # Select partner with highest score
                        combination[partner_component_idx] = np.nanargmax(selection_scores).item()

                    # Mark this combination as used
                    disallow_tracker.update(combination)
                    combinations.append(combination)

        return combinations
