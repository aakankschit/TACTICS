"""Enhanced warmup strategy with parallel pairing (legacy implementation)."""

import random
import numpy as np
from typing import List, TYPE_CHECKING
from .base import WarmupStrategy

if TYPE_CHECKING:
    from ..core.reagent import Reagent
    from ..legacy.disallow_tracker import DisallowTracker


class EnhancedWarmup(WarmupStrategy):
    """
    Enhanced warmup strategy using stochastic parallel pairing.

    This is the legacy warmup approach from EnhancedThompsonSampler. It shuffles
    all reagents and pairs them exhaustively in each trial, repeating for
    num_warmup_trials times.

    Key Characteristic: IMBALANCED sampling
    ----------------------------------------
    Small components get over-sampled relative to large components because they
    are repeated to match the size of the largest component.

    Example: 130 acids × 3844 amines
    - Each acid appears ~30 times per trial (repeated to match 3844)
    - Each amine appears 1 time per trial
    - Result: acids get 300 samples, amines get 10 samples (30x imbalance!)

    Characteristics:
    ----------------
    - Balance: ❌ Small component gets N×(rmax/rmin) samples, large gets N samples
    - Diversity: ✅ Guaranteed (shuffled pairing each trial)
    - Coverage: Excellent for small component (~7.8% of all partners)
    - Evaluations: max(component_sizes) × num_warmup_trials
    - Use case: When you WANT comprehensive small-component coverage

    Example:
    --------
    For 130 acids and 3844 amines with 10 trials:
    - Generates 3844 pairs per trial
    - Acids are repeated: [a1,a2,...a130] × 30 = 3844 slots
    - Amines used once: [m1,m2,...m3844]
    - Shuffle both, pair up: [(a_i, m_j) for all j in 1..3844]
    - Repeat 10 times with different shuffles
    - Total: 38,440 evaluations
    - Each acid tested with ~300 different amines
    - Each amine tested with 10 different acids

    Trade-off:
    ----------
    Provides comprehensive coverage of small component but creates imbalanced
    posteriors. Use this if you specifically want to thoroughly explore the
    small component's interactions with the large component.
    """

    def get_name(self) -> str:
        return "EnhancedWarmup"

    def get_expected_evaluations(
        self,
        reagent_lists: List[List['Reagent']],
        num_warmup_trials: int
    ) -> int:
        """
        Calculate expected evaluations for parallel pairing.

        Returns: max(component_sizes) × num_warmup_trials
        """
        reagent_counts = [len(rl) for rl in reagent_lists]
        rmax = max(reagent_counts)
        return rmax * num_warmup_trials

    def generate_warmup_combinations(
        self,
        reagent_lists: List[List['Reagent']],
        num_warmup_trials: int,
        disallow_tracker: 'DisallowTracker'
    ) -> List[List[int]]:
        """
        Generate warmup combinations using stochastic parallel pairing.

        Algorithm:
        ----------
        1. For each trial:
           a. Shuffle all reagent indices in each component
           b. Repeat small components to match largest component size
           c. Transpose to create pairs
           d. Add all pairs to combinations list

        This creates balanced trials (all reagents paired once per trial)
        but imbalanced sampling (small components over-sampled).
        """
        combinations = []
        n_components = len(reagent_lists)
        reagent_counts = [len(rl) for rl in reagent_lists]
        rmax = max(reagent_counts)
        idx_max = reagent_counts.index(rmax)

        # Calculate how many times each component needs to be repeated
        # to match the largest component
        packing = []
        for i, count in enumerate(reagent_counts):
            if count < rmax:
                repeats = rmax // count
                remainder = rmax % count
                packing.append((1, repeats, remainder))  # needs_repeat, full_repeats, extra
            else:
                packing.append((0,))  # no repeat needed

        # Generate trials
        for trial in range(num_warmup_trials):
            # Create matrix of reagent indices for this trial
            matrix = []

            for component_idx, n_reagents in enumerate(reagent_counts):
                # Get all reagent indices for this component
                indices = list(range(n_reagents))
                # Shuffle to randomize pairings
                random.shuffle(indices)

                # Repeat if necessary to match rmax
                if packing[component_idx][0]:  # needs repeat
                    repeats = packing[component_idx][1]
                    remainder = packing[component_idx][2]
                    # Repeat full list 'repeats' times, then add 'remainder' more
                    repeated_indices = indices * repeats + indices[:remainder]
                    matrix.append(repeated_indices)
                else:
                    matrix.append(indices)

            # Transpose matrix to get pairs
            # matrix[0] = [acid_indices...] (length rmax)
            # matrix[1] = [amine_indices...] (length rmax)
            # zip them to get [(acid_i, amine_i), ...]
            pairs = np.array(matrix).T

            # Add all pairs from this trial
            for pair in pairs:
                combination = list(pair)
                disallow_tracker.update(combination)
                combinations.append(combination)

        return combinations

    def get_description(self) -> str:
        return """
Enhanced warmup with stochastic parallel pairing.

This strategy provides comprehensive coverage of the small component by
pairing all reagents in a shuffled manner each trial. Small components
are repeated to match the largest component size.

WARNING: Creates imbalanced posteriors!
- Small components get many more samples than large components
- Use only if you specifically want to thoroughly explore small component
- For balanced posteriors, use StandardWarmup or StratifiedWarmup instead

Example imbalance for 130 acids × 3844 amines:
- Acids: ~300 samples per reagent
- Amines: ~10 samples per reagent
- 30x difference!
        """
