"""Base class for warmup strategies."""

from abc import ABC, abstractmethod
from typing import List, TYPE_CHECKING

if TYPE_CHECKING:
    from ..core.reagent import Reagent
    from ..legacy.disallow_tracker import DisallowTracker


class WarmupStrategy(ABC):
    """
    Abstract base class for warmup strategies.

    A warmup strategy determines how reagent combinations are generated
    during the warmup phase to initialize reagent posteriors.
    """

    @abstractmethod
    def generate_warmup_combinations(
        self,
        reagent_lists: List[List['Reagent']],
        num_warmup_trials: int,
        disallow_tracker: 'DisallowTracker'
    ) -> List[List[int]]:
        """
        Generate reagent combinations for warmup evaluation.

        Parameters:
        -----------
        reagent_lists : List[List[Reagent]]
            List of reagent lists, one for each component
        num_warmup_trials : int
            Number of trials per reagent
        disallow_tracker : DisallowTracker
            Tracker to prevent duplicate combinations

        Returns:
        --------
        List[List[int]]
            List of combinations, where each combination is [idx_comp1, idx_comp2, ...]
        """
        pass

    def get_expected_evaluations(
        self,
        reagent_lists: List[List['Reagent']],
        num_warmup_trials: int
    ) -> int:
        """
        Calculate expected number of evaluations for this strategy.

        Default implementation: sum of reagents × trials
        Override for strategies with different evaluation counts.

        Parameters:
        -----------
        reagent_lists : List[List[Reagent]]
            List of reagent lists
        num_warmup_trials : int
            Number of trials per reagent

        Returns:
        --------
        int
            Expected number of evaluations
        """
        return sum(len(rl) for rl in reagent_lists) * num_warmup_trials

    @abstractmethod
    def get_name(self) -> str:
        """
        Return strategy name for logging.

        Returns:
        --------
        str
            Human-readable strategy name
        """
        pass

    def get_description(self) -> str:
        """
        Return detailed strategy description.

        Returns:
        --------
        str
            Multi-line description of the strategy
        """
        return self.__doc__ or "No description available"
