from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Tuple
import numpy as np

class SelectionStrategy(ABC):
    """Abstract base class for reagent selection strategies"""

    def __init__(self, mode: str = "maximize"):
        self.mode = mode

    @abstractmethod
    def select_reagent(self,
                      reagent_list: List,
                      disallow_mask: set = None,
                      **kwargs) -> int:
        """
        Select a reagent index from the list

        Parameters:
        -----------
        reagent_list : List
            List of reagent objects with posterior distributions
        disallow_mask : set
            Indices to exclude from selection
        **kwargs : dict
            Strategy-specific parameters

        Returns:
        --------
        int : Selected reagent index
        """
        pass

    def select_batch(self,
                    reagent_list: List,
                    batch_size: int,
                    disallow_mask: set = None,
                    **kwargs) -> np.ndarray:
        """
        Select multiple reagent indices from the list (batch selection)

        Default implementation: call select_reagent multiple times.
        Strategies can override this for more efficient batch selection.

        Parameters:
        -----------
        reagent_list : List
            List of reagent objects with posterior distributions
        batch_size : int
            Number of reagents to select
        disallow_mask : set
            Indices to exclude from selection
        **kwargs : dict
            Strategy-specific parameters

        Returns:
        --------
        np.ndarray : Array of selected reagent indices
        """
        return np.array([
            self.select_reagent(reagent_list, disallow_mask, **kwargs)
            for _ in range(batch_size)
        ])

    def get_component_criticality(self, reagent_list: List) -> Optional[float]:
        """Return criticality score for a component, or None if not supported.

        Strategies with CATS (e.g., RouletteWheelSelection, BayesUCBSelection)
        override this to compute component criticality.

        Returns:
            Criticality score >= 0, or None if the strategy doesn't compute it.
        """
        return None

    def get_component_state(
        self,
        reagent_list: List,
        component_idx: int,
        current_cycle: int,
        total_cycles: int,
    ) -> Optional[Dict[str, Any]]:
        """Return full intermediate state for a component, or None if not supported.

        CATS-aware strategies (RouletteWheelSelection, BayesUCBSelection) override
        this to expose the complete criticality + temperature/percentile pipeline.

        Returns:
            Dict with all intermediate values, or None if the strategy
            doesn't compute component state.
        """
        return None

    def prepare_scores(self, reagent_list: List, rng: np.random.Generator) -> np.ndarray:
        """Sample scores from posterior distributions"""
        stds = np.array([r.std for r in reagent_list])
        mu = np.array([r.mean for r in reagent_list])
        return rng.normal(size=len(reagent_list)) * stds + mu