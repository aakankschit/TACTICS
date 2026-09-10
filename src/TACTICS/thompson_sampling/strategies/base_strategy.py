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
        """Select one reagent index from a component's reagent list.

        The sampler calls this once per component per cycle with the keyword
        context ``rng``, ``component_idx``, ``iteration``, ``current_cycle``
        and ``total_cycles``; strategies read what they need from ``kwargs``.

        Args:
            reagent_list: Reagent objects with posterior ``mean``/``std``/``n_samples``.
            disallow_mask: Indices that must not be selected (already sampled
                in combination with the other components' current picks).
            **kwargs: Per-cycle context from the sampler (see above).

        Returns:
            The selected index into ``reagent_list``.
        """
        pass

    def select_batch(self,
                    reagent_list: List,
                    batch_size: int,
                    disallow_mask: set = None,
                    **kwargs) -> np.ndarray:
        """Select ``batch_size`` reagent indices (with replacement).

        Default implementation calls :meth:`select_reagent` ``batch_size``
        times. Not used by :class:`~TACTICS.thompson_sampling.core.sampler.ThompsonSampler`,
        which builds batches itself; kept for strategies used standalone.

        Args:
            reagent_list: Reagent objects with posterior distributions.
            batch_size: Number of indices to return.
            disallow_mask: Indices to exclude from selection.
            **kwargs: Passed through to :meth:`select_reagent`.

        Returns:
            Array of selected indices, length ``batch_size``.
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