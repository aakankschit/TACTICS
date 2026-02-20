"""
Correlated UCB (C-UCB) posterior updater.

Propagates score observations to chemically similar reagents via
fingerprint similarity, speeding up posterior learning for all strategies.

Uses a hybrid approach:
- BulkTanimotoSimilarity builds the per-reagent neighbor graph
- iSIM complementary similarity weights propagation strength
  (typical reagents get more propagation, outliers get less)
"""

import logging
import warnings
from typing import List, Optional

import numpy as np
from rdkit import DataStructs

logger = logging.getLogger(__name__)


def _compute_complementary_weights(fp_list) -> np.ndarray:
    """
    Compute iSIM-based complementary similarity weights per reagent.

    Uses iSIM to measure how "typical" each reagent is within its component.
    Typical reagents (low complementary sim) get high propagation weight;
    outliers (high complementary sim) get low propagation weight.

    Falls back to uniform weights if iSIM is not installed.

    Parameters:
        fp_list: List of RDKit fingerprint bit vectors

    Returns:
        Array of weights in [0, 1], shape (n_reagents,)
    """
    n = len(fp_list)
    if n < 3:
        return np.ones(n)

    try:
        from iSIM import calculate_medoid, calculate_comp_sim
    except ImportError:
        warnings.warn(
            "iSIM not installed — using uniform complementary weights. "
            "Install with: conda install -c conda-forge isim",
            UserWarning,
            stacklevel=2,
        )
        return np.ones(n)

    # Convert fingerprints to binary numpy matrix (N x nbits)
    nbits = fp_list[0].GetNumBits()
    fp_matrix = np.zeros((n, nbits), dtype=np.int8)
    for i, fp in enumerate(fp_list):
        on_bits = list(fp.GetOnBits())
        if on_bits:
            fp_matrix[i, on_bits] = 1

    # Compute complementary similarity for each reagent
    # comp_sim[i] = how much similarity changes when reagent i is removed
    # High comp_sim → outlier, Low comp_sim → typical
    try:
        comp_sims = calculate_comp_sim(fp_matrix, n_ary="JT")
        comp_sims = np.array(comp_sims, dtype=float)

        # Normalize: low comp_sim → high weight (typical), high → low weight (outlier)
        if comp_sims.max() - comp_sims.min() > 1e-10:
            normalized = (comp_sims - comp_sims.min()) / (comp_sims.max() - comp_sims.min())
            weights = 1.0 - normalized  # Invert: typical = high weight
        else:
            weights = np.ones(n)

        return np.clip(weights, 0.3, 1.0)  # Floor at 0.3 to suppress outlier propagation
    except Exception as e:
        logger.warning(f"iSIM complementary similarity failed: {e}. Using uniform weights.")
        return np.ones(n)


class CorrelatedUpdater:
    """
    Propagates score observations to chemically similar reagents.

    For each reagent that receives a real score, nearby reagents (by
    Tanimoto similarity of Morgan fingerprints) receive a weighted
    pseudo-observation via Reagent.add_pseudo_score().

    Parameters:
        reagent_lists: List of lists of Reagent objects (one per component)
        similarity_threshold: Minimum Tanimoto similarity to be a neighbor (default: 0.3)
        max_neighbors: Maximum neighbors per reagent (default: 10)
    """

    def __init__(
        self,
        reagent_lists: List[list],
        similarity_threshold: float = 0.3,
        max_neighbors: int = 10,
    ):
        self.similarity_threshold = similarity_threshold
        self.max_neighbors = max_neighbors
        self.n_components = len(reagent_lists)

        # Pre-compute neighbor maps and complementary weights
        self._neighbor_maps: List[List[List[tuple]]] = []
        self._complementary_weights: List[np.ndarray] = []
        self._build(reagent_lists)

    def _build(self, reagent_lists: List[list]) -> None:
        """Build neighbor maps and complementary weights for all components."""
        total_neighbors = 0
        for comp_idx, reagent_list in enumerate(reagent_lists):
            n = len(reagent_list)

            # Collect fingerprints
            fps = []
            valid_mask = []
            for r in reagent_list:
                fp = r.fp  # Uses lazy-cached property
                fps.append(fp)
                valid_mask.append(fp is not None)

            # Compute complementary weights via iSIM
            valid_fps = [fp for fp, v in zip(fps, valid_mask) if v]
            if valid_fps:
                comp_weights = _compute_complementary_weights(valid_fps)
            else:
                comp_weights = np.ones(n)

            # Map valid_fps indices back to full indices
            full_weights = np.ones(n) * 0.1  # Default low weight for invalid fps
            valid_idx = 0
            for i in range(n):
                if valid_mask[i]:
                    full_weights[i] = comp_weights[valid_idx]
                    valid_idx += 1
            self._complementary_weights.append(full_weights)

            # Build neighbor lists using BulkTanimotoSimilarity
            neighbor_map = [[] for _ in range(n)]
            for i in range(n):
                if not valid_mask[i]:
                    continue

                # Compute similarities to all other reagents
                sims = DataStructs.BulkTanimotoSimilarity(
                    fps[i], [fp for fp in fps if fp is not None]
                )

                # Map back to full indices and filter
                candidates = []
                valid_j = 0
                for j in range(n):
                    if not valid_mask[j]:
                        continue
                    if j == i:
                        valid_j += 1
                        continue
                    sim = sims[valid_j]
                    valid_j += 1
                    if sim >= self.similarity_threshold:
                        # effective_weight = tanimoto_sim * complement_weight[j]
                        eff_weight = sim * full_weights[j]
                        candidates.append((j, eff_weight))

                # Keep top-K by effective weight
                candidates.sort(key=lambda x: x[1], reverse=True)
                neighbor_map[i] = candidates[: self.max_neighbors]
                total_neighbors += len(neighbor_map[i])

            self._neighbor_maps.append(neighbor_map)

        avg_neighbors = total_neighbors / max(
            sum(len(rl) for rl in reagent_lists), 1
        )
        logger.info(
            f"CorrelatedUpdater: {self.n_components} components, "
            f"avg {avg_neighbors:.1f} neighbors/reagent, "
            f"threshold={self.similarity_threshold}, max_k={self.max_neighbors}"
        )

    def propagate(
        self,
        component_idx: int,
        reagent_idx: int,
        score: float,
        reagent_lists: List[list],
    ) -> int:
        """
        Propagate a score observation to similar reagents.

        Parameters:
            component_idx: Which reaction component
            reagent_idx: Index of the reagent that received the real score
            score: The observed score
            reagent_lists: Current reagent lists (for calling add_pseudo_score)

        Returns:
            Number of pseudo-observations propagated
        """
        if component_idx >= self.n_components:
            return 0

        neighbors = self._neighbor_maps[component_idx][reagent_idx]
        count = 0
        for neighbor_idx, eff_weight in neighbors:
            reagent_lists[component_idx][neighbor_idx].add_pseudo_score(
                score, eff_weight
            )
            count += 1

        return count
