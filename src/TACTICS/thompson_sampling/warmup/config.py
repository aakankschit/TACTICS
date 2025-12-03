"""Pydantic configuration models for warmup strategies."""

from pydantic import BaseModel, Field
from typing import Literal, Optional


class StandardWarmupConfig(BaseModel):
    """
    Configuration for Standard warmup strategy.

    Standard warmup uses random partner selection with replacement.
    Each reagent is tested num_warmup_trials times with randomly selected partners.
    """

    warmup_type: Literal["standard"] = "standard"


class EnhancedWarmupConfig(BaseModel):
    """
    Configuration for Enhanced warmup strategy (legacy).

    Uses stochastic parallel pairing where all reagents are shuffled and paired
    exhaustively in each trial. Small components get over-sampled relative to
    large components.

    WARNING: Creates imbalanced posteriors. Use only if you specifically want
    comprehensive small-component coverage.
    """

    warmup_type: Literal["enhanced"] = "enhanced"


class BalancedWarmupConfig(BaseModel):
    """
    Configuration for Balanced warmup strategy (recommended).

    Guarantees exactly K observations per reagent using stratified partner selection.
    This ensures:
    1. Every reagent gets exactly `observations_per_reagent` observations
    2. Partners are selected from K different strata (no duplicates)
    3. Optional seeded RNG for reproducibility

    Total evaluations: sum(component_sizes) x observations_per_reagent

    Example:
        For 130 acids x 3844 amines with K=5:
        Total = (130 + 3844) x 5 = 19,870 evaluations
        Every reagent gets exactly 5 observations.
    """

    warmup_type: Literal["balanced"] = "balanced"

    observations_per_reagent: int = Field(
        default=5,
        ge=3,
        le=50,
        description="Number of observations guaranteed per reagent (K). "
                   "Must be >= 3 for reliable variance estimation."
    )

    seed: Optional[int] = Field(
        default=None,
        description="Random seed for reproducibility. None = random each run."
    )

    use_per_reagent_variance: bool = Field(
        default=True,
        description="If True, estimate variance per-reagent using warmup observations "
                   "with James-Stein shrinkage. If False, use global variance."
    )

    shrinkage_strength: float = Field(
        default=3.0,
        gt=0,
        description="Shrinkage parameter for per-reagent variance estimation. "
                   "Higher = more regularization toward global variance. "
                   "With n observations: weight = n / (n + shrinkage_strength)"
    )
