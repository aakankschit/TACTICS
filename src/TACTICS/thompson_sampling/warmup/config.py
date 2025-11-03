"""Pydantic configuration models for warmup strategies."""

from pydantic import BaseModel
from typing import Literal


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


class StratifiedWarmupConfig(BaseModel):
    """
    Configuration for Stratified warmup strategy (recommended).

    Uses random partner selection WITHOUT replacement, guaranteeing that each trial
    provides new information with no duplicate partners within a reagent's trials.

    This is the recommended default warmup strategy for most use cases.
    """

    warmup_type: Literal["stratified"] = "stratified"


class LatinHypercubeWarmupConfig(BaseModel):
    """
    Configuration for Latin Hypercube Sampling warmup strategy.

    Advanced strategy ensuring optimal space-filling properties. Divides partner
    space into strata and samples one partner from each stratum, guaranteeing
    coverage across the entire reagent range.

    Best for:
    - Reagents ordered by molecular properties
    - When diversity in warmup is critical
    - High-stakes applications requiring robust initialization
    """

    warmup_type: Literal["latin_hypercube"] = "latin_hypercube"
