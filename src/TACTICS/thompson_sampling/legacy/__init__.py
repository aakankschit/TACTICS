"""DEPRECATED: Legacy Thompson Sampling implementations.

These implementations are retained only for reproducing published results
from Klarich (2024) and Zhao et al. (2025). For new work, use
:class:`TACTICS.thompson_sampling.ThompsonSampler` with
:class:`~TACTICS.thompson_sampling.strategies.TopTwoSelection` or
:class:`~TACTICS.thompson_sampling.strategies.RouletteWheelSelection`.
"""

import warnings

warnings.warn(
    "TACTICS.thompson_sampling.legacy is deprecated. "
    "Use TACTICS.thompson_sampling.ThompsonSampler with TopTwoSelection or "
    "RouletteWheelSelection instead.",
    DeprecationWarning,
    stacklevel=2,
)
