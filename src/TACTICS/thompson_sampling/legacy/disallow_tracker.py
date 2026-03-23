"""Backward compatibility shim. DisallowTracker has moved to core/."""

import warnings

warnings.warn(
    "Import DisallowTracker from TACTICS.thompson_sampling.core.disallow_tracker "
    "instead of TACTICS.thompson_sampling.legacy.disallow_tracker. "
    "The legacy module is deprecated.",
    DeprecationWarning,
    stacklevel=2,
)

from ..core.disallow_tracker import DisallowTracker

__all__ = ["DisallowTracker"]
