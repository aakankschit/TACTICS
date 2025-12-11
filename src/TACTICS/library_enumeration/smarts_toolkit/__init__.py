"""
SMARTS Toolkit - Troubleshooting, validation, and multi-SMARTS support for reaction patterns
"""

from .pattern_validator import (
    SMARTSValidator,
    ValidationResult,
    ProtectingGroupInfo,
    DEFAULT_PROTECTING_GROUPS,
    DEFAULT_SALT_FRAGMENTS,
)
from .visualization import SMARTSVisualizer

# Multi-SMARTS and multi-step synthesis support
from .smarts_router import (
    SMARTSRouter,
    SMARTSPatternConfig,
    AlternativeSMARTSConfig
)
from .reaction_sequence import (
    ReactionSequence,
    ReactionInputSource,
    ReactionInputConfig,
    ReactionStepConfig,
    ProtectingGroupConfig,
    MultiStepSynthesisConfig
)

__all__ = [
    # Validation and analysis tools
    'SMARTSValidator',
    'ValidationResult',
    'ProtectingGroupInfo',
    'DEFAULT_PROTECTING_GROUPS',
    'DEFAULT_SALT_FRAGMENTS',
    'SMARTSVisualizer',
    # Multi-SMARTS routing
    'SMARTSRouter',
    'SMARTSPatternConfig',
    'AlternativeSMARTSConfig',
    # Multi-step synthesis
    'ReactionSequence',
    'ReactionInputSource',
    'ReactionInputConfig',
    'ReactionStepConfig',
    'ProtectingGroupConfig',
    'MultiStepSynthesisConfig'
]

__version__ = '0.2.0'