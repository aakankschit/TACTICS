from .generate_products import LibraryEnumerator
from .enumeration_utils import write_products_to_files
from .multiprocessing_utils import initializer

# Re-export SMARTS toolkit classes for convenience
from .smarts_toolkit import (
    SMARTSRouter,
    SMARTSPatternConfig,
    AlternativeSMARTSConfig,
    ReactionSequence,
    ReactionInputSource,
    ReactionInputConfig,
    ReactionStepConfig,
    ProtectingGroupConfig,
    MultiStepSynthesisConfig
)

__all__ = [
    'LibraryEnumerator',
    'initializer',
    'write_products_to_files',
    # SMARTS toolkit
    'SMARTSRouter',
    'SMARTSPatternConfig',
    'AlternativeSMARTSConfig',
    'ReactionSequence',
    'ReactionInputSource',
    'ReactionInputConfig',
    'ReactionStepConfig',
    'ProtectingGroupConfig',
    'MultiStepSynthesisConfig'
]