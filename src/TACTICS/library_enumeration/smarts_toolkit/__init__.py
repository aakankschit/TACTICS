"""
SMARTS Toolkit - Reaction definition and validation for drug discovery.

This toolkit provides:
- ReactionDef: Define and validate individual reactions
- ReactionConfig: Configure single or multi-step synthesis
- Validation utilities for reagent compatibility

Note: SynthesisPipeline has moved to TACTICS.library_enumeration.

Quick Start
-----------

**Define a reaction:**

    from TACTICS.library_enumeration.smarts_toolkit import ReactionDef, ReactionConfig

    rxn = ReactionDef(
        reaction_smarts="[C:1][NH2:2].[C:3](=O)Cl>>[C:1][NH:2][C:3]=O",
        step_index=0,
        description="Amide coupling"
    )

**Validate against reagents:**

    result = rxn.validate_reaction(reagent_files=["amines.smi", "acids.smi"])
    print(f"Coverage: {result.coverage_stats}")

**Create a configuration:**

    config = ReactionConfig(
        reactions=[rxn],
        reagent_file_list=["amines.smi", "acids.smi"]
    )

**Use with SynthesisPipeline:**

    from TACTICS.library_enumeration import SynthesisPipeline

    pipeline = SynthesisPipeline(config)
    result = pipeline.enumerate_single(reagent_mols, reagent_keys)

Classes
-------
ReactionDef
    Definition of a single reaction with built-in validation.

ReactionConfig
    Container for one or more ReactionDef objects.

DeprotectionSpec
    Specifies deprotection to apply before a reaction.

StepInput
    Configuration for input source to a reaction step.

InputSource
    Enum: REAGENT_FILE or PREVIOUS_STEP.

ValidationResult
    Results from SMARTS validation.

ProtectingGroupInfo
    Definition of a protecting group for detection/removal.
"""

__version__ = "0.5.0"

# Configuration classes
from .config import (
    ReactionDef,
    ReactionConfig,
    DeprotectionSpec,
    StepInput,
    InputSource,
    ProtectingGroupInfo,
)

# Validation result (public, from internal validator)
from ._validator import ValidationResult

# Constants
from .constants import (
    DEFAULT_PROTECTING_GROUPS,
    DEFAULT_SALT_FRAGMENTS,
    SALT_FRAGMENT_SMILES,
    PROTECTING_GROUP_MAP,
    get_protecting_group,
    get_all_protecting_group_names,
)

__all__ = [
    # Config
    "ReactionDef",
    "ReactionConfig",
    "DeprotectionSpec",
    "StepInput",
    "InputSource",
    "ProtectingGroupInfo",
    # Validation
    "ValidationResult",
    # Constants
    "DEFAULT_PROTECTING_GROUPS",
    "DEFAULT_SALT_FRAGMENTS",
    "SALT_FRAGMENT_SMILES",
    "PROTECTING_GROUP_MAP",
    "get_protecting_group",
    "get_all_protecting_group_names",
]
