"""
Reaction Sequence for orchestrating single-step and multi-step syntheses.

This module provides unified handling for:
1. Single-step reactions with a single SMARTS pattern
2. Single-step reactions with alternative SMARTS patterns
3. Multi-step synthesis with intermediate tracking
"""

from pydantic import BaseModel, field_validator, model_validator
from typing import List, Dict, Optional, Tuple, Set
from enum import Enum
from rdkit import Chem
from rdkit.Chem import AllChem
import logging

from .smarts_router import AlternativeSMARTSConfig, SMARTSRouter


class ReactionInputSource(str, Enum):
    """Source of input for a reaction step."""
    REAGENT_FILE = "reagent_file"      # From input .smi file
    PREVIOUS_STEP = "previous_step"     # From output of earlier step


class ReactionInputConfig(BaseModel):
    """
    Configuration for one input to a reaction step.

    Attributes:
        source: Where the input comes from
        component_idx: If source is REAGENT_FILE, which component (0-indexed)
        step_idx: If source is PREVIOUS_STEP, which step's output (0-indexed)

    Example (from reagent file):
        ReactionInputConfig(source="reagent_file", component_idx=0)

    Example (from previous step):
        ReactionInputConfig(source="previous_step", step_idx=1)
    """
    source: ReactionInputSource
    component_idx: Optional[int] = None
    step_idx: Optional[int] = None

    @model_validator(mode='after')
    def validate_source_fields(self) -> 'ReactionInputConfig':
        """Validate that correct fields are set based on source."""
        if self.source == ReactionInputSource.REAGENT_FILE:
            if self.component_idx is None:
                raise ValueError(
                    "component_idx required when source is 'reagent_file'"
                )
        elif self.source == ReactionInputSource.PREVIOUS_STEP:
            if self.step_idx is None:
                raise ValueError(
                    "step_idx required when source is 'previous_step'"
                )
        return self


class ReactionStepConfig(BaseModel):
    """
    Configuration for a single step in a multi-step synthesis.

    Attributes:
        step_idx: 0-indexed step number
        smarts_config: SMARTS configuration (can include alternatives)
        inputs: List of inputs to this reaction
        description: Optional human-readable description
        required_substructure: Optional SMARTS that must be present in inputs
        expected_substructure: Optional SMARTS that should appear in output

    Example:
        ReactionStepConfig(
            step_idx=0,
            smarts_config=AlternativeSMARTSConfig(
                primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"
            ),
            inputs=[
                ReactionInputConfig(source="reagent_file", component_idx=0),
                ReactionInputConfig(source="reagent_file", component_idx=1)
            ],
            description="Amide coupling"
        )
    """
    step_idx: int
    smarts_config: AlternativeSMARTSConfig
    inputs: List[ReactionInputConfig]
    description: Optional[str] = None
    required_substructure: Optional[str] = None
    expected_substructure: Optional[str] = None

    @field_validator('required_substructure', 'expected_substructure')
    @classmethod
    def validate_substructure_smarts(cls, v: Optional[str]) -> Optional[str]:
        """Validate substructure SMARTS if provided."""
        if v is not None:
            mol = Chem.MolFromSmarts(v)
            if mol is None:
                raise ValueError(f"Invalid substructure SMARTS: {v}")
        return v


class ProtectingGroupConfig(BaseModel):
    """
    Configuration for tracking a protecting group through synthesis.

    Attributes:
        name: Human-readable name (e.g., "Boc", "Fmoc")
        protected_smarts: SMARTS pattern to detect protected form
        deprotected_smarts: SMARTS pattern to detect deprotected form
        component_idx: Which reagent component has this protecting group
        protection_removed_at_step: Which step removes the protection

    Example:
        ProtectingGroupConfig(
            name="Boc",
            protected_smarts="[NH]C(=O)OC(C)(C)C",
            deprotected_smarts="[NH2]",
            component_idx=0,
            protection_removed_at_step=1
        )
    """
    name: str
    protected_smarts: str
    deprotected_smarts: str
    component_idx: int
    protection_removed_at_step: int

    @field_validator('protected_smarts', 'deprotected_smarts')
    @classmethod
    def validate_smarts(cls, v: str) -> str:
        """Validate SMARTS patterns."""
        mol = Chem.MolFromSmarts(v)
        if mol is None:
            raise ValueError(f"Invalid SMARTS: {v}")
        return v


class MultiStepSynthesisConfig(BaseModel):
    """
    Complete configuration for single or multi-step synthesis.

    Supports three modes:
    1. Single SMARTS (backwards compatible): Just provide reaction_smarts
    2. Single with alternatives: Provide alternative_smarts
    3. Multi-step: Provide steps list

    Attributes:
        reaction_smarts: Simple single SMARTS (mode 1)
        alternative_smarts: Single step with alternatives (mode 2)
        steps: Multi-step synthesis configuration (mode 3)
        protecting_groups: Optional protecting group tracking
        reagent_file_list: List of reagent file paths

    Example (single):
        MultiStepSynthesisConfig(
            reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            reagent_file_list=["acids.smi", "amines.smi"]
        )

    Example (multi-step):
        MultiStepSynthesisConfig(
            reagent_file_list=["diamines.smi", "acids.smi", "aldehydes.smi"],
            steps=[...],  # ReactionStepConfig list
            protecting_groups=[...]  # Optional
        )
    """
    # Mode 1: Simple single SMARTS
    reaction_smarts: Optional[str] = None

    # Mode 2: Single with alternatives
    alternative_smarts: Optional[AlternativeSMARTSConfig] = None

    # Mode 3: Multi-step
    steps: Optional[List[ReactionStepConfig]] = None
    protecting_groups: Optional[List[ProtectingGroupConfig]] = None

    # Common
    reagent_file_list: List[str]

    @model_validator(mode='after')
    def validate_mode(self) -> 'MultiStepSynthesisConfig':
        """Ensure exactly one mode is specified."""
        modes_set = sum([
            self.reaction_smarts is not None,
            self.alternative_smarts is not None,
            self.steps is not None and len(self.steps) > 0
        ])

        if modes_set == 0:
            raise ValueError(
                "Must specify one of: reaction_smarts, alternative_smarts, or steps"
            )
        if modes_set > 1:
            raise ValueError(
                "Cannot specify multiple modes. Use only one of: "
                "reaction_smarts, alternative_smarts, or steps"
            )

        return self

    def is_multi_step(self) -> bool:
        """Check if this is a multi-step synthesis."""
        return self.steps is not None and len(self.steps) > 0

    def get_smarts_config(self) -> AlternativeSMARTSConfig:
        """
        Get SMARTS config for single-step modes.

        Returns:
            AlternativeSMARTSConfig for mode 1 or 2

        Raises:
            ValueError: If called on multi-step config
        """
        if self.is_multi_step():
            raise ValueError(
                "Cannot get single SMARTS config for multi-step synthesis"
            )

        if self.alternative_smarts is not None:
            return self.alternative_smarts

        return AlternativeSMARTSConfig(primary_smarts=self.reaction_smarts)


class ReactionSequence:
    """
    Unified handler for single-step and multi-step synthesis.

    Uses SMARTSRouter at each step to handle alternative SMARTS patterns.
    Provides a consistent interface regardless of synthesis complexity.

    Attributes:
        config: MultiStepSynthesisConfig defining the synthesis
        logger: Logger instance for debugging

    Usage (single-step):
        config = MultiStepSynthesisConfig(
            reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            reagent_file_list=["acids.smi", "amines.smi"]
        )
        seq = ReactionSequence(config)
        product, patterns = seq.enumerate(reagent_mols, reagent_keys)

    Usage (multi-step):
        config = MultiStepSynthesisConfig(
            reagent_file_list=["diamines.smi", "acids.smi", "aldehydes.smi"],
            steps=[step1, step2, step3]
        )
        seq = ReactionSequence(config)
        product, patterns = seq.enumerate(reagent_mols, reagent_keys)
    """

    def __init__(self, config: MultiStepSynthesisConfig):
        """
        Initialize the reaction sequence.

        Args:
            config: Synthesis configuration
        """
        self.config = config
        self.logger = logging.getLogger(__name__)

        # One SMARTSRouter per step
        self._step_routers: List[SMARTSRouter] = []

        self._initialize()

    def _initialize(self) -> None:
        """Initialize SMARTS routers for each step."""
        if self.config.is_multi_step():
            for step in self.config.steps:
                router = SMARTSRouter(step.smarts_config)
                self._step_routers.append(router)
                self.logger.debug(
                    f"Initialized router for step {step.step_idx}: "
                    f"{step.description or 'no description'}"
                )
        else:
            # Single-step mode
            smarts_config = self.config.get_smarts_config()
            self._step_routers.append(SMARTSRouter(smarts_config))
            self.logger.debug("Initialized single-step router")

    def register_reagent_compatibility(
        self,
        component_idx: int,
        reagent_key: str,
        compatible_patterns: Set[str],
        step_idx: int = 0
    ) -> None:
        """
        Register reagent compatibility with a specific step's router.

        For single-step synthesis, step_idx should be 0.
        For multi-step, register with the appropriate step.

        Args:
            component_idx: Which reagent component (for logging/tracking)
            reagent_key: Reagent identifier
            compatible_patterns: Set of pattern IDs
            step_idx: Which step this applies to
        """
        if step_idx < len(self._step_routers):
            self._step_routers[step_idx].register_reagent_compatibility(
                reagent_key, compatible_patterns
            )

    def register_all_reagent_compatibilities(
        self,
        reagent_lists: List[List]  # List[List[Reagent]] but avoid circular import
    ) -> None:
        """
        Bulk register all reagent compatibilities from reagent lists.

        Reads compatibility from Reagent.compatible_smarts attribute.
        For single-step, registers all reagents with step 0.
        For multi-step, you may need custom logic based on which
        reagents participate in which steps.

        Args:
            reagent_lists: List of reagent lists (one per component)
        """
        if not self.config.is_multi_step():
            # Single-step: register all with step 0
            router = self._step_routers[0]
            for comp_idx, reagent_list in enumerate(reagent_lists):
                for reagent in reagent_list:
                    router.register_reagent_compatibility(
                        reagent.reagent_key,
                        reagent.compatible_smarts
                    )
        else:
            # Multi-step: register based on which steps use which components
            for step_idx, step_config in enumerate(self.config.steps):
                router = self._step_routers[step_idx]
                for input_config in step_config.inputs:
                    if input_config.source == ReactionInputSource.REAGENT_FILE:
                        comp_idx = input_config.component_idx
                        for reagent in reagent_lists[comp_idx]:
                            router.register_reagent_compatibility(
                                reagent.reagent_key,
                                reagent.compatible_smarts
                            )

    def enumerate(
        self,
        reagent_mols: List[Chem.Mol],
        reagent_keys: Optional[List[str]] = None
    ) -> Tuple[Optional[Chem.Mol], Dict[int, str]]:
        """
        Execute the full synthesis sequence.

        Args:
            reagent_mols: Reagent Mol objects, one per component
            reagent_keys: Optional reagent identifiers for SMARTS routing

        Returns:
            Tuple of:
            - product: Final product Mol, or None if synthesis failed
            - patterns_used: Dict mapping step_idx -> pattern_id used
        """
        if not self.config.is_multi_step():
            return self._enumerate_single_step(reagent_mols, reagent_keys)
        else:
            return self._enumerate_multi_step(reagent_mols, reagent_keys)

    def _enumerate_single_step(
        self,
        reagent_mols: List[Chem.Mol],
        reagent_keys: Optional[List[str]]
    ) -> Tuple[Optional[Chem.Mol], Dict[int, str]]:
        """Execute single-step synthesis."""
        product, pattern_id = self._step_routers[0].enumerate(
            reagent_mols, reagent_keys
        )
        patterns_used = {0: pattern_id} if pattern_id else {}
        return product, patterns_used

    def _enumerate_multi_step(
        self,
        reagent_mols: List[Chem.Mol],
        reagent_keys: Optional[List[str]]
    ) -> Tuple[Optional[Chem.Mol], Dict[int, str]]:
        """Execute multi-step synthesis with intermediate handling."""
        intermediates: Dict[int, Chem.Mol] = {}
        patterns_used: Dict[int, str] = {}

        for step_idx, step_config in enumerate(self.config.steps):
            # Gather inputs for this step
            inputs: List[Chem.Mol] = []
            input_keys: List[str] = []

            for input_config in step_config.inputs:
                if input_config.source == ReactionInputSource.REAGENT_FILE:
                    comp_idx = input_config.component_idx
                    inputs.append(reagent_mols[comp_idx])
                    if reagent_keys:
                        input_keys.append(reagent_keys[comp_idx])
                else:
                    # Input from previous step
                    prev_idx = input_config.step_idx
                    prev_mol = intermediates.get(prev_idx)
                    if prev_mol is None:
                        self.logger.debug(
                            f"Step {step_idx}: Missing intermediate from step {prev_idx}"
                        )
                        return None, patterns_used
                    inputs.append(prev_mol)
                    input_keys.append(f"intermediate_{prev_idx}")

            # Validate required substructure
            if step_config.required_substructure:
                if not self._validate_substructure(
                    inputs, step_config.required_substructure
                ):
                    self.logger.debug(
                        f"Step {step_idx}: Required substructure not found"
                    )
                    return None, patterns_used

            # Run this step's reaction
            router = self._step_routers[step_idx]
            product, pattern_id = router.enumerate(
                inputs,
                input_keys if input_keys else None
            )

            if product is None:
                self.logger.debug(f"Step {step_idx}: Enumeration failed")
                return None, patterns_used

            patterns_used[step_idx] = pattern_id
            intermediates[step_idx] = product

            # Validate expected output substructure
            if step_config.expected_substructure:
                if not self._has_substructure(
                    product, step_config.expected_substructure
                ):
                    self.logger.debug(
                        f"Step {step_idx}: Expected substructure not in product"
                    )
                    return None, patterns_used

        # Return final product
        final_idx = len(self.config.steps) - 1
        return intermediates.get(final_idx), patterns_used

    def _validate_substructure(
        self,
        mols: List[Chem.Mol],
        smarts: str
    ) -> bool:
        """Check if any molecule contains the required substructure."""
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            self.logger.warning(f"Invalid SMARTS pattern: {smarts}")
            return True  # Skip validation for invalid SMARTS
        return any(mol.HasSubstructMatch(pattern) for mol in mols)

    def _has_substructure(self, mol: Chem.Mol, smarts: str) -> bool:
        """Check if molecule contains the specified substructure."""
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            return True
        return mol.HasSubstructMatch(pattern)

    @property
    def num_components(self) -> int:
        """Number of reagent file components required."""
        return len(self.config.reagent_file_list)

    @property
    def num_steps(self) -> int:
        """Number of reaction steps."""
        return len(self._step_routers)

    def get_reaction(
        self,
        step_idx: int = 0,
        pattern_id: str = "primary"
    ) -> Optional[AllChem.ChemicalReaction]:
        """
        Get reaction object for backwards compatibility.

        Args:
            step_idx: Which step
            pattern_id: Which pattern within that step

        Returns:
            ChemicalReaction object or None
        """
        if step_idx < len(self._step_routers):
            return self._step_routers[step_idx].get_reaction(pattern_id)
        return None
