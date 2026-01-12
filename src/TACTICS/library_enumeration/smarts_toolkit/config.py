"""
Configuration classes for the SMARTS toolkit.

This module provides validated Pydantic configuration for:
- Single reaction definitions (ReactionDef)
- Multi-step synthesis pipelines (ReactionConfig)
- Deprotection specifications (DeprotectionSpec)
- Step input mappings (StepInput)

The configuration system supports three modes:
1. Single SMARTS: One ReactionDef with step_index=0
2. Alternative SMARTS: Multiple ReactionDef with same step_index, marked with step_modes
3. Multi-step: Multiple ReactionDef with different step_index values
"""

from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, List, Literal, Optional, Tuple, Union, TYPE_CHECKING

from pydantic import BaseModel, Field, PrivateAttr, field_validator, model_validator
from rdkit import Chem
from rdkit.Chem import AllChem

if TYPE_CHECKING:
    from ._validator import ValidationResult, _SMARTSValidator


class InputSource(str, Enum):
    """Source of input for a reaction step."""

    REAGENT_FILE = "reagent_file"  # From original reagent files
    PREVIOUS_STEP = "previous_step"  # From intermediate of a previous step


@dataclass
class ProtectingGroupInfo:
    """
    Definition of a protecting group for detection and optional removal.

    Attributes:
        name: Human-readable name (e.g., "Boc", "Fmoc")
        smarts: SMARTS pattern to detect the group
        deprotection_smarts: Optional reaction SMARTS for removal

    Example:
        >>> boc = ProtectingGroupInfo(
        ...     name="Boc",
        ...     smarts="[NX3][C](=O)OC(C)(C)C",
        ...     deprotection_smarts="[N:1][C](=O)OC(C)(C)C>>[N:1]"
        ... )
    """

    name: str
    smarts: str
    deprotection_smarts: Optional[str] = None

    def __post_init__(self):
        """Validate SMARTS patterns on creation."""
        # Validate detection SMARTS
        if Chem.MolFromSmarts(self.smarts) is None:
            raise ValueError(f"Invalid detection SMARTS for {self.name}: {self.smarts}")
        # Validate deprotection SMARTS if provided
        if self.deprotection_smarts:
            if AllChem.ReactionFromSmarts(self.deprotection_smarts) is None:
                raise ValueError(
                    f"Invalid deprotection SMARTS for {self.name}: {self.deprotection_smarts}"
                )


class StepInput(BaseModel):
    """
    Configuration for one input to a reaction step.

    Specifies where the input comes from:
    - REAGENT_FILE: From one of the original reagent files
    - PREVIOUS_STEP: From the output of a previous reaction step

    Attributes:
        source: Where input comes from (reagent_file or previous_step)
        file_index: Index into reagent_file_list (for REAGENT_FILE source)
        step_index: Index of previous step to use output from (for PREVIOUS_STEP source)

    Example:
        >>> # Input from first reagent file
        >>> input1 = StepInput(source=InputSource.REAGENT_FILE, file_index=0)
        >>> # Input from step 0's output
        >>> input2 = StepInput(source=InputSource.PREVIOUS_STEP, step_index=0)
    """

    source: InputSource
    file_index: Optional[int] = Field(
        None, description="Index into reagent_file_list (for REAGENT_FILE source)"
    )
    step_index: Optional[int] = Field(
        None,
        description="Index of previous step to use output from (for PREVIOUS_STEP source)",
    )

    @model_validator(mode="after")
    def validate_source_fields(self):
        """Ensure correct fields are provided based on source type."""
        if self.source == InputSource.REAGENT_FILE:
            if self.file_index is None:
                raise ValueError("file_index required for REAGENT_FILE source")
        elif self.source == InputSource.PREVIOUS_STEP:
            if self.step_index is None:
                raise ValueError("step_index required for PREVIOUS_STEP source")
        return self

    model_config = {"frozen": False}


class DeprotectionSpec(BaseModel):
    """
    Deprotection to apply to a reactant or product.

    Specifies which protecting group to remove and the target molecule.
    The target can be:
    - An integer (0, 1, 2, ...): Index of the reactant to deprotect BEFORE the reaction
    - The string "product": Apply deprotection to the product AFTER the reaction

    Multiple DeprotectionSpec can be applied to a single reaction.

    Attributes:
        group: Name of protecting group to remove (e.g., "Boc", "Fmoc")
        target: Either a reactant index (int >= 0) or "product" for post-reaction deprotection

    Examples:
        >>> # Remove Boc from the first reactant (before reaction)
        >>> deprot = DeprotectionSpec(group="Boc", target=0)
        >>>
        >>> # Remove Fmoc from the product (after reaction)
        >>> deprot = DeprotectionSpec(group="Fmoc", target="product")
    """

    group: str = Field(
        ..., description="Name of protecting group (e.g., 'Boc', 'Fmoc')"
    )
    target: Union[int, Literal["product"]] = Field(
        ...,
        description="Reactant index (int >= 0) for pre-reaction deprotection, "
        "or 'product' for post-reaction deprotection",
    )

    @field_validator("target")
    @classmethod
    def validate_target(cls, v):
        """Ensure target is valid: non-negative int or 'product'."""
        if isinstance(v, int) and v < 0:
            raise ValueError("target must be >= 0 when specifying a reactant index")
        return v

    @property
    def is_product_deprotection(self) -> bool:
        """Return True if this deprotection targets the product."""
        return self.target == "product"

    @property
    def reactant_index(self) -> Optional[int]:
        """Return reactant index if targeting a reactant, else None."""
        return self.target if isinstance(self.target, int) else None

    model_config = {"frozen": False}


class ReactionDef(BaseModel):
    """
    Definition of a single chemical reaction with built-in validation.

    This is the fundamental building block for all synthesis configurations.
    Every reaction must specify a step_index (0 for single-step reactions,
    0/1/2... for sequences).

    Validation is built-in: call validate() to check reagent compatibility.

    Attributes:
        reaction_smarts: The reaction SMARTS string (required)
        step_index: Which step in the sequence (0 = first step)
        pattern_id: Identifier for alternatives (auto-generated if not provided)
        description: Human-readable description
        deprotections: List of deprotections to apply before this reaction

    Example:
        >>> rxn = ReactionDef(
        ...     reaction_smarts="[C:1](=O)O.[N:2]H2>>[C:1](=O)[N:2]",
        ...     step_index=0,
        ...     description="Amide coupling"
        ... )
        >>> result = rxn.validate(reagent_files=["acids.smi", "amines.smi"])
        >>> print(rxn.coverage_stats)
    """

    reaction_smarts: str = Field(..., description="Reaction SMARTS string")
    step_index: int = Field(default=0, ge=0, description="Step index (0 = first step)")
    pattern_id: Optional[str] = Field(
        None, description="Identifier for this pattern (for alternatives)"
    )
    description: Optional[str] = Field(None, description="Human-readable description")
    deprotections: List[DeprotectionSpec] = Field(
        default_factory=list,
        description="Deprotections to apply before this reaction executes",
    )

    # Private attributes for internal state
    _reaction: Optional[Any] = PrivateAttr(default=None)
    _validator: Optional[Any] = PrivateAttr(default=None)
    _validation_result: Optional[Any] = PrivateAttr(default=None)

    @field_validator("reaction_smarts")
    @classmethod
    def validate_smarts(cls, v):
        """Validate that reaction SMARTS is parseable."""
        rxn = AllChem.ReactionFromSmarts(v)
        if rxn is None:
            raise ValueError(f"Invalid reaction SMARTS: {v}")
        return v

    def model_post_init(self, __context):
        """Initialize the RDKit reaction object."""
        self._reaction = AllChem.ReactionFromSmarts(self.reaction_smarts)

    # === Validation Methods ===

    def validate_reaction(
        self,
        reagent_files: Optional[List[str]] = None,
        reagent_smiles: Optional[List[List[Tuple[str, str]]]] = None,
        protecting_groups: Optional[List[ProtectingGroupInfo]] = None,
        deprotect: bool = False,
        desalt: bool = False,
        test_reactions: bool = False,
    ) -> "ValidationResult":
        """
        Validate this reaction against reagent files or SMILES lists.

        Args:
            reagent_files: Paths to reagent .smi files
            reagent_smiles: Direct SMILES lists as [(smiles, name), ...] per position
                           (useful for validating against intermediates)
            protecting_groups: Custom protecting group definitions
            deprotect: Remove protecting groups before checking compatibility
            desalt: Remove salt fragments before checking compatibility
            test_reactions: Actually run reaction on sample combinations

        Returns:
            ValidationResult with comprehensive compatibility information

        Example:
            >>> rxn = ReactionDef(reaction_smarts="[C:1](=O)O.[N:2]>>[C:1](=O)[N:2]")
            >>> result = rxn.validate_reaction(reagent_files=["acids.smi", "amines.smi"])
            >>> print(f"Coverage: {result.coverage_stats}")
        """
        from ._validator import _SMARTSValidator

        # Create validator
        self._validator = _SMARTSValidator(
            reaction_smarts=self.reaction_smarts,
            reagent_files=reagent_files,
            reagent_smiles_lists=reagent_smiles,
            protecting_groups=protecting_groups,
            deprotect_on_import=deprotect,
            desalt_on_import=desalt,
        )

        # Run validation
        self._validation_result = self._validator.validate(
            test_reactions=test_reactions
        )
        return self._validation_result

    def get_compatible_reagents(self, position: int) -> List[Tuple[str, str]]:
        """
        Get reagents compatible with template at position.

        Args:
            position: Reagent position (0-indexed)

        Returns:
            List of (smiles, name) tuples

        Raises:
            ValueError: If validate_reaction() hasn't been called
        """
        if self._validation_result is None:
            raise ValueError(
                "Must call validate_reaction() before get_compatible_reagents()"
            )
        return self._validation_result.compatible_reagents.get(position, [])

    def get_incompatible_reagents(self, position: int) -> List[Tuple[str, str]]:
        """
        Get reagents that don't match template at position.

        Args:
            position: Reagent position (0-indexed)

        Returns:
            List of (smiles, name) tuples

        Raises:
            ValueError: If validate_reaction() hasn't been called
        """
        if self._validation_result is None:
            raise ValueError(
                "Must call validate_reaction() before get_incompatible_reagents()"
            )
        return self._validation_result.incompatible_reagents.get(position, [])

    # === Visualization Methods ===

    def visualize_template_match(
        self,
        smiles: str,
        position: int,
        highlight_color: Tuple[float, float, float] = (0.0, 0.8, 0.0),
        size: Tuple[int, int] = (400, 300),
    ) -> Any:
        """
        Visualize which atoms in a molecule match the reaction template.

        This is the primary troubleshooting tool for understanding why a
        reagent or intermediate doesn't work with the reaction.

        Args:
            smiles: SMILES of molecule to visualize
            position: Which reactant position to check against
            highlight_color: RGB tuple for highlight color (default: green)
            size: Image dimensions (width, height)

        Returns:
            IPython Image for Jupyter display, or PIL Image
        """
        from rdkit.Chem import Draw

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        template = self.get_reactant_template(position)
        if template is None:
            return None

        matches = mol.GetSubstructMatches(template)

        if matches:
            highlight_atoms = list(matches[0])
            img = Draw.MolToImage(
                mol,
                size=size,
                highlightAtoms=highlight_atoms,
                highlightColor=highlight_color,
            )
        else:
            img = Draw.MolToImage(mol, size=size)

        return img

    def visualize_reaction(self, size: Tuple[int, int] = (800, 200)) -> Any:
        """
        Visualize the reaction scheme.

        Args:
            size: Image dimensions

        Returns:
            IPython Image or PIL Image
        """
        from rdkit.Chem import Draw

        if self._reaction is None:
            return None

        img = Draw.ReactionToImage(self._reaction, subImgSize=(size[0] // 3, size[1]))
        return img

    def get_reactant_template(self, position: int) -> Optional[Chem.Mol]:
        """
        Get the RDKit mol template for a reactant position.

        Args:
            position: Reactant position (0-indexed)

        Returns:
            RDKit Mol object representing the reactant pattern
        """
        if self._reaction is None:
            return None

        if position >= self._reaction.GetNumReactantTemplates():
            return None

        return self._reaction.GetReactantTemplate(position)

    # === Properties ===

    @property
    def num_reactants(self) -> int:
        """Number of reactants in this reaction."""
        if self._reaction is None:
            return 0
        return self._reaction.GetNumReactantTemplates()

    @property
    def is_validated(self) -> bool:
        """True if validate() has been called."""
        return self._validation_result is not None

    @property
    def coverage_stats(self) -> Dict[int, float]:
        """
        Coverage percentage per position.

        Returns:
            Dict mapping position to coverage percentage (0-100)

        Raises:
            ValueError: If validate_reaction() hasn't been called
        """
        if self._validation_result is None:
            raise ValueError(
                "Must call validate_reaction() before accessing coverage_stats"
            )
        return self._validation_result.coverage_stats

    @property
    def validation_result(self) -> Optional["ValidationResult"]:
        """Get the cached validation result."""
        return self._validation_result

    # === Utility ===

    def summary(self) -> str:
        """
        Human-readable validation summary.

        Returns:
            Multi-line string with validation summary
        """
        if self._validation_result is None:
            return f"ReactionDef(step={self.step_index}, not validated)"

        lines = [
            f"ReactionDef Summary (step={self.step_index})",
            "=" * 50,
            f"SMARTS: {self.reaction_smarts[:60]}{'...' if len(self.reaction_smarts) > 60 else ''}",
            "",
        ]

        for pos in sorted(self._validation_result.coverage_stats.keys()):
            coverage = self._validation_result.coverage_stats[pos]
            compatible = len(self._validation_result.compatible_reagents.get(pos, []))
            incompatible = len(
                self._validation_result.incompatible_reagents.get(pos, [])
            )
            lines.append(
                f"Position {pos}: {compatible} compatible, {incompatible} incompatible ({coverage:.1f}%)"
            )

        return "\n".join(lines)

    model_config = {"frozen": False, "arbitrary_types_allowed": True}


class ReactionConfig(BaseModel):
    """
    Container for synthesis configuration.

    Rules:
    - Every ReactionDef must have a step_index
    - step_inputs is REQUIRED when len(reactions) > 1
    - step_modes is only used to mark steps with alternative SMARTS patterns
    - pattern_ids are auto-generated for alternatives if not provided

    Attributes:
        reactions: List of ReactionDef objects (one or more)
        reagent_file_list: Paths to reagent files
        step_inputs: Mapping of step_index to list of StepInput (required if multiple reactions)
        step_modes: Mapping of step_index to "alternative" for steps with multiple SMARTS
        protecting_groups: Custom protecting group definitions

    Examples:
        # Single reaction (simplest case)
        >>> config = ReactionConfig(
        ...     reactions=[ReactionDef(reaction_smarts="...", step_index=0)],
        ...     reagent_file_list=["acids.smi", "amines.smi"]
        ... )

        # Alternatives at step 0
        >>> config = ReactionConfig(
        ...     reactions=[
        ...         ReactionDef(reaction_smarts="...", step_index=0, pattern_id="primary"),
        ...         ReactionDef(reaction_smarts="...", step_index=0, pattern_id="secondary"),
        ...     ],
        ...     reagent_file_list=["amines.smi", "acids.smi"],
        ...     step_inputs={0: [StepInput(...), StepInput(...)]},
        ...     step_modes={0: "alternative"}
        ... )

        # Multi-step sequence
        >>> config = ReactionConfig(
        ...     reactions=[
        ...         ReactionDef(reaction_smarts="...", step_index=0),
        ...         ReactionDef(reaction_smarts="...", step_index=1),
        ...     ],
        ...     reagent_file_list=["bb1.smi", "bb2.smi", "bb3.smi"],
        ...     step_inputs={
        ...         0: [StepInput(source=InputSource.REAGENT_FILE, file_index=0), ...],
        ...         1: [StepInput(source=InputSource.PREVIOUS_STEP, step_index=0), ...],
        ...     }
        ... )
    """

    reactions: List[ReactionDef] = Field(..., min_length=1)
    reagent_file_list: List[str] = Field(default_factory=list)
    step_inputs: Optional[Dict[int, List[StepInput]]] = Field(
        None,
        description="Mapping of step_index to input sources (required if multiple reactions)",
    )
    step_modes: Optional[Dict[int, Literal["alternative"]]] = Field(
        None,
        description="Mark steps with alternative SMARTS patterns",
    )
    protecting_groups: Optional[List[ProtectingGroupInfo]] = Field(
        None, description="Custom protecting group definitions"
    )

    @model_validator(mode="after")
    def validate_config(self):
        """
        Validate the configuration:
        1. Require step_inputs when multiple reactions
        2. Auto-generate pattern_ids for alternatives
        3. Validate step_modes references valid step indices
        4. Validate step_inputs references valid file/step indices
        """
        # Get unique step indices
        step_indices = set(r.step_index for r in self.reactions)

        # 1. Require step_inputs when multiple reactions
        if len(self.reactions) > 1 and self.step_inputs is None:
            raise ValueError(
                "step_inputs is required when multiple reactions are defined. "
                "Specify input sources for each step."
            )

        # 2. Auto-generate pattern_ids for steps with alternatives
        reactions_by_step: Dict[int, List[ReactionDef]] = {}
        for rxn in self.reactions:
            if rxn.step_index not in reactions_by_step:
                reactions_by_step[rxn.step_index] = []
            reactions_by_step[rxn.step_index].append(rxn)

        for step_idx, rxns in reactions_by_step.items():
            if len(rxns) > 1:
                # Multiple reactions at same step - ensure pattern_ids
                for i, rxn in enumerate(rxns):
                    if rxn.pattern_id is None:
                        rxn.pattern_id = f"alt_{i}"

        # 3. Validate step_modes references valid step indices
        if self.step_modes:
            for step_idx in self.step_modes.keys():
                if step_idx not in step_indices:
                    raise ValueError(
                        f"step_modes references non-existent step {step_idx}. "
                        f"Valid steps: {sorted(step_indices)}"
                    )
                # Check that step actually has alternatives
                if len(reactions_by_step.get(step_idx, [])) <= 1:
                    raise ValueError(
                        f"step_modes marks step {step_idx} as 'alternative' but "
                        f"it only has one reaction. Remove from step_modes or add alternatives."
                    )

        # 4. Validate step_inputs
        if self.step_inputs:
            for step_idx, inputs in self.step_inputs.items():
                if step_idx not in step_indices:
                    raise ValueError(
                        f"step_inputs references non-existent step {step_idx}. "
                        f"Valid steps: {sorted(step_indices)}"
                    )
                for inp in inputs:
                    if inp.source == InputSource.REAGENT_FILE:
                        if inp.file_index is not None and inp.file_index >= len(
                            self.reagent_file_list
                        ):
                            raise ValueError(
                                f"step_inputs for step {step_idx} references "
                                f"file_index {inp.file_index} but only "
                                f"{len(self.reagent_file_list)} reagent files provided"
                            )
                    elif inp.source == InputSource.PREVIOUS_STEP:
                        if inp.step_index is not None and inp.step_index >= step_idx:
                            raise ValueError(
                                f"step_inputs for step {step_idx} references "
                                f"step_index {inp.step_index} which is not a previous step"
                            )
                        if (
                            inp.step_index is not None
                            and inp.step_index not in step_indices
                        ):
                            raise ValueError(
                                f"step_inputs for step {step_idx} references "
                                f"non-existent step {inp.step_index}"
                            )

        return self

    # === Properties ===

    @property
    def num_steps(self) -> int:
        """Number of unique steps."""
        return len(set(r.step_index for r in self.reactions))

    @property
    def is_multi_step(self) -> bool:
        """True if more than one step."""
        return self.num_steps > 1

    @property
    def steps_with_alternatives(self) -> List[int]:
        """List of step indices that have alternative SMARTS."""
        if self.step_modes is None:
            return []
        return sorted(self.step_modes.keys())

    @property
    def step_indices(self) -> List[int]:
        """Sorted list of all step indices."""
        return sorted(set(r.step_index for r in self.reactions))

    # === Methods ===

    def get_reactions_for_step(self, step_index: int) -> List[ReactionDef]:
        """
        Get all ReactionDef objects for a step (including alternatives).

        Args:
            step_index: The step index

        Returns:
            List of ReactionDef objects for that step
        """
        return [r for r in self.reactions if r.step_index == step_index]

    def get_primary_reaction(self, step_index: int) -> Optional[ReactionDef]:
        """
        Get the primary reaction for a step.

        Returns the reaction with pattern_id='primary', or the first one if not found.

        Args:
            step_index: The step index

        Returns:
            The primary ReactionDef or None if step doesn't exist
        """
        rxns = self.get_reactions_for_step(step_index)
        if not rxns:
            return None

        # Look for pattern_id='primary'
        for rxn in rxns:
            if rxn.pattern_id == "primary":
                return rxn

        # Return first one
        return rxns[0]

    def get_inputs_for_step(self, step_index: int) -> List[StepInput]:
        """
        Get the input configuration for a step.

        For single-reaction configs without step_inputs, auto-generates
        inputs from reagent_file_list.

        Args:
            step_index: The step index

        Returns:
            List of StepInput objects for that step
        """
        if self.step_inputs is not None:
            return self.step_inputs.get(step_index, [])

        # Single reaction case - auto-generate from reagent files
        if len(self.reactions) == 1 and step_index == 0:
            return [
                StepInput(source=InputSource.REAGENT_FILE, file_index=i)
                for i in range(len(self.reagent_file_list))
            ]

        return []

    def has_alternatives_at_step(self, step_index: int) -> bool:
        """
        Check if a step has alternative SMARTS patterns.

        Args:
            step_index: The step index

        Returns:
            True if step has alternatives marked in step_modes
        """
        if self.step_modes is None:
            return False
        return self.step_modes.get(step_index) == "alternative"

    def validate_all(
        self,
        deprotect: bool = False,
        desalt: bool = False,
    ) -> Dict[int, Dict[str, "ValidationResult"]]:
        """
        Validate all reactions in the config.

        Args:
            deprotect: Apply deprotection during validation
            desalt: Apply desalting during validation

        Returns:
            Nested dict: {step_index: {pattern_id: ValidationResult}}
        """
        results: Dict[int, Dict[str, Any]] = {}

        for step_idx in self.step_indices:
            results[step_idx] = {}
            rxns = self.get_reactions_for_step(step_idx)

            # Get reagent files for this step
            inputs = self.get_inputs_for_step(step_idx)
            reagent_files = [
                self.reagent_file_list[inp.file_index]
                for inp in inputs
                if inp.source == InputSource.REAGENT_FILE and inp.file_index is not None
            ]

            for rxn in rxns:
                pattern_id = rxn.pattern_id or "primary"
                result = rxn.validate_reaction(
                    reagent_files=reagent_files,
                    protecting_groups=self.protecting_groups,
                    deprotect=deprotect,
                    desalt=desalt,
                )
                results[step_idx][pattern_id] = result

        return results

    model_config = {"frozen": False}
