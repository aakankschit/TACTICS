"""
SynthesisPipeline - Main class for chemical library enumeration.

The SynthesisPipeline is the primary interface for:
- Single compound generation (used by ThompsonSampler)
- Batch compound generation
- Full library enumeration
- Reagent validation and compatibility detection

It wraps ReactionConfig and provides execution capabilities for
single-step, multi-step, and alternative-SMARTS reactions.
"""

import logging
from typing import Any, Dict, List, Optional, Set, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem

from .smarts_toolkit import (
    ReactionConfig,
    ReactionDef,
    InputSource,
    DeprotectionSpec,
)
from .smarts_toolkit._validator import _SMARTSValidator, ValidationResult
from .smarts_toolkit.constants import DEFAULT_PROTECTING_GROUPS
from .enumeration_utils import (
    EnumerationResult,
    EnumerationError,
    AutoDetectionResult,
    read_reagent_file,
)
from .generate_products import (
    enumerate_products,
    run_single_reaction,
    apply_deprotections,
    apply_product_deprotections,
    generate_all_combinations,
)

logger = logging.getLogger(__name__)


class SynthesisPipeline:
    """
    Synthesis pipeline for single-step, multi-step, and alternative-SMARTS reactions.

    This is the main entry point for:
    - ThompsonSampler compound generation
    - Full library enumeration
    - Reagent validation and compatibility detection

    Attributes:
        config: ReactionConfig defining the synthesis
        reagent_file_list: List of reagent file paths

    Example:
        >>> from TACTICS.library_enumeration import SynthesisPipeline
        >>> from TACTICS.library_enumeration.smarts_toolkit import ReactionDef, ReactionConfig
        >>>
        >>> config = ReactionConfig(
        ...     reactions=[ReactionDef(reaction_smarts="...", step_index=0)],
        ...     reagent_file_list=["acids.smi", "amines.smi"]
        ... )
        >>> pipeline = SynthesisPipeline(config)
        >>> result = pipeline.enumerate_single(reagent_mols, reagent_keys)
    """

    def __init__(self, config: ReactionConfig):
        """
        Initialize pipeline from ReactionConfig.

        Args:
            config: ReactionConfig with reactions and reagent files
        """
        self.config = config

        # Organize reactions by step_index
        self._reactions_by_step: Dict[int, List[ReactionDef]] = {}
        for rxn in config.reactions:
            if rxn.step_index not in self._reactions_by_step:
                self._reactions_by_step[rxn.step_index] = []
            self._reactions_by_step[rxn.step_index].append(rxn)

        # Initialize RDKit reactions: {step_index: {pattern_id: ChemicalReaction}}
        self._rdkit_reactions: Dict[int, Dict[str, Any]] = {}
        self._init_reactions()

        # Compatibility cache: {(step_index, position, reagent_key): Set[pattern_ids]}
        self._compatibility_cache: Dict[Tuple[int, int, str], Set[str]] = {}

        # Validators: {step_index: {pattern_id: _SMARTSValidator}}
        self._validators: Dict[int, Dict[str, _SMARTSValidator]] = {}

        # Protecting groups for deprotection
        self._protecting_groups: Dict[str, Any] = {
            pg.name: pg for pg in DEFAULT_PROTECTING_GROUPS
        }
        if config.protecting_groups:
            for pg in config.protecting_groups:
                self._protecting_groups[pg.name] = pg

    def _init_reactions(self) -> None:
        """Build RDKit reaction objects from ReactionDef list."""
        for step_idx, rxns in self._reactions_by_step.items():
            self._rdkit_reactions[step_idx] = {}

            for rxn in rxns:
                pattern_id = rxn.pattern_id or "primary"
                rdkit_rxn = AllChem.ReactionFromSmarts(rxn.reaction_smarts)
                if rdkit_rxn is None:
                    logger.error(
                        f"Failed to parse SMARTS for step {step_idx}, "
                        f"pattern {pattern_id}"
                    )
                else:
                    self._rdkit_reactions[step_idx][pattern_id] = rdkit_rxn

    # === Properties ===

    @property
    def reagent_file_list(self) -> List[str]:
        """Reagent file paths from config."""
        return self.config.reagent_file_list

    @property
    def num_steps(self) -> int:
        """Number of reaction steps."""
        return len(self._reactions_by_step)

    @property
    def num_components(self) -> int:
        """Number of reagent components (files)."""
        return len(self.config.reagent_file_list)

    @property
    def has_alternatives(self) -> bool:
        """True if any step has alternative SMARTS patterns."""
        return self.config.step_modes is not None and len(self.config.step_modes) > 0

    @property
    def is_multi_step(self) -> bool:
        """True if this is a multi-step synthesis."""
        return len(self._reactions_by_step) > 1

    @property
    def pattern_ids(self) -> Dict[int, List[str]]:
        """Available pattern IDs at each step."""
        result = {}
        for step_idx, rxns in self._reactions_by_step.items():
            result[step_idx] = [r.pattern_id or "primary" for r in rxns]
        return result

    @property
    def reactions(self) -> List[ReactionDef]:
        """Get all ReactionDef objects."""
        return self.config.reactions

    # === Single Compound Enumeration ===

    def enumerate_single(
        self,
        reagent_mols: List[Chem.Mol],
        reagent_keys: Optional[List[str]] = None,
        store_intermediates: bool = False,
    ) -> EnumerationResult:
        """
        Enumerate a single product from reagent molecules.

        Used by ThompsonSampler for individual compound generation.

        Args:
            reagent_mols: List of RDKit Mol objects (one per reagent position)
            reagent_keys: Optional reagent identifiers for routing
            store_intermediates: If True, store intermediate products

        Returns:
            EnumerationResult with product or error details
        """
        patterns_used: Dict[int, str] = {}
        intermediates: Optional[Dict[int, Chem.Mol]] = (
            {} if store_intermediates else None
        )

        # Track intermediate products from previous steps
        step_outputs: Dict[int, Chem.Mol] = {}

        # Generate product name from reagent keys
        product_name = None
        if reagent_keys:
            product_name = "_".join(reagent_keys)

        for step_idx in sorted(self._reactions_by_step.keys()):
            rxns = self._reactions_by_step[step_idx]
            inputs = self.config.get_inputs_for_step(step_idx)

            # Collect inputs for this step
            step_inputs: List[Chem.Mol] = []
            step_input_keys: List[str] = []

            for i, inp in enumerate(inputs):
                if inp.source == InputSource.REAGENT_FILE:
                    if inp.file_index is not None and inp.file_index < len(
                        reagent_mols
                    ):
                        mol = reagent_mols[inp.file_index]
                        step_inputs.append(mol)
                        if reagent_keys and inp.file_index < len(reagent_keys):
                            step_input_keys.append(reagent_keys[inp.file_index])
                        else:
                            step_input_keys.append(f"mol_{inp.file_index}")
                elif inp.source == InputSource.PREVIOUS_STEP:
                    if inp.step_index is not None and inp.step_index in step_outputs:
                        mol = step_outputs[inp.step_index]
                        step_inputs.append(mol)
                        step_input_keys.append(f"step_{inp.step_index}_product")
                    else:
                        return EnumerationResult(
                            product_name=product_name,
                            error=EnumerationError(
                                step_index=step_idx,
                                pattern_id=None,
                                error_type="invalid_input",
                                message=f"Missing output from step {inp.step_index}",
                                reagent_smiles=[],
                                reagent_names=step_input_keys,
                            ),
                        )

            # Get the primary reaction for deprotection info
            primary_rxn = rxns[0]

            # Apply deprotection if configured
            if primary_rxn.deprotections:
                step_inputs = apply_deprotections(
                    step_inputs, primary_rxn.deprotections, self._protecting_groups
                )

            # Check if this step has alternatives
            step_has_alternatives = (
                self.config.step_modes is not None
                and step_idx in self.config.step_modes
            )

            # Get available patterns for this step
            available_patterns = list(
                self._rdkit_reactions.get(step_idx, {}).keys()
            )

            if not available_patterns:
                return EnumerationResult(
                    product_name=product_name,
                    error=EnumerationError(
                        step_index=step_idx,
                        pattern_id=None,
                        error_type="no_compatible_pattern",
                        message=f"No reactions available for step {step_idx}",
                        reagent_smiles=[Chem.MolToSmiles(m) for m in step_inputs if m],
                        reagent_names=step_input_keys,
                    ),
                )

            # Determine pattern order to try
            if step_has_alternatives:
                # Try cache first, then fall back to trying all patterns
                cached_pattern = self.get_compatible_patterns(step_input_keys, step_idx)
                if cached_pattern and cached_pattern in available_patterns:
                    # Put cached pattern first, then others
                    patterns_to_try = [cached_pattern] + [
                        p for p in available_patterns if p != cached_pattern
                    ]
                else:
                    # No cache hit - try all patterns starting with "primary" if available
                    if "primary" in available_patterns:
                        patterns_to_try = ["primary"] + [
                            p for p in available_patterns if p != "primary"
                        ]
                    else:
                        patterns_to_try = available_patterns
            else:
                # No alternatives - use cache or fallback to primary
                pattern_id = self.get_compatible_patterns(step_input_keys, step_idx)
                if pattern_id is None or pattern_id not in available_patterns:
                    pattern_id = "primary" if "primary" in available_patterns else available_patterns[0]
                patterns_to_try = [pattern_id]

            # Try each pattern until one succeeds
            product = None
            successful_pattern = None
            last_error_smarts = None

            for pattern_id in patterns_to_try:
                rxn = self._rdkit_reactions[step_idx][pattern_id]

                # Get reaction SMARTS for error reporting
                reaction_smarts = None
                for rxn_def in rxns:
                    if rxn_def.pattern_id == pattern_id or (
                        pattern_id == "primary" and rxn_def.pattern_id is None
                    ):
                        reaction_smarts = rxn_def.reaction_smarts
                        break

                # Run reaction
                product, product_smiles = run_single_reaction(rxn, step_inputs)

                if product is not None:
                    successful_pattern = pattern_id
                    break
                else:
                    last_error_smarts = reaction_smarts

            # If no pattern succeeded, return error
            if product is None:
                return EnumerationResult(
                    product_name=product_name,
                    error=EnumerationError(
                        step_index=step_idx,
                        pattern_id=patterns_to_try[-1] if patterns_to_try else None,
                        error_type="reaction_failed",
                        message=f"Reaction produced no products (tried {len(patterns_to_try)} pattern(s): {patterns_to_try})",
                        reagent_smiles=[Chem.MolToSmiles(m) for m in step_inputs if m],
                        reagent_names=step_input_keys,
                        reaction_smarts=last_error_smarts,
                    ),
                )

            pattern_id = successful_pattern

            # Apply product deprotections if configured
            if primary_rxn.deprotections:
                product = apply_product_deprotections(
                    product, primary_rxn.deprotections, self._protecting_groups
                )

            # Store result
            step_outputs[step_idx] = product
            patterns_used[step_idx] = pattern_id

            if intermediates is not None:
                intermediates[step_idx] = product

        # Get final product (from last step)
        final_step_idx = max(step_outputs.keys())
        final_product = step_outputs[final_step_idx]

        return EnumerationResult(
            product=final_product,
            product_smiles=Chem.MolToSmiles(final_product),
            product_name=product_name,
            patterns_used=patterns_used,
            intermediates=intermediates,
        )

    def enumerate_single_from_smiles(
        self,
        smiles_list: List[str],
        reagent_keys: Optional[List[str]] = None,
        store_intermediates: bool = False,
    ) -> EnumerationResult:
        """
        Enumerate a single product from SMILES strings.

        Args:
            smiles_list: List of SMILES strings
            reagent_keys: Optional reagent identifiers
            store_intermediates: If True, store intermediate products

        Returns:
            EnumerationResult
        """
        mols = []
        for i, smiles in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                reagent_name = (
                    reagent_keys[i]
                    if reagent_keys and i < len(reagent_keys)
                    else f"reagent_{i}"
                )
                return EnumerationResult(
                    product_name="_".join(reagent_keys) if reagent_keys else None,
                    error=EnumerationError(
                        step_index=0,
                        pattern_id=None,
                        error_type="invalid_input",
                        message=f"Invalid SMILES at position {i}: '{smiles}'",
                        reagent_smiles=smiles_list,
                        reagent_names=reagent_keys
                        or [f"reagent_{j}" for j in range(len(smiles_list))],
                    ),
                )
            mols.append(mol)

        return self.enumerate_single(mols, reagent_keys, store_intermediates)

    # === Alias Methods (for backward compatibility) ===

    def enumerate(
        self,
        reagent_mols: List[Chem.Mol],
        reagent_keys: Optional[List[str]] = None,
        store_intermediates: bool = False,
    ) -> EnumerationResult:
        """Alias for enumerate_single (backward compatibility)."""
        return self.enumerate_single(reagent_mols, reagent_keys, store_intermediates)

    def enumerate_from_smiles(
        self,
        smiles_list: List[str],
        reagent_keys: Optional[List[str]] = None,
        store_intermediates: bool = False,
    ) -> EnumerationResult:
        """Alias for enumerate_single_from_smiles (backward compatibility)."""
        return self.enumerate_single_from_smiles(
            smiles_list, reagent_keys, store_intermediates
        )

    # === Batch Enumeration ===

    def enumerate_batch(
        self,
        combinations: List[Tuple[List[Chem.Mol], Optional[List[str]]]],
        n_jobs: int = 1,
        show_progress: bool = False,
    ) -> List[EnumerationResult]:
        """
        Enumerate multiple products.

        Args:
            combinations: List of (reagent_mols, reagent_keys) tuples
            n_jobs: Parallel workers (1 = sequential)
            show_progress: Show progress bar

        Returns:
            List of EnumerationResult
        """
        return enumerate_products(self, combinations, n_jobs, show_progress)

    # === Full Library Enumeration ===

    def enumerate_library(
        self,
        n_jobs: int = 1,
        show_progress: bool = True,
    ) -> List[EnumerationResult]:
        """
        Enumerate entire combinatorial library.

        Generates all possible products from reagent files.

        Args:
            n_jobs: Number of parallel workers
            show_progress: Show progress bar

        Returns:
            List of EnumerationResult for all products
        """
        # Generate all combinations with separate reagent names for traceability
        mol_combinations, _, reagent_name_lists = generate_all_combinations(
            self.reagent_file_list,
            return_separate_names=True,
        )

        # Create input format for enumerate_batch: (mols, reagent_keys)
        combinations = [
            (list(mols), reagent_names)
            for mols, reagent_names in zip(mol_combinations, reagent_name_lists)
        ]

        return self.enumerate_batch(combinations, n_jobs, show_progress)

    # === Compatibility Detection ===

    def auto_detect_compatibility(
        self,
        reagent_lists: Optional[List[List[Any]]] = None,
        deprotect: bool = False,
        desalt: bool = False,
    ) -> AutoDetectionResult:
        """
        Automatically detect which reagents work with which patterns.

        This runs validation for each pattern and builds a compatibility
        map used during enumeration for routing.

        Args:
            reagent_lists: Optional list of Reagent lists from Thompson Sampling
            deprotect: Apply deprotection during detection
            desalt: Apply desalting during detection

        Returns:
            AutoDetectionResult with full compatibility information
        """
        result = AutoDetectionResult()

        for step_idx in sorted(self._reactions_by_step.keys()):
            result.pattern_results[step_idx] = {}
            result.coverage_by_pattern[step_idx] = {}

            rxns = self._reactions_by_step[step_idx]
            reagent_files = self._get_step_reagent_files(step_idx)

            for rxn in rxns:
                pattern_id = rxn.pattern_id or "primary"

                validator = _SMARTSValidator(
                    reaction_smarts=rxn.reaction_smarts,
                    reagent_files=reagent_files,
                    deprotect_on_import=deprotect,
                    desalt_on_import=desalt,
                )
                validation_result = validator.validate(test_reactions=False)

                result.pattern_results[step_idx][pattern_id] = validation_result

                # Store coverage by pattern
                result.coverage_by_pattern[step_idx][pattern_id] = dict(
                    validation_result.coverage_stats
                )

                # Cache validator
                if step_idx not in self._validators:
                    self._validators[step_idx] = {}
                self._validators[step_idx][pattern_id] = validator

                # Update compatibility cache
                for position, reagents in validation_result.compatible_reagents.items():
                    for smiles, name in reagents:
                        cache_key = (step_idx, position, name)
                        if cache_key not in self._compatibility_cache:
                            self._compatibility_cache[cache_key] = set()
                        self._compatibility_cache[cache_key].add(pattern_id)

        # Update Reagent objects if provided
        if reagent_lists is not None:
            self._update_reagent_compatibility(reagent_lists)

        return result

    def _get_step_reagent_files(self, step_idx: int) -> List[str]:
        """Get reagent files for a specific step."""
        inputs = self.config.get_inputs_for_step(step_idx)
        reagent_files = []

        for inp in inputs:
            if inp.source == InputSource.REAGENT_FILE:
                if inp.file_index is not None and inp.file_index < len(
                    self.config.reagent_file_list
                ):
                    reagent_files.append(self.config.reagent_file_list[inp.file_index])

        return reagent_files

    def _update_reagent_compatibility(self, reagent_lists: List[List[Any]]) -> None:
        """Update Reagent objects with detected compatibility."""
        for position, reagents in enumerate(reagent_lists):
            for reagent in reagents:
                cache_key = (0, position, reagent.reagent_name)
                if cache_key in self._compatibility_cache:
                    compatible = self._compatibility_cache[cache_key]
                    if hasattr(reagent, "set_compatible_smarts"):
                        reagent.set_compatible_smarts(compatible)

    def register_compatibility(
        self,
        position: int,
        reagent_key: str,
        compatible_patterns: Set[str],
        step_index: int = 0,
    ) -> None:
        """
        Manually register reagent-pattern compatibility.

        Args:
            position: Reagent position
            reagent_key: Unique reagent identifier
            compatible_patterns: Set of compatible pattern IDs
            step_index: Which step this applies to
        """
        cache_key = (step_index, position, reagent_key)
        self._compatibility_cache[cache_key] = compatible_patterns

    def get_compatible_patterns(
        self,
        reagent_keys: List[str],
        step_index: int = 0,
    ) -> Optional[str]:
        """
        Find a pattern compatible with all given reagents.

        Args:
            reagent_keys: Reagent identifiers for each position
            step_index: Which step

        Returns:
            First compatible pattern_id, or None if no pattern works
        """
        if not reagent_keys:
            return None

        # Get patterns available at this step
        rxns = self._reactions_by_step.get(step_index, [])
        if not rxns:
            return None

        available_patterns = [r.pattern_id or "primary" for r in rxns]

        # Find intersection of compatible patterns
        compatible_patterns: Optional[Set[str]] = None

        for position, key in enumerate(reagent_keys):
            cache_key = (step_index, position, key)
            if cache_key not in self._compatibility_cache:
                patterns = {"primary"}
            else:
                patterns = self._compatibility_cache[cache_key]

            if compatible_patterns is None:
                compatible_patterns = patterns.copy()
            else:
                compatible_patterns &= patterns

        if compatible_patterns is None or not compatible_patterns:
            return None

        # Return first available pattern in priority order
        for pattern_id in available_patterns:
            if pattern_id in compatible_patterns:
                return pattern_id

        return None

    def get_compatibility_map(self) -> Dict[Tuple[int, int, str], Set[str]]:
        """Get the full compatibility cache."""
        return self._compatibility_cache.copy()

    # === Validation ===

    def validate_all(
        self,
        deprotect: bool = False,
        desalt: bool = False,
    ) -> Dict[int, Dict[str, ValidationResult]]:
        """
        Validate all reactions against reagent files.

        Args:
            deprotect: Apply deprotection during validation
            desalt: Apply desalting during validation

        Returns:
            Nested dict: {step_index: {pattern_id: ValidationResult}}
        """
        result = self.auto_detect_compatibility(deprotect=deprotect, desalt=desalt)
        return result.pattern_results

    def get_validator(
        self,
        step_index: int,
        pattern_id: str = "primary",
    ) -> Optional[_SMARTSValidator]:
        """
        Get validator for a specific step/pattern.

        Args:
            step_index: Step index
            pattern_id: Pattern identifier

        Returns:
            _SMARTSValidator or None if not validated
        """
        return self._validators.get(step_index, {}).get(pattern_id)

    # === Serialization for Multiprocessing ===

    def prepare_worker_data(self) -> Dict:
        """
        Prepare serializable data for multiprocessing workers.

        Returns:
            Dict that workers can use to reconstruct the pipeline
        """
        config_dict = self.config.model_dump()

        # Convert compatibility cache to serializable format
        compat_cache = {
            f"{s}_{p}_{k}": list(v)
            for (s, p, k), v in self._compatibility_cache.items()
        }

        return {
            "config_dict": config_dict,
            "compatibility_cache": compat_cache,
        }

    @classmethod
    def from_worker_data(cls, data: Dict) -> "SynthesisPipeline":
        """
        Reconstruct pipeline in a worker process.

        Args:
            data: Dict from prepare_worker_data()

        Returns:
            Reconstructed SynthesisPipeline
        """
        config = ReactionConfig(**data["config_dict"])
        pipeline = cls(config)

        # Restore compatibility cache
        compat_cache = data["compatibility_cache"]
        for key_str, patterns in compat_cache.items():
            parts = key_str.split("_", 2)
            if len(parts) >= 3:
                step_idx = int(parts[0])
                position = int(parts[1])
                reagent_key = parts[2]
                pipeline._compatibility_cache[(step_idx, position, reagent_key)] = set(
                    patterns
                )

        return pipeline

    def __repr__(self) -> str:
        return (
            f"SynthesisPipeline(steps={self.num_steps}, "
            f"components={self.num_components}, "
            f"has_alternatives={self.has_alternatives})"
        )
