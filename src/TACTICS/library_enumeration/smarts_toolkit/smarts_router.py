"""
SMARTS Router for handling multiple reaction SMARTS patterns.

This module provides the foundation for routing reagent combinations to
appropriate SMARTS patterns based on compatibility information.
"""

from pydantic import BaseModel, field_validator
from typing import List, Dict, Optional, Tuple, Set
from rdkit import Chem
from rdkit.Chem import AllChem
import logging


class SMARTSPatternConfig(BaseModel):
    """
    Configuration for a single SMARTS pattern.

    Attributes:
        pattern_id: Unique identifier (e.g., "primary", "alt_1", "secondary_amine")
        reaction_smarts: The reaction SMARTS string
        description: Optional human-readable description
        component_compatibility: Optional dict mapping component_idx to list of
                                 substructure SMARTS that define compatibility.
                                 If None, compatibility is determined by reagent attributes.

    Example:
        SMARTSPatternConfig(
            pattern_id="secondary_amine",
            reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]",
            description="For secondary amines that lack a free NH2"
        )
    """
    pattern_id: str
    reaction_smarts: str
    description: Optional[str] = None
    component_compatibility: Optional[Dict[int, List[str]]] = None

    @field_validator('reaction_smarts')
    @classmethod
    def validate_smarts(cls, v: str) -> str:
        """Validate that the SMARTS string is parseable."""
        rxn = AllChem.ReactionFromSmarts(v)
        if rxn is None:
            raise ValueError(f"Invalid reaction SMARTS: {v}")
        return v


class AlternativeSMARTSConfig(BaseModel):
    """
    Configuration for libraries with multiple SMARTS patterns.

    Supports two modes:
    1. Single SMARTS (backwards compatible): Just provide primary_smarts
    2. Multiple SMARTS: Provide primary_smarts and alternatives list

    Attributes:
        primary_smarts: The default/primary reaction SMARTS
        alternatives: Optional list of alternative SMARTS patterns

    Example (single):
        AlternativeSMARTSConfig(
            primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"
        )

    Example (multiple):
        AlternativeSMARTSConfig(
            primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            alternatives=[
                SMARTSPatternConfig(
                    pattern_id="secondary_amine",
                    reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]"
                )
            ]
        )
    """
    primary_smarts: str
    alternatives: Optional[List[SMARTSPatternConfig]] = None

    @field_validator('primary_smarts')
    @classmethod
    def validate_primary_smarts(cls, v: str) -> str:
        """Validate primary SMARTS."""
        rxn = AllChem.ReactionFromSmarts(v)
        if rxn is None:
            raise ValueError(f"Invalid primary reaction SMARTS: {v}")
        return v

    def get_all_patterns(self) -> List[SMARTSPatternConfig]:
        """
        Get all patterns with primary first (priority order).

        Returns:
            List of SMARTSPatternConfig, primary pattern first
        """
        patterns = [
            SMARTSPatternConfig(
                pattern_id="primary",
                reaction_smarts=self.primary_smarts,
                description="Primary reaction SMARTS"
            )
        ]
        if self.alternatives:
            patterns.extend(self.alternatives)
        return patterns

    def has_alternatives(self) -> bool:
        """Check if alternative SMARTS patterns are defined."""
        return self.alternatives is not None and len(self.alternatives) > 0


class SMARTSRouter:
    """
    Routes reagent combinations to appropriate SMARTS patterns.

    This class is responsible for:
    1. Managing multiple SMARTS patterns for a single reaction step
    2. Determining which SMARTS pattern to use based on reagent compatibility
    3. Executing reactions with the selected SMARTS
    4. Handling enumeration failures gracefully

    The router serves as the foundation for both:
    - Single-step reactions with alternative SMARTS
    - Individual steps within multi-step syntheses

    Attributes:
        config: AlternativeSMARTSConfig defining available patterns
        logger: Logger instance for debugging

    Usage:
        # Initialize router
        config = AlternativeSMARTSConfig(
            primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            alternatives=[
                SMARTSPatternConfig(
                    pattern_id="alt_1",
                    reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]"
                )
            ]
        )
        router = SMARTSRouter(config)

        # Register reagent compatibility (done during initialization)
        router.register_reagent_compatibility("reagent_A", {"primary", "alt_1"})
        router.register_reagent_compatibility("reagent_B", {"primary"})
        router.register_reagent_compatibility("reagent_C", {"alt_1"})

        # Enumerate products
        product, pattern_used = router.enumerate(reagent_mols, reagent_keys)
    """

    def __init__(self, config: AlternativeSMARTSConfig):
        """
        Initialize the SMARTS router.

        Args:
            config: Configuration containing primary and alternative SMARTS patterns
        """
        self.config = config
        self.logger = logging.getLogger(__name__)

        # Parse all patterns
        self._patterns: List[SMARTSPatternConfig] = config.get_all_patterns()

        # Pre-compiled reactions: pattern_id -> ChemicalReaction
        self._reactions: Dict[str, AllChem.ChemicalReaction] = {}

        # Reagent compatibility cache: reagent_key -> Set[pattern_id]
        self._compatibility_cache: Dict[str, Set[str]] = {}

        # Initialize reactions
        self._initialize_reactions()

    def _initialize_reactions(self) -> None:
        """Parse and validate all reaction SMARTS patterns."""
        for pattern in self._patterns:
            rxn = AllChem.ReactionFromSmarts(pattern.reaction_smarts)
            if rxn is None:
                raise ValueError(
                    f"Failed to parse reaction SMARTS for pattern "
                    f"'{pattern.pattern_id}': {pattern.reaction_smarts}"
                )
            self._reactions[pattern.pattern_id] = rxn
            self.logger.debug(
                f"Initialized reaction pattern '{pattern.pattern_id}'"
            )

    def register_reagent_compatibility(
        self,
        reagent_key: str,
        compatible_patterns: Set[str]
    ) -> None:
        """
        Register which SMARTS patterns a reagent is compatible with.

        This should be called during initialization after SMARTS toolkit
        validation has determined compatibility.

        Args:
            reagent_key: Unique identifier for the reagent (name or SMILES)
            compatible_patterns: Set of pattern_ids this reagent works with

        Raises:
            ValueError: If any pattern_id is not recognized
        """
        valid_pattern_ids = {p.pattern_id for p in self._patterns}
        invalid = compatible_patterns - valid_pattern_ids
        if invalid:
            raise ValueError(
                f"Unknown pattern IDs for reagent '{reagent_key}': {invalid}. "
                f"Valid patterns: {valid_pattern_ids}"
            )

        self._compatibility_cache[reagent_key] = compatible_patterns
        self.logger.debug(
            f"Registered reagent '{reagent_key}' with patterns: {compatible_patterns}"
        )

    def find_compatible_smarts(
        self,
        reagent_keys: List[str]
    ) -> Optional[str]:
        """
        Find a SMARTS pattern compatible with all given reagents.

        Uses the compatibility cache to find the intersection of compatible
        patterns across all reagents. Returns the first match in priority
        order (primary first).

        Args:
            reagent_keys: List of reagent identifiers, one per reaction component

        Returns:
            pattern_id of first compatible SMARTS in priority order,
            or None if no pattern is compatible with all reagents
        """
        # If no compatibility info registered, default to primary
        if not self._compatibility_cache:
            self.logger.debug("No compatibility cache, defaulting to primary")
            return "primary"

        # Gather compatible pattern sets for each reagent
        compatible_sets: List[Set[str]] = []
        for key in reagent_keys:
            if key in self._compatibility_cache:
                compatible_sets.append(self._compatibility_cache[key])
            else:
                # Unknown reagent - assume compatible with primary only
                self.logger.warning(
                    f"Reagent '{key}' not in compatibility cache, "
                    f"assuming primary only"
                )
                compatible_sets.append({"primary"})

        # Find intersection of all sets
        if not compatible_sets:
            return "primary"

        common_patterns = set.intersection(*compatible_sets)

        if not common_patterns:
            self.logger.debug(
                f"No common patterns for reagents {reagent_keys}"
            )
            return None

        # Return first match in priority order
        for pattern in self._patterns:
            if pattern.pattern_id in common_patterns:
                return pattern.pattern_id

        return None

    def enumerate(
        self,
        reagent_mols: List[Chem.Mol],
        reagent_keys: Optional[List[str]] = None,
        pattern_id: Optional[str] = None
    ) -> Tuple[Optional[Chem.Mol], Optional[str]]:
        """
        Enumerate product from reagents using appropriate SMARTS.

        Determines the correct SMARTS pattern (via routing or explicit
        specification), runs the reaction, and returns the product.

        Args:
            reagent_mols: List of RDKit Mol objects, one per reaction component
            reagent_keys: Optional list of reagent identifiers for SMARTS selection.
                          If None and pattern_id is None, uses primary pattern.
            pattern_id: Optional explicit pattern to use (bypasses routing)

        Returns:
            Tuple of (product_mol, pattern_id_used):
            - product_mol: RDKit Mol of product, or None if enumeration failed
            - pattern_id_used: The pattern that was used, or None if routing failed
        """
        # Determine which SMARTS to use
        if pattern_id is None:
            if reagent_keys is not None:
                pattern_id = self.find_compatible_smarts(reagent_keys)
                if pattern_id is None:
                    self.logger.debug(
                        f"No compatible SMARTS found for {reagent_keys}"
                    )
                    return None, None
            else:
                pattern_id = "primary"

        # Get the reaction
        rxn = self._reactions.get(pattern_id)
        if rxn is None:
            self.logger.error(f"No reaction found for pattern '{pattern_id}'")
            return None, None

        # Run reaction
        try:
            products = rxn.RunReactants(tuple(reagent_mols))
            if products and len(products) > 0 and len(products[0]) > 0:
                prod_mol = products[0][0]
                Chem.SanitizeMol(prod_mol)
                return prod_mol, pattern_id
            else:
                self.logger.debug(
                    f"No products from pattern '{pattern_id}'"
                )
                return None, pattern_id
        except Exception as e:
            self.logger.debug(
                f"Enumeration failed with pattern '{pattern_id}': {e}"
            )
            return None, pattern_id

    def get_reaction(
        self,
        pattern_id: str = "primary"
    ) -> Optional[AllChem.ChemicalReaction]:
        """
        Get the reaction object for a specific pattern.

        Useful for backwards compatibility or direct reaction access.

        Args:
            pattern_id: The pattern identifier

        Returns:
            The ChemicalReaction object, or None if pattern not found
        """
        return self._reactions.get(pattern_id)

    @property
    def pattern_ids(self) -> List[str]:
        """Get list of all pattern IDs in priority order."""
        return [p.pattern_id for p in self._patterns]

    @property
    def num_patterns(self) -> int:
        """Get the number of available SMARTS patterns."""
        return len(self._patterns)

    def get_compatibility_summary(self) -> Dict[str, int]:
        """
        Get summary of reagent compatibility counts per pattern.

        Returns:
            Dict mapping pattern_id -> number of compatible reagents
        """
        summary = {p.pattern_id: 0 for p in self._patterns}
        for patterns in self._compatibility_cache.values():
            for pid in patterns:
                if pid in summary:
                    summary[pid] += 1
        return summary
