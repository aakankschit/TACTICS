"""
Internal SMARTS pattern validation module.

This module provides the internal _SMARTSValidator class used by ReactionDef
for validation. Users should interact with validation through ReactionDef.validate_reaction()
rather than using this class directly.

The ValidationResult dataclass is public and returned by validation methods.
"""

import logging
import random
from dataclasses import dataclass, field
from itertools import product as itertools_product
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import polars as pl
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

from .config import ProtectingGroupInfo
from .constants import DEFAULT_PROTECTING_GROUPS, DEFAULT_SALT_FRAGMENTS

logger = logging.getLogger(__name__)


@dataclass
class ValidationResult:
    """
    Comprehensive results from SMARTS validation.

    Attributes:
        compatible_reagents: {position: [(smiles, name), ...]} - Reagents that match template
        incompatible_reagents: {position: [(smiles, name), ...]} - Reagents that don't match
        invalid_smiles: {position: [(smiles, name), ...]} - Unparseable SMILES
        duplicate_smiles: {position: [(smiles, name), ...]} - Duplicate entries
        protected_reagents: {position: [(smiles, name, [groups]), ...]} - With protecting groups
        multi_fragment_reagents: {position: [(smiles, name, [fragments]), ...]} - With salts
        coverage_stats: {position: float} - Percent compatible per position (0-100)
        reaction_success_rate: float - Percent of test reactions that succeeded
        error_messages: List[str] - Critical errors
        warnings: List[str] - Non-critical warnings
    """

    compatible_reagents: Dict[int, List[Tuple[str, str]]] = field(default_factory=dict)
    incompatible_reagents: Dict[int, List[Tuple[str, str]]] = field(
        default_factory=dict
    )
    invalid_smiles: Dict[int, List[Tuple[str, str]]] = field(default_factory=dict)
    duplicate_smiles: Dict[int, List[Tuple[str, str]]] = field(default_factory=dict)
    protected_reagents: Dict[int, List[Tuple[str, str, List[str]]]] = field(
        default_factory=dict
    )
    multi_fragment_reagents: Dict[int, List[Tuple[str, str, List[str]]]] = field(
        default_factory=dict
    )
    coverage_stats: Dict[int, float] = field(default_factory=dict)
    reaction_success_rate: float = 0.0
    error_messages: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    def is_valid(self) -> bool:
        """True if all positions have >0% coverage and no critical errors."""
        if self.error_messages:
            return False
        return all(cov > 0 for cov in self.coverage_stats.values())

    @property
    def total_compatible(self) -> int:
        """Total number of compatible reagents across all positions."""
        return sum(len(v) for v in self.compatible_reagents.values())

    @property
    def total_incompatible(self) -> int:
        """Total number of incompatible reagents across all positions."""
        return sum(len(v) for v in self.incompatible_reagents.values())


class _SMARTSValidator:
    """
    Internal validator class used by ReactionDef.

    Users should use ReactionDef.validate_reaction() instead of this class directly.
    """

    def __init__(
        self,
        reaction_smarts: str,
        reagent_files: Optional[List[str]] = None,
        reagent_smiles_lists: Optional[List[List[Tuple[str, str]]]] = None,
        reagent_dataframes: Optional[List[pl.DataFrame]] = None,
        protecting_groups: Optional[List[ProtectingGroupInfo]] = None,
        salt_fragments: Optional[List[Tuple[str, str]]] = None,
        deprotect_on_import: bool = False,
        deprotect_groups: Optional[List[str]] = None,
        desalt_on_import: bool = False,
    ):
        """
        Initialize validator.

        Args:
            reaction_smarts: Reaction SMARTS pattern to validate against
            reagent_files: List of paths to reagent .smi files
            reagent_smiles_lists: Direct SMILES lists as [(smiles, name), ...] per position
                                  (useful for validating against intermediates)
            reagent_dataframes: Alternative: pre-loaded DataFrames with 'smiles', 'name' columns
            protecting_groups: Custom protecting groups (replaces defaults if provided)
            salt_fragments: Custom salt fragments (replaces defaults if provided)
            deprotect_on_import: If True, remove protecting groups when loading reagents
            deprotect_groups: Specific groups to remove (None = all)
            desalt_on_import: If True, remove salts when loading reagents
        """
        self.reaction_smarts = reaction_smarts
        self.reagent_files = reagent_files or []

        # Import-time processing options
        self.deprotect_on_import = deprotect_on_import
        self.deprotect_groups = deprotect_groups
        self.desalt_on_import = desalt_on_import

        # Protecting groups configuration
        self.protecting_groups: List[ProtectingGroupInfo] = (
            protecting_groups
            if protecting_groups is not None
            else list(DEFAULT_PROTECTING_GROUPS)
        )

        # Salt fragments configuration
        self.salt_fragments: List[Tuple[str, str]] = (
            salt_fragments
            if salt_fragments is not None
            else list(DEFAULT_SALT_FRAGMENTS)
        )

        # Internal state
        self.reaction: Optional[Any] = None  # RDKit ChemicalReaction
        self.reagents: List[List[Tuple[str, str]]] = []
        self.errors: List[str] = []
        self.warnings: List[str] = []

        # Preprocessing results
        self.has_original_names: Dict[int, bool] = {}
        self.invalid_smiles: Dict[int, List[Tuple[str, str]]] = {}
        self.duplicate_smiles: Dict[int, List[Tuple[str, str]]] = {}
        self.protected_reagents: Dict[int, List[Tuple[str, str, List[str]]]] = {}
        self.multi_fragment_reagents: Dict[int, List[Tuple[str, str, List[str]]]] = {}
        self.deprotected_on_import: Dict[
            int, List[Tuple[str, str, str, List[str]]]
        ] = {}
        self.desalted_on_import: Dict[int, List[Tuple[str, str, str, List[str]]]] = {}

        # Cached validation result
        self._cached_result: Optional[ValidationResult] = None

        # Load reaction
        self._load_reaction()

        # Load reagents from various sources (priority order)
        if reagent_smiles_lists is not None:
            self._load_reagents_from_smiles_lists(reagent_smiles_lists)
        elif reagent_dataframes is not None:
            self._load_reagents_from_dataframes(reagent_dataframes)
        elif reagent_files:
            self._load_reagents()

    def _load_reaction(self) -> None:
        """Load and validate the reaction SMARTS."""
        try:
            self.reaction = AllChem.ReactionFromSmarts(self.reaction_smarts)
            if self.reaction is None:
                raise ValueError("Invalid reaction SMARTS")

            num_reactants = self.reaction.GetNumReactantTemplates()
            num_products = self.reaction.GetNumProductTemplates()

            if num_reactants == 0:
                self.errors.append("No reactant templates found in reaction SMARTS")
            if num_products == 0:
                self.errors.append("No product templates found in reaction SMARTS")

            logger.debug(
                f"Loaded reaction with {num_reactants} reactants and {num_products} products"
            )

        except Exception as e:
            self.errors.append(f"Failed to parse reaction SMARTS: {str(e)}")
            logger.error(f"Reaction loading failed: {e}")

    def _load_reagents_from_smiles_lists(
        self, smiles_lists: List[List[Tuple[str, str]]]
    ) -> None:
        """Load reagents from direct SMILES lists."""
        self.reagents = []

        for i, smiles_list in enumerate(smiles_lists):
            self.has_original_names[i] = True
            self.reagents.append(list(smiles_list))
            logger.debug(f"Loaded {len(smiles_list)} reagents from SMILES list {i}")

        self._preprocess_reagents()

    def _load_reagents(self) -> None:
        """Load reagents from files and preprocess them."""
        self.reagents = []

        for i, filepath in enumerate(self.reagent_files):
            path = Path(filepath)
            if not path.exists():
                self.errors.append(f"Reagent file {filepath} not found")
                self.reagents.append([])
                self.has_original_names[i] = False
                continue

            reagent_list = []
            has_names = False

            if path.suffix == ".smi":
                reagent_list, has_names = self._load_smi_file(path)
            elif path.suffix == ".csv":
                reagent_list, has_names = self._load_csv_file(path)
            else:
                self.warnings.append(f"Unknown file format: {path.suffix}")

            self.has_original_names[i] = has_names
            self.reagents.append(reagent_list)
            logger.debug(f"Loaded {len(reagent_list)} reagents from {filepath}")

        self._preprocess_reagents()

    def _load_reagents_from_dataframes(self, dataframes: List[pl.DataFrame]) -> None:
        """Load reagents from pre-loaded Polars DataFrames."""
        self.reagents = []

        for i, df in enumerate(dataframes):
            # Find SMILES and name columns
            smiles_col = None
            name_col = None

            for col in df.columns:
                col_upper = col.upper()
                if "SMILES" in col_upper:
                    smiles_col = col
                elif "NAME" in col_upper or "ID" in col_upper:
                    name_col = col

            if smiles_col is None:
                self.errors.append(f"No SMILES column found in DataFrame {i}")
                self.reagents.append([])
                self.has_original_names[i] = False
                continue

            has_names = name_col is not None
            self.has_original_names[i] = has_names

            reagent_list = []
            for row in df.iter_rows(named=True):
                smiles = row[smiles_col]
                name = row[name_col] if name_col else f"R{len(reagent_list)}"
                reagent_list.append((smiles, name))

            self.reagents.append(reagent_list)
            logger.debug(f"Loaded {len(reagent_list)} reagents from DataFrame {i}")

        self._preprocess_reagents()

    def _load_smi_file(self, path: Path) -> Tuple[List[Tuple[str, str]], bool]:
        """Load SMILES from .smi file."""
        reagents = []
        has_original_names = None

        with open(path, "r") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 1:
                    smiles = parts[0]
                    if len(parts) > 1:
                        name = parts[1]
                        if has_original_names is None:
                            has_original_names = True
                    else:
                        name = f"R{len(reagents)}"
                        if has_original_names is None:
                            has_original_names = False
                    reagents.append((smiles, name))

        return reagents, has_original_names if has_original_names is not None else False

    def _load_csv_file(self, path: Path) -> Tuple[List[Tuple[str, str]], bool]:
        """Load SMILES from CSV file."""
        import pandas as pd

        df = pd.read_csv(path)

        smiles_col = None
        name_col = None

        for col in df.columns:
            if "SMILES" in col.upper():
                smiles_col = col
            elif "NAME" in col.upper() or "ID" in col.upper():
                name_col = col

        if smiles_col is None:
            self.errors.append(f"No SMILES column found in {path}")
            return [], False

        has_original_names = name_col is not None

        reagents = []
        for idx, row in df.iterrows():
            smiles = row[smiles_col]
            name = row[name_col] if name_col else f"R{idx}"
            reagents.append((smiles, name))

        return reagents, has_original_names

    def _preprocess_reagents(self) -> None:
        """Validate SMILES, apply transformations, and remove duplicates."""
        for position, reagent_list in enumerate(self.reagents):
            valid_reagents = []
            seen_smiles: set = set()
            invalid_list: List[Tuple[str, str]] = []
            duplicate_list: List[Tuple[str, str]] = []
            desalted_list: List[Tuple[str, str, str, List[str]]] = []
            deprotected_list: List[Tuple[str, str, str, List[str]]] = []

            for smiles, name in reagent_list:
                original_smiles = smiles

                # Check for valid SMILES
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    invalid_list.append((smiles, name))
                    self.warnings.append(
                        f"Invalid SMILES at position {position}: {smiles} ({name})"
                    )
                    continue

                # Apply desalting if requested
                if self.desalt_on_import and "." in smiles:
                    desalted_smiles, removed_salts = self.desalt_smiles(smiles)
                    if removed_salts:
                        desalted_list.append(
                            (smiles, desalted_smiles, name, removed_salts)
                        )
                        smiles = desalted_smiles
                        logger.debug(
                            f"Position {position}: Desalted {name}: {original_smiles} -> {smiles}"
                        )

                # Apply deprotection if requested
                if self.deprotect_on_import:
                    deprotected_smiles, removed_groups = self.deprotect_smiles(
                        smiles, self.deprotect_groups
                    )
                    if removed_groups:
                        deprotected_list.append(
                            (smiles, deprotected_smiles, name, removed_groups)
                        )
                        smiles = deprotected_smiles
                        logger.debug(
                            f"Position {position}: Deprotected {name}: removed {removed_groups}"
                        )

                # Check for duplicates
                if smiles in seen_smiles:
                    duplicate_list.append((smiles, name))
                    continue

                seen_smiles.add(smiles)
                valid_reagents.append((smiles, name))

            # Store preprocessing results
            self.invalid_smiles[position] = invalid_list
            self.duplicate_smiles[position] = duplicate_list
            self.desalted_on_import[position] = desalted_list
            self.deprotected_on_import[position] = deprotected_list
            self.reagents[position] = valid_reagents

            # Log summary
            if invalid_list:
                logger.warning(
                    f"Position {position}: Removed {len(invalid_list)} invalid SMILES"
                )
            if duplicate_list:
                logger.debug(
                    f"Position {position}: Removed {len(duplicate_list)} duplicate SMILES"
                )

        # Detect protecting groups and multi-fragment reagents
        self._detect_protecting_groups()
        self._detect_multi_fragment_reagents()

    def _detect_protecting_groups(self) -> None:
        """Detect protecting groups in all reagents."""
        for position, reagent_list in enumerate(self.reagents):
            protected_list: List[Tuple[str, str, List[str]]] = []

            for smiles, name in reagent_list:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue

                found_groups: List[str] = []
                for pg in self.protecting_groups:
                    pattern = Chem.MolFromSmarts(pg.smarts)
                    if pattern and mol.HasSubstructMatch(pattern):
                        found_groups.append(pg.name)

                if found_groups:
                    protected_list.append((smiles, name, found_groups))

            self.protected_reagents[position] = protected_list

            if protected_list:
                logger.debug(
                    f"Position {position}: Found {len(protected_list)} reagents "
                    f"with protecting groups"
                )

    def _detect_multi_fragment_reagents(self) -> None:
        """Detect reagents containing multiple fragments (salts, counterions)."""
        # Build canonical lookup for salt fragments
        salt_canonical_lookup: Dict[str, str] = {}
        for salt_smiles, salt_name in self.salt_fragments:
            mol = Chem.MolFromSmiles(salt_smiles)
            if mol:
                canonical = Chem.MolToSmiles(mol)
                salt_canonical_lookup[canonical] = salt_name

        for position, reagent_list in enumerate(self.reagents):
            multi_frag_list: List[Tuple[str, str, List[str]]] = []

            for smiles, name in reagent_list:
                if "." not in smiles:
                    continue

                fragments = smiles.split(".")
                detected_salts: List[str] = []

                for frag in fragments:
                    frag_mol = Chem.MolFromSmiles(frag)
                    if frag_mol:
                        frag_canonical = Chem.MolToSmiles(frag_mol)
                        if frag_canonical in salt_canonical_lookup:
                            detected_salts.append(salt_canonical_lookup[frag_canonical])

                multi_frag_list.append((smiles, name, detected_salts))

            self.multi_fragment_reagents[position] = multi_frag_list

            if multi_frag_list:
                logger.debug(
                    f"Position {position}: Found {len(multi_frag_list)} reagents "
                    f"with multiple fragments"
                )

    # ===== Core Validation =====

    def validate(
        self, test_reactions: bool = False, sample_size: Optional[int] = None
    ) -> ValidationResult:
        """
        Run full validation of reagents against reaction template.

        Args:
            test_reactions: If True, test actual reaction execution on sample combinations
            sample_size: Number of random combinations to test (default: 100)

        Returns:
            ValidationResult with comprehensive compatibility information
        """
        if self.errors:
            return ValidationResult(
                compatible_reagents={},
                incompatible_reagents={},
                invalid_smiles=self.invalid_smiles,
                duplicate_smiles=self.duplicate_smiles,
                protected_reagents=self.protected_reagents,
                multi_fragment_reagents=self.multi_fragment_reagents,
                coverage_stats={},
                reaction_success_rate=0.0,
                error_messages=self.errors,
                warnings=self.warnings,
            )

        # Check template matching
        compatible, incompatible = self._check_template_matching()

        # Calculate coverage statistics
        coverage = self._calculate_coverage(compatible)

        # Test reactions if requested
        success_rate = 0.0
        if test_reactions:
            success_rate = self._test_reaction_combinations(compatible, sample_size)

        result = ValidationResult(
            compatible_reagents=compatible,
            incompatible_reagents=incompatible,
            invalid_smiles=self.invalid_smiles,
            duplicate_smiles=self.duplicate_smiles,
            protected_reagents=self.protected_reagents,
            multi_fragment_reagents=self.multi_fragment_reagents,
            coverage_stats=coverage,
            reaction_success_rate=success_rate,
            error_messages=self.errors,
            warnings=self.warnings,
        )

        self._cached_result = result
        return result

    def _check_template_matching(self) -> Tuple[Dict, Dict]:
        """Check which reagents match reaction templates."""
        compatible: Dict[int, List[Tuple[str, str]]] = {}
        incompatible: Dict[int, List[Tuple[str, str]]] = {}

        if self.reaction is None:
            return compatible, incompatible

        num_templates = self.reaction.GetNumReactantTemplates()

        for i in range(num_templates):
            template = self.reaction.GetReactantTemplate(i)
            compatible[i] = []
            incompatible[i] = []

            if i >= len(self.reagents):
                self.warnings.append(f"No reagents provided for position {i}")
                continue

            for smiles, name in self.reagents[i]:
                mol = Chem.MolFromSmiles(smiles)
                if mol and mol.HasSubstructMatch(template):
                    compatible[i].append((smiles, name))
                else:
                    incompatible[i].append((smiles, name))

        return compatible, incompatible

    def _calculate_coverage(self, compatible: Dict) -> Dict[int, float]:
        """Calculate coverage statistics."""
        coverage = {}

        for i, reagent_list in enumerate(self.reagents):
            if i in compatible:
                num_compatible = len(compatible[i])
                total = len(reagent_list)
                coverage[i] = (num_compatible / total * 100) if total > 0 else 0
            else:
                coverage[i] = 0.0

        return coverage

    def _test_reaction_combinations(
        self, compatible: Dict, sample_size: Optional[int] = None
    ) -> float:
        """Test actual reaction combinations."""
        if self.reaction is None:
            return 0.0

        if all(i in compatible for i in range(len(self.reagents))):
            reagent_lists = [compatible[i] for i in range(len(self.reagents))]

            # Calculate total combinations
            total_combinations = 1
            for rlist in reagent_lists:
                total_combinations *= len(rlist)

            # Sample if needed
            if sample_size and sample_size < total_combinations:
                combinations = []
                for _ in range(sample_size):
                    combo = [random.choice(rlist) for rlist in reagent_lists]
                    combinations.append(combo)
            else:
                combinations = list(itertools_product(*reagent_lists))

            # Test reactions
            successful = 0
            total = len(combinations)

            for combo in combinations:
                reagent_mols = []
                for smiles, name in combo:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        reagent_mols.append(mol)

                if len(reagent_mols) == len(combo):
                    products = self.reaction.RunReactants(tuple(reagent_mols))
                    if products:
                        successful += 1

            return (successful / total * 100) if total > 0 else 0

        return 0.0

    def get_compatible_reagents(
        self, position: int, deprotect: bool = False, desalt: bool = False
    ) -> List[Tuple[str, str]]:
        """
        Get reagents compatible with template at a position.

        Args:
            position: Reagent position (0-indexed)
            deprotect: If True, return deprotected SMILES
            desalt: If True, return desalted SMILES

        Returns:
            List of (smiles, name) tuples
        """
        if self._cached_result is None:
            self.validate(test_reactions=False)

        if self._cached_result is None:
            return []

        compatible = self._cached_result.compatible_reagents.get(position, [])

        if not deprotect and not desalt:
            return compatible

        # Apply transformations
        result = []
        for smiles, name in compatible:
            output_smiles = smiles

            if desalt and "." in output_smiles:
                output_smiles, _ = self.desalt_smiles(output_smiles)

            if deprotect:
                output_smiles, _ = self.deprotect_smiles(output_smiles)

            result.append((output_smiles, name))

        return result

    def get_incompatible_reagents(self, position: int) -> List[Tuple[str, str]]:
        """Get reagents that don't match template at a position."""
        if self._cached_result is None:
            self.validate(test_reactions=False)

        if self._cached_result is None:
            return []

        return self._cached_result.incompatible_reagents.get(position, [])

    def get_reactant_template(self, position: int) -> Optional[Chem.Mol]:
        """
        Get the reactant template molecule for a position.

        Returns:
            RDKit Mol object representing the reactant pattern
        """
        if self.reaction is None:
            return None

        if position >= self.reaction.GetNumReactantTemplates():
            return None

        return self.reaction.GetReactantTemplate(position)

    # ===== Deprotection / Desalting =====

    def deprotect_smiles(
        self, smiles: str, groups_to_remove: Optional[List[str]] = None
    ) -> Tuple[str, List[str]]:
        """
        Remove protecting groups from a SMILES.

        Args:
            smiles: Input SMILES
            groups_to_remove: Specific groups to remove (default: all detected)

        Returns:
            (cleaned_smiles, list_of_removed_group_names)
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles, []

        removed_groups: List[str] = []

        for pg in self.protecting_groups:
            # Skip if not in the list to remove
            if groups_to_remove is not None and pg.name not in groups_to_remove:
                continue

            # Skip if no deprotection SMARTS available
            if pg.deprotection_smarts is None:
                continue

            # Check if this protecting group is present
            pattern = Chem.MolFromSmarts(pg.smarts)
            if pattern is None or not mol.HasSubstructMatch(pattern):
                continue

            # Apply deprotection reaction
            rxn = AllChem.ReactionFromSmarts(pg.deprotection_smarts)
            if rxn is None:
                continue

            # Run reaction (may need multiple iterations)
            max_iterations = 10
            for _ in range(max_iterations):
                if not mol.HasSubstructMatch(pattern):
                    break

                products = rxn.RunReactants((mol,))
                if products and len(products) > 0 and len(products[0]) > 0:
                    mol = products[0][0]
                    try:
                        Chem.SanitizeMol(mol)
                    except Exception:
                        break
                    if pg.name not in removed_groups:
                        removed_groups.append(pg.name)
                else:
                    break

        return Chem.MolToSmiles(mol), removed_groups

    def desalt_smiles(self, smiles: str) -> Tuple[str, List[str]]:
        """
        Remove salt fragments from a multi-fragment SMILES.

        Returns the largest organic fragment.

        Args:
            smiles: Input SMILES (may contain "." separators)

        Returns:
            (cleaned_smiles, list_of_removed_fragments)
        """
        if "." not in smiles:
            return smiles, []

        # Build canonical lookup for salt fragments
        salt_canonical_lookup: Dict[str, str] = {}
        for salt_smiles, salt_name in self.salt_fragments:
            mol = Chem.MolFromSmiles(salt_smiles)
            if mol:
                canonical = Chem.MolToSmiles(mol)
                salt_canonical_lookup[canonical] = salt_name

        fragments = smiles.split(".")
        kept_fragments: List[str] = []
        removed_names: List[str] = []

        for frag in fragments:
            frag_mol = Chem.MolFromSmiles(frag)
            if frag_mol is None:
                continue

            frag_canonical = Chem.MolToSmiles(frag_mol)

            if frag_canonical in salt_canonical_lookup:
                removed_names.append(salt_canonical_lookup[frag_canonical])
            else:
                kept_fragments.append(frag)

        if not kept_fragments:
            return smiles, []

        # Keep largest fragment by heavy atom count
        if len(kept_fragments) > 1:
            frag_sizes = []
            for frag in kept_fragments:
                mol = Chem.MolFromSmiles(frag)
                if mol:
                    frag_sizes.append((frag, mol.GetNumHeavyAtoms()))
                else:
                    frag_sizes.append((frag, 0))
            frag_sizes.sort(key=lambda x: x[1], reverse=True)
            kept_fragments = [frag_sizes[0][0]]
            for frag, _ in frag_sizes[1:]:
                removed_names.append(f"fragment ({frag})")

        return kept_fragments[0], removed_names

    def detect_protecting_groups(self, smiles: str) -> List[Tuple[str, str]]:
        """
        Detect protecting groups in a molecule.

        Args:
            smiles: SMILES string

        Returns:
            List of (group_name, matched_smarts) tuples
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return []

        found = []
        for pg in self.protecting_groups:
            pattern = Chem.MolFromSmarts(pg.smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                found.append((pg.name, pg.smarts))

        return found

    # ===== Reports =====

    @property
    def compatibility_report(self) -> pl.DataFrame:
        """
        Tabular report of validation results.

        Columns: Position, Total_Reagents, Compatible, Incompatible, Coverage_Percent
        """
        if self._cached_result is None:
            self.validate(test_reactions=False)

        if self._cached_result is None:
            return pl.DataFrame()

        report_data = []

        for position in range(len(self.reagents)):
            compatible = self._cached_result.compatible_reagents.get(position, [])
            incompatible = self._cached_result.incompatible_reagents.get(position, [])

            report_data.append(
                {
                    "Position": position,
                    "Total_Reagents": len(self.reagents[position]),
                    "Compatible": len(compatible),
                    "Incompatible": len(incompatible),
                    "Coverage_Percent": self._cached_result.coverage_stats.get(
                        position, 0.0
                    ),
                }
            )

        return pl.DataFrame(report_data)

    def summary(self) -> str:
        """
        Human-readable summary of validation results.

        Returns:
            Multi-line string with validation summary
        """
        if self._cached_result is None:
            self.validate(test_reactions=False)

        if self._cached_result is None:
            return "Validation failed - check errors"

        lines = [
            "SMARTS Validation Summary",
            "=" * 50,
            f"Reaction: {self.reaction_smarts[:60]}{'...' if len(self.reaction_smarts) > 60 else ''}",
            "",
        ]

        for position in range(len(self.reagents)):
            compatible = self._cached_result.compatible_reagents.get(position, [])
            incompatible = self._cached_result.incompatible_reagents.get(position, [])
            total = len(self.reagents[position])
            coverage = self._cached_result.coverage_stats.get(position, 0.0)

            file_name = (
                Path(self.reagent_files[position]).name
                if position < len(self.reagent_files)
                else f"Position {position}"
            )

            lines.append(f"Position {position} ({file_name}):")
            lines.append(f"  Compatible: {len(compatible)}/{total} ({coverage:.1f}%)")
            lines.append(f"  Incompatible: {len(incompatible)}")

            # Show protecting groups if any
            protected = self.protected_reagents.get(position, [])
            if protected:
                lines.append(f"  Protected: {len(protected)}")

            lines.append("")

        # Overall stats
        total_compatible = self._cached_result.total_compatible
        total_incompatible = self._cached_result.total_incompatible
        total = total_compatible + total_incompatible
        overall_pct = (total_compatible / total * 100) if total > 0 else 0

        lines.append(f"Overall: {overall_pct:.1f}% compatible")

        if self._cached_result.reaction_success_rate > 0:
            lines.append(
                f"Reaction success rate: {self._cached_result.reaction_success_rate:.1f}%"
            )

        if self._cached_result.error_messages:
            lines.append("")
            lines.append("Errors:")
            for err in self._cached_result.error_messages:
                lines.append(f"  - {err}")

        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"_SMARTSValidator(reaction='{self.reaction_smarts[:40]}...', "
            f"positions={len(self.reagents)})"
        )
