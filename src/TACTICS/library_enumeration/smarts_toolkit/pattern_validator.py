"""
Core SMARTS pattern validation and testing functionality
"""

from typing import List, Dict, Tuple, Optional, Union
from dataclasses import dataclass
import pandas as pd
import polars as pl
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


@dataclass
class ProtectingGroupInfo:
    """Information about a protecting group for detection and deprotection.

    Attributes
    ----------
    name : str
        Human-readable name (e.g., "Boc", "Fmoc")
    smarts : str
        SMARTS pattern to detect the protecting group
    deprotection_smarts : Optional[str]
        Reaction SMARTS for deprotection (None if manual deprotection required)
    """

    name: str
    smarts: str
    deprotection_smarts: Optional[str] = None

    def __post_init__(self):
        """Validate SMARTS patterns."""
        if Chem.MolFromSmarts(self.smarts) is None:
            raise ValueError(f"Invalid detection SMARTS for {self.name}: {self.smarts}")
        if self.deprotection_smarts is not None:
            rxn = AllChem.ReactionFromSmarts(self.deprotection_smarts)
            if rxn is None:
                raise ValueError(
                    f"Invalid deprotection SMARTS for {self.name}: {self.deprotection_smarts}"
                )


# Common salts and counterions found in reagent storage
# These are typically separated by "." in SMILES and should be removed
# Format: (SMILES, name) - exact SMILES match for the fragment
DEFAULT_SALT_FRAGMENTS: List[Tuple[str, str]] = [
    # Halide counterions
    ("[F-]", "fluoride"),
    ("[Cl-]", "chloride"),
    ("[Br-]", "bromide"),
    ("[I-]", "iodide"),
    # Common acid counterions
    ("O=S(=O)([O-])C(F)(F)F", "triflate"),
    ("O=S(=O)([O-])C", "mesylate"),
    ("O=S(=O)([O-])c1ccc(C)cc1", "tosylate"),
    ("[O-][N+](=O)=O", "nitrate"),
    # Phosphorus-containing (stability agents)
    ("O=P(O)(O)O", "phosphoric acid"),
    ("O=P([O-])([O-])[O-]", "phosphate"),
    ("[O-]P([O-])([O-])=O", "phosphate (alt)"),
    # Boron-containing
    ("O=B([O-])[O-]", "borate"),
    ("FB(F)F", "trifluoroborate"),
    ("[B-](F)(F)(F)F", "tetrafluoroborate"),
    # Carboxylate counterions
    ("CC(=O)[O-]", "acetate"),
    ("OC(=O)C(O)C(=O)O", "tartrate"),
    ("O=C([O-])C([O-])=O", "oxalate"),
    ("OC(=O)CC(O)(CC(=O)O)C(=O)O", "citrate"),
    ("O=C([O-])C(F)(F)F", "trifluoroacetate (TFA)"),
    # Sodium/potassium (common cations)
    ("[Na+]", "sodium"),
    ("[K+]", "potassium"),
    ("[Li+]", "lithium"),
    ("[Cs+]", "cesium"),
    # Ammonium salts
    ("[NH4+]", "ammonium"),
    # Water and solvents (sometimes included)
    ("O", "water"),
    ("CO", "methanol"),
    ("CCO", "ethanol"),
    # Hydrochloride/hydrobromide (protonated forms)
    ("Cl", "HCl"),
    ("Br", "HBr"),
]


# Default protecting groups with detection and deprotection SMARTS
# Note: SMARTS patterns are ordered from most specific to most general
DEFAULT_PROTECTING_GROUPS: List[ProtectingGroupInfo] = [
    # Amine protecting groups (carbamates)
    ProtectingGroupInfo(
        name="Boc",
        smarts="[NX3;H1,H2]C(=O)OC(C)(C)C",  # N-Boc (carbamate)
        deprotection_smarts="[N:1]C(=O)OC(C)(C)C>>[N:1]",
    ),
    ProtectingGroupInfo(
        name="Fmoc",
        smarts="[NX3;H1,H2]C(=O)OCC1c2ccccc2-c3ccccc13",  # N-Fmoc
        deprotection_smarts="[N:1]C(=O)OCC1c2ccccc2-c3ccccc13>>[N:1]",
    ),
    ProtectingGroupInfo(
        name="Cbz",
        smarts="[NX3;H1,H2]C(=O)OCc1ccccc1",  # N-Cbz
        deprotection_smarts="[N:1]C(=O)OCc1ccccc1>>[N:1]",
    ),
    ProtectingGroupInfo(
        name="Acetamide",
        smarts="[NX3;H1,H2]C(=O)[CH3]",  # N-acetyl (acetamide, not generic amide)
        deprotection_smarts="[N:1]C(=O)[CH3]>>[N:1]",
    ),
    # Hydroxyl protecting groups
    ProtectingGroupInfo(
        name="TBS",
        smarts="[OX2][Si](C)(C)C(C)(C)C",  # O-TBS/TBDMS
        deprotection_smarts="[O:1][Si](C)(C)C(C)(C)C>>[O:1]",
    ),
    ProtectingGroupInfo(
        name="O-Bn",
        smarts="[OX2;H0]Cc1ccccc1",  # O-benzyl (not phenol)
        deprotection_smarts="[O:1]Cc1ccccc1>>[O:1]",
    ),
    ProtectingGroupInfo(
        name="Trityl",
        smarts="[OX2,NX3]C(c1ccccc1)(c2ccccc2)c3ccccc3",  # Trityl
        deprotection_smarts=None,  # Manual deprotection required
    ),
    # Carboxylic acid protecting groups (esters, not carbamates)
    ProtectingGroupInfo(
        name="tBu-ester",
        smarts="[CX3;!$(C-N)](=O)OC(C)(C)C",  # tBu ester (exclude carbamates)
        deprotection_smarts="[C:1](=O)OC(C)(C)C>>[C:1](=O)O",
    ),
    ProtectingGroupInfo(
        name="Me-ester",
        smarts="[CX3;!$(C-N)](=O)O[CH3]",  # Methyl ester (exclude carbamates)
        deprotection_smarts="[C:1](=O)O[CH3]>>[C:1](=O)O",
    ),
    ProtectingGroupInfo(
        name="Et-ester",
        smarts="[CX3;!$(C-N)](=O)OCC",  # Ethyl ester
        deprotection_smarts="[C:1](=O)OCC>>[C:1](=O)O",
    ),
]


@dataclass
class ValidationResult:
    """Container for validation results"""

    compatible_reagents: Dict[int, List[Tuple[str, str]]]  # (SMILES, name)
    incompatible_reagents: Dict[int, List[Tuple[str, str]]]  # (SMILES, name)
    invalid_smiles: Dict[int, List[Tuple[str, str]]]  # Failed RDKit parsing
    duplicate_smiles: Dict[int, List[Tuple[str, str]]]  # Removed duplicates
    protected_reagents: Dict[
        int, List[Tuple[str, str, List[str]]]
    ]  # (SMILES, name, [protecting_group_names])
    multi_fragment_reagents: Dict[
        int, List[Tuple[str, str, List[str]]]
    ]  # (SMILES, name, [detected_salt_names])
    coverage_stats: Dict[int, float]
    reaction_success_rate: float
    error_messages: List[str]
    warnings: List[str]

class SMARTSValidator:
    """
    Validate SMARTS patterns against reagent libraries
    """

    def __init__(
        self,
        reaction_smarts: str,
        reagent_files: List[str],
        additional_protecting_groups: Optional[List[ProtectingGroupInfo]] = None,
        additional_salt_fragments: Optional[List[Tuple[str, str]]] = None,
        deprotect_on_import: bool = False,
        deprotect_groups: Optional[List[str]] = None,
        desalt_on_import: bool = False,
    ):
        """
        Initialize validator with reaction SMARTS and reagent files

        Parameters
        ----------
        reaction_smarts : str
            Reaction SMARTS pattern to validate
        reagent_files : List[str]
            List of reagent file paths (.smi, .csv, etc.)
        additional_protecting_groups : Optional[List[ProtectingGroupInfo]]
            Additional protecting groups to detect beyond the defaults.
            These are added to DEFAULT_PROTECTING_GROUPS.
        additional_salt_fragments : Optional[List[Tuple[str, str]]]
            Additional salt/counterion fragments to detect beyond the defaults.
            Format: [(SMILES, name), ...]. These are added to DEFAULT_SALT_FRAGMENTS.
        deprotect_on_import : bool
            If True, apply deprotection reactions to all reagents during import.
            This can help match reagents that would otherwise fail substructure
            matching due to protecting groups. Default is False.
        deprotect_groups : Optional[List[str]]
            List of protecting group names to remove during import (e.g., ["Boc", "Fmoc"]).
            Only used if deprotect_on_import=True. If None, removes all detectable
            protecting groups.
        desalt_on_import : bool
            If True, remove salt/counterion fragments from multi-fragment SMILES
            during import. This keeps the largest organic fragment. Default is False.
        """
        self.reaction_smarts = reaction_smarts
        self.reagent_files = reagent_files
        self.reaction = None
        self.reagents: List[List[Tuple[str, str]]] = []
        self.errors: List[str] = []
        self.warnings: List[str] = []

        # Import-time processing options
        self.deprotect_on_import = deprotect_on_import
        self.deprotect_groups = deprotect_groups
        self.desalt_on_import = desalt_on_import

        # Protecting groups configuration
        self.protecting_groups: List[ProtectingGroupInfo] = list(DEFAULT_PROTECTING_GROUPS)
        if additional_protecting_groups:
            self.protecting_groups.extend(additional_protecting_groups)

        # Salt fragments configuration
        self.salt_fragments: List[Tuple[str, str]] = list(DEFAULT_SALT_FRAGMENTS)
        if additional_salt_fragments:
            self.salt_fragments.extend(additional_salt_fragments)

        # Preprocessing results
        self.has_original_names: Dict[int, bool] = {}
        self.invalid_smiles: Dict[int, List[Tuple[str, str]]] = {}
        self.duplicate_smiles: Dict[int, List[Tuple[str, str]]] = {}
        self.protected_reagents: Dict[int, List[Tuple[str, str, List[str]]]] = {}
        self.multi_fragment_reagents: Dict[int, List[Tuple[str, str, List[str]]]] = {}
        # Track modifications made during import
        self.deprotected_on_import: Dict[int, List[Tuple[str, str, str, List[str]]]] = {}
        self.desalted_on_import: Dict[int, List[Tuple[str, str, str, List[str]]]] = {}

        # Load reaction
        self._load_reaction()

        # Load and preprocess reagents
        self._load_reagents()
        
    def _load_reaction(self):
        """Load and validate the reaction SMARTS"""
        try:
            self.reaction = AllChem.ReactionFromSmarts(self.reaction_smarts)
            if self.reaction is None:
                raise ValueError("Invalid reaction SMARTS")
            
            # Check reaction validity
            num_reactants = self.reaction.GetNumReactantTemplates()
            num_products = self.reaction.GetNumProductTemplates()
            
            if num_reactants == 0:
                self.errors.append("No reactant templates found in reaction SMARTS")
            if num_products == 0:
                self.errors.append("No product templates found in reaction SMARTS")
                
            logger.info(f"Loaded reaction with {num_reactants} reactants and {num_products} products")
            
        except Exception as e:
            self.errors.append(f"Failed to parse reaction SMARTS: {str(e)}")
            logger.error(f"Reaction loading failed: {e}")
    
    def _load_reagents(self):
        """Load reagents from files and preprocess them"""
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

            # Handle different file formats
            if path.suffix == ".smi":
                reagent_list, has_names = self._load_smi_file(path)
            elif path.suffix == ".csv":
                reagent_list, has_names = self._load_csv_file(path)
            else:
                self.warnings.append(f"Unknown file format: {path.suffix}")

            self.has_original_names[i] = has_names
            self.reagents.append(reagent_list)
            logger.info(f"Loaded {len(reagent_list)} reagents from {filepath}")

        # Preprocess all reagents (validate SMILES, remove duplicates)
        self._preprocess_reagents()

    def _preprocess_reagents(self):
        """Validate SMILES, apply transformations, and remove duplicates from loaded reagents.

        This method:
        1. Validates each SMILES using RDKit (invalid ones are removed and tracked)
        2. Optionally applies desalting (if desalt_on_import=True)
        3. Optionally applies deprotection (if deprotect_on_import=True)
        4. Removes duplicate SMILES within each position (first occurrence kept)

        Results are stored in:
        - self.invalid_smiles: Dict mapping position -> list of invalid (SMILES, name)
        - self.duplicate_smiles: Dict mapping position -> list of duplicate (SMILES, name)
        - self.desalted_on_import: Dict mapping position -> list of (original, new, name, removed)
        - self.deprotected_on_import: Dict mapping position -> list of (original, new, name, removed)
        - self.reagents: Updated to contain only valid, unique reagents
        """
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
                    self.warnings.append(f"Invalid SMILES at position {position}: {smiles} ({name})")
                    continue

                # Apply desalting if requested (before deprotection)
                if self.desalt_on_import and "." in smiles:
                    desalted_smiles, removed_salts = self.desalt_smiles(smiles)
                    if removed_salts:
                        desalted_list.append((smiles, desalted_smiles, name, removed_salts))
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
                        deprotected_list.append((smiles, deprotected_smiles, name, removed_groups))
                        smiles = deprotected_smiles
                        logger.debug(
                            f"Position {position}: Deprotected {name}: removed {removed_groups}"
                        )

                # Check for duplicates (within this position only, using final SMILES)
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

            # Update reagents list with only valid, unique entries
            self.reagents[position] = valid_reagents

            # Log summary
            if invalid_list:
                logger.warning(
                    f"Position {position}: Removed {len(invalid_list)} invalid SMILES"
                )
            if duplicate_list:
                logger.info(
                    f"Position {position}: Removed {len(duplicate_list)} duplicate SMILES"
                )
            if desalted_list:
                logger.info(
                    f"Position {position}: Desalted {len(desalted_list)} reagents on import"
                )
            if deprotected_list:
                logger.info(
                    f"Position {position}: Deprotected {len(deprotected_list)} reagents on import"
                )

        # Detect protecting groups and multi-fragment reagents after preprocessing
        self._detect_protecting_groups()
        self._detect_multi_fragment_reagents()

    def _detect_protecting_groups(self):
        """Detect protecting groups in all reagents.

        Scans each reagent for known protecting group substructures and stores
        the results in self.protected_reagents.

        Results are stored in:
        - self.protected_reagents: Dict mapping position -> list of
          (SMILES, name, [protecting_group_names])
        """
        for position, reagent_list in enumerate(self.reagents):
            protected_list: List[Tuple[str, str, List[str]]] = []

            for smiles, name in reagent_list:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue  # Already filtered in preprocessing, but be safe

                # Check for each protecting group
                found_groups: List[str] = []
                for pg in self.protecting_groups:
                    pattern = Chem.MolFromSmarts(pg.smarts)
                    if pattern and mol.HasSubstructMatch(pattern):
                        found_groups.append(pg.name)

                if found_groups:
                    protected_list.append((smiles, name, found_groups))

            self.protected_reagents[position] = protected_list

            # Log summary
            if protected_list:
                logger.info(
                    f"Position {position}: Found {len(protected_list)} reagents "
                    f"with protecting groups"
                )

    def _detect_multi_fragment_reagents(self):
        """Detect reagents containing multiple fragments (salts, counterions, solvents).

        Scans each reagent SMILES for "." separators indicating multiple fragments,
        then identifies known salt/counterion fragments.

        Results are stored in:
        - self.multi_fragment_reagents: Dict mapping position -> list of
          (SMILES, name, [detected_fragment_names])
        """
        # Build a lookup dict for canonical SMILES of salt fragments
        salt_canonical_lookup: Dict[str, str] = {}
        for salt_smiles, salt_name in self.salt_fragments:
            mol = Chem.MolFromSmiles(salt_smiles)
            if mol:
                canonical = Chem.MolToSmiles(mol)
                salt_canonical_lookup[canonical] = salt_name

        for position, reagent_list in enumerate(self.reagents):
            multi_frag_list: List[Tuple[str, str, List[str]]] = []

            for smiles, name in reagent_list:
                # Check if SMILES contains multiple fragments
                if "." not in smiles:
                    continue

                # Split into fragments and identify each
                fragments = smiles.split(".")
                detected_salts: List[str] = []

                for frag in fragments:
                    frag_mol = Chem.MolFromSmiles(frag)
                    if frag_mol:
                        frag_canonical = Chem.MolToSmiles(frag_mol)
                        if frag_canonical in salt_canonical_lookup:
                            detected_salts.append(salt_canonical_lookup[frag_canonical])

                # Record this reagent as multi-fragment
                multi_frag_list.append((smiles, name, detected_salts))

            self.multi_fragment_reagents[position] = multi_frag_list

            # Log summary
            if multi_frag_list:
                logger.info(
                    f"Position {position}: Found {len(multi_frag_list)} reagents "
                    f"with multiple fragments (salts/counterions)"
                )

    def _load_smi_file(self, path: Path) -> Tuple[List[Tuple[str, str]], bool]:
        """Load SMILES from .smi file

        Returns:
            Tuple of (reagent_list, has_original_names)
            - reagent_list: List of (SMILES, name) tuples
            - has_original_names: True if file contained name column
        """
        reagents = []
        has_original_names = None  # Will be determined from first valid line

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

        # Default to False if file was empty
        return reagents, has_original_names if has_original_names is not None else False
    
    def _load_csv_file(self, path: Path) -> Tuple[List[Tuple[str, str]], bool]:
        """Load SMILES from CSV file

        Returns:
            Tuple of (reagent_list, has_original_names)
            - reagent_list: List of (SMILES, name) tuples
            - has_original_names: True if file contained name/ID column
        """
        df = pd.read_csv(path)

        # Find SMILES column
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
    
    def validate(self, test_reactions: bool = True, 
                 sample_size: Optional[int] = None) -> ValidationResult:
        """
        Perform comprehensive validation
        
        Parameters:
        -----------
        test_reactions : bool
            Whether to test actual reactions (can be slow for large libraries)
        sample_size : Optional[int]
            If set, only test this many random combinations
            
        Returns:
        --------
        ValidationResult containing all validation information
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

        return ValidationResult(
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
    
    def _check_template_matching(self) -> Tuple[Dict, Dict]:
        """Check which reagents match reaction templates.

        Note: Invalid SMILES are already filtered out during preprocessing,
        so this method only checks template compatibility.
        """
        compatible: Dict[int, List[Tuple[str, str]]] = {}
        incompatible: Dict[int, List[Tuple[str, str]]] = {}

        num_templates = self.reaction.GetNumReactantTemplates()

        for i in range(num_templates):
            template = self.reaction.GetReactantTemplate(i)
            compatible[i] = []
            incompatible[i] = []

            if i >= len(self.reagents):
                self.warnings.append(f"No reagents provided for position {i}")
                continue

            for smiles, name in self.reagents[i]:
                # SMILES already validated during preprocessing
                mol = Chem.MolFromSmiles(smiles)

                if mol.HasSubstructMatch(template):
                    compatible[i].append((smiles, name))
                else:
                    incompatible[i].append((smiles, name))

        return compatible, incompatible
    
    def _calculate_coverage(self, compatible: Dict) -> Dict[int, float]:
        """Calculate coverage statistics"""
        coverage = {}
        
        for i, reagent_list in enumerate(self.reagents):
            if i in compatible:
                num_compatible = len(compatible[i])
                total = len(reagent_list)
                coverage[i] = (num_compatible / total * 100) if total > 0 else 0
            else:
                coverage[i] = 0.0
                
        return coverage
    
    def _test_reaction_combinations(self, compatible: Dict, 
                                   sample_size: Optional[int] = None) -> float:
        """Test actual reaction combinations"""
        import random
        import itertools
        
        # Get all possible combinations
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
                # Test all combinations (be careful with large libraries!)
                combinations = list(itertools.product(*reagent_lists))
                
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

    def deprotect_smiles(
        self, smiles: str, groups_to_remove: Optional[List[str]] = None
    ) -> Tuple[str, List[str]]:
        """Apply deprotection reactions to remove protecting groups.

        This is a standalone method that can be used independently to deprotect
        any SMILES string, either during import (via constructor options) or
        on-demand for individual molecules.

        Parameters
        ----------
        smiles : str
            Input SMILES string
        groups_to_remove : Optional[List[str]]
            List of protecting group names to remove (e.g., ["Boc", "Fmoc"]).
            If None, removes all detectable protecting groups.

        Returns
        -------
        Tuple[str, List[str]]
            Tuple of (deprotected_smiles, list_of_removed_groups)

        Examples
        --------
        >>> validator = SMARTSValidator(reaction_smarts, reagent_files)
        >>> # Deprotect a single SMILES
        >>> deprotected, removed = validator.deprotect_smiles("CC(C)(C)OC(=O)NCc1ccccc1")
        >>> print(f"Deprotected: {deprotected}, removed: {removed}")
        >>> # Deprotect only specific groups
        >>> deprotected, removed = validator.deprotect_smiles(smiles, groups_to_remove=["Boc"])
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

            # Run reaction (may need multiple iterations for multiple groups)
            max_iterations = 10  # Safety limit
            for _ in range(max_iterations):
                if not mol.HasSubstructMatch(pattern):
                    break

                products = rxn.RunReactants((mol,))
                if products and len(products) > 0 and len(products[0]) > 0:
                    mol = products[0][0]
                    Chem.SanitizeMol(mol)
                    if pg.name not in removed_groups:
                        removed_groups.append(pg.name)
                else:
                    break

        return Chem.MolToSmiles(mol), removed_groups

    def desalt_smiles(self, smiles: str) -> Tuple[str, List[str]]:
        """Remove salt fragments from a multi-fragment SMILES.

        This is a standalone method that can be used independently to desalt
        any SMILES string, either during import (via constructor options) or
        on-demand for individual molecules.

        Identifies and removes known salt/counterion fragments, keeping the
        largest organic fragment as the main molecule.

        Parameters
        ----------
        smiles : str
            Input SMILES string (may contain "." separators)

        Returns
        -------
        Tuple[str, List[str]]
            Tuple of (desalted_smiles, list_of_removed_fragment_names)

        Examples
        --------
        >>> validator = SMARTSValidator(reaction_smarts, reagent_files)
        >>> # Desalt a single SMILES
        >>> desalted, removed = validator.desalt_smiles("CCN.[Cl-]")
        >>> print(f"Desalted: {desalted}, removed: {removed}")
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
            # If all fragments were salts, return the original (shouldn't happen normally)
            return smiles, []

        # If multiple non-salt fragments remain, keep the largest by heavy atom count
        if len(kept_fragments) > 1:
            # Sort by heavy atom count (descending) and keep largest
            frag_sizes = []
            for frag in kept_fragments:
                mol = Chem.MolFromSmiles(frag)
                if mol:
                    frag_sizes.append((frag, mol.GetNumHeavyAtoms()))
                else:
                    frag_sizes.append((frag, 0))
            frag_sizes.sort(key=lambda x: x[1], reverse=True)
            # Keep only the largest fragment
            kept_fragments = [frag_sizes[0][0]]
            # Add the smaller non-salt fragments to removed (as "unknown fragment")
            for frag, _ in frag_sizes[1:]:
                removed_names.append(f"fragment ({frag})")

        return kept_fragments[0], removed_names

    def export_compatible_reagents(
        self,
        output_dir: str = ".",
        prefix: Optional[Union[str, List[str]]] = None,
        deprotect: bool = False,
        deprotect_groups: Optional[List[str]] = None,
        desalt: bool = False,
    ) -> List[Path]:
        """Export compatible reagents to TACTICS-compatible .smi files.

        Parameters
        ----------
        output_dir : str
            Directory to write output files.
        prefix : Optional[Union[str, List[str]]]
            Prefix for reagent naming. Can be:
            - A single string (applied to all positions, e.g., "reagent")
            - A list of strings (one per position, e.g., ["amine", "acid"])
            Used only when original names are not available in the input file.
        deprotect : bool
            If True, apply deprotection reactions to remove protecting groups.
        deprotect_groups : Optional[List[str]]
            List of protecting group names to remove (e.g., ["Boc", "Fmoc"]).
            If None and deprotect=True, removes all detectable protecting groups.
        desalt : bool
            If True, remove salt/counterion fragments from multi-fragment SMILES.
            Keeps the largest organic fragment.

        Returns
        -------
        List[Path]
            List of paths to created files.

        Raises
        ------
        ValueError
            If a position has no original names and no prefix is provided.
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        result = self.validate(test_reactions=False)
        created_files: List[Path] = []

        # Build lookup for protected reagents
        protected_lookup: Dict[int, Dict[str, List[str]]] = {}
        for position, protected_list in result.protected_reagents.items():
            protected_lookup[position] = {
                smiles: groups for smiles, name, groups in protected_list
            }

        # Build lookup for multi-fragment reagents
        multi_frag_lookup: Dict[int, set] = {}
        for position, multi_frag_list in result.multi_fragment_reagents.items():
            multi_frag_lookup[position] = {smiles for smiles, name, _ in multi_frag_list}

        for position, reagents in result.compatible_reagents.items():
            # Determine naming strategy for this position
            use_original = self.has_original_names.get(position, False)
            pos_prefix = None

            if not use_original:
                # Need a prefix to generate names
                if prefix is None:
                    raise ValueError(
                        f"Position {position} has no original names in input file. "
                        f"Provide a prefix for naming."
                    )

                # Get prefix for this position
                if isinstance(prefix, list):
                    if position >= len(prefix):
                        raise ValueError(
                            f"Prefix list has {len(prefix)} items but position {position} "
                            f"requires a prefix. Provide a prefix for all positions."
                        )
                    pos_prefix = prefix[position]
                else:
                    pos_prefix = prefix

            # Write file
            filename = output_path / f"compatible_position_{position}.smi"
            deprotected_count = 0
            desalted_count = 0

            with open(filename, "w") as f:
                for idx, (smiles, original_name) in enumerate(reagents):
                    # Determine name
                    if use_original:
                        name = original_name
                    else:
                        name = f"{pos_prefix}_{idx}"

                    output_smiles = smiles

                    # Apply desalting first (if requested)
                    if desalt:
                        is_multi_frag = smiles in multi_frag_lookup.get(position, set())
                        if is_multi_frag:
                            output_smiles, removed_salts = self.desalt_smiles(output_smiles)
                            if removed_salts:
                                desalted_count += 1
                                logger.debug(
                                    f"Desalted {name}: removed {removed_salts}"
                                )

                    # Apply deprotection if requested
                    if deprotect:
                        # Check if this reagent has protecting groups
                        has_protection = smiles in protected_lookup.get(position, {})
                        if has_protection:
                            output_smiles, removed = self.deprotect_smiles(
                                output_smiles, deprotect_groups
                            )
                            if removed:
                                deprotected_count += 1
                                logger.debug(
                                    f"Deprotected {name}: removed {removed}"
                                )

                    f.write(f"{output_smiles}\t{name}\n")

            created_files.append(filename)
            log_msg = f"Exported {len(reagents)} compatible reagents to {filename}"
            if desalt and desalted_count > 0:
                log_msg += f" ({desalted_count} desalted)"
            if deprotect and deprotected_count > 0:
                log_msg += f" ({deprotected_count} deprotected)"
            logger.info(log_msg)

        return created_files

    def generate_compatibility_report(self) -> pl.DataFrame:
        """Generate a compatibility report as a Polars DataFrame.

        This provides a tabular summary of reagent compatibility at each position,
        including coverage statistics and failure reasons.

        Returns
        -------
        pl.DataFrame
            DataFrame with columns:
            - Position: Reagent position index
            - Total_Reagents: Total number of reagents at this position
            - Compatible: Number of compatible reagents
            - Incompatible: Number of incompatible reagents
            - Coverage_Percent: Percentage of compatible reagents (0-100)
            - Top_Failure_Reason: Most common reason for incompatibility
            - Failure_Count: Count of the top failure reason
        """
        result = self.validate(test_reactions=False)

        report_data = []

        for position in range(len(self.reagents)):
            compatible = result.compatible_reagents.get(position, [])
            incompatible = result.incompatible_reagents.get(position, [])

            failure_count = len(incompatible)

            report_data.append({
                'Position': position,
                'Total_Reagents': len(self.reagents[position]),
                'Compatible': len(compatible),
                'Incompatible': len(incompatible),
                'Coverage_Percent': result.coverage_stats.get(position, 0.0),
                'Top_Failure_Reason': 'Template mismatch' if failure_count > 0 else 'N/A',
                'Failure_Count': failure_count
            })

        return pl.DataFrame(report_data)

    @property
    def compatibility_report(self) -> pl.DataFrame:
        """Compatibility report as a class attribute (Polars DataFrame).

        This is a convenience property that calls generate_compatibility_report().
        For repeated access, consider caching the result.

        Returns
        -------
        pl.DataFrame
            Compatibility report DataFrame.
        """
        return self.generate_compatibility_report()