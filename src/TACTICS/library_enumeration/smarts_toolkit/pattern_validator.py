"""
Core SMARTS pattern validation and testing functionality
"""

from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

@dataclass
class ValidationResult:
    """Container for validation results"""
    compatible_reagents: Dict[int, List[str]]
    incompatible_reagents: Dict[int, List[str]]
    coverage_stats: Dict[int, float]
    reaction_success_rate: float
    error_messages: List[str]
    warnings: List[str]

class SMARTSValidator:
    """
    Validate SMARTS patterns against reagent libraries
    """
    
    def __init__(self, reaction_smarts: str, reagent_files: List[str]):
        """
        Initialize validator with reaction SMARTS and reagent files
        
        Parameters:
        -----------
        reaction_smarts : str
            Reaction SMARTS pattern to validate
        reagent_files : List[str]
            List of reagent file paths (.smi, .csv, etc.)
        """
        self.reaction_smarts = reaction_smarts
        self.reagent_files = reagent_files
        self.reaction = None
        self.reagents = []
        self.errors = []
        self.warnings = []
        
        # Load reaction
        self._load_reaction()
        
        # Load reagents
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
        """Load reagents from files"""
        self.reagents = []
        
        for i, filepath in enumerate(self.reagent_files):
            path = Path(filepath)
            if not path.exists():
                self.errors.append(f"Reagent file {filepath} not found")
                continue
                
            reagent_list = []
            
            # Handle different file formats
            if path.suffix == '.smi':
                reagent_list = self._load_smi_file(path)
            elif path.suffix == '.csv':
                reagent_list = self._load_csv_file(path)
            else:
                self.warnings.append(f"Unknown file format: {path.suffix}")
                
            self.reagents.append(reagent_list)
            logger.info(f"Loaded {len(reagent_list)} reagents from {filepath}")
    
    def _load_smi_file(self, path: Path) -> List[Tuple[str, str]]:
        """Load SMILES from .smi file"""
        reagents = []
        with open(path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 1:
                    smiles = parts[0]
                    name = parts[1] if len(parts) > 1 else f"R{len(reagents)}"
                    reagents.append((smiles, name))
        return reagents
    
    def _load_csv_file(self, path: Path) -> List[Tuple[str, str]]:
        """Load SMILES from CSV file"""
        df = pd.read_csv(path)
        
        # Find SMILES column
        smiles_col = None
        name_col = None
        
        for col in df.columns:
            if 'SMILES' in col.upper():
                smiles_col = col
            elif 'NAME' in col.upper() or 'ID' in col.upper():
                name_col = col
                
        if smiles_col is None:
            self.errors.append(f"No SMILES column found in {path}")
            return []
            
        reagents = []
        for idx, row in df.iterrows():
            smiles = row[smiles_col]
            name = row[name_col] if name_col else f"R{idx}"
            reagents.append((smiles, name))
            
        return reagents
    
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
                coverage_stats={},
                reaction_success_rate=0.0,
                error_messages=self.errors,
                warnings=self.warnings
            )
        
        # Check template matching
        compatible, incompatible = self._check_template_matching()
        
        # Calculate coverage statistics
        coverage = self._calculate_coverage(compatible)
        
        # Test reactions if requested
        success_rate = 0.0
        if test_reactions:
            success_rate = self._test_reaction_combinations(
                compatible, 
                sample_size
            )
        
        return ValidationResult(
            compatible_reagents=compatible,
            incompatible_reagents=incompatible,
            coverage_stats=coverage,
            reaction_success_rate=success_rate,
            error_messages=self.errors,
            warnings=self.warnings
        )
    
    def _check_template_matching(self) -> Tuple[Dict, Dict]:
        """Check which reagents match reaction templates"""
        compatible = {}
        incompatible = {}
        
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
                if mol is None:
                    self.warnings.append(f"Invalid SMILES: {smiles} ({name})")
                    incompatible[i].append((smiles, name, "Invalid SMILES"))
                    continue
                
                if mol.HasSubstructMatch(template):
                    compatible[i].append((smiles, name))
                else:
                    # Try to identify why it doesn't match
                    reason = self._diagnose_mismatch(mol, template)
                    incompatible[i].append((smiles, name, reason))
                    
        return compatible, incompatible
    
    def _diagnose_mismatch(self, mol: Chem.Mol, template: Chem.Mol) -> str:
        """Diagnose why a molecule doesn't match a template"""
        # Get template SMARTS
        template_smarts = Chem.MolToSmarts(template)
        
        # Check for specific functional groups
        reasons = []
        
        # Check atom counts
        template_atoms = template.GetNumAtoms()
        mol_atoms = mol.GetNumAtoms()
        
        if mol_atoms < template_atoms:
            reasons.append("Too few atoms")
        
        # Check for required functional groups
        # This is simplified - could be expanded
        functional_groups = {
            '[CX3](=O)[OX2H1]': 'carboxylic acid',
            '[NX3;H2,H1]': 'primary/secondary amine',
            '[CX3H1]=O': 'aldehyde',
            '[CX3]=O': 'carbonyl',
            'c[Br,I,Cl]': 'aryl halide'
        }
        
        for fg_smarts, fg_name in functional_groups.items():
            if fg_smarts in template_smarts:
                fg_mol = Chem.MolFromSmarts(fg_smarts)
                if not mol.HasSubstructMatch(fg_mol):
                    reasons.append(f"Missing {fg_name}")
                    
        return "; ".join(reasons) if reasons else "Template mismatch"
    
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
    
    def export_compatible_reagents(self, output_dir: str = "."):
        """Export compatible reagents to new files"""
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        result = self.validate(test_reactions=False)
        
        for position, reagents in result.compatible_reagents.items():
            filename = output_path / f"compatible_position_{position}.smi"
            with open(filename, 'w') as f:
                for smiles, name in reagents:
                    f.write(f"{smiles}\t{name}\n")
            logger.info(f"Exported {len(reagents)} compatible reagents to {filename}")