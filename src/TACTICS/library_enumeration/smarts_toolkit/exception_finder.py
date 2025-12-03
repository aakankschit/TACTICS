"""
Automatically discover exception rules for SMARTS patterns
"""

from typing import List, Tuple, Dict, Optional
from collections import defaultdict, Counter
from rdkit import Chem
from rdkit.Chem import AllChem
import logging

logger = logging.getLogger(__name__)

class ExceptionFinder:
    """
    Find and suggest exception rules for reaction SMARTS
    """
    
    def __init__(self, validator):
        """
        Initialize with a SMARTSValidator instance
        """
        self.validator = validator
        self.exceptions = []
        
    def find_exceptions(self, min_occurrence: int = 3) -> List[Tuple[str, int, str]]:
        """
        Find patterns that require exception rules
        
        Parameters:
        -----------
        min_occurrence : int
            Minimum number of failures with same pattern to suggest exception
            
        Returns:
        --------
        List of (substructure_smarts, position, exception_smarts) tuples
        """
        result = self.validator.validate(test_reactions=False)
        
        exceptions = []
        
        for position, incompatible_list in result.incompatible_reagents.items():
            if len(incompatible_list) < min_occurrence:
                continue
                
            # Group incompatible reagents by structural features
            pattern_groups = self._group_by_patterns(incompatible_list)
            
            # For each pattern group, try to create exception rule
            for pattern, reagents in pattern_groups.items():
                if len(reagents) >= min_occurrence:
                    exception_smarts = self._create_exception_smarts(
                        pattern, 
                        position,
                        reagents
                    )
                    if exception_smarts:
                        exceptions.append((pattern, position, exception_smarts))
                        logger.info(f"Found exception for position {position}: {pattern}")
        
        self.exceptions = exceptions
        return exceptions
    
    def _group_by_patterns(self, incompatible_list: List[Tuple]) -> Dict[str, List]:
        """Group incompatible reagents by common structural patterns"""
        
        # Patterns to check for
        test_patterns = [
            # Protecting groups
            ('Boc protected', '[C](=O)OC(C)(C)C'),
            ('Fmoc protected', 'c1ccc2c(c1)cc1ccccc1c2C'),
            ('Cbz protected', '[C](=O)OCc1ccccc1'),
            
            # Functional groups that might interfere
            ('Phenol', 'c[OH]'),
            ('Thiol', '[SH]'),
            ('Secondary amine', '[NH]([#6])[#6]'),
            ('Tertiary amine', '[N]([#6])([#6])[#6]'),
            ('Ester', 'C(=O)O[#6]'),
            ('Anhydride', 'C(=O)OC(=O)'),
            ('Acyl halide', 'C(=O)[Cl,Br,I]'),
            
            # Sterically hindered
            ('Tertiary carbon', '[C]([#6])([#6])([#6])'),
            ('Quaternary carbon', '[C]([#6])([#6])([#6])([#6])'),
            
            # Electronic effects
            ('Electron withdrawing', '[N+](=O)[O-],C(F)(F)F,C#N'),
            ('Electron donating', 'O[CH3],[N]([CH3])[CH3]'),
        ]
        
        pattern_groups = defaultdict(list)
        
        for smiles, name, reason in incompatible_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
                
            for pattern_name, pattern_smarts in test_patterns:
                pattern_mol = Chem.MolFromSmarts(pattern_smarts)
                if pattern_mol and mol.HasSubstructMatch(pattern_mol):
                    pattern_groups[pattern_smarts].append((smiles, name))
                    break
            else:
                # No specific pattern found, group as "other"
                pattern_groups['other'].append((smiles, name))
        
        return pattern_groups
    
    def _create_exception_smarts(self, pattern: str, position: int, 
                                reagents: List[Tuple]) -> Optional[str]:
        """
        Create an exception SMARTS pattern for problematic substructures
        """
        if pattern == 'other':
            return None
            
        # Get the original template
        template = self.validator.reaction.GetReactantTemplate(position)
        template_smarts = Chem.MolToSmarts(template)
        
        # Try to create a modified SMARTS that excludes the problematic pattern
        # This is simplified - in practice would need more sophisticated logic
        
        # For protecting groups, might need different reaction conditions
        protecting_group_alternatives = {
            '[C](=O)OC(C)(C)C': f'{template_smarts}.{{Boc-deprotection required}}',
            'c1ccc2c(c1)cc1ccccc1c2C': f'{template_smarts}.{{Fmoc-deprotection required}}',
            '[C](=O)OCc1ccccc1': f'{template_smarts}.{{Cbz-deprotection required}}',
        }
        
        if pattern in protecting_group_alternatives:
            return protecting_group_alternatives[pattern]
        
        # For other patterns, try to create an exclusion
        # This would need to be customized based on the reaction type
        modified_smarts = f'{template_smarts}.[!{pattern}]'
        
        return modified_smarts
    
    def validate_exceptions(self) -> Dict:
        """
        Validate that exception rules improve compatibility
        """
        if not self.exceptions:
            self.find_exceptions()
            
        validation_results = {}
        
        for pattern, position, exception_smarts in self.exceptions:
            # Test the exception rule
            try:
                exception_rxn = AllChem.ReactionFromSmarts(exception_smarts)
                
                # Count how many previously incompatible reagents now work
                result = self.validator.validate(test_reactions=False)
                incompatible = result.incompatible_reagents.get(position, [])
                
                fixed_count = 0
                still_broken = 0
                
                for smiles, name, reason in incompatible:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        pattern_mol = Chem.MolFromSmarts(pattern)
                        if mol.HasSubstructMatch(pattern_mol):
                            # This reagent has the problematic pattern
                            # Check if exception rule helps
                            # (simplified - would need actual testing)
                            fixed_count += 1
                
                validation_results[pattern] = {
                    'position': position,
                    'exception_smarts': exception_smarts,
                    'potentially_fixed': fixed_count,
                    'pattern': pattern
                }
                
            except Exception as e:
                logger.error(f"Invalid exception SMARTS: {exception_smarts}")
                validation_results[pattern] = {
                    'error': str(e)
                }
        
        return validation_results
    
    def export_exceptions(self, filename: str = "exception_rules.txt"):
        """Export exception rules to a file"""
        if not self.exceptions:
            self.find_exceptions()
            
        with open(filename, 'w') as f:
            f.write("# Exception Rules for Reaction SMARTS\n")
            f.write(f"# Original SMARTS: {self.validator.reaction_smarts}\n\n")
            
            for pattern, position, exception_smarts in self.exceptions:
                f.write(f"# Position {position} - Pattern: {pattern}\n")
                f.write(f"({pattern!r}, {position}, {exception_smarts!r}),\n\n")
                
        logger.info(f"Exported {len(self.exceptions)} exception rules to {filename}")