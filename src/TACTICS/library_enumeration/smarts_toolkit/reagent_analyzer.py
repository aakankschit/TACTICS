"""
Analyze reagent compatibility and suggest modifications
"""

from typing import List, Dict, Set, Tuple, Optional
from collections import Counter, defaultdict
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, Draw
from rdkit.Chem.Draw import IPythonConsole
import pandas as pd
import numpy as np

class ReagentCompatibilityAnalyzer:
    """
    Detailed analysis of reagent compatibility with reaction SMARTS
    """
    
    def __init__(self, validator):
        """
        Initialize with a SMARTSValidator instance
        
        Parameters:
        -----------
        validator : SMARTSValidator
            Initialized validator with reaction and reagents loaded
        """
        self.validator = validator
        self.analysis_results = {}
        
    def analyze_incompatible_reagents(self) -> Dict:
        """
        Analyze patterns in incompatible reagents
        """
        result = self.validator.validate(test_reactions=False)
        
        analysis = {}
        
        for position, incompatible_list in result.incompatible_reagents.items():
            if not incompatible_list:
                continue
                
            # Group by failure reason
            reason_groups = defaultdict(list)
            for smiles, name, reason in incompatible_list:
                reason_groups[reason].append((smiles, name))
            
            # Analyze common substructures in each group
            pattern_analysis = {}
            for reason, reagents in reason_groups.items():
                smiles_list = [r[0] for r in reagents]
                patterns = self._find_common_substructures(smiles_list)
                pattern_analysis[reason] = {
                    'count': len(reagents),
                    'examples': reagents[:5],  # First 5 examples
                    'common_patterns': patterns
                }
            
            analysis[position] = pattern_analysis
            
        self.analysis_results = analysis
        return analysis
    
    def _find_common_substructures(self, smiles_list: List[str]) -> List[str]:
        """Find common substructures in a list of SMILES"""
        if len(smiles_list) < 2:
            return []
            
        mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]
        mols = [m for m in mols if m is not None]
        
        if len(mols) < 2:
            return []
        
        # Common functional groups to check
        functional_groups = {
            '[OH]': 'hydroxyl',
            '[NH2]': 'primary amine',
            '[NH]': 'secondary amine',
            'C(=O)O': 'carboxylic acid',
            'C(=O)N': 'amide',
            'C(=O)': 'carbonyl',
            '[N+](=O)[O-]': 'nitro',
            'S(=O)(=O)': 'sulfonyl',
            'C#N': 'nitrile',
            'c': 'aromatic',
            'C=C': 'alkene',
            'C#C': 'alkyne',
            '[F,Cl,Br,I]': 'halogen',
            'O[CH3]': 'methoxy',
            '[Si]': 'silicon',
            '[P]': 'phosphorus',
            '[B]': 'boron'
        }
        
        found_patterns = []
        pattern_counts = Counter()
        
        for smarts, name in functional_groups.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern:
                count = sum(1 for mol in mols if mol.HasSubstructMatch(pattern))
                if count > len(mols) * 0.5:  # Present in >50% of molecules
                    pattern_counts[name] = count
                    
        # Sort by frequency
        found_patterns = [f"{name} ({count}/{len(mols)})" 
                          for name, count in pattern_counts.most_common(5)]
        
        return found_patterns
    
    def suggest_template_modifications(self) -> Dict[int, List[str]]:
        """
        Suggest modifications to reaction templates to increase compatibility
        """
        suggestions = {}
        
        result = self.validator.validate(test_reactions=False)
        
        for position in range(self.validator.reaction.GetNumReactantTemplates()):
            template = self.validator.reaction.GetReactantTemplate(position)
            template_smarts = Chem.MolToSmarts(template)
            
            compatible_count = len(result.compatible_reagents.get(position, []))
            incompatible_count = len(result.incompatible_reagents.get(position, []))
            
            if incompatible_count == 0:
                continue
                
            position_suggestions = []
            
            # Analyze if template is too restrictive
            if compatible_count < incompatible_count:
                position_suggestions.append(
                    f"Template may be too restrictive. Only {compatible_count} of "
                    f"{compatible_count + incompatible_count} reagents match."
                )
                
                # Suggest relaxing constraints
                relaxed_patterns = self._suggest_relaxed_patterns(template_smarts)
                for pattern in relaxed_patterns:
                    position_suggestions.append(f"Try relaxed pattern: {pattern}")
            
            # Check for common incompatible patterns
            if position in self.analysis_results:
                for reason, data in self.analysis_results[position].items():
                    if data['count'] > 5:  # Significant number of failures
                        position_suggestions.append(
                            f"{data['count']} reagents fail due to: {reason}"
                        )
                        
            suggestions[position] = position_suggestions
            
        return suggestions
    
    def _suggest_relaxed_patterns(self, template_smarts: str) -> List[str]:
        """Suggest relaxed versions of a SMARTS pattern"""
        suggestions = []
        
        # Common relaxations
        relaxations = {
            '[CX3]': '[C]',  # Any carbon instead of sp2
            '[NX3]': '[N]',  # Any nitrogen
            '[OX2]': '[O]',  # Any oxygen
            'H1': '',        # Remove H count requirement
            'H2': 'H1,H2',   # Allow more flexibility in H count
            '!': '',         # Remove NOT conditions
        }
        
        relaxed = template_smarts
        for strict, relaxed_version in relaxations.items():
            if strict in template_smarts:
                new_pattern = template_smarts.replace(strict, relaxed_version)
                if new_pattern != template_smarts:
                    suggestions.append(new_pattern)
                    
        return suggestions[:3]  # Return top 3 suggestions
    
    def generate_compatibility_report(self) -> pd.DataFrame:
        """
        Generate a detailed compatibility report as a DataFrame
        """
        result = self.validator.validate(test_reactions=False)
        
        report_data = []
        
        for position in range(len(self.validator.reagents)):
            compatible = result.compatible_reagents.get(position, [])
            incompatible = result.incompatible_reagents.get(position, [])
            
            # Aggregate failure reasons
            failure_reasons = Counter()
            for _, _, reason in incompatible:
                failure_reasons[reason] += 1
            
            report_data.append({
                'Position': position,
                'Total_Reagents': len(self.validator.reagents[position]),
                'Compatible': len(compatible),
                'Incompatible': len(incompatible),
                'Coverage_Percent': result.coverage_stats.get(position, 0),
                'Top_Failure_Reason': failure_reasons.most_common(1)[0][0] if failure_reasons else 'N/A',
                'Failure_Count': failure_reasons.most_common(1)[0][1] if failure_reasons else 0
            })
        
        return pd.DataFrame(report_data)
    
    def find_universal_reagents(self) -> Dict[int, List[Tuple[str, str]]]:
        """
        Find reagents that work with many partners (good for initial testing)
        """
        # This would require testing combinations - placeholder for now
        # Could be expanded to test compatibility across partners
        universal = {}
        
        result = self.validator.validate(test_reactions=False)
        
        for position, compatible in result.compatible_reagents.items():
            # For now, just return the first few compatible reagents
            # In a full implementation, would test cross-compatibility
            universal[position] = compatible[:5] if compatible else []
            
        return universal