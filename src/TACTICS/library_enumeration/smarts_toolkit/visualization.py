"""
Visualization tools for SMARTS pattern debugging
"""

from typing import List, Optional, Dict, Tuple
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import IPythonConsole
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from IPython.display import display, HTML
import io
import base64

class SMARTSVisualizer:
    """
    Visualize SMARTS patterns and their matches
    """
    
    def __init__(self, validator):
        """Initialize with a SMARTSValidator instance"""
        self.validator = validator
        
    def visualize_reaction(self) -> None:
        """Visualize the reaction SMARTS pattern"""
        if not self.validator.reaction:
            print("No reaction loaded")
            return
            
        # Draw reaction
        rxn_image = Draw.ReactionToImage(self.validator.reaction, subImgSize=(200, 200))
        display(rxn_image)
        
    def visualize_template_matches(self, position: int, num_examples: int = 6):
        """
        Visualize template matching for a specific position
        
        Parameters:
        -----------
        position : int
            Reaction component position
        num_examples : int
            Number of examples to show
        """
        result = self.validator.validate(test_reactions=False)
        
        template = self.validator.reaction.GetReactantTemplate(position)
        
        # Get examples
        compatible = result.compatible_reagents.get(position, [])[:num_examples//2]
        incompatible = result.incompatible_reagents.get(position, [])[:num_examples//2]
        
        # Create grid visualization
        all_examples = []
        labels = []
        colors = []
        
        for smiles, name in compatible:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                all_examples.append(mol)
                labels.append(f"✓ {name[:15]}")
                colors.append('green')
                
        for smiles, name, reason in incompatible:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                all_examples.append(mol)
                labels.append(f"✗ {name[:15]}\n{reason[:20]}")
                colors.append('red')
        
        if all_examples:
            # Highlight matching substructure in compatible molecules
            for i, mol in enumerate(all_examples[:len(compatible)]):
                matches = mol.GetSubstructMatches(template)
                if matches:
                    all_examples[i].__sssAtoms = matches[0]
            
            img = Draw.MolsToGridImage(
                all_examples,
                molsPerRow=3,
                subImgSize=(250, 250),
                legends=labels,
                legendFontSize=10
            )
            
            display(img)
            print(f"\nTemplate SMARTS: {Chem.MolToSmarts(template)}")
            print(f"Compatible: {len(compatible)}, Incompatible: {len(incompatible)}")
    
    def visualize_exception_patterns(self):
        """Visualize patterns that need exception rules"""
        from .exception_finder import ExceptionFinder
        
        finder = ExceptionFinder(self.validator)
        exceptions = finder.find_exceptions()
        
        if not exceptions:
            print("No exception patterns found")
            return
        
        for pattern_smarts, position, exception_smarts in exceptions:
            print(f"\n=== Position {position} ===")
            print(f"Problem pattern: {pattern_smarts}")
            
            # Show examples of molecules with this pattern
            result = self.validator.validate(test_reactions=False)
            incompatible = result.incompatible_reagents.get(position, [])
            
            examples = []
            pattern_mol = Chem.MolFromSmarts(pattern_smarts)
            
            for smiles, name, reason in incompatible[:6]:
                mol = Chem.MolFromSmiles(smiles)
                if mol and mol.HasSubstructMatch(pattern_mol):
                    matches = mol.GetSubstructMatches(pattern_mol)
                    if matches:
                        mol.__sssAtoms = matches[0]
                    examples.append(mol)
                    
            if examples:
                img = Draw.MolsToGridImage(
                    examples,
                    molsPerRow=3,
                    subImgSize=(200, 200),
                    legends=[f"Has {pattern_smarts}" for _ in examples]
                )
                display(img)
    
    def generate_summary_plot(self):
        """Generate a summary plot of compatibility statistics"""
        result = self.validator.validate(test_reactions=False)
        
        positions = list(range(len(self.validator.reagents)))
        compatible_counts = [len(result.compatible_reagents.get(i, [])) for i in positions]
        incompatible_counts = [len(result.incompatible_reagents.get(i, [])) for i in positions]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Bar chart of compatible vs incompatible
        x = range(len(positions))
        width = 0.35
        
        ax1.bar([i - width/2 for i in x], compatible_counts, width, label='Compatible', color='green', alpha=0.7)
        ax1.bar([i + width/2 for i in x], incompatible_counts, width, label='Incompatible', color='red', alpha=0.7)
        ax1.set_xlabel('Reaction Position')
        ax1.set_ylabel('Number of Reagents')
        ax1.set_title('Reagent Compatibility by Position')
        ax1.set_xticks(x)
        ax1.set_xticklabels([f'Pos {i}' for i in positions])
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Coverage pie chart
        coverage_values = [result.coverage_stats.get(i, 0) for i in positions]
        colors_pie = ['green' if v > 75 else 'orange' if v > 50 else 'red' for v in coverage_values]
        
        ax2.bar(x, coverage_values, color=colors_pie, alpha=0.7)
        ax2.set_xlabel('Reaction Position')
        ax2.set_ylabel('Coverage (%)')
        ax2.set_title('Template Coverage by Position')
        ax2.set_xticks(x)
        ax2.set_xticklabels([f'Pos {i}' for i in positions])
        ax2.axhline(y=75, color='g', linestyle='--', alpha=0.5, label='Good (>75%)')
        ax2.axhline(y=50, color='orange', linestyle='--', alpha=0.5, label='Fair (>50%)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
        
        # Print summary statistics
        print("\n=== Summary Statistics ===")
        for i in positions:
            print(f"Position {i}:")
            print(f"  Coverage: {result.coverage_stats.get(i, 0):.1f}%")
            print(f"  Compatible: {compatible_counts[i]}")
            print(f"  Incompatible: {incompatible_counts[i]}")