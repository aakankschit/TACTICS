"""
Visualization tools for SMARTS pattern debugging
"""

from typing import Tuple
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
from IPython.display import display

class SMARTSVisualizer:
    """
    Visualize SMARTS patterns and their matches
    """

    def __init__(
        self,
        validator,
        img_size: Tuple[int, int] = (250, 250),
        mols_per_row: int = 5
    ):
        """Initialize with a SMARTSValidator instance.

        Parameters
        ----------
        validator : SMARTSValidator
            Initialized validator with reaction and reagents loaded
        img_size : Tuple[int, int]
            Size of each molecule image in grid (width, height)
        mols_per_row : int
            Number of molecules per row in grid visualizations
        """
        self.validator = validator
        self.img_size = img_size
        self.mols_per_row = mols_per_row
        
    def visualize_reaction(self) -> None:
        """Visualize the reaction SMARTS pattern"""
        if not self.validator.reaction:
            print("No reaction loaded")
            return
            
        # Draw reaction
        rxn_image = Draw.ReactionToImage(self.validator.reaction, subImgSize=(200, 200))
        display(rxn_image)
        
    def visualize_compatible_reagents(
        self, position: int, max_molecules: int = 20
    ) -> None:
        """Display grid of compatible reagents at a position with template highlights.

        Parameters
        ----------
        position : int
            Reaction component position (0-indexed)
        max_molecules : int
            Maximum number of molecules to display
        """
        result = self.validator.validate(test_reactions=False)
        template = self.validator.reaction.GetReactantTemplate(position)

        compatible = result.compatible_reagents.get(position, [])[:max_molecules]

        if not compatible:
            print(f"No compatible reagents found at position {position}")
            return

        mols = []
        labels = []

        for smiles, name in compatible:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # Highlight matching substructure
                matches = mol.GetSubstructMatches(template)
                if matches:
                    mol.__sssAtoms = matches[0]
                mols.append(mol)
                labels.append(name[:20])

        if mols:
            img = Draw.MolsToGridImage(
                mols,
                molsPerRow=self.mols_per_row,
                subImgSize=self.img_size,
                legends=labels,
                legendFontSize=10
            )
            display(img)

            print(f"\nPosition {position} - Compatible Reagents")
            print(f"Template SMARTS: {Chem.MolToSmarts(template)}")
            print(f"Showing {len(mols)} of {len(result.compatible_reagents.get(position, []))} compatible reagents")

    def visualize_incompatible_reagents(
        self, position: int, max_molecules: int = 20
    ) -> None:
        """Display grid of incompatible reagents at a position.

        Parameters
        ----------
        position : int
            Reaction component position (0-indexed)
        max_molecules : int
            Maximum number of molecules to display
        """
        result = self.validator.validate(test_reactions=False)
        template = self.validator.reaction.GetReactantTemplate(position)

        incompatible = result.incompatible_reagents.get(position, [])[:max_molecules]

        if not incompatible:
            print(f"No incompatible reagents found at position {position}")
            return

        mols = []
        labels = []

        for smiles, name in incompatible:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mols.append(mol)
                labels.append(name[:20])

        if mols:
            img = Draw.MolsToGridImage(
                mols,
                molsPerRow=self.mols_per_row,
                subImgSize=self.img_size,
                legends=labels,
                legendFontSize=10
            )
            display(img)

            print(f"\nPosition {position} - Incompatible Reagents")
            print(f"Template SMARTS: {Chem.MolToSmarts(template)}")
            print(f"Showing {len(mols)} of {len(result.incompatible_reagents.get(position, []))} incompatible reagents")
    
    def generate_summary_plot(self) -> None:
        """Generate a summary plot of compatibility statistics.

        Uses the validator's compatibility_report for statistics and displays:
        - Bar chart of compatible vs incompatible counts by position
        - Pie charts showing compatibility percentage per position
        """
        report = self.validator.compatibility_report

        positions = report["Position"].to_list()
        compatible_counts = report["Compatible"].to_list()
        incompatible_counts = report["Incompatible"].to_list()

        num_positions = len(positions)

        # Create figure with bar chart and pie charts
        fig = plt.figure(figsize=(14, 5))

        # Left: Bar chart of compatible vs incompatible
        ax1 = fig.add_subplot(1, 2, 1)
        x = range(num_positions)
        width = 0.35

        ax1.bar([i - width/2 for i in x], compatible_counts, width,
                label='Compatible', color='green', alpha=0.7)
        ax1.bar([i + width/2 for i in x], incompatible_counts, width,
                label='Incompatible', color='red', alpha=0.7)
        ax1.set_xlabel('Reaction Position')
        ax1.set_ylabel('Number of Reagents')
        ax1.set_title('Reagent Compatibility by Position')
        ax1.set_xticks(x)
        ax1.set_xticklabels([f'Pos {i}' for i in positions])
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Right: Pie charts for each position
        ax2 = fig.add_subplot(1, 2, 2)

        # Create pie chart showing overall or per-position breakdown
        if num_positions == 1:
            # Single position: show one pie chart
            sizes = [compatible_counts[0], incompatible_counts[0]]
            labels = ['Compatible', 'Incompatible']
            colors = ['green', 'red']
            ax2.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%',
                    startangle=90, explode=(0.05, 0))
            ax2.set_title(f'Position 0 Compatibility')
        else:
            # Multiple positions: show stacked or individual pie charts
            # Use a grid of small pie charts
            ax2.axis('off')
            for i, pos in enumerate(positions):
                # Create small subplot for each position
                sub_ax = fig.add_axes([0.55 + (i % 3) * 0.15,
                                       0.55 - (i // 3) * 0.45,
                                       0.12, 0.35])
                sizes = [compatible_counts[i], incompatible_counts[i]]
                if sum(sizes) == 0:
                    sizes = [1, 0]  # Avoid empty pie
                colors = ['green', 'red']
                sub_ax.pie(sizes, colors=colors, autopct='%1.0f%%',
                           startangle=90, textprops={'fontsize': 8})
                sub_ax.set_title(f'Pos {pos}', fontsize=10)

        plt.tight_layout()
        plt.show()

        # Print summary using the compatibility report
        print("\n=== Compatibility Report ===")
        print(report)