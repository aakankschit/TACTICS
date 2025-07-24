import polars as pl
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
from multiprocessing import Pool, cpu_count, get_context
from itertools import product
import dill as pickle
from collections import Counter
from typing import List, Tuple, Optional, Dict, Any
from .multiprocessing_utils import initializer

class LibraryEnumerator:
    def __init__(self, building_block_files: List[str]):
        """
        Initialize the LibraryEnumerator with building block files.
        
        Parameters:
        -----------
        building_block_files : List[str]
            List of paths to building block files containing SMILES and names
        """
        self.building_block_files = building_block_files
        self.building_blocks = None
        self.enumerated_products = None
        self._load_building_blocks()
    
    def _load_building_blocks(self) -> None:
        """Load building block data from files into polars DataFrames."""
        data_frames = []
        for file_path in self.building_block_files:
            data = []
            with open(file_path, 'r') as f:
                for line in f:
                    smiles, name = line.split(' ')
                    name = name.strip()
                    data.append({'SMILES': smiles, 'Name': name})
            df = pl.DataFrame(data)
            assert not df.is_empty(), f"Building block file {file_path} resulted in an empty DataFrame"
            data_frames.append(df)
        self.building_blocks = data_frames
    
    def _generate_reactant_combinations_parallel(self, smiles_column: str = 'SMILES', 
                                              name_column: str = 'Name') -> Tuple[List[Tuple], List[str]]:
        """Generate reactant combinations in parallel for library enumeration."""
        columns = [df[smiles_column].to_list() for df in self.building_blocks]
        names = [df[name_column].to_list() for df in self.building_blocks]

        for i, col in enumerate(columns):
            assert col, f"DataFrame {i} has no data in column '{smiles_column}'"

        num_cores = cpu_count()
        num_chunks = num_cores
        chunk_size = len(columns[0]) // num_chunks
        
        chunks = [columns[0][i:i + chunk_size] for i in range(0, len(columns[0]), chunk_size)]
        name_chunks = [names[0][i:i + chunk_size] for i in range(0, len(names[0]), chunk_size)]

        assert chunks, "Reactant combinations could not be generated"

        with get_context("spawn").Pool(num_cores) as pool:
            results = pool.starmap(product, [(chunk, *columns[1:]) for chunk in chunks])
            name_results = pool.starmap(product, [(name_chunk, *names[1:]) for name_chunk in name_chunks])

        reactant_combinations = [item for sublist in results for item in sublist]
        name_combinations = ['_'.join([''.join(name.split()) for name in item]) 
                           for sublist in name_results for item in sublist]

        assert reactant_combinations, "No reactant combinations were generated"
        print(f"Generated {len(reactant_combinations)} reactant combinations")

        return reactant_combinations, name_combinations
    
    def _apply_reaction(self, smarts: str, reactants: Tuple[str, ...], 
                       exception_rules: Optional[List[Tuple[str, int, str]]] = None) -> Optional[str]:
        """Apply reaction SMARTS pattern to reactants for library enumeration."""
        if exception_rules is None:
            exception_rules = []

        reactant_mols = tuple(Chem.MolFromSmiles(reactant) for reactant in reactants)

        for i, mol in enumerate(reactant_mols):
            assert mol is not None, f"Reactant {reactants[i]} could not be converted to an RDKit molecule"

        for i, reactant_mol in enumerate(reactant_mols):
            for substructure_smarts, position, exception_smarts in exception_rules:
                substructure = Chem.MolFromSmarts(substructure_smarts)
                if reactant_mol.HasSubstructMatch(substructure) and i == position:
                    smarts = exception_smarts
                    break

        rxn = AllChem.ReactionFromSmarts(smarts)
        
        if len(reactant_mols) == 1:
            products = rxn.RunReactants((reactant_mols[0],))
        else:
            products = rxn.RunReactants(reactant_mols)
        
        if products:
            return Chem.MolToSmiles(products[0][0], isomericSmiles=True)
        return None
    
    def enumerate_library(self, smarts: str, exception_rules: Optional[List[Tuple[str, int, str]]] = None) -> None:
        """
        Enumerate the chemical library using the reaction SMARTS pattern.
        
        Parameters:
        -----------
        smarts : str
            Reaction SMARTS pattern
        exception_rules : Optional[List[Tuple[str, int, str]]]
            List of exception rules (substructure_smarts, position, exception_smarts)
        """
        reactant_combinations, name_combinations = self._generate_reactant_combinations_parallel()
        
        theoretical_product_count = 1
        for df in self.building_blocks:
            theoretical_product_count *= len(df)
        
        assert theoretical_product_count > 0, "Theoretical product count is zero"
        
        args = [(smarts, reactants, exception_rules) for reactants in reactant_combinations]
        
        with get_context("spawn").Pool(cpu_count(), initializer=initializer) as pool:
            results = list(tqdm(pool.starmap(self._apply_reaction, args), total=len(args)))
        
        combined_data = [{'Product_SMILES': product, 'Product_Name': ' '.join(names)} 
                        for product, names in zip(results, name_combinations) if product]
        
        failed_reactants = [reactants for reactants, product in zip(reactant_combinations, results) if not product]
        
        print(f"Enumerated {len(combined_data)} products")
        assert len(combined_data) == theoretical_product_count, (
            f"Mismatch in product count: Expected {theoretical_product_count}, but got {len(combined_data)}"
        )
        
        if failed_reactants:
            distinct_failed_reactants = set(tuple(reactants) for reactants in failed_reactants)
            print(f"Failed to enumerate products for {len(distinct_failed_reactants)} distinct reactant combinations")
            print(f"Failed reactant combinations: {repr(list(distinct_failed_reactants))}")
            
            reactant_counter = Counter(reactant for reactants in distinct_failed_reactants 
                                     for reactant in reactants)
            repeated_failed_reactants = [reactant for reactant, count in reactant_counter.items() 
                                       if count > 1]
            print(f"failed reactants: {repr(repeated_failed_reactants)}")
        else:
            print("All reactants were successfully enumerated")
        
        self.enumerated_products = pl.DataFrame(combined_data)
        self.enumerated_products = self.enumerated_products.with_columns([
            pl.col('Product_Name').str.replace_all(' ', '')
        ])
    
    def get_product_smiles(self, product_name: str) -> Optional[str]:
        """
        Get the SMILES string for a specific product by name.
        
        Parameters:
        -----------
        product_name : str
            Name of the product to look up
            
        Returns:
        --------
        Optional[str]
            SMILES string if found, None otherwise
        """
        if self.enumerated_products is None:
            raise ValueError("No products have been enumerated yet")
            
        result = self.enumerated_products.filter(pl.col('Product_Name') == product_name)
        if result.is_empty():
            return None
        return result['Product_SMILES'][0]