from collections import Counter
from typing import Union, List, Tuple, Dict, Literal, Optional
import pandas as pd
import polars as pl
from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import display
import os

class LibraryAnalysis:
    def __init__(self, df: Union[pd.DataFrame, pl.DataFrame],
                 smiles_files: Union[str, List[str]],
                 product_code_column: str = "Product_Code",
                 score_column: str = "Scores"):
        """
        Initialize the LibraryAnalysis class.
        
        Parameters:
        -----------
        df : Union[pandas.DataFrame, polars.DataFrame]
            DataFrame containing product codes and scores
        smiles_files : Union[str, List[str]]
            Path(s) to .smi file(s) containing SMILES and building block codes.
            Each file should have columns for SMILES and building block codes.
        product_code_column : str
            Name of the column containing product codes
        score_column : str
            Name of the column containing scores to sort by
            
        Raises:
        -------
        ValueError
            If smiles_files is not provided
            If smiles_files is provided but files don't exist
            If required columns are missing in the files
        """
        if smiles_files is None:
            raise ValueError("smiles_files parameter is required")
            
        self.df = pl.from_pandas(df) if isinstance(df, pd.DataFrame) else df
        self.product_code_column = product_code_column
        self.score_column = score_column
        
        # Read and combine SMILES from files
        self.smiles_dict = self._build_smiles_dictionary(smiles_files)
        
        self.position_counters = None
        self.total_molecules = None
        self.current_cutoff = None
        self.overlap_results = None

    def _build_smiles_dictionary(self, files: Union[str, List[str]]) -> Dict[str, str]:
        """
        Build a combined dictionary mapping building block codes to SMILES strings from .smi files.
        Assumes first column is SMILES and second column is building block code. 
        
        This is used to find the appropriate SMILES for each building block code.
        TODO: This wont handle isomeric SMILES since the product code is the same for such molecules. 
        
        Parameters:
        -----------
        files : Union[str, List[str]]
            Path(s) to .smi file(s) containing SMILES and building block codes
            
        Returns:
        --------
        Dict[str, str]
            Combined dictionary mapping building block codes to SMILES strings
            
        Raises:
        -------
        ValueError
            If any of the files don't exist
            If files are not properly formatted
            If duplicate building block codes are found
        """
        if isinstance(files, str):
            files = [files]
            
        # Check if all files exist
        for file in files:
            if not os.path.exists(file):
                raise ValueError(f"File not found: {file}")
                
        combined_dict = {}
        for file in files:
            # Read .smi file
            with open(file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:  # Skip empty lines
                        continue
                        
                    # Split line into SMILES and code
                    parts = line.split()
                    if len(parts) < 2:
                        raise ValueError(f"Invalid format in {file} at line {line_num}. Expected SMILES and code.")
                        
                    smiles = parts[0]
                    code = parts[1]
                    
                    # Check for duplicate codes
                    if code in combined_dict:
                        raise ValueError(f"Duplicate building block code '{code}' found in {file} at line {line_num}")
                        
                    # Add to combined dictionary
                    combined_dict[code] = smiles
            
        # Verify all building blocks in product codes have corresponding SMILES
        all_product_codes = set()
        for code in list(combined_dict.keys()):
            all_product_codes.update(code.split('_'))
            
        missing_codes = all_product_codes - set(combined_dict.keys())
        if missing_codes:
            raise ValueError(f"Missing SMILES for building block codes: {missing_codes}")
            
        return combined_dict

    def find_top_building_blocks(self,
                               cutoff: int = None,
                               sort_scores_by: Literal["maximize", "minimize"] = "maximize",
                               top_n: int = 20) -> Tuple[List[Counter], int]:
        """
        Find the top N most commonly occurring building blocks for each position in a list of product codes.
        
        Parameters:
        -----------
        cutoff : int, optional
            Number of top compounds to consider based on score. If None, uses all compounds.
        sort_scores_by : Literal["maximize", "minimize"]
            Whether to sort scores in descending ("maximize") or ascending ("minimize") order
        top_n : int
            Number of top building blocks to return for each position
            
        Returns:
        --------
        tuple
            A tuple containing:
            - List of Counters (one for each position)
            - Total number of molecules analyzed
        """
        # Sort the dataframe based on the score column
        if sort_scores_by == "maximize":
            df_sorted = self.df.sort(self.score_column, descending=True)
        else:
            df_sorted = self.df.sort(self.score_column, descending=False)
        
        # Apply cutoff if specified
        if cutoff is not None:
            df_sorted = df_sorted.head(cutoff)
            self.current_cutoff = cutoff
        
        # Extract the product codes as a list
        product_codes = df_sorted[self.product_code_column].to_list()
        self.total_molecules = len(product_codes)
        
        # Initialize counters for each position
        self.position_counters = []
        
        # Iterate through the product codes
        for product_code in product_codes:
            building_blocks = product_code.split("_")  # Split the product code by "_"
            # Ensure the position_counters list is large enough to handle all positions
            while len(self.position_counters) < len(building_blocks):
                self.position_counters.append(Counter())
            # Update the counters for each position
            for i, block in enumerate(building_blocks):
                self.position_counters[i][block] += 1
                
        return self.position_counters, self.total_molecules

    def _calculate_building_blocks_for_cutoff(self, cutoff: int, sort_scores_by: Literal["maximize", "minimize"] = "maximize") -> List[Counter]:
        """
        Calculate building blocks for a specific cutoff without updating the instance state.
        This is based on the top n building blocks, for example if the cutoff is 100, it will return the top building blocks for each position.

        Ensure that you use the appropriate sorting method when calling this function. For example if you want the top 10 building blocks for each position sorted by the maximum score, you should call this function with sort_scores_by="maximize".
        
        Parameters:
        -----------
        cutoff : int
            Number of top compounds to consider based on score
        sort_scores_by : Literal["maximize", "minimize"]
            Whether to sort scores in descending ("maximize") or ascending ("minimize") order
            
        Returns:
        --------
        List[Counter]
            List of Counters (one for each position) for the specified cutoff
        """
        # Sort the dataframe based on the score column
        if sort_scores_by == "maximize":
            df_sorted = self.df.sort(self.score_column, descending=True)
        else:
            df_sorted = self.df.sort(self.score_column, descending=False)
        
        # Apply cutoff
        df_sorted = df_sorted.head(cutoff)
        
        # Extract the product codes as a list
        product_codes = df_sorted[self.product_code_column].to_list()
        
        # Initialize counters for each position
        position_counters = []
        
        # Iterate through the product codes
        for product_code in product_codes:
            building_blocks = product_code.split("_")  # Split the product code by "_"
            # Ensure the position_counters list is large enough to handle all positions
            while len(position_counters) < len(building_blocks):
                position_counters.append(Counter())
            # Update the counters for each position
            for i, block in enumerate(building_blocks):
                position_counters[i][block] += 1
                
        return position_counters

    def check_overlap(self, new_cutoff: int, sort_scores_by: Literal["maximize", "minimize"] = "maximize") -> List[Tuple[int, int, set]]:
        """
        The goal here is to check if the top building blocks change when the cutoff is expanded. For example, if the top 10 building blocks for a cutoff of 100 are the same as the top 10 building blocks for a cutoff of 200, then the top building blocks have not changed. This function prints the number of overlapping building blocks at each position. And the total number of overlapping building blocks across all positions.
        
        Compare building blocks between the current cutoff and a new cutoff.
        
        Parameters:
        -----------
        new_cutoff : int
            New cutoff value to compare against the current cutoff
        sort_scores_by : Literal["maximize", "minimize"]
            Whether to sort scores in descending ("maximize") or ascending ("minimize") order
            
        Returns:
        --------
        List[Tuple[int, int, set]]
            List of tuples containing:
            - Position number (1-based)
            - Number of overlapping building blocks
            - Set of overlapping building blocks
        """
        if self.position_counters is None:
            raise ValueError("Must call find_top_building_blocks before checking overlap")
            
        if self.current_cutoff is None:
            raise ValueError("Current cutoff not set. Please run find_top_building_blocks with a cutoff first")
            
        # Get building blocks for the new cutoff without updating instance state
        new_counters = self._calculate_building_blocks_for_cutoff(
            cutoff=new_cutoff,
            sort_scores_by=sort_scores_by
        )
        
        # Get the top 20 building blocks for both cutoffs
        current_top_20 = [set(block for block, _ in counter.most_common(20)) 
                         for counter in self.position_counters]
        new_top_20 = [set(block for block, _ in counter.most_common(20)) 
                     for counter in new_counters]
        
        # Calculate overlap
        overlap_results = []
        total_overlap = 0
        
        for position, (current_blocks, new_blocks) in enumerate(zip(current_top_20, new_top_20)):
            overlap = current_blocks.intersection(new_blocks)
            overlap_results.append((position + 1, len(overlap), overlap))
            total_overlap += len(overlap)
            print(f"Position {position + 1}: {len(overlap)} building blocks overlap")
            
        print(f"Total number of overlapping compounds across all positions: {total_overlap}")
        self.overlap_results = overlap_results
        return overlap_results

    def compare_analysis_overlap(self, other_analysis: 'LibraryAnalysis', top_n: int = 20) -> List[Tuple[int, int, set]]:
        """
        Compare building blocks between two different LibraryAnalysis instances.
        
        Parameters:
        -----------
        other_analysis : LibraryAnalysis
            Another LibraryAnalysis instance to compare with
        top_n : int
            Number of top building blocks to consider for each position
            
        Returns:
        --------
        List[Tuple[int, int, set]]
            List of tuples containing:
            - Position number (1-based)
            - Number of overlapping building blocks
            - Set of overlapping building blocks
            
        Raises:
        -------
        ValueError
            If either analysis instance hasn't run find_top_building_blocks
        """
        if self.position_counters is None:
            raise ValueError("Current analysis instance must run find_top_building_blocks before comparison")
            
        if other_analysis.position_counters is None:
            raise ValueError("Other analysis instance must run find_top_building_blocks before comparison")
            
        # Get the top N building blocks for both analyses
        current_top_n = [set(block for block, _ in counter.most_common(top_n)) 
                        for counter in self.position_counters]
        other_top_n = [set(block for block, _ in counter.most_common(top_n)) 
                      for counter in other_analysis.position_counters]
        
        # Calculate overlap
        overlap_results = []
        total_overlap = 0
        
        for position, (current_blocks, other_blocks) in enumerate(zip(current_top_n, other_top_n)):
            overlap = current_blocks.intersection(other_blocks)
            overlap_results.append((position + 1, len(overlap), overlap))
            total_overlap += len(overlap)
            print(f"Position {position + 1}: {len(overlap)} building blocks overlap")
            
        print(f"Total number of overlapping compounds across all positions: {total_overlap}")
        return overlap_results 