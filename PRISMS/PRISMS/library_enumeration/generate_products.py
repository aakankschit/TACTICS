from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
from multiprocessing import cpu_count, get_context
import polars as pl
from itertools import product

def get_bb_data(file_path_list):
    """
    Load building block data from a list of file paths.
    Each file should contain a list of SMILES strings and the names of the building blocks.
    Returns a list of polars data frames, each corresponding to a building block file.
    """
    data_frames = []
    for file_path in file_path_list:
        data = []
        with open(file_path, 'r') as f:
            for line in f:
                smiles, name = line.split(' ')
                data.append({'SMILES': smiles, 'Name': name})
        df = pl.DataFrame(data)
        data_frames.append(df)
    
    return data_frames

def generate_reactant_combinations_parallel(data_frames, column_name='SMILES'):
    """
    Generate reactant combinations from a list of building block data frames in parallel.
    Returns a list of tuples, each containing a combination of values from the specified column.
    """
    # Extract the values from the specified column in each DataFrame
    columns = [df[column_name].to_list() for df in data_frames]

    # Determine the number of CPU cores
    num_cores = cpu_count()

    # Set the number of chunks to the number of CPU cores
    num_chunks = num_cores

    # Split the first column into chunks
    chunk_size = len(columns[0]) // num_chunks
    chunks = [columns[0][i:i + chunk_size] for i in range(0, len(columns[0]), chunk_size)]

    # Generate combinations in parallel
    with get_context("fork").Pool(num_cores) as pool:
        results = pool.starmap(product, [(chunk, *columns[1:]) for chunk in chunks])

    # Flatten the list of results
    reactant_combinations = [item for sublist in results for item in sublist]
    return reactant_combinations

def apply_reaction(smarts, *reactants):
    """
    Apply the reaction SMARTS pattern to a set of reactants and return the product SMILES string.
    """
    rxn = AllChem.ReactionFromSmarts(smarts)
    reactant_mols = tuple(Chem.MolFromSmiles(reactant) for reactant in reactants)
    products = rxn.RunReactants(reactant_mols)
    Chem.SanitizeMol(products)
    if products:
        return Chem.MolToSmiles(products[0][0])
    return None

def apply_reaction_wrapper(smarts, reactants):
    """
    Wrapper function to apply the reaction.
    """
    return apply_reaction(smarts, *reactants)

def combine_building_blocks(data_frames, smarts, column_name='SMILES'):
    """
    Combine building blocks using a reaction SMARTS pattern.
    Uses multiprocessing to parallelize the combination process.
    """
    from .data_processing import generate_reactant_combinations_parallel

    # Generate reactant combinations in parallel
    reactant_combinations = generate_reactant_combinations_parallel(data_frames, column_name)

    # Prepare arguments for apply_reaction_wrapper
    args = [(smarts, reactants) for reactants in reactant_combinations]

    # Initialize a multiprocessing pool with the number of available CPU cores
    with get_context("fork").Pool(cpu_count()) as pool:
        # Use tqdm to display a progress bar and pool.starmap to apply the reaction in parallel
        results = list(tqdm(pool.starmap(apply_reaction_wrapper, args), total=len(args)))

    # Collect the results into a list of dictionaries
    combined_data = [{'Product_SMILES': product} for product in results if product]

    # Convert the list of dictionaries into a polars DataFrame
    combined_df = pl.DataFrame(combined_data)
    return combined_df