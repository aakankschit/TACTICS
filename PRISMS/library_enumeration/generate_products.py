import polars as pl
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
from multiprocessing import Pool, cpu_count, get_context
from itertools import product
import dill as pickle
from collections import Counter

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
                name = name.strip()
                data.append({'SMILES': smiles, 'Name': name})
        df = pl.DataFrame(data)
        data_frames.append(df)
    
    return data_frames

def generate_reactant_combinations_parallel(data_frames, smiles_column='SMILES', name_column='Name'):
    """
    Generate reactant combinations from a list of building block data frames in parallel.
    Returns a list of tuples, each containing a combination of values from the specified column.
    """
    # Extract the values from the specified column in each DataFrame
    columns = [df[smiles_column].to_list() for df in data_frames]
    names = [df[name_column].to_list() for df in data_frames]

    # Determine the number of CPU cores
    num_cores = cpu_count()

    # Set the number of chunks to the number of CPU cores
    num_chunks = num_cores

    # Split the first column into chunks
    chunk_size = len(columns[0]) // num_chunks
    chunks = [columns[0][i:i + chunk_size] for i in range(0, len(columns[0]), chunk_size)]
    name_chunks = [names[0][i:i + chunk_size] for i in range(0, len(names[0]), chunk_size)]

    # Generate combinations in parallel
    with get_context("spawn").Pool(num_cores) as pool:
        results = pool.starmap(product, [(chunk, *columns[1:]) for chunk in chunks])
        name_results = pool.starmap(product, [(name_chunk, *names[1:]) for name_chunk in name_chunks])

    # Flatten the list of results
    reactant_combinations = [item for sublist in results for item in sublist]
    name_combinations = ['_'.join([''.join(name.split()) for name in item]) for sublist in name_results for item in sublist]

    # Print the number of reactant combinations
    print(f"Generated {len(reactant_combinations)} reactant combinations")

    return reactant_combinations, name_combinations

def apply_reaction(smarts, reactants, exception_rules=None):
    """
    Apply the reaction SMARTS pattern to a set of reactants and return the product SMILES string.
    Handles exception molecules by using specific reaction SMARTS patterns based on the position of the reactant.
    """
    if exception_rules is None:
        exception_rules = []

    reactant_mols = tuple(Chem.MolFromSmiles(reactant) for reactant in reactants)

    # Check for exception molecules and use the corresponding reaction SMARTS pattern
    for i, reactant_mol in enumerate(reactant_mols):
        for substructure_smarts, position, exception_smarts in exception_rules:
            substructure = Chem.MolFromSmarts(substructure_smarts)
            if reactant_mol.HasSubstructMatch(substructure) and i == position:
                smarts = exception_smarts
                break

    rxn = AllChem.ReactionFromSmarts(smarts)
    
    # Apply the main reaction
    if len(reactant_mols) == 1:
        products = rxn.RunReactants((reactant_mols[0],))
    else:
        products = rxn.RunReactants(reactant_mols)
    
    if products:
        return Chem.MolToSmiles(products[0][0])
    return None

def apply_reaction_wrapper(smarts, reactants, exception_rules=None):
    """
    Wrapper function to apply the reaction.
    """
    return apply_reaction(smarts, reactants, exception_rules)

def initializer():
    """
    Initializer function for the multiprocessing pool.
    """
    pickle.settings.update({'recurse': True})

def combine_building_blocks(data_frames, smarts, smiles_column='SMILES', name_column='Name', exception_rules=None):
    """
    Combine building blocks using a reaction SMARTS pattern.
    Uses multiprocessing to parallelize the combination process.
    Handles exception molecules by using specific reaction SMARTS patterns based on the position of the reactant.
    """
    # Generate reactant combinations in parallel
    reactant_combinations, name_combinations = generate_reactant_combinations_parallel(data_frames, smiles_column, name_column)

    # Prepare arguments for apply_reaction_wrapper
    args = [(smarts, reactants, exception_rules) for reactants in reactant_combinations]

    # Initialize a multiprocessing pool with the number of available CPU cores
    with get_context("spawn").Pool(cpu_count(), initializer=initializer) as pool:
        # Use tqdm to display a progress bar and pool.starmap to apply the reaction in parallel
        results = list(tqdm(pool.starmap(apply_reaction_wrapper, args), total=len(args)))

    # Collect the results into a list of dictionaries
    combined_data = [{'Product_SMILES': product, 'Product_Name': ' '.join(names)} for product, names in zip(results, name_combinations) if product]

    # Collect the reactants that could not be used to create a product
    failed_reactants = [reactants for reactants, product in zip(reactant_combinations, results) if not product]

    # Print the number of products generated
    print(f"Generated {len(combined_data)} products")

    # Print the distinct reactants that could not be used to create a product
    if failed_reactants:
        distinct_failed_reactants = set(tuple(reactants) for reactants in failed_reactants)
        print(f"Failed to generate products for {len(distinct_failed_reactants)} distinct reactant combinations")
        print(f"Failed reactant combinations: {repr(list(distinct_failed_reactants))}")
        
        # Find and print reactants that occur more than once in the failed combinations
        reactant_counter = Counter(reactant for reactants in distinct_failed_reactants for reactant in reactants)
        repeated_failed_reactants = [reactant for reactant, count in reactant_counter.items() if count > 1]
        print(f"failed reactants: {repr(repeated_failed_reactants)}")
    else:
        print("All reactants were used to generate products")

    # Convert the list of dictionaries into a polars DataFrame
    combined_df = pl.DataFrame(combined_data)

    # Reformat the 'Product_Name' column to remove spaces
    combined_df = combined_df.with_columns([
        pl.col('Product_Name').str.replace_all(' ', '')
    ])

    return combined_df