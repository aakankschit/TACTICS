import os
import polars as pl
from multiprocessing import Pool, cpu_count, get_context

def find_reactants_from_product_code(product_df, reactant_df, product_smiles_dict):
    """
    Find all reactants from a given Polars DataFrame based on the product code.
    Generates a new DataFrame with columns for product SMILES, score, and subsequent columns for each reactant used in the product.
    
    :param product_df: A Polars DataFrame containing 'Product_Code' and 'Score' columns.
    :param reactant_df: A Polars DataFrame containing reactant data.
    :param product_smiles_dict: A dictionary mapping product codes to SMILES.
    :return: A Polars DataFrame with columns for product SMILES, score, and subsequent columns for each reactant used in the product.
    """
    reactant_columns = [f'Reactant_{i+1}' for i in range(len(product_df[0, 'Product_Code'].split('_')))]
    
    def process_product_code(row):
        product_code = row['Product_Code']
        product_smiles = product_smiles_dict.get(product_code, None)
        score = row['Score']
        reactant_codes = product_code.split('_')
        reactant_data = {col: [] for col in reactant_columns}
        reactant_data['Product_SMILES'] = [product_smiles]
        reactant_data['Score'] = [score]

        for code in reactant_codes:
            reactant_smiles = reactant_df.filter(pl.col('Name') == code)['SMILES'].to_list()
            if reactant_smiles:
                reactant_data[f'Reactant_{reactant_codes.index(code)+1}'].append(reactant_smiles[0])
            else:
                reactant_data[f'Reactant_{reactant_codes.index(code)+1}'].append(None)

        return reactant_data

    # Use multiprocessing to process product codes in parallel
    num_cores = cpu_count()
    with get_context("spawn").Pool(num_cores) as pool:
        results = pool.map(process_product_code, product_df.rows())

    # Combine results into a single dictionary
    combined_data = {col: [] for col in reactant_columns}
    combined_data['Product_SMILES'] = []
    combined_data['Score'] = []
    for result in results:
        for key, value in result.items():
            combined_data[key].extend(value)

    # Convert the dictionary to a Polars DataFrame
    result_df = pl.DataFrame(combined_data)
    
    return result_df

def write_products_to_files(df, output_dir, products_per_file=5000):
    """
    Write the combined DataFrame of products to text files with a specified number of products per file.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Determine the number of files needed
    num_files = (len(df) + products_per_file - 1) // products_per_file

    # Split the DataFrame into chunks and write each chunk to a separate file
    for i in range(num_files):
        start_idx = i * products_per_file
        end_idx = min((i + 1) * products_per_file, len(df))
        chunk = df[start_idx:end_idx]

        # Write the chunk to a temporary .csv file without the header
        temp_output_file = os.path.join(output_dir, f'temp_products_{i + 1}.csv')
        chunk.write_csv(temp_output_file, include_header=False)

        # Read the temporary .csv file and replace commas with tabs
        with open(temp_output_file, 'r') as temp_file:
            content = temp_file.read().replace(',', ' ')

        # Write the modified content to the final .smi file
        output_file = os.path.join(output_dir, f'products_{i + 1}.smi')
        with open(output_file, 'w') as final_file:
            final_file.write(content)

        # Remove the temporary .csv file
        os.remove(temp_output_file)