import os
import polars as pl

def compile_product_scores(directory):
    """
    Read multiple text files that have "products" in their names followed by a number.
    Each file contains two columns: "Product_Code" and "Scores".
    Compile the data from all these files into a Polars DataFrame and return it.
    
    :param directory: The directory containing the product files.
    :return: A Polars DataFrame containing the compiled data.
    """
    # Initialize an empty list to store the data frames
    data_frames = []

    # Iterate over all files in the directory
    for filename in os.listdir(directory):
        if filename.startswith("products") and filename.endswith(".txt"):
            file_path = os.path.join(directory, filename)
            # Read the file into a Polars DataFrame with specified column names
            df = pl.read_csv(file_path, separator="\t", has_header=True, new_columns=["Product_Code", "Scores"])
            # Append the DataFrame to the list
            data_frames.append(df)

    # Concatenate all the DataFrames into a single DataFrame
    compiled_df = pl.concat(data_frames, how='vertical')
    
    return compiled_df

def compile_product_smiles(directory):
    """
    Read multiple .smi files from a directory that have "products" in their names.
    Each file contains product codes and product SMILES without a header.
    Return a dictionary where the product codes are the keys and the product SMILES are the values.
    
    :param directory: The directory containing the .smi files.
    :return: A dictionary with product codes as keys and product SMILES as values.
    """
    product_dict = {}

    # Iterate over all files in the directory
    for filename in os.listdir(directory):
        if "products" in filename and filename.endswith(".smi"):
            file_path = os.path.join(directory, filename)
            # Read the file line by line
            with open(file_path, 'r') as file:
                for line in file:
                    product_smiles, product_code = line.strip().split()
                    product_dict[product_code] = product_smiles
    
    return product_dict