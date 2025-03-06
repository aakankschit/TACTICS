import os

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