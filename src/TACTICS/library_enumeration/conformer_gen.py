import os
import multiprocessing
from openeye import oechem, oeomega


def read_smi_file(file_path):
    """
    Reads a .smi file containing SMILES and product codes.

    Parameters:
        file_path (str): Path to the .smi file.

    Returns:
        list: A list of tuples (SMILES, product_code).
    """
    smiles_data = []
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 2:  # Ensure the line has both SMILES and product code
                smiles, product_code = parts
                smiles_data.append((smiles, product_code))
            else:
                print(f"Skipping invalid line: {line.strip()}")
    return smiles_data


def generate_conformers_for_molecule(smiles_product):
    """
    Generate conformers for a single molecule using OpenEye Omega.

    Parameters:
        smiles_product (tuple): A tuple (SMILES, product_code).

    Returns:
        tuple: (OEMol object or None, product_code, random_stereo_flag).
               - OEMol object if conformers are successfully generated.
               - None if conformer generation fails.
               - random_stereo_flag indicates whether random stereoisomer selection was used.
    """
    smiles, product_code = smiles_product
    mol = oechem.OEMol()
    if not oechem.OESmilesToMol(mol, smiles):
        print(f"Failed to parse SMILES: {smiles} for product code: {product_code}")
        return None, product_code, False

    mol.SetTitle(product_code)  # Set the molecule title to the product code

    omega = oeomega.OEOmega()
    omega.SetMaxConfs(2000)  # Set max conformers
    omega.SetStrictStereo(True)  # Enable strict stereo

    mol_copy = oechem.OEMol(mol)  # Create a copy of the molecule

    # Attempt to generate conformers with strict stereo
    if not omega(mol_copy):
        # If failed, allow random stereoisomer selection
        print(f"Stereo issue detected for product code: {product_code}. Generating random stereoisomer.")
        omega.SetStrictStereo(False)  # Disable strict stereo for this molecule
        if omega(mol_copy):
            return mol_copy, product_code, True  # Random stereoisomer used
        else:
            print(f"Failed to generate conformers even with relaxed stereo for: {product_code}.")
            return None, product_code, False

    return mol_copy, product_code, False


def generate_conformers_with_multiprocessing(smiles_data):
    """
    Generate conformers for products using multiprocessing.

    Parameters:
        smiles_data (list): List of tuples (SMILES, product_code).

    Returns:
        dict: A dictionary with keys:
              - "processed_products": List of tuples (OEMol object, product_code).
              - "random_stereo_products": List of tuples (SMILES, product_code) for which stereoisomers were randomly selected.
              - "failed_products": List of tuples (SMILES, product_code) that failed completely.
    """
    processed_products = []
    random_stereo_products = []
    failed_products = []

    # Use multiprocessing to parallelize conformer generation
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())  # Use all available CPU cores
    results = pool.map(generate_conformers_for_molecule, smiles_data)

    pool.close()
    pool.join()

    for mol_copy, product_code, random_stereo_flag in results:
        smiles_product = next((sp for sp in smiles_data if sp[1] == product_code), None)
        if mol_copy is not None:
            processed_products.append((mol_copy, product_code))
            if random_stereo_flag and smiles_product:
                random_stereo_products.append(smiles_product)
        else:
            if smiles_product:
                failed_products.append(smiles_product)

    return {
        "processed_products": processed_products,
        "random_stereo_products": random_stereo_products,
        "failed_products": failed_products,
    }


def write_conformers_to_oeb(output_file_path, processed_products):
    """
    Write generated conformers to a .oeb.gz file.

    Parameters:
        output_file_path (str): Path to the output .oeb.gz file.
        processed_products (list): List of tuples (OEMol object, product_code).
                                   Each OEMol object contains generated conformers.
    """
    with oechem.oemolostream(output_file_path) as ofs:
        for mol, _ in processed_products:
            if not oechem.OEWriteMolecule(ofs, mol):
                print(f"Failed to write molecule with title: {mol.GetTitle()}")


def write_random_stereo_to_smi(input_file_path, random_stereo_products):
    """
    Write SMILES and product codes for molecules where random stereoisomers were chosen to a .smi file.

    Parameters:
        input_file_path (str): Path to the original input .smi file.
        random_stereo_products (list): List of tuples (SMILES, product_code).
                                       Each tuple corresponds to a molecule where random stereoisomers were chosen.
    
    Returns:
        str: Path to the generated .smi file.
    """
    output_file_path = os.path.splitext(input_file_path)[0] + "_random_stereo.smi"
    
    with open(output_file_path, 'w') as ofs:
        for smiles, product_code in random_stereo_products:
            ofs.write(f"{smiles} {product_code}\n")
    
    return output_file_path


def process_all_smi_files(directory):
    """
    Process all .smi files in the given directory.

    Parameters:
        directory (str): Path to the directory containing .smi files.
    
    Returns:
        None
    """
    smi_files = [f for f in os.listdir(directory) if f.endswith(".smi")]

    for smi_file in smi_files:
        smi_file_path = os.path.join(directory, smi_file)
        
        # Read SMILES and product codes from the .smi file
        smiles_data = read_smi_file(smi_file_path)

        # Generate conformers using multiprocessing
        result = generate_conformers_with_multiprocessing(smiles_data)

        # Create output filenames based on input filename
        conformer_output_path = os.path.splitext(smi_file_path)[0] + "_conformers.oeb.gz"
        
        # Write generated conformers to an output .oeb.gz file
        write_conformers_to_oeb(conformer_output_path, result["processed_products"])
        
        print(f"Conformers written to: {conformer_output_path}")

        # Write SMILES and codes for molecules with random stereoisomers to a .smi file
        random_stereo_output_path = write_random_stereo_to_smi(smi_file_path, result["random_stereo_products"])
        
        print(f"Random stereoisomers written to: {random_stereo_output_path}")

        # Print summary statistics
        print(f"Processed {len(result['processed_products'])} molecules.")
        print(f"{len(result['random_stereo_products'])} molecules had random stereoisomers.")
        print(f"{len(result['failed_products'])} molecules failed.")


# Example usage
if __name__ == "__main__":
    directory_path = "."  # Specify your directory path here
    process_all_smi_files(directory_path)
