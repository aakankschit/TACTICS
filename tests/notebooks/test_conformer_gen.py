import os
import argparse
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


def generate_conformers_for_molecule(smiles_product, max_conformers=200):
    """
    Generate conformers for a single molecule using OpenEye Omega.

    Parameters:
        smiles_product (tuple): A tuple (SMILES, product_code).
        max_conformers (int): Maximum number of conformers to generate. Default is 200.

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
    oechem.OESetSDData(mol, "product_code", product_code)  # Set product code as SD data

    omega = oeomega.OEOmega()
    omega.SetMaxConfs(max_conformers)  # Set maximum number of conformers
    omega.SetStrictStereo(True)  # Enable strict stereo

    # Attempt to generate conformers with strict stereo
    if not omega(mol):
        print(f"Stereo issue detected for product code: {product_code}. Generating random stereoisomer.")
        omega.SetStrictStereo(False)  # Disable strict stereo for this molecule
        if omega(mol):
            print(f"Conformers successfully generated with relaxed stereo for product code: {product_code}.")
            return mol, product_code, True  # Random stereoisomer used
        else:
            print(f"Failed to generate conformers even with relaxed stereo for: {product_code}.")
            return None, product_code, False

    print(f"Conformers successfully generated with strict stereo for product code: {product_code}.")
    return mol, product_code, False


def write_conformers_to_oeb_gz(output_file_path, smiles_data, max_conformers=200):
    """
    Write generated conformers directly to a .oeb.gz file.

    Parameters:
        output_file_path (str): Path to the output .oeb.gz file.
        smiles_data (list): List of tuples (SMILES, product_code).
                            Each tuple corresponds to a molecule whose conformers need to be generated.
        max_conformers (int): Maximum number of conformers to generate per molecule. Default is 200.
    
    Returns:
        None
    """
    processed_count = 0
    failed_count = 0

    # Open an oemolostream with .oeb.gz format
    ofs = oechem.oemolostream(output_file_path)
    
    if not ofs.IsValid():
        print(f"Error: Could not open output file '{output_file_path}' for writing.")
        exit(1)

    # Process each molecule sequentially and write directly to the output file
    for smiles_product in smiles_data:
        mol, product_code, random_stereo_flag = generate_conformers_for_molecule(smiles_product, max_conformers)
        if mol is not None:
            if oechem.OEWriteMolecule(ofs, mol):
                processed_count += 1
                if random_stereo_flag:
                    print(f"Random stereoisomer was chosen for product code: {product_code}.")
            else:
                print(f"Failed to write molecule with title: {mol.GetTitle()}")
        else:
            failed_count += 1

    ofs.close()

    print(f"\nSummary:")
    print(f" - Processed molecules: {processed_count}")
    print(f" - Failed molecules: {failed_count}")
    print(f"\nConformers written directly to compressed file: {output_file_path}")


if __name__ == "__main__":
    
    # parser = argparse.ArgumentParser(description="Process a .smi file to generate conformers.")
    
    # parser.add_argument("smi_file", help="Path to the .smi file containing SMILES and Product Codes.")
    
    # parser.add_argument(
    #     "--max_conformers",
    #     type=int,
    #     default=2000,
    #     help="Maximum number of conformers to generate per molecule. Default is 2000.",
    # )
    
    # args = parser.parse_args()
    
    smi_file_path = '/Users/aakankschitnandkeolyar/Desktop/TS_Chem_Space/Thrombin/Linear_amide/test.smi'
    max_conformers = 10

    # Validate input file exists
    if not os.path.isfile(smi_file_path):
        print(f"Error: File '{smi_file_path}' does not exist.")
        exit(1)

    # Read SMILES and Product Codes from the .smi file
    smiles_data = read_smi_file(smi_file_path)

    # Generate Output File Name based on input filename
    oeb_gz_output_path = os.path.splitext(smi_file_path)[0] + "_conformers.oeb.gz"

    # Write conformers directly to .oeb.gz file
    write_conformers_to_oeb_gz(oeb_gz_output_path, smiles_data, max_conformers)
