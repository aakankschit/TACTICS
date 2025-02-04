from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
import os
from multiprocessing import Pool, cpu_count

"""
This script is used to generate the library of molecules for the M_pro dataset.
The idea is to split the enumerated library into a set of smaller files containing 5000 molecules each.
In total it will generate 18561 files, each of which will become a job on hpc3.
TODO: This number will need to be modified to after running a filter to remove any building blocks that are not feasible.
"""

def get_bb_data(aldehydes_file_path, acids_file_path, coupled_aa_file_path):
    # Load aldehyde building blocks
    with open(aldehydes_file_path, 'r') as f1:
        aldehydes_data = f1.read().splitlines()
        aldehydes = [line.split(' ')[0] for line in aldehydes_data]
        aldehydes_names = [line.split(' ')[1] for line in aldehydes_data]
    
    # Load carboxylic acid building blocks
    with open(acids_file_path, 'r') as f2:
        acids_data = f2.read().splitlines()
        acids = [line.split(' ')[0] for line in acids_data]
        acids_names = [line.split(' ')[1] for line in acids_data]
    
    # Load coupled amino acid building blocks
    with open(coupled_aa_file_path, 'r') as f3:
        coupled_aa_data = f3.read().splitlines()
        coupled_aa = [line.split(' ')[0] for line in coupled_aa_data]
        coupled_aa_names = [line.split(' ')[1] for line in coupled_aa_data]

    return aldehydes, aldehydes_names, acids, acids_names, coupled_aa, coupled_aa_names

def process_combination(args):
    """Modified to accept single tuple argument for multiprocessing"""
    r1, r2, r3, n1, n2, n3, rxn = args
    try:
        products = rxn.RunReactants((r1, r2, r3))
        if products:
            smiles = Chem.MolToSmiles(products[0][0])
            reagent_no = f"{n1}_{n2}_{n3}"
            return smiles, reagent_no
    except:
        return None

def generate_combinations(coupled_aa_mols, aldehyde_mols, acid_mols, coupled_aa_names, aldehydes_names, acids_names, rxn):
    for n1, r1 in zip(coupled_aa_names, coupled_aa_mols):
        for n2, r2 in zip(aldehydes_names, aldehyde_mols):
            for n3, r3 in zip(acids_names, acid_mols):
                yield (r1, r2, r3, n1, n2, n3, rxn)

def main():
    # Get building block data
    aldehydes, aldehydes_names, acids, acids_names, coupled_aa, coupled_aa_names = get_bb_data(
        '/Users/aakankschitnandkeolyar/Desktop/TS_library_search/Data/M_pro/input_files/aldehydes_input_filtered.smi',
        '/Users/aakankschitnandkeolyar/Desktop/TS_library_search/Data/M_pro/input_files/acids_input_filtered.smi',
        '/Users/aakankschitnandkeolyar/Desktop/TS_library_search/Data/M_pro/input_files/coupled_aa_filtered.smi'
    )
    
    # Convert SMILES to molecules
    acid_mols = [Chem.MolFromSmiles(smiles) for smiles in acids]
    coupled_aa_mols = [Chem.MolFromSmiles(smiles) for smiles in coupled_aa]
    aldehyde_mols = [Chem.MolFromSmiles(smiles) for smiles in aldehydes]

    # Define reaction
    rxn = AllChem.ReactionFromSmarts('[#7H2;!$(NC=O):1].[#6:2](=O).[#6:3](=O)[O]>>[#6:2][#7:1][#6:3](=O)')
    
    # Setup output
    output_dir = './Data/M_pro/enumerated_library_files'
    os.makedirs(output_dir, exist_ok=True)
    file_index = 0
    molecule_count = 0
    
    # Process combinations using multiprocessing
    try:
        with Pool(processes=cpu_count()) as pool:
            # Open first output file
            current_file = open(os.path.join(output_dir, f'full_enum_lib_1.smi'), 'w')
            
            # Process combinations
            for result in tqdm(
                pool.imap_unordered(process_combination, generate_combinations(coupled_aa_mols, aldehyde_mols, acid_mols, coupled_aa_names, aldehydes_names, acids_names, rxn)),
                total=len(coupled_aa_mols) * len(aldehyde_mols) * len(acid_mols),
                desc="Processing combinations"
            ):
                if result:
                    smiles, reagent_no = result
                    
                    # Handle file writing
                    if molecule_count > 0 and molecule_count % 5000 == 0:
                        current_file.close()
                        file_index += 1
                        current_file = open(os.path.join(output_dir, f'full_enum_lib_{file_index+1}.smi'), 'w')
                    
                    current_file.write(f"{smiles} {reagent_no}\n")
                    molecule_count += 1
            
            # Close final file
            current_file.close()
    except KeyboardInterrupt:
        print("Process interrupted by user.")
        if not current_file.closed:
            current_file.close()

if __name__ == '__main__':
    main()