import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdDeprotect
from rdkit.Chem.rdDeprotect import DeprotectData

sys.path.append('/Users/aakankschitnandkeolyar/Desktop/PRISMS')  # Adjust the path to point to the project root directory
from PRISMS.library_enumeration.generate_products import *
from PRISMS.library_enumeration.enumeration_utils import *

input_files = [['../Data/Thrombin/input_files/acids.smi',
                '../Data/Thrombin/input_files/coupled_aa_sub.smi']]

output_dir = '../Data/Thrombin/output_files'

def tbc_deprotect(mol):
    """
    Deprotect a molecule using the tert-butyloxycarbonyl (tBoc) protecting group.
    :param mol: A molecule object
    :return: A molecule object
    """
    reaction_class = "amine"
    reaction_smarts = "[#6:1][#7X3;H0,H1:2][C](=O)OC([C])([C])[C]>>[#6:1][#7:2]"
    abbreviation = "tbc"
    full_name = "tert-butyloxycarbonyl"

    deprotect_data = DeprotectData(reaction_class, reaction_smarts, abbreviation, full_name)
    deprotected_mol = rdDeprotect.Deprotect(mol, [deprotect_data])
    return deprotected_mol


def pbf_deprotect(mol):
    """
    Deprotect a molecule using the pentamethyldihydrobenzofuran-5-sulfonyl (Pbf) protecting group.
    :param mol: A molecule object
    :return: A molecule object
    """
    reaction_class = "amine"
    reaction_smarts = "[#6:1](=N)[#7X3;H0,H1:2]S(=O)(=O)c1c(c2c(c(c1[C])[C])OC([C]2)([C])[C])[C]>>[#6:1](=N)[#7X3:2]"
    abbreviation = "pbf"
    full_name = "pentamethyldihydrobenzofuran-5-sulfonyl"

    deprotect_data = DeprotectData(reaction_class, reaction_smarts, abbreviation, full_name)
    deprotected_mol = rdDeprotect.Deprotect(mol, [deprotect_data])
    return deprotected_mol

def tpg_deprotect(mol):
    """
    Deprotect a molecule using the trityl protecting group.
    :param mol: A molecule object
    :return: A molecule object
    """
    reaction_class = "amine"
    reaction_smarts = "[#6:1][#7X3;H0,H1:2]C(c1[c][c][c][c][c]1)(c2[c][c][c][c][c]2)c3[c][c][c][c][c]3>>[#6:1][#7:2][H]"
    abbreviation = "tpg"
    full_name = "trityl"

    deprotect_data = DeprotectData(reaction_class, reaction_smarts, abbreviation, full_name)
    deprotected_mol = rdDeprotect.Deprotect(mol, [deprotect_data])
    return deprotected_mol

def tb_deprotect(mol):
    """
    Deprotect a molecule using the tert-butyl protecting group.
    :param mol: A molecule object
    :return: A molecule object
    """
    reaction_class = "hydroxyl"
    reaction_smarts = "[#8:1]C([C])([C])[C]>>[#8:1]"
    abbreviation = "tb"
    full_name = "tert-Butyl"

    deprotect_data = DeprotectData(reaction_class, reaction_smarts, abbreviation, full_name)
    deprotected_mol = rdDeprotect.Deprotect(mol, [deprotect_data])
    return deprotected_mol

def main():
    # Load Data
    df_input_data = get_bb_data(input_files)

    # Reaction SMARTS
    amidation = '[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]>>[#6:1](=[O:2])[#7:3]'

    # Enumerate Library
    enum_library = combine_building_blocks(df_input_data, smarts=amidation)

    # Deprotect Library
    deprotected_prods = []
    for prod_smiles in enum_library['Product_SMILES']:
        mol = Chem.MolFromSmiles(prod_smiles)
        # Check if it has a pbf protecting group
        if mol.HasSubstructMatch(Chem.MolFromSmarts('S(=O)(=O)c1c(c2c(c(c1[C])[C])OC([C]2)([C])[C])[C]')):
            deprotected_mol = pbf_deprotect(mol)
            deprotected_prods.append(Chem.MolToSmiles(deprotected_mol))
        # Check if it has the tbc protecting group
        elif mol.HasSubstructMatch(Chem.MolFromSmarts('[#7][C](=O)OC([C])([C])[C]')):
            deprotected_mol = tbc_deprotect(mol)
            deprotected_prods.append(Chem.MolToSmiles(deprotected_mol))
        # Check if it has the tpg protecting group
        elif mol.HasSubstructMatch(Chem.MolFromSmarts('C(c1[c][c][c][c][c]1)(c2[c][c][c][c][c]2)c3[c][c][c][c][c]3')):
            deprotected_mol = tpg_deprotect(mol)
            deprotected_prods.append(Chem.MolToSmiles(deprotected_mol))
        # Check if it has the tb protecting group
        elif mol.HasSubstructMatch(Chem.MolFromSmarts('[#8]C([C])([C])[C]')):
            deprotected_mol = tb_deprotect(mol)
            deprotected_prods.append(Chem.MolToSmiles(deprotected_mol))
        else:
            deprotected_prods.append(Chem.MolToSmiles(mol))

    # Update enumerated library with deprotected products
    enum_deprotect_df = enum_library.clone()
    enum_deprotect_df = enum_deprotect_df.with_columns([pl.Series("Product_SMILES", deprotected_prods)])

    # Write enumerated library to files
    write_products_to_files(enum_deprotect_df, output_dir)

