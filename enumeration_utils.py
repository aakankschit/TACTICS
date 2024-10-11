from rdkit import Chem

def check_substructure(mol_list, substructure_SMARTS):
	"""
	Check whether a SMARTS pattern of a substructure is present in a set of molecules and identify where the substructure is present and where it is not present.
	:param mol_list: A list of SMILES
	:param substructure_SMARTS: SMARTS pattern for recognizing a specific substructure
	:return: selected_mol_list: A list of SMILES of molecules containing the substructure.
	"""
	molecule_list = [Chem.MolFromSmiles(x) for x in mol_list]
	sub_structure = Chem.MolFromSmarts(substructure_SMARTS)
	selected_molecules_list = []
	unselected_molecules_list = []
	for molecule in molecule_list:
		if molecule.HasSubstructMatch(sub_structure):
			selected_molecules_list.append(Chem.MolToSmiles(molecule))
		else:
			unselected_molecules_list.append(Chem.MolToSmiles(molecule))
	return selected_molecules_list, unselected_molecules_list

