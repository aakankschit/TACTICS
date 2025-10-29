from typing import List, Optional

from ..core.reagent import Reagent


def create_reagents(filename: str, num_to_select: Optional[int] = None) -> List[Reagent]:
    """
    Creates a list of Reagents from a file.

    Parameters:
        filename: Path to a SMILES file containing reagents
        num_to_select: Optional limit on number of reagents to read

    Returns:
        List of Reagent objects
    """
    reagent_list = []
    with open(filename, 'r') as f:
        for line in f.readlines():
            smiles, reagent_name = line.split()
            reagent = Reagent(reagent_name=reagent_name, smiles=smiles)
            reagent_list.append(reagent)

    if num_to_select is not None and len(reagent_list) > num_to_select:
        reagent_list = reagent_list[:num_to_select]

    return reagent_list


def read_reagents(reagent_file_list: List[str], num_to_select: Optional[int] = None) -> List[List[Reagent]]:
    """
    Read reagents from multiple SMILES files.

    Parameters:
        reagent_file_list: List of file paths containing reagents for each reaction component
        num_to_select: Optional limit on number of reagents to read per file

    Returns:
        List of lists of Reagent objects (one list per reaction component)
    """
    reagents = []
    for reagent_filename in reagent_file_list:
        reagent_list = create_reagents(filename=reagent_filename, num_to_select=num_to_select)
        reagents.append(reagent_list)
    return reagents
