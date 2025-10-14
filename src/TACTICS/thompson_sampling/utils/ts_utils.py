from typing import List, Optional

from ..legacy.reagent import Standard_TS_Reagent, Enhanced_TS_Reagent


def create_reagents(filename: str, ts_mode: str, num_to_select: Optional[int] = None) -> List[Standard_TS_Reagent | Enhanced_TS_Reagent]:
    """
    Creates a list of Reagents from a file
    :param filename: a smiles file containing the reagents
    :param num_to_select: For dev purposes; the number of molecules to return
    :return: List of Reagents
    """
    reagent_list = []
    with open(filename, 'r') as f:
        for line in f.readlines():
            smiles, reagent_name = line.split()
            if ts_mode == "standard":
                reagent = Standard_TS_Reagent(reagent_name=reagent_name, smiles=smiles)
            elif ts_mode == "enhanced":
                reagent = Enhanced_TS_Reagent(reagent_name=reagent_name, smiles=smiles)
            else:
                raise ValueError(f"Invalid ts_mode: {ts_mode}")
            reagent_list.append(reagent)
    if num_to_select is not None and len(reagent_list) > num_to_select:
        reagent_list = reagent_list[:num_to_select]
    return reagent_list


def read_reagents(reagent_file_list, ts_mode: str, num_to_select: Optional[int]) -> List[Standard_TS_Reagent | Enhanced_TS_Reagent]:
    """
    Read the reagents SMILES files
    :param reagent_file_list: a list of filenames containing reagents for the reaction. Each file list contains smiles
    strings for a single component of the reaction.
    :param ts_mode: "standard" or "enhanced"
    :param num_to_select: select how many reagents to read, mostly a development function
    :return: List of Reagents
    """
    reagents = []
    for reagent_filename in reagent_file_list:
        reagent_list = create_reagents(filename=reagent_filename, ts_mode=ts_mode, num_to_select=num_to_select)
        reagents.append(reagent_list)
    return reagents
