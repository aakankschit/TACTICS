import json
from typing import List, Optional, Union

from pydantic import BaseModel

from .config import TSConfig, RWSConfig
from .reagent import Reagent, Standard_TS_Reagent, Enhanced_TS_Reagent


def load_config(json_filename: str) -> Union[TSConfig, RWSConfig]:
    """
    Loads, validates, and returns the configuration from a JSON file.

    :param json_filename: Path to the input JSON file.
    :return: A validated Pydantic config object (TSConfig or RWSConfig).
    """
    with open(json_filename, 'r') as f:
        data = json.load(f)

    if "num_ts_iterations" in data:
        return TSConfig(**data)
    else:
        return RWSConfig(**data)


def create_reagents(filename: str, ts_mode: str = "standard", num_to_select: Optional[int] = None) -> List:
    """
    Creates a list of Reagents from a file
    :param filename: a smiles file containing the reagents
    :param ts_mode: "standard" for Standard_TS_Reagent, "enhanced" for Enhanced_TS_Reagent
    :param num_to_select: For dev purposes; the number of molecules to return
    :return: List of Reagents
    """
    # Select reagent class based on mode
    if ts_mode == "standard":
        reagent_class = Standard_TS_Reagent
    elif ts_mode == "enhanced":
        reagent_class = Enhanced_TS_Reagent
    else:
        reagent_class = Reagent

    reagent_list = []
    with open(filename, "r") as f:
        for line in f.readlines():
            smiles, reagent_name = line.split()
            reagent = reagent_class(reagent_name=reagent_name, smiles=smiles)
            reagent_list.append(reagent)
    if num_to_select is not None and len(reagent_list) > num_to_select:
        reagent_list = reagent_list[:num_to_select]
    return reagent_list


def read_reagents(reagent_file_list, ts_mode: str = "standard", num_to_select: Optional[int] = None) -> List:
    """
    Read the reagents SMILES files
    :param reagent_file_list: a list of filenames containing reagents for the reaction. Each file list contains smiles
    strings for a single component of the reaction.
    :param ts_mode: "standard" for Standard_TS_Reagent, "enhanced" for Enhanced_TS_Reagent
    :param num_to_select: select how many reagents to read, mostly a development function
    :return: List of Reagents
    """
    reagents = []
    for reagent_filename in reagent_file_list:
        reagent_list = create_reagents(filename=reagent_filename, ts_mode=ts_mode, num_to_select=num_to_select)
        reagents.append(reagent_list)
    return reagents
