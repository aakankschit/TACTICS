from typing import List, Optional

from ..core.reagent import Reagent


def create_reagents(filename: str, num_to_select: Optional[int] = None, use_boltzmann_weighting: bool = False, mode: str = "maximize") -> List[Reagent]:
    """
    Creates a list of Reagents from a file.

    Parameters:
        filename: Path to a SMILES file containing reagents
        num_to_select: Optional limit on number of reagents to read
        use_boltzmann_weighting: If True, reagents will use Boltzmann-weighted Bayesian updates (legacy RWS).
                                If False, reagents will use standard uniform-weighted Bayesian updates (default).
        mode: "maximize" or "minimize" - affects Boltzmann weighting direction

    Returns:
        List of Reagent objects
    """
    reagent_list = []
    with open(filename, 'r') as f:
        for line in f.readlines():
            smiles, reagent_name = line.split()
            reagent = Reagent(reagent_name=reagent_name, smiles=smiles, use_boltzmann_weighting=use_boltzmann_weighting, mode=mode)
            reagent_list.append(reagent)

    if num_to_select is not None and len(reagent_list) > num_to_select:
        reagent_list = reagent_list[:num_to_select]

    return reagent_list


def read_reagents(
    reagent_file_list: List[str],
    num_to_select: Optional[int] = None,
    use_boltzmann_weighting: bool = False,
    mode: str = "maximize"
) -> List[List[Reagent]]:
    """
    Read reagents from multiple SMILES files.

    Parameters:
        reagent_file_list: List of file paths containing reagents for each reaction component
        num_to_select: Optional limit on number of reagents to read per file
        use_boltzmann_weighting: If True, reagents will use Boltzmann-weighted Bayesian updates (legacy RWS).
                                If False, reagents will use standard uniform-weighted Bayesian updates (default).
        mode: "maximize" or "minimize" - affects Boltzmann weighting direction

    Returns:
        List of lists of Reagent objects (one list per reaction component)

    """
    reagents = []
    for reagent_filename in reagent_file_list:
        reagent_list = create_reagents(
            filename=reagent_filename,
            num_to_select=num_to_select,
            use_boltzmann_weighting=use_boltzmann_weighting,
            mode=mode
        )
        reagents.append(reagent_list)

    return reagents
