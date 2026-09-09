"""
Product generation utilities for chemical library enumeration.

This module provides functions for:
- Single product enumeration from reagents
- Batch enumeration (sequential and parallel)
- Reaction execution helpers
- Deprotection application

These functions are used by SynthesisPipeline for library enumeration.
"""

import logging
from itertools import product as itertools_product
from multiprocessing import get_context
from typing import Any, Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem

from .enumeration_utils import EnumerationError, EnumerationResult

logger = logging.getLogger(__name__)


# === Single Reaction Execution ===


def run_single_reaction(
    rdkit_reaction: Any,
    reagent_mols: List[Chem.Mol],
) -> Tuple[Optional[Chem.Mol], Optional[str]]:
    """
    Run a single RDKit reaction on reagent molecules.

    Args:
        rdkit_reaction: RDKit ChemicalReaction object
        reagent_mols: List of reagent Mol objects

    Returns:
        Tuple of (product_mol, product_smiles) or (None, None) if failed
    """
    try:
        products = rdkit_reaction.RunReactants(tuple(reagent_mols))

        if not products or len(products) == 0 or len(products[0]) == 0:
            return None, None

        product = products[0][0]
        Chem.SanitizeMol(product)
        smiles = Chem.MolToSmiles(product, isomericSmiles=True)

        return product, smiles

    except Exception as e:
        logger.debug(f"Reaction failed: {e}")
        return None, None


def apply_deprotection(
    mol: Chem.Mol,
    group_name: str,
    protecting_groups: Dict[str, Any],
    max_iterations: int = 10,
) -> Optional[Chem.Mol]:
    """
    Apply deprotection to remove a protecting group from a molecule.

    Args:
        mol: RDKit Mol to deprotect
        group_name: Name of protecting group (e.g., "Boc", "Fmoc")
        protecting_groups: Dict of {name: ProtectingGroupInfo}
        max_iterations: Max deprotection iterations

    Returns:
        Deprotected molecule or None if failed
    """
    if group_name not in protecting_groups:
        logger.warning(f"Unknown protecting group: {group_name}")
        return mol

    pg = protecting_groups[group_name]
    if pg.deprotection_smarts is None:
        logger.warning(f"No deprotection SMARTS for {group_name}")
        return mol

    rxn = AllChem.ReactionFromSmarts(pg.deprotection_smarts)
    pattern = Chem.MolFromSmarts(pg.smarts)

    if rxn is None or pattern is None:
        return mol

    result = mol
    for _ in range(max_iterations):
        if not result.HasSubstructMatch(pattern):
            break

        products = rxn.RunReactants((result,))
        if products and len(products) > 0 and len(products[0]) > 0:
            result = products[0][0]
            try:
                Chem.SanitizeMol(result)
            except Exception:
                break
        else:
            break

    return result


def apply_deprotections(
    mols: List[Chem.Mol],
    deprotection_specs: List[Any],
    protecting_groups: Dict[str, Any],
    reactant_only: bool = True,
) -> List[Chem.Mol]:
    """
    Apply reactant deprotections to a list of molecules.

    This function only processes deprotections that target reactants (integer targets).
    Product deprotections (target="product") are handled separately after the reaction.

    Args:
        mols: List of reagent molecules
        deprotection_specs: List of DeprotectionSpec objects
        protecting_groups: Dict of protecting group definitions
        reactant_only: If True (default), only apply reactant deprotections (skip product targets)

    Returns:
        List of molecules with deprotections applied
    """
    result = list(mols)

    for deprot in deprotection_specs:
        # Skip product deprotections if reactant_only=True
        if reactant_only and deprot.is_product_deprotection:
            continue

        position = deprot.reactant_index
        if position is None or position >= len(result):
            continue

        mol = result[position]
        if mol is None:
            continue

        result[position] = apply_deprotection(mol, deprot.group, protecting_groups)

    return result


def apply_product_deprotections(
    product: Chem.Mol,
    deprotection_specs: List[Any],
    protecting_groups: Dict[str, Any],
) -> Chem.Mol:
    """
    Apply product deprotections to a molecule.

    This function only processes deprotections that target the product (target="product").

    Args:
        product: The product molecule to deprotect
        deprotection_specs: List of DeprotectionSpec objects
        protecting_groups: Dict of protecting group definitions

    Returns:
        Deprotected product molecule
    """
    result = product

    for deprot in deprotection_specs:
        if not deprot.is_product_deprotection:
            continue

        if result is None:
            break

        result = apply_deprotection(result, deprot.group, protecting_groups)

    return result


# === Batch Enumeration ===


def enumerate_products(
    pipeline: "SynthesisPipeline",
    combinations: List[Tuple[List[Chem.Mol], Optional[List[str]]]],
    n_jobs: int = 1,
    show_progress: bool = False,
) -> List[EnumerationResult]:
    """
    Enumerate products from reagent combinations.

    Args:
        pipeline: SynthesisPipeline instance (for reaction config and state)
        combinations: List of (reagent_mols, reagent_keys) tuples
        n_jobs: Number of parallel workers (1 = sequential)
        show_progress: Show progress bar

    Returns:
        List of EnumerationResult objects
    """
    if n_jobs == 1:
        return _enumerate_sequential(pipeline, combinations, show_progress)
    else:
        return _enumerate_parallel(pipeline, combinations, n_jobs, show_progress)


def _enumerate_sequential(
    pipeline: "SynthesisPipeline",
    combinations: List[Tuple[List[Chem.Mol], Optional[List[str]]]],
    show_progress: bool = False,
) -> List[EnumerationResult]:
    """Sequential enumeration of products."""
    results = []
    iterator = combinations

    if show_progress:
        try:
            from tqdm import tqdm

            iterator = tqdm(combinations, desc="Enumerating")
        except ImportError:
            pass

    for mols, keys in iterator:
        result = pipeline.enumerate_single(mols, keys)
        results.append(result)

    return results


def _enumerate_parallel(
    pipeline: "SynthesisPipeline",
    combinations: List[Tuple[List[Chem.Mol], Optional[List[str]]]],
    n_jobs: int,
    show_progress: bool = False,
) -> List[EnumerationResult]:
    """Parallel enumeration using multiprocessing."""
    from functools import partial

    # Prepare worker data (serializable)
    worker_data = pipeline.prepare_worker_data()

    # Convert molecules to SMILES for serialization
    smiles_combinations = []
    for mols, keys in combinations:
        smiles_list = [Chem.MolToSmiles(m) if m else "" for m in mols]
        smiles_combinations.append((smiles_list, keys))

    # Run in parallel
    ctx = get_context("spawn")
    with ctx.Pool(n_jobs) as pool:
        func = partial(
            _worker_enumerate,
            worker_data=worker_data,
        )

        if show_progress:
            try:
                from tqdm import tqdm

                results = list(
                    tqdm(
                        pool.imap(func, smiles_combinations),
                        total=len(smiles_combinations),
                        desc="Enumerating (parallel)",
                    )
                )
            except ImportError:
                results = pool.map(func, smiles_combinations)
        else:
            results = pool.map(func, smiles_combinations)

    return results


def _worker_enumerate(
    combo: Tuple[List[str], Optional[List[str]]],
    worker_data: Dict,
) -> EnumerationResult:
    """Worker function for parallel enumeration."""
    # Import here to avoid circular imports in worker processes
    from .synthesis_pipeline import SynthesisPipeline

    pipeline = SynthesisPipeline.from_worker_data(worker_data)
    smiles_list, keys = combo

    # Convert SMILES back to molecules
    mols = []
    for smiles in smiles_list:
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            mols.append(mol)
        else:
            mols.append(None)

    return pipeline.enumerate_single(mols, keys)


# === Reagent Combination Generation ===


def generate_all_combinations(
    reagent_files: List[str],
    return_separate_names: bool = False,
) -> Tuple[List[Tuple[Chem.Mol, ...]], List[str], Optional[List[List[str]]]]:
    """
    Generate all reagent combinations from files.

    Args:
        reagent_files: List of reagent file paths
        return_separate_names: If True, also return separate reagent names per combination

    Returns:
        Tuple of (mol_combinations, name_combinations, reagent_name_lists)
        - mol_combinations: List of tuples of Mol objects
        - name_combinations: List of product names (reagent names joined by _)
        - reagent_name_lists: List of [reagent_name, ...] lists (only if return_separate_names=True)
    """
    from .enumeration_utils import read_reagent_file

    # Load reagents from each file
    all_reagents = []
    all_names = []

    for file_path in reagent_files:
        reagents = read_reagent_file(file_path)
        smiles_list = [r[0] for r in reagents]
        names_list = [r[1] for r in reagents]
        all_reagents.append(smiles_list)
        all_names.append(names_list)

    # Generate all combinations
    smiles_combinations = list(itertools_product(*all_reagents))
    all_name_combos = list(itertools_product(*all_names))
    name_combinations = ["_".join(names) for names in all_name_combos]

    # Convert SMILES to Mol objects
    mol_combinations = []
    valid_names = []
    valid_name_lists = [] if return_separate_names else None

    for i, smiles_combo in enumerate(smiles_combinations):
        mols = []
        valid = True
        for smiles in smiles_combo:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                valid = False
                break
            mols.append(mol)

        if valid:
            mol_combinations.append(tuple(mols))
            valid_names.append(name_combinations[i])
            if return_separate_names:
                valid_name_lists.append(list(all_name_combos[i]))

    return mol_combinations, valid_names, valid_name_lists
