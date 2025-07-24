#!/usr/bin/env python

import heapq
import math
import json
from itertools import product
from typing import Optional

import numpy as np
import pandas as pd
import polars as pl
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm.auto import tqdm

from .config import RandomBaselineConfig
from .utils import get_logger
from .utils.ts_utils import create_reagents


def keep_largest(items, n):
    """Keeps the n largest items in a list, designed to work with a list of [score,SMILES]
    :param items: the list of items to keep
    :param n: the number of items to keep
    :return: list of the n largest items
    """
    heap = []
    for item in items:
        if len(heap) < n:
            heapq.heappush(heap, item)
        else:
            if item[0] > heap[0][0]:
                heapq.heapreplace(heap, item)
    return heap


def setup_baseline_from_config(config: RandomBaselineConfig):
    """Setup baseline components from Pydantic config.
    
    Args:
        config: RandomBaselineConfig object
        
    Returns:
        tuple: (evaluator, reaction, reagent_lists)
    """
    import importlib
    
    # Create evaluator
    module = importlib.import_module("PRISMS.thompson_sampling.core.evaluators")
    class_ = getattr(module, config.evaluator_class_name)
    
    # Handle evaluator_arg - convert dict to JSON string if needed
    evaluator_arg = config.evaluator_arg
    if isinstance(evaluator_arg, dict):
        evaluator_arg = json.dumps(evaluator_arg)
    
    evaluator = class_(evaluator_arg)
    
    # Create reaction
    rxn = AllChem.ReactionFromSmarts(config.reaction_smarts)
    
    # Create reagent lists
    reagent_lists = []
    for reagent_filename in config.reagent_file_list:
        reagent_list = create_reagents(filename=reagent_filename, num_to_select=None, ts_mode="standard")
        reagent_lists.append(reagent_list)
    
    return evaluator, rxn, reagent_lists


def run_random_baseline(config: RandomBaselineConfig, hide_progress: bool = False) -> pl.DataFrame:
    """Run random baseline sampling using Pydantic configuration.
    
    Args:
        config: RandomBaselineConfig object
        hide_progress: hide the progress bar
        
    Returns:
        pl.DataFrame: Results dataframe with scores, SMILES, and names
    """
    # Setup components
    evaluator, rxn, reagent_lists = setup_baseline_from_config(config)
    
    # Setup logging
    log_filename = getattr(config, "log_filename", None)
    logger = get_logger(__name__, filename=log_filename)
    
    # Calculate total products
    num_reagents = len(reagent_lists)
    len_list = [len(x) for x in reagent_lists]
    total_prods = math.prod(len_list)
    
    if not hide_progress:
        print(f"{total_prods:.2e} products")
    
    # Run random sampling
    score_list = []
    for i in tqdm(range(0, config.num_trials), disable=hide_progress):
        reagent_mol_list = []
        reagent_name_list = []
        
        for j in range(0, num_reagents):
            reagent_idx = np.random.randint(0, len_list[j] - 1)
            reagent_mol_list.append(reagent_lists[j][reagent_idx].mol)
            reagent_name_list.append(reagent_lists[j][reagent_idx].reagent_name)
        
        prod = rxn.RunReactants(reagent_mol_list)
        if len(prod):
            product_mol = prod[0][0]
            Chem.SanitizeMol(product_mol)
            product_smiles = Chem.MolToSmiles(product_mol)
            product_name = "_".join(reagent_name_list)
            score = evaluator.evaluate(product_name)  # Using the name here for the lookup evaluator
            score_list = keep_largest(score_list + [[score, product_smiles, product_name]], config.num_to_save)
    
    # Create results dataframe
    score_df = pl.DataFrame(score_list, schema=["score", "SMILES", "Name"], orient="row")
    
    # Sort results
    if config.ascending_output:
        score_df = score_df.sort("score", descending=False)
    else:
        score_df = score_df.sort("score", descending=True)
    
    # Save results if filename provided
    if config.outfile_name is not None:
        score_df.write_csv(config.outfile_name)
        logger.info(f"Saved results to: {config.outfile_name}")
    
    # Log summary
    total_evaluations = evaluator.counter
    percent_searched = total_evaluations / total_prods * 100
    logger.info(f"{total_evaluations} evaluations | {percent_searched:.3f}% of total")
    
    if not hide_progress:
        print(score_df.head(10))
    
    return score_df


def run_exhaustive_baseline(config: RandomBaselineConfig, num_to_select: Optional[int] = None, hide_progress: bool = False) -> pl.DataFrame:
    """Run exhaustive baseline sampling using Pydantic configuration.
    
    Args:
        config: RandomBaselineConfig object
        num_to_select: Number of reagents to use per file (None for all)
        hide_progress: hide the progress bar
        
    Returns:
        pl.DataFrame: Results dataframe with scores, SMILES, and names
    """
    # Setup components
    evaluator, rxn, reagent_lists = setup_baseline_from_config(config)
    
    # Setup logging
    log_filename = getattr(config, "log_filename", None)
    logger = get_logger(__name__, filename=log_filename)
    
    # Create reagent lists with selection if specified
    if num_to_select is not None:
        reagent_lists = []
        for reagent_filename in config.reagent_file_list:
            reagent_list = create_reagents(filename=reagent_filename, num_to_select=num_to_select, ts_mode="standard")
            reagent_lists.append(reagent_list)
    
    # Calculate total products
    len_list = [len(x) for x in reagent_lists]
    total_prods = math.prod(len_list)
    
    if not hide_progress:
        print(f"{total_prods:.2e} products")
    
    # Run exhaustive sampling
    score_list = []
    for reagents in tqdm(product(*reagent_lists), total=total_prods, disable=hide_progress):
        reagent_mol_list = [x.mol for x in reagents]
        prod = rxn.RunReactants(reagent_mol_list)
        if len(prod):
            product_mol = prod[0][0]
            Chem.SanitizeMol(product_mol)
            product_smiles = Chem.MolToSmiles(product_mol)
            product_name = "_".join([x.reagent_name for x in reagents])
            score = evaluator.evaluate(product_name)
            score_list = keep_largest(score_list + [[score, product_smiles, product_name]], config.num_to_save)
    
    # Create results dataframe
    score_df = pl.DataFrame(score_list, schema=["score", "SMILES", "Name"], orient="row")
    
    # Sort results
    if config.ascending_output:
        score_df = score_df.sort("score", descending=False)
    else:
        score_df = score_df.sort("score", descending=True)
    
    # Save results if filename provided
    if config.outfile_name is not None:
        score_df.write_csv(config.outfile_name)
        logger.info(f"Saved results to: {config.outfile_name}")
    
    # Log summary
    total_evaluations = evaluator.counter
    percent_searched = total_evaluations / total_prods * 100
    logger.info(f"{total_evaluations} evaluations | {percent_searched:.3f}% of total")
    
    if not hide_progress:
        print(score_df.head(10))
    
    return score_df


def main():
    """CLI entry point for Random Baseline."""
    import sys
    
    # For CLI usage: load config from JSON file and parse with Pydantic
    json_filename = sys.argv[1]
    
    with open(json_filename, "r") as f:
        data = json.load(f)
    
    config = RandomBaselineConfig.model_validate(data)
    run_random_baseline(config)


if __name__ == "__main__":
    main() 