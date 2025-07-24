#!/usr/bin/env python

import importlib
import sys
import json
from datetime import timedelta
from timeit import default_timer as timer

import polars as pl

from .config import StandardSamplerConfig, EnhancedSamplerConfig, TSConfig
from .core import StandardThompsonSampler, EnhancedThompsonSampler
from .utils import get_logger


def run_ts(config: TSConfig, hide_progress: bool = False) -> pl.DataFrame:
    """
    Perform a Thompson sampling run using a Pydantic config object.
    
    Args:
        config: Pydantic config object (StandardSamplerConfig or EnhancedSamplerConfig)
        hide_progress: hide the progress bar
        
    Returns:
        pl.DataFrame: Results dataframe with scores, SMILES, and names
    """
    # Dynamically import evaluator class
    module = importlib.import_module("PRISMS.thompson_sampling.core.evaluators")
    class_ = getattr(module, config.evaluator_class_name)
    
    # Handle evaluator_arg - convert dict to JSON string if needed
    evaluator_arg = config.evaluator_arg
    if isinstance(evaluator_arg, dict):
        evaluator_arg = json.dumps(evaluator_arg)
    
    evaluator = class_(evaluator_arg)
    result_filename = getattr(config, "results_filename", None)
    log_filename = getattr(config, "log_filename", None)
    logger = get_logger(__name__, filename=log_filename)

    # Switch between modes
    if config.sampler_type == "standard":
        ts = StandardThompsonSampler(mode=config.ts_mode)
    elif config.sampler_type == "enhanced":
        ts = EnhancedThompsonSampler(
            processes=config.processes,
            scaling=config.scaling,
            percent_lib=config.percent_of_library,
            search_stop=config.stopping_criteria,
            min_cpds_per_core=config.minimum_no_of_compounds_per_core
        )
    else:
        raise ValueError(f"Unknown sampler_type: {config.sampler_type}")

    ts.set_hide_progress(hide_progress)
    ts.set_evaluator(evaluator)
    ts.read_reagents(reagent_file_list=config.reagent_file_list, num_to_select=None)
    ts.set_reaction(config.reaction_smarts)
    ts.warm_up(num_warmup_trials=config.num_warmup_trials)
    
    if config.sampler_type == "standard":
        out_list = ts.search(num_cycles=config.num_ts_iterations)
    elif config.sampler_type == "enhanced":
        out_list = ts.search(
            percent_of_library=config.percent_of_library,
            stop=config.stopping_criteria,
            min_cpds_per_core=config.minimum_no_of_compounds_per_core,
            results_filename=result_filename
        )
    
    total_evaluations = evaluator.counter
    percent_searched = total_evaluations / ts.get_num_prods() * 100
    logger.info(f"{total_evaluations} evaluations | {percent_searched:.3f}% of total")
    
    # write the results to disk
    out_df = pl.DataFrame(out_list, schema=["score", "SMILES", "Name"], orient="row")
    if result_filename is not None:
        out_df.write_csv(result_filename)
        logger.info(f"Saved results to: {result_filename}")
    
    if not hide_progress:
        # Use sampler_type to determine sorting since ts_mode only exists for standard sampler
        if config.sampler_type == "standard" and getattr(config, "ts_mode", None) in ["maximize", "boltzmann_maximize"]:
            print(out_df.sort("score", descending=False).unique(subset="SMILES").head(10))
        else:
            print(out_df.sort("score", descending=True).unique(subset="SMILES").head(10))
    
    return out_df


def main():
    """CLI entry point for Thompson Sampling."""
    start = timer()
    
    # For CLI usage: load config from JSON file and parse with Pydantic
    json_filename = sys.argv[1]
    
    with open(json_filename, "r") as f:
        data = json.load(f)
    
    # Determine config type
    if data.get("sampler_type") == "standard":
        config = StandardSamplerConfig.model_validate(data)
    elif data.get("sampler_type") == "enhanced":
        config = EnhancedSamplerConfig.model_validate(data)
    else:
        raise ValueError(f"Unknown sampler_type: {data.get('sampler_type')}")
    
    run_ts(config)
    end = timer()
    print("Elapsed time", timedelta(seconds=end - start))


if __name__ == "__main__":
    main() 