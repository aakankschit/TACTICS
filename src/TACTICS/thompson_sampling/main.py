#!/usr/bin/env python

import importlib
import sys
import json
from datetime import timedelta
from timeit import default_timer as timer

import polars as pl

from .config import ThompsonSamplingConfig
from .core import ThompsonSampler
from .strategies import GreedySelection, RouletteWheelSelection
from .utils import get_logger


def run_ts(config: ThompsonSamplingConfig, hide_progress: bool = False) -> pl.DataFrame:
    """
    Perform a Thompson sampling run using a Pydantic config object.

    Args:
        config: Pydantic config object (ThompsonSamplingConfig)
        hide_progress: hide the progress bar

    Returns:
        pl.DataFrame: Results dataframe with scores, SMILES, and names
    """
    # Dynamically import evaluator class
    module = importlib.import_module("TACTICS.thompson_sampling.core.evaluators")
    class_ = getattr(module, config.evaluator_class_name)

    # Handle evaluator_arg - convert dict to JSON string if needed
    evaluator_arg = config.evaluator_arg
    if isinstance(evaluator_arg, dict):
        evaluator_arg = json.dumps(evaluator_arg)

    evaluator = class_(evaluator_arg)
    result_filename = getattr(config, "results_filename", None)
    log_filename = getattr(config, "log_filename", None)
    logger = get_logger(__name__, filename=log_filename)

    # Create selection strategy based on config
    if config.selection_strategy == "greedy":
        strategy = GreedySelection(mode=config.mode)
    elif config.selection_strategy == "roulette_wheel":
        strategy = RouletteWheelSelection(mode=config.mode)
    else:
        raise ValueError(f"Unknown selection_strategy: {config.selection_strategy}")

    # Create unified Thompson sampler with strategy injection
    ts = ThompsonSampler(
        selection_strategy=strategy,
        warmup_strategy=GreedySelection(mode=config.mode),
        log_filename=log_filename,
        batch_size=getattr(config, "batch_size", 1),
        max_resamples=getattr(config, "max_resamples", None),
        processes=getattr(config, "processes", 1),
        min_cpds_per_core=getattr(config, "min_cpds_per_core", 10)
    )

    ts.set_hide_progress(hide_progress)
    ts.set_evaluator(evaluator)
    ts.read_reagents(reagent_file_list=config.reagent_file_list, num_to_select=None)
    ts.set_reaction(config.reaction_smarts)
    ts.warm_up(num_warmup_trials=config.num_warmup_trials)

    # Run search
    out_list = ts.search(num_cycles=config.num_ts_iterations)

    total_evaluations = evaluator.counter
    percent_searched = total_evaluations / ts.get_num_prods() * 100
    logger.info(f"{total_evaluations} evaluations | {percent_searched:.3f}% of total")

    # write the results to disk
    out_df = pl.DataFrame(out_list, schema=["score", "SMILES", "Name"], orient="row")
    if result_filename is not None:
        out_df.write_csv(result_filename)
        logger.info(f"Saved results to: {result_filename}")

    if not hide_progress:
        # Sort based on mode
        if config.mode in ["maximize", "maximize_boltzmann"]:
            print(out_df.sort("score", descending=True).unique(subset="SMILES").head(10))
        else:
            print(out_df.sort("score", descending=False).unique(subset="SMILES").head(10))

    return out_df


def main():
    """CLI entry point for Thompson Sampling."""
    start = timer()

    # For CLI usage: load config from JSON file and parse with Pydantic
    json_filename = sys.argv[1]

    with open(json_filename, "r") as f:
        data = json.load(f)

    config = ThompsonSamplingConfig.model_validate(data)

    run_ts(config)
    end = timer()
    print("Elapsed time", timedelta(seconds=end - start))


if __name__ == "__main__":
    main() 