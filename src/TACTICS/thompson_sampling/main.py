#!/usr/bin/env python

import sys
import json
from datetime import timedelta
from timeit import default_timer as timer
from typing import Tuple, Optional

import polars as pl

from .config import ThompsonSamplingConfig
from .core import ThompsonSampler
from .utils import get_logger


def run_ts(
    config: ThompsonSamplingConfig,
    hide_progress: bool = False,
    return_warmup: bool = False
) -> pl.DataFrame | Tuple[pl.DataFrame, pl.DataFrame]:
    """
    Convenience wrapper for running Thompson Sampling with automatic logging and file saving.

    This function provides an "out-of-the-box" interface that handles all setup automatically.
    It works seamlessly with preset configurations and supports all selection strategies,
    warmup strategies, and evaluators via the modern Pydantic config system.

    For maximum flexibility, use ThompsonSampler.from_config() directly instead.

    Args:
        config: ThompsonSamplingConfig (works with both modern and legacy format)
            - Modern (recommended): Use nested configs (strategy_config, warmup_config, evaluator_config)
            - Legacy: Use string-based configs (selection_strategy, evaluator_class_name, etc.)
        hide_progress: Hide progress bars during execution
        return_warmup: If True, return tuple of (search_results, warmup_results)

    Returns:
        pl.DataFrame: Search results with columns ["score", "SMILES", "Name"]
        OR Tuple[pl.DataFrame, pl.DataFrame]: (search_results, warmup_results) if return_warmup=True

    Example with Presets (Recommended):
        >>> from TACTICS.thompson_sampling import run_ts, get_preset
        >>> from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig
        >>>
        >>> # 1. Create evaluator config
        >>> evaluator = LookupEvaluatorConfig(ref_filename="scores.csv")
        >>>
        >>> # 2. Get preset configuration
        >>> config = get_preset(
        ...     "fast_exploration",
        ...     reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        ...     reagent_file_list=["acids.smi", "amines.smi"],
        ...     evaluator_config=evaluator,
        ...     mode="minimize"  # For docking scores
        ... )
        >>>
        >>> # 3. Run and get results
        >>> results = run_ts(config)

    Example with Custom Config:
        >>> from TACTICS.thompson_sampling import ThompsonSamplingConfig, run_ts
        >>> from TACTICS.thompson_sampling.strategies.config import EpsilonGreedyConfig
        >>> from TACTICS.thompson_sampling.warmup.config import StratifiedWarmupConfig
        >>> from TACTICS.thompson_sampling.core.evaluator_config import DBEvaluatorConfig
        >>>
        >>> config = ThompsonSamplingConfig(
        ...     reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        ...     reagent_file_list=["acids.smi", "amines.smi"],
        ...     num_ts_iterations=1000,
        ...     strategy_config=EpsilonGreedyConfig(mode="maximize", epsilon=0.2),
        ...     warmup_config=StratifiedWarmupConfig(),
        ...     evaluator_config=DBEvaluatorConfig(db_filename="scores.db"),
        ...     results_filename="results.csv"
        ... )
        >>> results = run_ts(config)

    Available Presets:
        - "fast_exploration": Epsilon-greedy strategy, quick screening
        - "parallel_batch": Roulette wheel with batch processing and multiprocessing
        - "conservative_exploit": Greedy strategy, focus on best reagents
        - "balanced_sampling": UCB strategy, theoretical guarantees
        - "diverse_coverage": Roulette wheel with maximum diversity
    """
    # Override hide_progress if specified in config
    if hasattr(config, 'hide_progress') and config.hide_progress:
        hide_progress = True

    # Set up logging
    log_filename = getattr(config, "log_filename", None)
    logger = get_logger(__name__, filename=log_filename)

    # Create sampler from config (handles both modern and legacy formats)
    sampler = ThompsonSampler.from_config(config)
    sampler.set_hide_progress(hide_progress)

    # Set up reaction
    sampler.set_reaction(config.reaction_smarts)

    # Run warmup
    logger.info("Starting warmup phase...")
    warmup_df = sampler.warm_up(num_warmup_trials=config.num_warmup_trials)

    # Run search
    logger.info("Starting search phase...")
    search_df = sampler.search(num_cycles=config.num_ts_iterations)

    # Log statistics
    evaluator = sampler.evaluator
    total_evaluations = evaluator.counter
    percent_searched = total_evaluations / sampler.get_num_prods() * 100
    logger.info(f"{total_evaluations} evaluations | {percent_searched:.3f}% of total")

    # Save results to disk if filename provided
    result_filename = getattr(config, "results_filename", None)
    if result_filename is not None:
        search_df.write_csv(result_filename)
        logger.info(f"Saved results to: {result_filename}")

    # Display top results (unless hidden)
    if not hide_progress:
        # Determine sort direction from strategy config or legacy mode
        if hasattr(config, 'strategy_config') and config.strategy_config is not None:
            mode = config.strategy_config.mode
        else:
            mode = getattr(config, 'mode', 'maximize')

        # Sort and display
        if mode in ["maximize", "maximize_boltzmann"]:
            print("\nTop 10 results (highest scores):")
            print(search_df.sort("score", descending=True).unique(subset="SMILES").head(10))
        else:
            print("\nTop 10 results (lowest scores):")
            print(search_df.sort("score", descending=False).unique(subset="SMILES").head(10))

    # Clean up multiprocessing resources
    sampler.close()

    # Return results
    if return_warmup:
        return search_df, warmup_df
    else:
        return search_df


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