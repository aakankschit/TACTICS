"""
Experiment Template for TACTICS

Copy this file to create a new experiment:
    cp experiments/_templates/experiment_template.py experiments/{category}/{name}.py

Then modify the EXPERIMENT CONFIGURATION section below.

Usage:
    python experiments/{category}/{name}.py
    
Outputs will be saved to:
    outputs/runs/{YYYY-MM-DD}_{EXPERIMENT_NAME}/
"""

# =============================================================================
# EXPERIMENT CONFIGURATION - Modify this section
# =============================================================================

EXPERIMENT_NAME = "template_experiment"  # Change this!
DATASET = "thrombin"  # Dataset to use: thrombin, seh, adenine, quinazoline
DESCRIPTION = """
Brief description of what this experiment tests.
"""

# Experiment parameters
PARAMS = {
    "iterations": 500,
    "warmup_samples": 50,
    "strategy": "rws",
    "seeds": [42, 123, 456],  # Multiple seeds for statistical validity
}

# =============================================================================
# IMPORTS
# =============================================================================

import sys
from pathlib import Path

# Add project root to path for imports
project_root = Path(__file__).resolve().parents[2]
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))
    sys.path.insert(0, str(project_root / "src"))

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, Any, List

# Research infrastructure (not part of package)
from experiment_utils import ExperimentRun, load_dataset

# TACTICS package
from TACTICS.thompson_sampling import ThompsonSamplingConfig, run_ts
from TACTICS.thompson_sampling.presets import get_preset


# =============================================================================
# EXPERIMENT SETUP
# =============================================================================

def setup_experiment() -> ExperimentRun:
    """Initialize experiment run with logging and config."""
    run = ExperimentRun(EXPERIMENT_NAME, dataset=DATASET)
    
    # Log experiment description
    run.log.info(f"Description: {DESCRIPTION.strip()}")
    
    # Save configuration
    config = {
        "description": DESCRIPTION.strip(),
        **PARAMS,
    }
    run.save_config(config)
    
    return run


# =============================================================================
# EXPERIMENT LOGIC
# =============================================================================

def run_single_trial(
    run: ExperimentRun,
    seed: int,
    **kwargs
) -> Dict[str, Any]:
    """
    Run a single trial of the experiment.
    
    Args:
        run: ExperimentRun instance for logging
        seed: Random seed for reproducibility
        **kwargs: Additional parameters
    
    Returns:
        Dictionary of results for this trial
    """
    run.log.info(f"Running trial with seed={seed}")
    
    # Get base configuration
    config = get_preset("balanced")
    config.iterations = PARAMS["iterations"]
    config.warmup.num_samples = PARAMS["warmup_samples"]
    config.seed = seed
    
    # Configure evaluator for dataset
    dataset = run.dataset
    if dataset.has_scores:
        config.evaluator.type = "lookup"
        config.evaluator.score_file = str(dataset.score_file)
    
    # Run Thompson Sampling
    results = run_ts(config)
    
    # Extract metrics
    trial_results = {
        "seed": seed,
        "best_score": float(results.best_scores[-1]),
        "iterations_to_best": int(np.argmax(results.best_scores)),
        # Add more metrics as needed
    }
    
    run.log.info(f"Trial complete: best_score={trial_results['best_score']:.4f}")
    
    return trial_results


def run_experiment(run: ExperimentRun) -> List[Dict[str, Any]]:
    """
    Run all trials of the experiment.
    
    Returns:
        List of trial results
    """
    all_results = []
    
    for seed in PARAMS["seeds"]:
        result = run_single_trial(run, seed)
        all_results.append(result)
    
    return all_results


# =============================================================================
# ANALYSIS & VISUALIZATION
# =============================================================================

def analyze_results(run: ExperimentRun, results: List[Dict[str, Any]]):
    """Compute summary statistics from trial results."""
    best_scores = [r["best_score"] for r in results]
    
    summary = {
        "n_trials": len(results),
        "best_score_mean": float(np.mean(best_scores)),
        "best_score_std": float(np.std(best_scores)),
        "best_score_min": float(np.min(best_scores)),
        "best_score_max": float(np.max(best_scores)),
    }
    
    run.log.info(f"Summary: mean={summary['best_score_mean']:.4f} ± {summary['best_score_std']:.4f}")
    
    return summary


def create_figures(run: ExperimentRun, results: List[Dict[str, Any]]):
    """Generate and save visualization figures."""
    
    # Example: Box plot of best scores
    fig, ax = plt.subplots(figsize=(8, 6))
    
    best_scores = [r["best_score"] for r in results]
    ax.boxplot(best_scores)
    ax.set_ylabel("Best Score")
    ax.set_title(f"{EXPERIMENT_NAME}: Score Distribution")
    
    fig.savefig(run.figure_path("score_distribution"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    run.log.info("Saved figure: score_distribution.png")
    
    # Add more figures as needed...


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Main experiment entry point."""
    # Setup
    run = setup_experiment()
    
    try:
        # Run experiment
        results = run_experiment(run)
        
        # Analyze
        summary = analyze_results(run, results)
        
        # Save results
        run.save_results({
            "trials": results,
            "summary": summary,
        })
        
        # Create figures
        create_figures(run, results)
        
        run.log.info(f"Experiment complete. Outputs: {run.run_dir}")
        
    except Exception as e:
        run.log.error(f"Experiment failed: {e}")
        raise


if __name__ == "__main__":
    main()