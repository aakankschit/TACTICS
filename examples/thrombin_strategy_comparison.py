"""
Thrombin Dataset: Selection Strategy Comparison

This script compares different Thompson Sampling selection strategies
on the Thrombin inhibitor screening dataset.
"""

import sys
import json
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Add TACTICS to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from TACTICS.thompson_sampling import (
    ThompsonSampler,
    GreedySelection,
    RouletteWheelSelection,
    UCBSelection,
    EpsilonGreedySelection,
    LookupEvaluator
)

# ============================================================================
# Configuration
# ============================================================================

# Paths
REAGENT_FILES = [
    "../examples/input_data/acids.smi",
    "../examples/input_data/coupled_aa_sub.smi"
]

REACTION_SMARTS = "[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]>>[#6:1](=[O:2])[#7:3]"

# Reference scores (adjust path to your data)
EVALUATOR_CONFIG = {
    "ref_filename": "/Users/aakankschitnandkeolyar/Desktop/TS_Chem_Space/Thrombin/Linear_amide/docking_scores/product_scores.csv"
}

# Experiment parameters
NUM_WARMUP_TRIALS = 3
NUM_TS_ITERATIONS = 1000  # Number of compounds to screen
NUM_TO_SELECT = None  # Use all reagents (or set a number for testing)

# ============================================================================
# Strategy Configurations
# ============================================================================

STRATEGIES = {
    "Greedy": {
        "strategy": GreedySelection(mode="minimize"),
        "batch_size": 1,
        "color": "#2ecc71"
    },
    "Roulette_Wheel": {
        "strategy": RouletteWheelSelection(
            mode="minimize",
            alpha=0.1,
            beta=0.1,
            alpha_increment=0.01,
            beta_increment=0.001,
            efficiency_threshold=0.1
        ),
        "batch_size": 1,
        "color": "#e74c3c"
    },
    "Roulette_Wheel_Batch": {
        "strategy": RouletteWheelSelection(
            mode="minimize",
            alpha=0.1,
            beta=0.1,
            alpha_increment=0.01,
            beta_increment=0.001,
            efficiency_threshold=0.1
        ),
        "batch_size": 100,
        "max_resamples": 1000,
        "color": "#e67e22"
    },
    "UCB": {
        "strategy": UCBSelection(mode="minimize", c=2.0),
        "batch_size": 1,
        "color": "#3498db"
    },
    "Epsilon_Greedy": {
        "strategy": EpsilonGreedySelection(
            mode="minimize",
            epsilon=0.1,
            decay=0.995
        ),
        "batch_size": 1,
        "color": "#9b59b6"
    }
}

# ============================================================================
# Helper Functions
# ============================================================================

def run_strategy(strategy_name, strategy_config, verbose=True):
    """
    Run a single Thompson Sampling strategy.

    Returns:
        dict: Results including scores, timing, and statistics
    """
    if verbose:
        print(f"\n{'='*80}")
        print(f"Running: {strategy_name}")
        print(f"{'='*80}")

    # Create sampler
    sampler = ThompsonSampler(
        selection_strategy=strategy_config["strategy"],
        batch_size=strategy_config.get("batch_size", 1),
        max_resamples=strategy_config.get("max_resamples", None)
    )

    # Setup
    evaluator = LookupEvaluator(json.dumps(EVALUATOR_CONFIG))
    sampler.set_evaluator(evaluator)
    sampler.read_reagents(REAGENT_FILES, num_to_select=NUM_TO_SELECT)
    sampler.set_reaction(REACTION_SMARTS)

    if verbose:
        print(f"Total possible products: {sampler.get_num_prods():.2e}")

    # Warmup
    if verbose:
        print(f"\nRunning warmup ({NUM_WARMUP_TRIALS} trials per reagent)...")
    sampler.warm_up(num_warmup_trials=NUM_WARMUP_TRIALS)

    # Search
    if verbose:
        print(f"\nRunning search ({NUM_TS_ITERATIONS} iterations)...")

    import time
    start_time = time.time()
    results = sampler.search(num_cycles=NUM_TS_ITERATIONS)
    elapsed_time = time.time() - start_time

    # Convert to DataFrame
    df = pd.DataFrame(results, columns=["score", "SMILES", "Name"])

    # Calculate statistics
    stats = {
        "strategy": strategy_name,
        "n_compounds": len(df),
        "n_unique": df["SMILES"].nunique(),
        "best_score": df["score"].min(),
        "mean_score": df["score"].mean(),
        "median_score": df["score"].median(),
        "top10_mean": df.nsmallest(10, "score")["score"].mean(),
        "top50_mean": df.nsmallest(50, "score")["score"].mean(),
        "elapsed_time": elapsed_time,
        "batch_size": strategy_config.get("batch_size", 1)
    }

    if verbose:
        print(f"\n{'='*80}")
        print(f"Results for {strategy_name}:")
        print(f"{'='*80}")
        print(f"  Compounds evaluated: {stats['n_compounds']}")
        print(f"  Unique compounds: {stats['n_unique']}")
        print(f"  Best score: {stats['best_score']:.3f}")
        print(f"  Mean score: {stats['mean_score']:.3f}")
        print(f"  Median score: {stats['median_score']:.3f}")
        print(f"  Top 10 mean: {stats['top10_mean']:.3f}")
        print(f"  Top 50 mean: {stats['top50_mean']:.3f}")
        print(f"  Time: {stats['elapsed_time']:.2f}s")
        print(f"\nTop 10 compounds:")
        print(df.nsmallest(10, "score")[["score", "Name"]])

    return {
        "results": df,
        "stats": stats,
        "config": strategy_config
    }


def plot_comparison(all_results):
    """
    Create comparison plots for all strategies.
    """
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle("Thompson Sampling Strategy Comparison - Thrombin Dataset", fontsize=16, fontweight='bold')

    # Extract data
    strategies = list(all_results.keys())
    stats = [all_results[s]["stats"] for s in strategies]
    colors = [all_results[s]["config"]["color"] for s in strategies]

    # Plot 1: Best scores found
    ax = axes[0, 0]
    best_scores = [s["best_score"] for s in stats]
    bars = ax.bar(strategies, best_scores, color=colors, alpha=0.7, edgecolor='black')
    ax.set_ylabel("Best Score (lower is better)", fontsize=12)
    ax.set_title("Best Score Found", fontsize=14, fontweight='bold')
    ax.tick_params(axis='x', rotation=45)
    for i, (bar, score) in enumerate(zip(bars, best_scores)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                f'{score:.3f}', ha='center', va='bottom', fontsize=10)

    # Plot 2: Top 10 mean scores
    ax = axes[0, 1]
    top10_means = [s["top10_mean"] for s in stats]
    bars = ax.bar(strategies, top10_means, color=colors, alpha=0.7, edgecolor='black')
    ax.set_ylabel("Mean Score (lower is better)", fontsize=12)
    ax.set_title("Top 10 Mean Score", fontsize=14, fontweight='bold')
    ax.tick_params(axis='x', rotation=45)
    for i, (bar, score) in enumerate(zip(bars, top10_means)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                f'{score:.3f}', ha='center', va='bottom', fontsize=10)

    # Plot 3: Cumulative best score over iterations
    ax = axes[1, 0]
    for strategy_name, result in all_results.items():
        df = result["results"]
        cumulative_best = df["score"].cummin()
        ax.plot(range(len(cumulative_best)), cumulative_best,
                label=strategy_name, color=result["config"]["color"], linewidth=2)
    ax.set_xlabel("Iteration", fontsize=12)
    ax.set_ylabel("Best Score Found (lower is better)", fontsize=12)
    ax.set_title("Cumulative Best Score", fontsize=14, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    # Plot 4: Score distributions (violin plot)
    ax = axes[1, 1]
    data_for_violin = []
    labels_for_violin = []
    colors_for_violin = []
    for strategy_name, result in all_results.items():
        df = result["results"]
        data_for_violin.append(df["score"].values)
        labels_for_violin.append(strategy_name)
        colors_for_violin.append(result["config"]["color"])

    parts = ax.violinplot(data_for_violin, positions=range(len(data_for_violin)),
                          showmeans=True, showmedians=True)
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors_for_violin[i])
        pc.set_alpha(0.7)
    ax.set_xticks(range(len(labels_for_violin)))
    ax.set_xticklabels(labels_for_violin, rotation=45)
    ax.set_ylabel("Score Distribution", fontsize=12)
    ax.set_title("Score Distributions", fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    return fig


def create_summary_table(all_results):
    """
    Create a summary table comparing all strategies.
    """
    summary_data = []
    for strategy_name, result in all_results.items():
        stats = result["stats"]
        summary_data.append({
            "Strategy": strategy_name,
            "Batch Size": stats["batch_size"],
            "N Compounds": stats["n_compounds"],
            "N Unique": stats["n_unique"],
            "Best Score": f"{stats['best_score']:.3f}",
            "Top 10 Mean": f"{stats['top10_mean']:.3f}",
            "Top 50 Mean": f"{stats['top50_mean']:.3f}",
            "Median": f"{stats['median_score']:.3f}",
            "Time (s)": f"{stats['elapsed_time']:.2f}"
        })

    df = pd.DataFrame(summary_data)
    return df


# ============================================================================
# Main Execution
# ============================================================================

def main():
    """
    Run all strategies and compare results.
    """
    print("="*80)
    print("THROMBIN DATASET: SELECTION STRATEGY COMPARISON")
    print("="*80)
    print(f"\nConfiguration:")
    print(f"  Reagent files: {REAGENT_FILES}")
    print(f"  Warmup trials: {NUM_WARMUP_TRIALS}")
    print(f"  Search iterations: {NUM_TS_ITERATIONS}")
    print(f"  Strategies to test: {', '.join(STRATEGIES.keys())}")

    # Run all strategies
    all_results = {}
    for strategy_name, strategy_config in STRATEGIES.items():
        try:
            result = run_strategy(strategy_name, strategy_config, verbose=True)
            all_results[strategy_name] = result
        except Exception as e:
            print(f"\n❌ Error running {strategy_name}: {e}")
            import traceback
            traceback.print_exc()
            continue

    if not all_results:
        print("\n❌ No strategies completed successfully!")
        return

    # Create summary
    print("\n" + "="*80)
    print("SUMMARY COMPARISON")
    print("="*80)
    summary_df = create_summary_table(all_results)
    print(summary_df.to_string(index=False))

    # Save summary
    summary_df.to_csv("thrombin_strategy_comparison.csv", index=False)
    print("\n✅ Summary saved to: thrombin_strategy_comparison.csv")

    # Create plots
    print("\n📊 Generating comparison plots...")
    fig = plot_comparison(all_results)
    fig.savefig("thrombin_strategy_comparison.png", dpi=300, bbox_inches='tight')
    print("✅ Plots saved to: thrombin_strategy_comparison.png")

    # Save detailed results
    for strategy_name, result in all_results.items():
        filename = f"thrombin_results_{strategy_name}.csv"
        result["results"].to_csv(filename, index=False)
        print(f"✅ Detailed results saved to: {filename}")

    print("\n" + "="*80)
    print("✅ ALL COMPARISONS COMPLETE!")
    print("="*80)

    # Show best strategy
    best_strategy = min(all_results.items(),
                       key=lambda x: x[1]["stats"]["best_score"])
    print(f"\n🏆 Best Strategy: {best_strategy[0]}")
    print(f"   Best Score: {best_strategy[1]['stats']['best_score']:.3f}")
    print(f"   Top 10 Mean: {best_strategy[1]['stats']['top10_mean']:.3f}")

    plt.show()


if __name__ == "__main__":
    main()
