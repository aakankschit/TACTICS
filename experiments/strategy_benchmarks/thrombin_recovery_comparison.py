"""
Thrombin Dataset: Recovery of Top 100 Actives - Strategy Comparison

This script replicates the legacy notebook analysis to compare how well
different selection strategies recover the top 100 compounds from brute force docking.

Matches the parameters from legacy_thrombin_results.ipynb.
"""

import sys
import json
from pathlib import Path
import pandas as pd
import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

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
# Configuration - Matching Legacy Notebook
# ============================================================================

REAGENT_FILES = [
    "/Users/aakankschitnandkeolyar/Desktop/TACTICS/examples/input_data/acids.smi",
    "/Users/aakankschitnandkeolyar/Desktop/TACTICS/examples/input_data/coupled_aa_sub.smi"
]

REACTION_SMARTS = "[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]>>[#6:1](=[O:2])[#7:3]"

EVALUATOR_CONFIG = {
    "ref_filename": "data/scores/thrombin/product_scores.csv"
}

# Match legacy parameters
NUM_WARMUP_TRIALS = 10      # Legacy uses 10
NUM_TS_ITERATIONS = 5000    # Legacy uses 5000
NUM_CYCLES = 10             # Run 10 replicates like legacy
BATCH_SIZE = 1              # Single compound mode

# Load reference (brute force) scores
REFERENCE_FILE = "data/scores/thrombin/product_scores.csv"

# ============================================================================
# Strategy Configurations
# ============================================================================

STRATEGIES = {
    "TS_Greedy": {
        "strategy": GreedySelection(mode="minimize"),
        "color": "#e41a1c"
    },
    "TS_RouletteWheel": {
        "strategy": RouletteWheelSelection(
            mode="minimize",
            alpha=0.1,
            beta=0.1
        ),
        "color": "#377eb8"
    },
    "TS_UCB": {
        "strategy": UCBSelection(mode="minimize", c=2.0),
        "color": "#4daf4a"
    },
    "TS_EpsilonGreedy": {
        "strategy": EpsilonGreedySelection(
            mode="minimize",
            epsilon=0.1,
            decay=0.995
        ),
        "color": "#984ea3"
    }
}

# ============================================================================
# Helper Functions
# ============================================================================

def load_reference_data():
    """Load brute force reference data and get top 100."""
    print("Loading reference (brute force) data...")

    # Load reference scores
    ref_df = pd.read_csv(REFERENCE_FILE)

    # Standardize column names
    if "Product_Code" in ref_df.columns:
        ref_df = ref_df.rename(columns={"Product_Code": "Name", "Scores": "score"})

    # Sort and get top 100
    ref_df = ref_df.sort_values("score")
    top_100_ref = ref_df.head(100).copy()

    print(f"  Total products in reference: {len(ref_df):,}")
    print(f"  Best score in reference: {ref_df['score'].min():.3f}")
    print(f"  Top 100 score range: {top_100_ref['score'].min():.3f} to {top_100_ref['score'].max():.3f}")

    return ref_df, top_100_ref


def run_ts_replicate(strategy, strategy_name, replicate_num, verbose=False):
    """
    Run a single Thompson Sampling replicate.

    Returns:
        pd.DataFrame: Results with scores, SMILES, and names
    """
    if verbose:
        print(f"\n{'='*80}")
        print(f"Running: {strategy_name} - Replicate {replicate_num + 1}/{NUM_CYCLES}")
        print(f"{'='*80}")

    # Create sampler
    sampler = ThompsonSampler(
        selection_strategy=strategy,
        batch_size=BATCH_SIZE
    )

    # Setup
    evaluator = LookupEvaluator(json.dumps(EVALUATOR_CONFIG))
    sampler.set_evaluator(evaluator)
    sampler.read_reagents(REAGENT_FILES)
    sampler.set_reaction(REACTION_SMARTS)
    sampler.set_hide_progress(not verbose)

    if verbose:
        print(f"Total possible products: {sampler.get_num_prods():.2e}")

    # Warmup
    if verbose:
        print(f"Running warmup ({NUM_WARMUP_TRIALS} trials)...")
    sampler.warm_up(num_warmup_trials=NUM_WARMUP_TRIALS)

    # Search
    if verbose:
        print(f"Running search ({NUM_TS_ITERATIONS} iterations)...")
    results = sampler.search(num_cycles=NUM_TS_ITERATIONS)

    # Convert to DataFrame
    df = pd.DataFrame(results, columns=["score", "SMILES", "Name"])

    if verbose:
        print(f"\nResults:")
        print(f"  Compounds evaluated: {len(df)}")
        print(f"  Unique compounds: {df['SMILES'].nunique()}")
        print(f"  Best score: {df['score'].min():.3f}")
        print(f"  Top 10 mean: {df.nsmallest(10, 'score')['score'].mean():.3f}")

    return df


def calculate_recovery(ts_results, reference_top_100, top_n=100):
    """
    Calculate how many of the top N reference compounds were recovered.

    Returns:
        int: Number of reference compounds recovered
        float: Fraction recovered
    """
    # Get top N from TS results
    ts_top_n = ts_results.nsmallest(top_n, "score")

    # Get reference names
    ref_names = set(reference_top_100["Name"].values)
    ts_names = set(ts_top_n["Name"].values)

    # Count matches
    matches = len(ref_names.intersection(ts_names))
    fraction = matches / len(ref_names)

    return matches, fraction


def run_all_replicates():
    """Run all strategies for all replicates."""

    # Load reference data
    ref_df, top_100_ref = load_reference_data()

    print(f"\n{'='*80}")
    print(f"RUNNING ALL STRATEGIES AND REPLICATES")
    print(f"{'='*80}")
    print(f"Strategies: {', '.join(STRATEGIES.keys())}")
    print(f"Replicates per strategy: {NUM_CYCLES}")
    print(f"Iterations per replicate: {NUM_TS_ITERATIONS}")
    print(f"Warmup trials: {NUM_WARMUP_TRIALS}")

    # Store all results
    all_results = []
    recovery_stats = []

    for strategy_name, config in STRATEGIES.items():
        print(f"\n{'='*80}")
        print(f"Strategy: {strategy_name}")
        print(f"{'='*80}")

        strategy_replicates = []

        for rep in tqdm(range(NUM_CYCLES), desc=f"{strategy_name}"):
            # Run replicate
            df = run_ts_replicate(
                strategy=config["strategy"],
                strategy_name=strategy_name,
                replicate_num=rep,
                verbose=False
            )

            # Add metadata
            df["method"] = strategy_name
            df["cycle"] = str(rep)

            strategy_replicates.append(df)

            # Calculate recovery
            matches, fraction = calculate_recovery(df, top_100_ref, top_n=100)
            recovery_stats.append({
                "method": strategy_name,
                "cycle": str(rep),
                "recovered": matches,
                "fraction": fraction,
                "best_score": df["score"].min(),
                "top_100_mean": df.nsmallest(100, "score")["score"].mean()
            })

        # Combine replicates
        all_results.extend(strategy_replicates)

        # Summary for this strategy
        mean_recovery = np.mean([s["recovered"] for s in recovery_stats if s["method"] == strategy_name])
        std_recovery = np.std([s["recovered"] for s in recovery_stats if s["method"] == strategy_name])
        print(f"\n  Recovery: {mean_recovery:.1f} ± {std_recovery:.1f} out of 100")

    return all_results, recovery_stats, top_100_ref


def create_recovery_plots(all_results, recovery_stats, top_100_ref):
    """Create grouped barplot showing recovery per cycle like legacy notebook."""

    # Get reference product names
    reference_products = set(top_100_ref["Name"].values)

    # Get numbered cycles (0-9)
    numbered_cycles = sorted([str(i) for i in range(NUM_CYCLES)])
    all_cycles = numbered_cycles + ["concat"]

    # Calculate match counts for each cycle and method
    match_counts = []

    for cycle in all_cycles:
        for strategy_name in STRATEGIES.keys():
            if cycle == "concat":
                # For concat, get unique recovered compounds across all cycles
                cycle_dfs = [df for df in all_results if df["method"].iloc[0] == strategy_name]
                combined_df = pd.concat(cycle_dfs, ignore_index=True)
                unique_names = combined_df["Name"].unique()
                match_count = len(set(unique_names).intersection(reference_products))
            else:
                # For individual cycles, get matches for that specific replicate
                cycle_df = [df for df in all_results
                           if df["method"].iloc[0] == strategy_name and df["cycle"].iloc[0] == cycle]
                if cycle_df:
                    names = cycle_df[0]["Name"].values
                    match_count = len(set(names).intersection(reference_products))
                else:
                    match_count = 0

            match_counts.append({
                "cycle": cycle,
                "method": strategy_name,
                "match_count": match_count
            })

    match_counts_df = pd.DataFrame(match_counts)

    # Create grouped bar plot
    palette_colors = [STRATEGIES[method]["color"] for method in STRATEGIES.keys()]
    fig, ax = plt.subplots(figsize=(16, 8))

    sns.barplot(data=match_counts_df, x="cycle", y="match_count", hue="method",
                palette=palette_colors, ax=ax, order=all_cycles)

    # Add value labels on top of bars
    for container in ax.containers:
        ax.bar_label(container, fontsize=12, padding=3)

    ax.set_ylabel("Count of Recovered Actives (out of 100)", fontsize=18)
    ax.set_xlabel("Cycle", fontsize=18)
    ax.set_title("Recovery of Top 100 Actives per Cycle", fontsize=20, fontweight='bold', pad=20)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=14, title="Strategy")
    ax.set_ylim([0, 105])
    ax.axhline(y=100, color='red', linestyle='--', alpha=0.3, linewidth=1)

    plt.tight_layout()
    return fig


def create_summary_table(recovery_stats):
    """Create summary table of results."""

    recovery_df = pd.DataFrame(recovery_stats)

    summary = recovery_df.groupby("method").agg({
        "recovered": ["mean", "std", "min", "max"],
        "fraction": "mean",
        "best_score": ["mean", "std"],
        "top_100_mean": ["mean", "std"]
    }).round(2)

    # Flatten column names
    summary.columns = ['_'.join(col).strip() for col in summary.columns.values]
    summary = summary.reset_index()

    # Format for display
    display_summary = pd.DataFrame({
        "Strategy": summary["method"],
        "Recovered (mean±std)": summary.apply(
            lambda x: f"{x['recovered_mean']:.1f}±{x['recovered_std']:.1f}", axis=1
        ),
        "Recovered (min-max)": summary.apply(
            lambda x: f"{x['recovered_min']:.0f}-{x['recovered_max']:.0f}", axis=1
        ),
        "Recovery %": (summary["fraction_mean"] * 100).round(1),
        "Best Score": summary.apply(
            lambda x: f"{x['best_score_mean']:.2f}±{x['best_score_std']:.2f}", axis=1
        ),
        "Top 100 Mean": summary.apply(
            lambda x: f"{x['top_100_mean_mean']:.2f}±{x['top_100_mean_std']:.2f}", axis=1
        )
    })

    return display_summary


# ============================================================================
# Main Execution
# ============================================================================

def main():
    """Run the full comparison."""

    print("="*80)
    print("THROMBIN DATASET: TOP 100 ACTIVES RECOVERY COMPARISON")
    print("="*80)
    print("\nThis script replicates the legacy notebook analysis.")
    print("It compares how well different strategies recover the")
    print("top 100 compounds from brute force docking.")

    # Run all strategies and replicates
    all_results, recovery_stats, top_100_ref = run_all_replicates()

    # Create summary table
    print("\n" + "="*80)
    print("SUMMARY RESULTS")
    print("="*80)
    summary = create_summary_table(recovery_stats)
    print(summary.to_string(index=False))

    # Save summary
    summary.to_csv("thrombin_recovery_summary.csv", index=False)
    print("\n✅ Summary saved to: thrombin_recovery_summary.csv")

    # Save detailed results
    detailed_df = pd.DataFrame(recovery_stats)
    detailed_df.to_csv("thrombin_recovery_detailed.csv", index=False)
    print("✅ Detailed results saved to: thrombin_recovery_detailed.csv")

    # Create plots
    print("\n📊 Generating comparison plots...")
    fig = create_recovery_plots(all_results, recovery_stats, top_100_ref)
    fig.savefig("thrombin_recovery_comparison.png", dpi=300, bbox_inches='tight')
    print("✅ Plots saved to: thrombin_recovery_comparison.png")

    # Find best strategy
    mean_recovery = detailed_df.groupby("method")["recovered"].mean()
    best_strategy = mean_recovery.idxmax()
    best_recovery = mean_recovery.max()

    print("\n" + "="*80)
    print("✅ COMPARISON COMPLETE!")
    print("="*80)
    print(f"\n🏆 Best Strategy: {best_strategy}")
    print(f"   Mean Recovery: {best_recovery:.1f}/100 compounds")
    print(f"   Recovery Rate: {(best_recovery/100)*100:.1f}%")

    plt.show()


if __name__ == "__main__":
    main()
