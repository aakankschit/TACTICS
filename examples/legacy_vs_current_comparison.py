"""
Legacy vs Current Implementation Comparison
Compares TS Greedy and Roulette Wheel selection between legacy and current implementations.
"""

import sys
import json
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm



# Add TACTICS to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

# Import current implementation
from TACTICS.thompson_sampling import (
    ThompsonSampler,
    GreedySelection,
    RouletteWheelSelection,
    LookupEvaluator  # NEW - warmup without replacement (recommended)
)
from TACTICS.thompson_sampling.warmup.latin_hypercube import LatinHypercubeWarmup
from TACTICS.thompson_sampling.warmup.stratified import StratifiedWarmup

# Import legacy implementation
from TACTICS.thompson_sampling.legacy.standard_thompson_sampling import StandardThompsonSampler
from TACTICS.thompson_sampling.legacy.enhanced_thompson_sampling import EnhancedThompsonSampler
from TACTICS.thompson_sampling.legacy.evaluators import LookupEvaluator as LegacyLookupEvaluator
import polars as pl
import tempfile
import os

# ============================================================================
# Configuration
# ============================================================================

REAGENT_FILES = [
    "/Users/aakankschitnandkeolyar/Desktop/TACTICS/examples/input_data/acids.smi",
    "/Users/aakankschitnandkeolyar/Desktop/TACTICS/examples/input_data/coupled_aa_sub.smi"
]

REACTION_SMARTS = "[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]>>[#6:1](=[O:2])[#7:3]"

EVALUATOR_CONFIG = {
    "ref_filename": "/Users/aakankschitnandkeolyar/Desktop/TS_Chem_Space/Thrombin/Linear_amide/docking_scores/product_scores.csv"
}

# Reference data for recovery calculation
REFERENCE_FILE = "/Users/aakankschitnandkeolyar/Desktop/TS_Chem_Space/Thrombin/Linear_amide/docking_scores/product_scores.csv"

# Experiment parameters - matching legacy
NUM_WARMUP_TRIALS = 10
NUM_TS_ITERATIONS = 5000
NUM_CYCLES = 10  # Number of replicates

# ============================================================================
# Helper Functions
# ============================================================================

def load_reference_data():
    """Load brute force reference data and get top 100."""
    print("Loading reference (brute force) data...")

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


def run_current_implementation(strategy_name, selection_strategy, replicate_num, verbose=False):
    """Run current implementation."""

    if verbose:
        print(f"Running current {strategy_name} - Replicate {replicate_num + 1}")

    # Create sampler with current implementation
    # Using StratifiedWarmup to ensure each trial provides new information (no duplicates)
    sampler = ThompsonSampler(
        selection_strategy=selection_strategy,
        warmup_strategy=LatinHypercubeWarmup(),  # NEW - eliminates within-reagent redundancy
        batch_size=1
    )

    # Setup
    evaluator = LookupEvaluator(json.dumps(EVALUATOR_CONFIG))
    sampler.set_evaluator(evaluator)
    sampler.read_reagents(REAGENT_FILES)
    sampler.set_reaction(REACTION_SMARTS)
    sampler.set_hide_progress(not verbose)

    # Warmup
    sampler.warm_up(num_warmup_trials=NUM_WARMUP_TRIALS)

    # Search
    results = sampler.search(num_cycles=NUM_TS_ITERATIONS)

    # Convert to DataFrame
    df = pd.DataFrame(results, columns=["score", "SMILES", "Name"])

    return df


def run_legacy_implementation(strategy_name, replicate_num, verbose=False):
    """Run legacy implementation."""

    if verbose:
        print(f"Running legacy {strategy_name} - Replicate {replicate_num + 1}")

    if strategy_name == "Greedy":
        # Use StandardThompsonSampler for Greedy
        sampler = StandardThompsonSampler(mode="minimize")

        # Set evaluator
        evaluator = LegacyLookupEvaluator(json.dumps(EVALUATOR_CONFIG))
        sampler.set_evaluator(evaluator)

        # Read reagents
        sampler.read_reagents(REAGENT_FILES)

        # Set reaction
        sampler.set_reaction(REACTION_SMARTS)

        # Set progress visibility
        sampler.set_hide_progress(not verbose)

        # Run warmup
        sampler.warm_up(num_warmup_trials=NUM_WARMUP_TRIALS)

        # Run search
        results = sampler.search(num_cycles=NUM_TS_ITERATIONS)

        # Convert to DataFrame
        df = pd.DataFrame(results, columns=["score", "SMILES", "Name"])

    elif strategy_name == "RouletteWheel":
        # Use EnhancedThompsonSampler for Roulette Wheel (thermal cycling + batch sampling)
        # This is the legacy "TS_enhanced" / roulette wheel method
        #
        # To match StandardTS evaluation budget:
        # - StandardTS: 39,740 warmup + 5,000 search = 44,740 total
        # - EnhancedTS: 38,440 warmup + 6,300 search = 44,740 total
        # - percent_lib = 6300 / 499720 ≈ 0.01261

        temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv')
        temp_filename = temp_file.name
        temp_file.close()

        try:
            sampler = EnhancedThompsonSampler(
                processes=1,
                scaling=-1,  # For minimization
                percent_lib=6300/499720,  # Match StandardTS evaluation budget
                search_stop=1000,
                min_cpds_per_core=1
            )

            # Set evaluator
            evaluator = LegacyLookupEvaluator(json.dumps(EVALUATOR_CONFIG))
            sampler.set_evaluator(evaluator)

            # Read reagents
            sampler.read_reagents(REAGENT_FILES)

            # Set reaction
            sampler.set_reaction(REACTION_SMARTS)

            # Set progress visibility
            sampler.set_hide_progress(not verbose)

            # Run warmup
            sampler.warm_up(num_warmup_trials=NUM_WARMUP_TRIALS, results_filename=temp_filename)

            # Run search
            results = sampler.search(
                percent_of_library=6300/499720,
                stop=1000,
                min_cpds_per_core=1,
                results_filename=temp_filename
            )

            # Read all results from file
            df_pl = pl.read_csv(temp_filename)
            df = df_pl.to_pandas()

            # Sort by score
            df = df.sort_values("score").reset_index(drop=True)

        finally:
            # Clean up temp file
            if os.path.exists(temp_filename):
                os.remove(temp_filename)

    else:
        raise ValueError(f"Unknown strategy: {strategy_name}")

    return df


def calculate_recovery(ts_results, reference_top_100, top_n=100):
    """Calculate how many of the top N reference compounds were recovered."""

    # Get top N from TS results
    ts_top_n = ts_results.nsmallest(top_n, "score")

    # Get reference names
    ref_names = set(reference_top_100["Name"].values)
    ts_names = set(ts_top_n["Name"].values)

    # Count matches
    matches = len(ref_names.intersection(ts_names))
    fraction = matches / len(ref_names)

    return matches, fraction


def run_all_experiments():
    """Run all experiments for both implementations and strategies."""

    # Load reference data
    ref_df, top_100_ref = load_reference_data()

    print(f"\n{'='*80}")
    print(f"RUNNING LEGACY VS CURRENT COMPARISON")
    print(f"{'='*80}")
    print(f"Strategies: TS Greedy, Roulette Wheel")
    print(f"Implementations: Legacy, Current")
    print(f"Replicates per configuration: {NUM_CYCLES}")
    print(f"Iterations per replicate: {NUM_TS_ITERATIONS}")
    print(f"Warmup trials: {NUM_WARMUP_TRIALS}")

    # Results storage - store full dataframes for recovery calculation
    all_dataframes = []

    # Strategy configurations
    strategies = {
        "TS_Greedy": {
            "current_strategy": GreedySelection(mode="minimize"),
            "legacy_name": "Greedy"
        },
        "TS_RouletteWheel": {
            "current_strategy": RouletteWheelSelection(mode="minimize", alpha=0.1, beta=0.1),
            "legacy_name": "RouletteWheel"
        }
    }

    for strategy_name, config in strategies.items():
        print(f"\n{'='*80}")
        print(f"Strategy: {strategy_name}")
        print(f"{'='*80}")

        # Run current implementation
        print(f"\n  Running CURRENT implementation...")
        for rep in tqdm(range(NUM_CYCLES), desc=f"  {strategy_name} Current"):
            df = run_current_implementation(
                strategy_name=strategy_name,
                selection_strategy=config["current_strategy"],
                replicate_num=rep,
                verbose=False
            )

            # Add metadata
            df["strategy"] = strategy_name
            df["implementation"] = "Current"
            df["cycle"] = rep

            all_dataframes.append(df)

        # Run legacy implementation
        print(f"\n  Running LEGACY implementation...")
        for rep in tqdm(range(NUM_CYCLES), desc=f"  {strategy_name} Legacy"):
            df = run_legacy_implementation(
                strategy_name=config["legacy_name"],
                replicate_num=rep,
                verbose=False
            )

            # Add metadata
            df["strategy"] = strategy_name
            df["implementation"] = "Legacy"
            df["cycle"] = rep

            all_dataframes.append(df)

    # Combine all dataframes
    combined_df = pd.concat(all_dataframes, ignore_index=True)

    return combined_df, top_100_ref


def create_grouped_barplot(results_df, top_100_ref):
    """Create two-panel grouped bar plots matching the legacy notebook style."""

    # Get reference products
    reference_products = set(top_100_ref["Name"].values)

    # Get numbered cycles (0-9)
    numbered_cycles = sorted([c for c in results_df["cycle"].unique()])
    all_cycles = [str(c) for c in numbered_cycles] + ["concat"]

    # Calculate match counts for each cycle, implementation, and strategy
    match_counts = []

    for implementation in ["Legacy", "Current"]:
        for strategy in ["TS_Greedy", "TS_RouletteWheel"]:
            for cycle in all_cycles:
                if cycle == "concat":
                    # For concat, get unique recovered compounds across all cycles
                    cycle_data = results_df[
                        (results_df["implementation"] == implementation) &
                        (results_df["strategy"] == strategy)
                    ]
                    # Get top 100 from each cycle, then find unique
                    top_100_per_cycle = cycle_data.groupby("cycle").apply(
                        lambda x: x.nsmallest(100, "score")
                    ).reset_index(drop=True)
                    unique_names = top_100_per_cycle["Name"].unique()
                    match_count = len(set(unique_names).intersection(reference_products))
                else:
                    # For individual cycles, get top 100
                    cycle_data = results_df[
                        (results_df["implementation"] == implementation) &
                        (results_df["strategy"] == strategy) &
                        (results_df["cycle"] == int(cycle))
                    ]
                    if len(cycle_data) > 0:
                        top_100 = cycle_data.nsmallest(100, "score")
                        names = top_100["Name"].values
                        match_count = len(set(names).intersection(reference_products))
                    else:
                        match_count = 0

                match_counts.append({
                    "cycle": cycle,
                    "strategy": strategy,
                    "implementation": implementation,
                    "match_count": match_count
                })

    match_counts_df = pd.DataFrame(match_counts)

    # Create two-panel plot
    fig, axes = plt.subplots(2, 1, figsize=(16, 12))
    fig.suptitle("Legacy vs Current Implementation: Recovery of Top 100 Actives\nTS Greedy and Roulette Wheel Selection",
                 fontsize=18, fontweight='bold', y=0.995)

    # Colors for strategies
    strategy_colors = {
        "TS_Greedy": "#e41a1c",
        "TS_RouletteWheel": "#377eb8"
    }

    # Top panel: Legacy
    ax = axes[0]
    legacy_data = match_counts_df[match_counts_df["implementation"] == "Legacy"]

    sns.barplot(data=legacy_data, x="cycle", y="match_count", hue="strategy",
                palette=strategy_colors, ax=ax, order=all_cycles)

    # Add value labels on top of bars
    for container in ax.containers:
        ax.bar_label(container, fontsize=12, padding=3)

    ax.set_ylabel("Count of Recovered Actives (out of 100)", fontsize=16)
    ax.set_xlabel("Cycle", fontsize=16)
    ax.set_title("LEGACY Implementation", fontsize=16, fontweight='bold', pad=15)
    ax.tick_params(axis='both', which='major', labelsize=14)

    # Get legend handles and labels, then update labels
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, ["TS Greedy", "TS Roulette Wheel"],
              title="Strategy", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=14)

    ax.set_ylim([0, 105])
    ax.axhline(y=100, color='red', linestyle='--', alpha=0.3, linewidth=1)

    # Bottom panel: Current
    ax = axes[1]
    current_data = match_counts_df[match_counts_df["implementation"] == "Current"]

    sns.barplot(data=current_data, x="cycle", y="match_count", hue="strategy",
                palette=strategy_colors, ax=ax, order=all_cycles)

    # Add value labels on top of bars
    for container in ax.containers:
        ax.bar_label(container, fontsize=12, padding=3)

    ax.set_ylabel("Count of Recovered Actives (out of 100)", fontsize=16)
    ax.set_xlabel("Cycle", fontsize=16)
    ax.set_title("CURRENT Implementation", fontsize=16, fontweight='bold', pad=15)
    ax.tick_params(axis='both', which='major', labelsize=14)

    # Get legend handles and labels, then update labels
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, ["TS Greedy", "TS Roulette Wheel"],
              title="Strategy", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=14)

    ax.set_ylim([0, 105])
    ax.axhline(y=100, color='red', linestyle='--', alpha=0.3, linewidth=1)

    plt.tight_layout()
    return fig


def create_summary_table(results_df, top_100_ref):
    """Create summary statistics table."""

    reference_products = set(top_100_ref["Name"].values)

    summary_data = []

    for implementation in ["Legacy", "Current"]:
        for strategy in ["TS_Greedy", "TS_RouletteWheel"]:
            # Get all cycles for this combination
            strategy_data = results_df[
                (results_df["implementation"] == implementation) &
                (results_df["strategy"] == strategy)
            ]

            # Calculate recovery per cycle
            recoveries = []
            best_scores = []
            top_100_means = []

            for cycle in strategy_data["cycle"].unique():
                cycle_data = strategy_data[strategy_data["cycle"] == cycle]
                top_100 = cycle_data.nsmallest(100, "score")

                # Recovery
                matches = len(set(top_100["Name"].values).intersection(reference_products))
                recoveries.append(matches)

                # Best score
                best_scores.append(cycle_data["score"].min())

                # Top 100 mean
                top_100_means.append(top_100["score"].mean())

            summary_data.append({
                "Implementation": implementation,
                "Strategy": strategy.replace("TS_", "TS "),
                "Recovered (mean±std)": f"{np.mean(recoveries):.1f}±{np.std(recoveries):.1f}",
                "Recovered (min-max)": f"{int(np.min(recoveries))}-{int(np.max(recoveries))}",
                "Recovery %": f"{(np.mean(recoveries)/100)*100:.1f}%",
                "Best Score": f"{np.mean(best_scores):.2f}±{np.std(best_scores):.2f}",
                "Top 100 Mean": f"{np.mean(top_100_means):.2f}±{np.std(top_100_means):.2f}"
            })

    return pd.DataFrame(summary_data)


# ============================================================================
# Main Execution
# ============================================================================

def main():
    """Run the full comparison."""

    print("="*80)
    print("LEGACY VS CURRENT IMPLEMENTATION COMPARISON")
    print("="*80)
    print("\nThis script compares TS Greedy and Roulette Wheel selection")
    print("between the legacy code and the current refactored implementation.")
    print("\nFocus: Recovery of top 100 actives from brute force docking.")

    # Run all experiments
    results_df, top_100_ref = run_all_experiments()

    # Create summary table
    print("\n" + "="*80)
    print("SUMMARY TABLE")
    print("="*80)
    summary = create_summary_table(results_df, top_100_ref)
    print(summary.to_string(index=False))

    # Save summary
    summary.to_csv("legacy_vs_current_summary.csv", index=False)
    print("\n✅ Summary saved to: legacy_vs_current_summary.csv")

    # Save detailed results
    results_df.to_csv("legacy_vs_current_detailed.csv", index=False)
    print("✅ Detailed results saved to: legacy_vs_current_detailed.csv")

    # Create plots
    print("\n📊 Generating comparison plots...")
    fig = create_grouped_barplot(results_df, top_100_ref)
    fig.savefig("legacy_vs_current_comparison.png", dpi=300, bbox_inches='tight')
    print("✅ Plots saved to: legacy_vs_current_comparison.png")

    print("\n" + "="*80)
    print("✅ COMPARISON COMPLETE!")
    print("="*80)

    plt.show()


if __name__ == "__main__":
    main()
