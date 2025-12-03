"""
Thompson Sampling Strategy Comparison Benchmark

Compares multiple selection strategies for Thompson Sampling:
- Greedy (Legacy & Current)
- Roulette Wheel with Thermal Cycling (Legacy & Current)
- UCB (Current only)
- ε-Greedy with decay (Current only)
- Bayes-UCB with Thermal Cycling (Current only - NEW!)

This benchmark evaluates recovery of top 100 actives from brute force docking,
allowing direct comparison of different exploration/exploitation strategies.
"""

import sys
import json
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

# Import current implementation
from TACTICS.thompson_sampling import (

# Add project root to path to import experiment_utils
import sys
from pathlib import Path
project_root = Path(__file__).resolve().parents[2]
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

    ThompsonSampler,
    GreedySelection,
    RouletteWheelSelection,
    LookupEvaluator
)
from TACTICS.thompson_sampling.warmup import StandardWarmup, EnhancedWarmup

# Import legacy implementation
from TACTICS.thompson_sampling.legacy.thompson_sampling import ThompsonSampler as LegacyThompsonSampler
from TACTICS.thompson_sampling.legacy.rws_sampling import RWSSampler as LegacyRWSSampler
from TACTICS.thompson_sampling.legacy.evaluators import LookupEvaluator as LegacyLookupEvaluator
import polars as pl
import tempfile
import os

# ============================================================================
# Configuration
# ============================================================================

# 2-component Thrombin library (acids + coupled amino acids)
REAGENT_FILES = [
    "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/reagents/thrombin/acids.smi",
    "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/reagents/thrombin/coupled_aa_sub.smi"
]

# NOTE: No reaction SMARTS needed for LookupEvaluator (uses product codes)

# Current TACTICS evaluator config (uses compound_col and score_col)
EVALUATOR_CONFIG = {
    "ref_filename": "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/scores/thrombin/product_scores.csv",
    "compound_col": "Product_Code",  # Column name for compound identifiers
    "score_col": "Scores"            # Column name for scores
}

# Legacy evaluator config (uses ref_colname)
LEGACY_EVALUATOR_CONFIG = {
    "ref_filename": "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/scores/thrombin/product_scores.csv",
    "ref_colname": "Scores"  # Legacy expects ref_colname for score column
}

# Reference data for recovery calculation
REFERENCE_FILE = "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/scores/thrombin/product_scores.csv"

# Experiment parameters
NUM_WARMUP_TRIALS = 10
NUM_TS_ITERATIONS = 5000
NUM_CYCLES = 5  # Number of replicates (reduced for testing)

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
    # Use EnhancedWarmup and Boltzmann weighting for RouletteWheel to match legacy RWS
    if strategy_name == "TS_RouletteWheel":
        warmup_strategy = EnhancedWarmup()
        use_boltzmann = True
        batch_size = 100  # Match legacy RWS batch size for proper adaptive temperature control
    else:
        warmup_strategy = StandardWarmup()
        use_boltzmann = False
        batch_size = 1  # Greedy uses single-compound selection

    sampler = ThompsonSampler(
        selection_strategy=selection_strategy,
        warmup_strategy=warmup_strategy,
        batch_size=batch_size,
        use_boltzmann_weighting=use_boltzmann
    )

    # Setup - LookupEvaluator now expects dict, not JSON string
    evaluator = LookupEvaluator(EVALUATOR_CONFIG)
    sampler.set_evaluator(evaluator)
    sampler.read_reagents(REAGENT_FILES)
    # NOTE: No reaction SMARTS needed for LookupEvaluator
    sampler.set_hide_progress(not verbose)

    # Warmup
    sampler.warm_up(num_warmup_trials=NUM_WARMUP_TRIALS)

    # Search
    # Use max_evaluations parameter to ensure fair comparison
    # Both Greedy and RWS will evaluate exactly NUM_TS_ITERATIONS unique compounds
    # DisallowTracker ensures no duplicates, so max_evaluations guarantees exact count
    if batch_size > 1:
        # RWS (batch_size=100): Run with large num_cycles, stop after NUM_TS_ITERATIONS unique evals
        # Each cycle generates batch_size combinations, so need fewer cycles
        results = sampler.search(num_cycles=100000, max_evaluations=NUM_TS_ITERATIONS)
        if verbose:
            actual_evals = len(results) if not isinstance(results, pl.DataFrame) else len(results.to_pandas())
            print(f"  Evaluated {actual_evals} unique compounds (target: {NUM_TS_ITERATIONS})")
    else:
        # Greedy (batch_size=1): 1 combination per cycle = 1 eval per cycle
        # With DisallowTracker preventing duplicates, num_cycles = num_evaluations
        # Using max_evaluations for consistency and safety
        results = sampler.search(num_cycles=NUM_TS_ITERATIONS, max_evaluations=NUM_TS_ITERATIONS)
        if verbose:
            actual_evals = len(results) if not isinstance(results, pl.DataFrame) else len(results.to_pandas())
            print(f"  Evaluated {actual_evals} unique compounds (target: {NUM_TS_ITERATIONS})")

    # Convert polars DataFrame to pandas DataFrame
    if isinstance(results, pl.DataFrame):
        df = results.to_pandas()
    else:
        df = pd.DataFrame(results, columns=["score", "SMILES", "Name"])

    return df


def run_legacy_implementation(strategy_name, replicate_num, verbose=False):
    """Run legacy implementation."""

    if verbose:
        print(f"Running legacy {strategy_name} - Replicate {replicate_num + 1}")

    if strategy_name == "Greedy":
        # Use legacy ThompsonSampler for Greedy
        sampler = LegacyThompsonSampler(mode="minimize")

        # Set evaluator (legacy uses json.dumps and ref_colname)
        # Set evaluator - convert config to legacy format
        legacy_config = {
            "ref_filename": EVALUATOR_CONFIG["ref_filename"],
            "compound_col": EVALUATOR_CONFIG["compound_col"],
            "ref_colname": EVALUATOR_CONFIG["score_col"]  # Legacy expects "ref_colname"
        }
        evaluator = LegacyLookupEvaluator(legacy_config)
        sampler.set_evaluator(evaluator)

        # Read reagents
        sampler.read_reagents(REAGENT_FILES)

        # NOTE: No reaction SMARTS needed for LookupEvaluator

        # Set progress visibility
        sampler.set_hide_progress(not verbose)

        # Run warmup
        sampler.warm_up(num_warmup_trials=NUM_WARMUP_TRIALS)

        # Run search
        results = sampler.search(num_cycles=NUM_TS_ITERATIONS)

        # Convert to DataFrame
        df = pd.DataFrame(results, columns=["score", "SMILES", "Name"])

    elif strategy_name == "RouletteWheel":
        # Use legacy RWSSampler for Roulette Wheel
        sampler = LegacyRWSSampler(processes=1, scaling=-1)  # scaling=-1 for minimization

        # Set evaluator (legacy uses json.dumps and ref_colname)
        legacy_config = {
            "ref_filename": EVALUATOR_CONFIG["ref_filename"],
            "compound_col": EVALUATOR_CONFIG["compound_col"],
            "ref_colname": EVALUATOR_CONFIG["score_col"]  # Legacy expects "ref_colname"
        }
        evaluator = LegacyLookupEvaluator(legacy_config)
        sampler.set_evaluator(evaluator)

        # Read reagents
        sampler.read_reagents(REAGENT_FILES)

        # NOTE: No reaction SMARTS needed for LookupEvaluator

        # Set progress visibility
        sampler.set_hide_progress(not verbose)

        # Create temp file for results
        temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv')
        temp_filename = temp_file.name
        temp_file.close()

        try:
            # Run warmup (writes to results file)
            num_warmup_evals = sampler.warm_up(num_warmup_trials=NUM_WARMUP_TRIALS, results_filename=temp_filename)

            # Run search (appends to results file)
            # The search method subtracts num_warmup from the target internally
            # So we need: total_evaluations = num_warmup + NUM_TS_ITERATIONS
            total_target = num_warmup_evals + NUM_TS_ITERATIONS
            percent = total_target / sampler.get_num_prods()

            sampler.search(
                percent_of_library=percent,
                stop=6000,  # Match original RWS config (was 1000, causing premature stopping)
                min_cpds_per_core=50,  # Match original RWS config (was 1, causing inefficiency)
                results_filename=temp_filename
            )

            # Read results from file
            df = pd.read_csv(temp_filename)

            # Log actual number of evaluations performed
            total_evals = len(df)
            search_evals = total_evals - num_warmup_evals
            if verbose:
                print(f"  RWS evaluations: warmup={num_warmup_evals}, search={search_evals}, total={total_evals}")
                print(f"  Target was {NUM_TS_ITERATIONS} search evaluations")

            # Legacy RWS includes warmup in results, but we want to match the evaluation budget
            # Keep only the search results (after warmup)
            if len(df) > num_warmup_evals:
                df = df.iloc[num_warmup_evals:num_warmup_evals + NUM_TS_ITERATIONS].reset_index(drop=True)
                # Warn if we got fewer evaluations than expected
                if search_evals < NUM_TS_ITERATIONS and verbose:
                    print(f"  WARNING: RWS stopped early! Got {search_evals}/{NUM_TS_ITERATIONS} evaluations")

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
    print(f"RUNNING STRATEGY COMPARISON BENCHMARK")
    print(f"{'='*80}")
    print(f"Strategies: Greedy, Roulette Wheel, UCB, ε-Greedy, Bayes-UCB")
    print(f"Implementations: Legacy (Greedy, RW only), Current (all)")
    print(f"Replicates per configuration: {NUM_CYCLES}")
    print(f"Iterations per replicate: {NUM_TS_ITERATIONS}")
    print(f"Warmup trials: {NUM_WARMUP_TRIALS}")

    # Results storage - store full dataframes for recovery calculation
    all_dataframes = []

    # Strategy configurations - Focus on Greedy and RouletteWheel for legacy comparison
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

        # Run legacy implementation (only if available)
        if config["legacy_name"] is not None:
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
        else:
            print(f"\n  LEGACY implementation not available for {strategy_name} (Current only)")

    # Combine all dataframes
    combined_df = pd.concat(all_dataframes, ignore_index=True)

    return combined_df, top_100_ref


def create_grouped_barplot(results_df, top_100_ref):
    """Create multi-panel grouped bar plots for all strategies."""

    # Get reference products
    reference_products = set(top_100_ref["Name"].values)

    # Get numbered cycles (0-9)
    numbered_cycles = sorted([c for c in results_df["cycle"].unique()])
    all_cycles = [str(c) for c in numbered_cycles] + ["concat"]

    # Get all strategies
    all_strategies = sorted(results_df["strategy"].unique())

    # Calculate match counts for each cycle, implementation, and strategy
    match_counts = []

    for implementation in results_df["implementation"].unique():
        for strategy in all_strategies:
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

    # Create strategy comparison plot (Current implementation only)
    current_data = match_counts_df[match_counts_df["implementation"] == "Current"]

    # Create single panel plot with all strategies
    fig, ax = plt.subplots(1, 1, figsize=(18, 10))
    fig.suptitle("Strategy Comparison: Recovery of Top 100 Actives (Current Implementation)",
                 fontsize=18, fontweight='bold', y=0.98)

    # Colors for strategies (expanded palette)
    strategy_colors = {
        "TS_Greedy": "#e41a1c",
        "TS_RouletteWheel": "#377eb8",
        "TS_UCB": "#4daf4a",
        "TS_EpsilonGreedy": "#984ea3",
        "TS_BayesUCB": "#ff7f00"
    }

    sns.barplot(data=current_data, x="cycle", y="match_count", hue="strategy",
                palette=strategy_colors, ax=ax, order=all_cycles)

    # Add value labels on top of bars
    for container in ax.containers:
        ax.bar_label(container, fontsize=9, padding=3)

    ax.set_ylabel("Count of Recovered Actives (out of 100)", fontsize=16)
    ax.set_xlabel("Cycle", fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)

    # Update legend labels to be more readable
    handles, labels = ax.get_legend_handles_labels()
    nice_labels = [l.replace("TS_", "").replace("RouletteWheel", "Roulette Wheel").replace("BayesUCB", "Bayes-UCB").replace("EpsilonGreedy", "ε-Greedy") for l in labels]
    ax.legend(handles, nice_labels, title="Strategy", bbox_to_anchor=(1.02, 1),
              loc='upper left', fontsize=12, title_fontsize=14)

    ax.set_ylim([0, 105])
    ax.axhline(y=100, color='red', linestyle='--', alpha=0.3, linewidth=1)

    plt.tight_layout()

    # Also create legacy comparison plot if legacy data exists
    if "Legacy" in match_counts_df["implementation"].values:
        legacy_strategies = ["TS_Greedy", "TS_RouletteWheel"]

        # Add mean±std as another "cycle" category
        numeric_cycles_legacy = match_counts_df[
            (match_counts_df["cycle"] != "concat") &
            (match_counts_df["strategy"].isin(legacy_strategies))
        ]
        mean_recovery_legacy = numeric_cycles_legacy.groupby(["implementation", "strategy"])["match_count"].agg(["mean", "std"]).reset_index()

        # Add mean values as "mean" cycle to the match_counts list
        legacy_mean_data = []
        for _, row in mean_recovery_legacy.iterrows():
            legacy_mean_data.append({
                "cycle": "mean",
                "strategy": row["strategy"],
                "implementation": row["implementation"],
                "match_count": row["mean"],
                "match_std": row["std"]
            })

        # Recreate match_counts_df with mean data
        match_counts_with_mean_df = pd.concat([
            match_counts_df,
            pd.DataFrame(legacy_mean_data)
        ], ignore_index=True)

        # Update cycle order to include mean
        all_cycles_with_mean = all_cycles + ["mean"]

        fig_legacy, axes = plt.subplots(2, 1, figsize=(20, 12))
        fig_legacy.suptitle("Legacy vs Current Implementation: Recovery of Top 100 Actives\n(Greedy and Roulette Wheel only)",
                     fontsize=18, fontweight='bold', y=0.995)

        # Top panel: Legacy
        ax = axes[0]
        legacy_data = match_counts_with_mean_df[
            (match_counts_with_mean_df["implementation"] == "Legacy") &
            (match_counts_with_mean_df["strategy"].isin(legacy_strategies))
        ]

        sns.barplot(data=legacy_data, x="cycle", y="match_count", hue="strategy",
                    palette=strategy_colors, ax=ax, order=all_cycles_with_mean)

        for container in ax.containers:
            ax.bar_label(container, fontsize=12, padding=3)

        ax.set_ylabel("Count of Recovered Actives (out of 100)", fontsize=16)
        ax.set_xlabel("Cycle", fontsize=16)
        ax.set_title("LEGACY Implementation", fontsize=16, fontweight='bold', pad=15)
        ax.tick_params(axis='both', which='major', labelsize=14)

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, ["Greedy", "Roulette Wheel"],
                  title="Strategy", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=14)

        ax.set_ylim([0, 105])
        ax.axhline(y=100, color='red', linestyle='--', alpha=0.3, linewidth=1)
        ax.grid(axis='y', alpha=0.3)

        # Add error bars for the "mean" bars
        mean_data_legacy = legacy_data[legacy_data["cycle"] == "mean"]
        for i, strategy in enumerate(legacy_strategies):
            strategy_data = mean_data_legacy[mean_data_legacy["strategy"] == strategy]
            if len(strategy_data) > 0 and "match_std" in strategy_data.columns:
                mean_idx = all_cycles_with_mean.index("mean")
                offset = (i - 0.5) * 0.8 / len(legacy_strategies)
                std_val = strategy_data["match_std"].values[0]
                mean_val = strategy_data["match_count"].values[0]
                ax.errorbar(mean_idx + offset, mean_val, yerr=std_val, fmt='none',
                           ecolor='black', capsize=5, capthick=2, linewidth=2)

        # Bottom panel: Current
        ax = axes[1]
        current_legacy_comp = match_counts_with_mean_df[
            (match_counts_with_mean_df["implementation"] == "Current") &
            (match_counts_with_mean_df["strategy"].isin(legacy_strategies))
        ]

        sns.barplot(data=current_legacy_comp, x="cycle", y="match_count", hue="strategy",
                    palette=strategy_colors, ax=ax, order=all_cycles_with_mean)

        for container in ax.containers:
            ax.bar_label(container, fontsize=12, padding=3)

        ax.set_ylabel("Count of Recovered Actives (out of 100)", fontsize=16)
        ax.set_xlabel("Cycle", fontsize=16)
        ax.set_title("CURRENT Implementation", fontsize=16, fontweight='bold', pad=15)
        ax.tick_params(axis='both', which='major', labelsize=14)

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, ["Greedy", "Roulette Wheel"],
                  title="Strategy", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=14)

        ax.set_ylim([0, 105])
        ax.axhline(y=100, color='red', linestyle='--', alpha=0.3, linewidth=1)
        ax.grid(axis='y', alpha=0.3)

        # Add error bars for the "mean" bars
        mean_data_current = current_legacy_comp[current_legacy_comp["cycle"] == "mean"]
        for i, strategy in enumerate(legacy_strategies):
            strategy_data = mean_data_current[mean_data_current["strategy"] == strategy]
            if len(strategy_data) > 0 and "match_std" in strategy_data.columns:
                mean_idx = all_cycles_with_mean.index("mean")
                offset = (i - 0.5) * 0.8 / len(legacy_strategies)
                std_val = strategy_data["match_std"].values[0]
                mean_val = strategy_data["match_count"].values[0]
                ax.errorbar(mean_idx + offset, mean_val, yerr=std_val, fmt='none',
                           ecolor='black', capsize=5, capthick=2, linewidth=2)

        plt.tight_layout()

        return fig, fig_legacy

    return fig, None


def create_summary_table(results_df, top_100_ref):
    """Create summary statistics table."""

    reference_products = set(top_100_ref["Name"].values)

    summary_data = []

    # Get all implementations and strategies
    for implementation in sorted(results_df["implementation"].unique()):
        for strategy in sorted(results_df["strategy"].unique()):
            # Get all cycles for this combination
            strategy_data = results_df[
                (results_df["implementation"] == implementation) &
                (results_df["strategy"] == strategy)
            ]

            if len(strategy_data) == 0:
                continue

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

            # Nice strategy name
            nice_name = (strategy.replace("TS_", "")
                        .replace("RouletteWheel", "Roulette Wheel")
                        .replace("BayesUCB", "Bayes-UCB")
                        .replace("EpsilonGreedy", "ε-Greedy"))

            summary_data.append({
                "Implementation": implementation,
                "Strategy": nice_name,
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
    print("STRATEGY COMPARISON BENCHMARK")
    print("="*80)
    print("\nThis script compares multiple Thompson Sampling selection strategies:")
    print("  • Greedy (Legacy & Current)")
    print("  • Roulette Wheel with Thermal Cycling (Legacy & Current)")
    print("  • UCB (Current only)")
    print("  • ε-Greedy (Current only)")
    print("  • Bayes-UCB with Thermal Cycling (Current only)")
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
    summary.to_csv("strategy_comparison_summary.csv", index=False)
    print("\n✅ Summary saved to: strategy_comparison_summary.csv")

    # Save detailed results
    results_df.to_csv("strategy_comparison_detailed.csv", index=False)
    print("✅ Detailed results saved to: strategy_comparison_detailed.csv")

    # Create plots
    print("\n📊 Generating comparison plots...")
    fig_main, fig_legacy = create_grouped_barplot(results_df, top_100_ref)

    # Save main strategy comparison plot
    fig_main.savefig("strategy_comparison.pdf", dpi=300, bbox_inches='tight')
    print("✅ Main comparison plot saved to: strategy_comparison.png")

    # Save legacy comparison plot if it exists
    if fig_legacy is not None:
        fig_legacy.savefig("legacy_vs_current_comparison.png", dpi=300, bbox_inches='tight')
        print("✅ Legacy comparison plot saved to: legacy_vs_current_comparison.png")

    print("\n" + "="*80)
    print("✅ BENCHMARK COMPLETE!")
    print("="*80)

    plt.show()


if __name__ == "__main__":
    main()
