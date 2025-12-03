"""
Legacy vs Current Thompson Sampling Comparison - ROCS Similarity

Compares legacy vs current Thompson Sampling implementation on ROCS similarity scores:
- Greedy (Legacy & Current)
- Roulette Wheel with Thermal Cycling (Legacy & Current)

This benchmark evaluates recovery of top 100 most similar compounds from ROCS,
using maximization mode (higher scores = more similar).
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

# 3-component Quinazoline library
REAGENT_FILES = [
    "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/reagents/thrombin/quinazoline/niementowski_reagent_0.smi",
    "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/reagents/thrombin/quinazoline/niementowski_reagent_1.smi",
    "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/reagents/thrombin/quinazoline/niementowski_reagent_2.smi"
]

# NOTE: No reaction SMARTS needed (uses product codes)

EVALUATOR_CONFIG = {
    "ref_filename": "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/reagents/thrombin/quinazoline/niementowski.parquet",
    "compound_col": "Name",
    "score_col": "query_001"
}

# Reference data for recovery calculation
REFERENCE_FILE = "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/reagents/thrombin/quinazoline/niementowski.parquet"

# Experiment parameters
NUM_WARMUP_TRIALS = 5
NUM_TS_ITERATIONS = 18500
NUM_CYCLES = 10  # Number of replicates

# ============================================================================
# Helper Functions
# ============================================================================

def load_reference_data():
    """Load ROCS reference data and get top 100 most similar."""
    print("Loading reference (ROCS) data...")

    # Read parquet file
    ref_df = pd.read_parquet(REFERENCE_FILE)

    # Standardize column names
    score_col = EVALUATOR_CONFIG["score_col"]
    compound_col = EVALUATOR_CONFIG["compound_col"]

    if score_col != "score":
        ref_df = ref_df.rename(columns={score_col: "score"})
    if compound_col != "Name":
        ref_df = ref_df.rename(columns={compound_col: "Name"})

    # Sort in descending order (maximize) and get top 100
    ref_df = ref_df.sort_values("score", ascending=False)
    top_100_ref = ref_df.head(100).copy()

    print(f"  Total products in reference: {len(ref_df):,}")
    print(f"  Best score in reference: {ref_df['score'].max():.3f}")
    print(f"  Top 100 score range: {top_100_ref['score'].max():.3f} to {top_100_ref['score'].min():.3f}")

    return ref_df, top_100_ref


def run_current_implementation(strategy_name, selection_strategy, replicate_num, verbose=False):
    """Run current implementation."""

    if verbose:
        print(f"Running current {strategy_name} - Replicate {replicate_num + 1}")

    # Choose warmup strategy and Boltzmann weighting based on selection strategy
    if strategy_name == "RouletteWheel":
        warmup_strategy = EnhancedWarmup()
        use_boltzmann = True  # Enable Boltzmann weighting for RWS
    else:
        warmup_strategy = EnhancedWarmup()  # Also use EnhancedWarmup for Greedy for fair comparison
        use_boltzmann = False

    # Create sampler with current implementation
    sampler = ThompsonSampler(
        selection_strategy=selection_strategy,
        warmup_strategy=warmup_strategy,
        batch_size=1,
        use_boltzmann_weighting=use_boltzmann
    )

    # Setup
    evaluator = LookupEvaluator(EVALUATOR_CONFIG)
    sampler.set_evaluator(evaluator)
    sampler.read_reagents(REAGENT_FILES)
    # NOTE: No reaction SMARTS needed for LookupEvaluator
    sampler.set_hide_progress(not verbose)

    # Warmup
    sampler.warm_up(num_warmup_trials=NUM_WARMUP_TRIALS)

    # Search
    results = sampler.search(num_cycles=NUM_TS_ITERATIONS)

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
        sampler = LegacyThompsonSampler(mode="maximize")

        # Set evaluator - convert config to legacy format
        legacy_config = {
            "ref_filename": EVALUATOR_CONFIG["ref_filename"],
            "ref_colname": EVALUATOR_CONFIG["score_col"],  # Legacy expects "ref_colname"
            "compound_col": EVALUATOR_CONFIG["compound_col"]
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
        sampler = LegacyRWSSampler(processes=1, scaling=1)  # scaling=1 for maximization

        # Set evaluator - convert config to legacy format
        legacy_config = {
            "ref_filename": EVALUATOR_CONFIG["ref_filename"],
            "ref_colname": EVALUATOR_CONFIG["score_col"],  # Legacy expects "ref_colname"
            "compound_col": EVALUATOR_CONFIG["compound_col"]
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
            num_prods = sampler.get_num_prods()
            total_target = num_warmup_evals + NUM_TS_ITERATIONS
            percent = total_target / num_prods

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

    # Get top N from TS results (maximize = nlargest)
    ts_top_n = ts_results.nlargest(top_n, "score")

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
    print(f"RUNNING LEGACY VS CURRENT COMPARISON - ROCS SIMILARITY")
    print(f"{'='*80}")
    print(f"Strategies: Greedy, Roulette Wheel")
    print(f"Implementations: Legacy vs Current")
    print(f"Replicates per configuration: {NUM_CYCLES}")
    print(f"Iterations per replicate: {NUM_TS_ITERATIONS}")
    print(f"Warmup trials: {NUM_WARMUP_TRIALS}")
    print(f"Mode: MAXIMIZE (higher ROCS scores = more similar)")

    # Results storage
    all_dataframes = []

    # Strategy configurations
    strategies = {
        "Greedy": {
            "current_strategy": GreedySelection(mode="maximize"),
            "legacy_name": "Greedy"
        },
        "RouletteWheel": {
            "current_strategy": RouletteWheelSelection(mode="maximize", alpha=0.1, beta=0.1),
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


def create_comparison_plots(results_df, top_100_ref):
    """Create comparison plots for legacy vs current."""

    # Get reference products
    reference_products = set(top_100_ref["Name"].values)

    # Get numbered cycles
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
                    # Unique recovered compounds across all cycles
                    cycle_data = results_df[
                        (results_df["implementation"] == implementation) &
                        (results_df["strategy"] == strategy)
                    ]
                    # Get top 100 from each cycle, then find unique
                    top_100_per_cycle = cycle_data.groupby("cycle").apply(
                        lambda x: x.nlargest(100, "score")
                    ).reset_index(drop=True)
                    unique_names = top_100_per_cycle["Name"].unique()
                    match_count = len(set(unique_names).intersection(reference_products))
                else:
                    # Individual cycle
                    cycle_data = results_df[
                        (results_df["implementation"] == implementation) &
                        (results_df["strategy"] == strategy) &
                        (results_df["cycle"] == int(cycle))
                    ]
                    if len(cycle_data) > 0:
                        top_100 = cycle_data.nlargest(100, "score")
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

    # Add mean±std as another "cycle" category
    numeric_cycles = match_counts_df[match_counts_df["cycle"] != "concat"]
    mean_recovery = numeric_cycles.groupby(["implementation", "strategy"])["match_count"].agg(["mean", "std"]).reset_index()

    # Add mean values as "mean" cycle
    for _, row in mean_recovery.iterrows():
        match_counts.append({
            "cycle": "mean",
            "strategy": row["strategy"],
            "implementation": row["implementation"],
            "match_count": row["mean"],
            "match_std": row["std"]
        })

    match_counts_df = pd.DataFrame(match_counts)

    # Update cycle order to include mean
    all_cycles_with_mean = [str(c) for c in numbered_cycles] + ["concat", "mean"]

    # Create legacy comparison plot
    fig, axes = plt.subplots(2, 1, figsize=(20, 12))
    fig.suptitle("Legacy vs Current Implementation: ROCS Similarity Recovery\n(Top 100 Most Similar Compounds)",
                 fontsize=18, fontweight='bold', y=0.995)

    # Use Set1 color palette
    colors = sns.color_palette("Set1", n_colors=8)
    strategy_colors = {
        "Greedy": colors[0],
        "RouletteWheel": colors[1]
    }

    all_strategies = ["Greedy", "RouletteWheel"]

    # Top panel: Legacy
    ax = axes[0]
    legacy_data = match_counts_df[match_counts_df["implementation"] == "Legacy"]

    sns.barplot(data=legacy_data, x="cycle", y="match_count", hue="strategy",
                palette=strategy_colors, ax=ax, order=all_cycles_with_mean)

    for container in ax.containers:
        ax.bar_label(container, fontsize=10, padding=3)

    ax.set_ylabel("Count of Recovered Similar Compounds (out of 100)", fontsize=14)
    ax.set_xlabel("Cycle (Replicate Number)", fontsize=14)
    ax.set_title("LEGACY Implementation", fontsize=16, fontweight='bold', pad=15)
    ax.tick_params(axis='both', which='major', labelsize=12)

    handles, labels = ax.get_legend_handles_labels()
    nice_labels = [l.replace("RouletteWheel", "Roulette Wheel") for l in labels]
    ax.legend(handles, nice_labels, title="Strategy", bbox_to_anchor=(1.02, 1),
              loc='upper left', fontsize=12, title_fontsize=13)

    ax.set_ylim([0, 105])
    ax.axhline(y=100, color='red', linestyle='--', alpha=0.3, linewidth=1)
    ax.grid(axis='y', alpha=0.3)

    # Add error bars for the "mean" bars
    mean_data = legacy_data[legacy_data["cycle"] == "mean"]
    for i, strategy in enumerate(all_strategies):
        strategy_data = mean_data[mean_data["strategy"] == strategy]
        if len(strategy_data) > 0 and "match_std" in strategy_data.columns:
            # Find the x position of the mean bar
            mean_idx = all_cycles_with_mean.index("mean")
            # seaborn uses different x positions for grouped bars
            # Need to calculate offset based on number of strategies
            offset = (i - 0.5) * 0.8 / len(all_strategies)
            std_val = strategy_data["match_std"].values[0]
            mean_val = strategy_data["match_count"].values[0]
            ax.errorbar(mean_idx + offset, mean_val, yerr=std_val, fmt='none',
                       ecolor='black', capsize=5, capthick=2, linewidth=2)

    # Bottom panel: Current
    ax = axes[1]
    current_data = match_counts_df[match_counts_df["implementation"] == "Current"]

    sns.barplot(data=current_data, x="cycle", y="match_count", hue="strategy",
                palette=strategy_colors, ax=ax, order=all_cycles_with_mean)

    for container in ax.containers:
        ax.bar_label(container, fontsize=10, padding=3)

    ax.set_ylabel("Count of Recovered Similar Compounds (out of 100)", fontsize=14)
    ax.set_xlabel("Cycle (Replicate Number)", fontsize=14)
    ax.set_title("CURRENT Implementation", fontsize=16, fontweight='bold', pad=15)
    ax.tick_params(axis='both', which='major', labelsize=12)

    handles, labels = ax.get_legend_handles_labels()
    nice_labels = [l.replace("RouletteWheel", "Roulette Wheel") for l in labels]
    ax.legend(handles, nice_labels, title="Strategy", bbox_to_anchor=(1.02, 1),
              loc='upper left', fontsize=12, title_fontsize=13)

    ax.set_ylim([0, 105])
    ax.axhline(y=100, color='red', linestyle='--', alpha=0.3, linewidth=1)
    ax.grid(axis='y', alpha=0.3)

    # Add error bars for the "mean" bars
    mean_data_current = current_data[current_data["cycle"] == "mean"]
    for i, strategy in enumerate(all_strategies):
        strategy_data = mean_data_current[mean_data_current["strategy"] == strategy]
        if len(strategy_data) > 0 and "match_std" in strategy_data.columns:
            mean_idx = all_cycles_with_mean.index("mean")
            offset = (i - 0.5) * 0.8 / len(all_strategies)
            std_val = strategy_data["match_std"].values[0]
            mean_val = strategy_data["match_count"].values[0]
            ax.errorbar(mean_idx + offset, mean_val, yerr=std_val, fmt='none',
                       ecolor='black', capsize=5, capthick=2, linewidth=2)

    plt.tight_layout()

    return fig


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
                top_100 = cycle_data.nlargest(100, "score")

                # Recovery
                matches = len(set(top_100["Name"].values).intersection(reference_products))
                recoveries.append(matches)

                # Best score (maximize = max)
                best_scores.append(cycle_data["score"].max())

                # Top 100 mean
                top_100_means.append(top_100["score"].mean())

            # Nice strategy name
            nice_name = strategy.replace("RouletteWheel", "Roulette Wheel")

            summary_data.append({
                "Implementation": implementation,
                "Strategy": nice_name,
                "Recovered (mean±std)": f"{np.mean(recoveries):.1f}±{np.std(recoveries):.1f}",
                "Recovered (min-max)": f"{int(np.min(recoveries))}-{int(np.max(recoveries))}",
                "Recovery %": f"{(np.mean(recoveries)/100)*100:.1f}%",
                "Best Score": f"{np.mean(best_scores):.3f}±{np.std(best_scores):.3f}",
                "Top 100 Mean": f"{np.mean(top_100_means):.3f}±{np.std(top_100_means):.3f}"
            })

    return pd.DataFrame(summary_data)


# ============================================================================
# Main Execution
# ============================================================================

def main():
    """Run the full comparison."""

    print("="*80)
    print("LEGACY VS CURRENT COMPARISON - ROCS SIMILARITY")
    print("="*80)
    print("\nThis script compares Legacy vs Current Thompson Sampling")
    print("implementations for recovering the top 100 most similar compounds")
    print("based on ROCS similarity scores.")
    print("\nStrategies tested:")
    print("  • Greedy (Legacy & Current)")
    print("  • Roulette Wheel with Thermal Cycling (Legacy & Current)")

    # Run all experiments
    results_df, top_100_ref = run_all_experiments()

    # Create summary table
    print("\n" + "="*80)
    print("SUMMARY TABLE")
    print("="*80)
    summary = create_summary_table(results_df, top_100_ref)
    print(summary.to_string(index=False))

    # Save summary
    summary.to_csv("legacy_vs_current_rocs_summary.csv", index=False)
    print("\n✅ Summary saved to: legacy_vs_current_rocs_summary.csv")

    # Save detailed results
    results_df.to_csv("legacy_vs_current_rocs_detailed.csv", index=False)
    print("✅ Detailed results saved to: legacy_vs_current_rocs_detailed.csv")

    # Create plots
    print("\n📊 Generating comparison plots...")
    fig = create_comparison_plots(results_df, top_100_ref)
    fig.savefig("legacy_vs_current_rocs_comparison.pdf", dpi=300, bbox_inches='tight')
    print("✅ Comparison plot saved to: legacy_vs_current_rocs_comparison.png")

    # Print key findings
    print("\n" + "="*80)
    print("KEY FINDINGS")
    print("="*80)

    for strategy in ["Greedy", "RouletteWheel"]:
        print(f"\n{strategy.replace('RouletteWheel', 'Roulette Wheel')}:")
        legacy_recovery = results_df[
            (results_df["strategy"] == strategy) &
            (results_df["implementation"] == "Legacy")
        ].groupby("cycle").apply(
            lambda x: len(set(x.nlargest(100, "score")["Name"].values).intersection(set(top_100_ref["Name"].values)))
        )
        current_recovery = results_df[
            (results_df["strategy"] == strategy) &
            (results_df["implementation"] == "Current")
        ].groupby("cycle").apply(
            lambda x: len(set(x.nlargest(100, "score")["Name"].values).intersection(set(top_100_ref["Name"].values)))
        )

        print(f"  Legacy:  {legacy_recovery.mean():.1f} ± {legacy_recovery.std():.1f}")
        print(f"  Current: {current_recovery.mean():.1f} ± {current_recovery.std():.1f}")
        diff = current_recovery.mean() - legacy_recovery.mean()
        print(f"  Difference: {diff:+.1f} compounds ({(diff/legacy_recovery.mean())*100:+.1f}%)")

    print("\n" + "="*80)
    print("✅ COMPARISON COMPLETE!")
    print("="*80)

    plt.show()


if __name__ == "__main__":
    main()
