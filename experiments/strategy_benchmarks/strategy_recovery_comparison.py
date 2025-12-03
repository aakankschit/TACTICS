"""
Selection Strategy Recovery Comparison

This script compares how well different selection strategies recover the top 100
most similar compounds, with and without CATS enhancement.

Experimental Design:
- 5 strategies: Greedy, RouletteWheel, UCB, EpsilonGreedy, BayesUCB
- 2 conditions per strategy: ±CATS
- 10 trials per condition (no triplicates)
- 5000 iterations per trial

Visualizations:
- Per-strategy comparison with two plots:
  1. Concatenated recovery: Individual cycles + cumulative unique compounds
  2. Mean recovery ± std with percentage improvement annotation
"""

import sys
import json
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

# Add TACTICS to path if needed
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from TACTICS.thompson_sampling import (
    ThompsonSampler,
    GreedySelection,
    RouletteWheelSelection,
    UCBSelection,
    EpsilonGreedySelection,
    BayesUCBSelection,
    LookupEvaluator
)
from TACTICS.thompson_sampling.cats import CATSConfig, CATSManager
from TACTICS.thompson_sampling.warmup import StandardWarmup

# ============================================================================
# Configuration
# ============================================================================

REAGENT_FILES = [
    "/Users/aakankschitnandkeolyar/Desktop/TACTICS/examples/input_files/quinazoline/niementowski_reagent_0.smi",
    "/Users/aakankschitnandkeolyar/Desktop/TACTICS/examples/input_files/quinazoline/niementowski_reagent_1.smi",
    "/Users/aakankschitnandkeolyar/Desktop/TACTICS/examples/input_files/quinazoline/niementowski_reagent_2.smi"
]

EVALUATOR_CONFIG = {
    "ref_filename": "/Users/aakankschitnandkeolyar/Desktop/TACTICS/examples/input_files/quinazoline/niementowski.parquet",
    "compound_col": "Name",
    "score_col": "query_001"
}

REFERENCE_FILE = "/Users/aakankschitnandkeolyar/Desktop/TACTICS/examples/input_files/quinazoline/niementowski.parquet"

# Experiment parameters
NUM_WARMUP_TRIALS = 5
NUM_TS_ITERATIONS = 18500
NUM_TRIALS = 10  # No triplicates, just 10 trials
BATCH_SIZE = 1

# CATS configuration
CATS_CONFIG = CATSConfig(
    enabled=True,
    phase_boundaries=(0.2, 0.6),
    update_frequency=100,
    component_names=["AminoBenzoic","Amines","Acids"],
    verbose=False
)

# Use seaborn Set1 color palette
COLORS = sns.color_palette("Set1", n_colors=8)

# ============================================================================
# Helper Functions
# ============================================================================

def load_reference_data():
    """Load brute force reference data and get top 100."""
    print("Loading reference (brute force) data...")

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


def run_experiment(strategy_name, selection_strategy, use_cats, trial_num, verbose=False):
    """
    Run a single experiment.

    Args:
        strategy_name: Name of strategy
        selection_strategy: SelectionStrategy instance
        use_cats: Whether to enable CATS
        trial_num: Trial number
        verbose: Print progress

    Returns:
        pd.DataFrame: Results with scores, SMILES, and names
    """
    if verbose:
        cats_str = "WITH CATS" if use_cats else "NO CATS"
        print(f"\nRunning {strategy_name} {cats_str} - Trial {trial_num + 1}")

    # Initialize CATS manager if enabled
    cats_manager = None
    if use_cats:
        cats_manager = CATSManager(
            n_components=3,
            target_compounds=NUM_TS_ITERATIONS,
            mode='maximize',
            phase_boundaries=CATS_CONFIG.phase_boundaries,
            update_frequency=CATS_CONFIG.update_frequency,
            component_names=CATS_CONFIG.component_names,
            verbose=CATS_CONFIG.verbose
        )

    # Create sampler with optional CATS integration
    sampler = ThompsonSampler(
        selection_strategy=selection_strategy,
        warmup_strategy=StandardWarmup(),
        batch_size=BATCH_SIZE,
        cats_manager=cats_manager
    )

    # Setup evaluator and reagents
    evaluator = LookupEvaluator(EVALUATOR_CONFIG)
    sampler.set_evaluator(evaluator)
    sampler.read_reagents(REAGENT_FILES)
    sampler.set_hide_progress(not verbose)

    # Warmup
    sampler.warm_up(num_warmup_trials=NUM_WARMUP_TRIALS)

    # Search
    results = sampler.search(num_cycles=NUM_TS_ITERATIONS)

    # Convert to DataFrame
    if hasattr(results, 'to_pandas'):
        df = results.to_pandas()
    else:
        df = pd.DataFrame(results, columns=["score", "SMILES", "Name"])

    return df


def calculate_recovery(ts_results, reference_top_100, top_n=100):
    """Calculate how many of the top N reference compounds were recovered."""
    ts_top_n = ts_results.nlargest(top_n, "score")
    ref_names = set(reference_top_100["Name"].values)
    ts_names = set(ts_top_n["Name"].values)
    matches = len(ref_names.intersection(ts_names))
    fraction = matches / len(ref_names)
    return matches, fraction


# ============================================================================
# Main Experiment Loop
# ============================================================================

def run_all_experiments():
    """Run all experimental conditions."""

    # Load reference data
    ref_df, top_100_ref = load_reference_data()

    print(f"\n{'='*80}")
    print(f"SELECTION STRATEGY RECOVERY COMPARISON")
    print(f"{'='*80}")
    print(f"Strategies: Greedy, RouletteWheel, UCB, EpsilonGreedy, BayesUCB")
    print(f"Conditions: ±CATS (2 per strategy)")
    print(f"Trials per condition: {NUM_TRIALS}")
    print(f"Iterations per trial: {NUM_TS_ITERATIONS}")
    print(f"Warmup trials: {NUM_WARMUP_TRIALS}")
    print(f"Total runs: {NUM_TRIALS} × 10 conditions = {NUM_TRIALS * 10}")

    # Define experimental conditions with strategy factories
    conditions = [
        # Greedy
        {
            "name": "Greedy",
            "strategy_factory": lambda: GreedySelection(mode="maximize"),
            "use_cats": False,
            "cats_label": "NoCATS",
            "color_idx": 0
        },
        {
            "name": "Greedy",
            "strategy_factory": lambda: GreedySelection(mode="maximize"),
            "use_cats": True,
            "cats_label": "CATS",
            "color_idx": 1
        },
        # Roulette Wheel
        {
            "name": "RouletteWheel",
            "strategy_factory": lambda: RouletteWheelSelection(mode="maximize", alpha=0.1, beta=0.1),
            "use_cats": False,
            "cats_label": "NoCATS",
            "color_idx": 0
        },
        {
            "name": "RouletteWheel",
            "strategy_factory": lambda: RouletteWheelSelection(mode="maximize", alpha=0.1, beta=0.1),
            "use_cats": True,
            "cats_label": "CATS",
            "color_idx": 1
        },
        # UCB
        {
            "name": "UCB",
            "strategy_factory": lambda: UCBSelection(mode="maximize", c=2.0),
            "use_cats": False,
            "cats_label": "NoCATS",
            "color_idx": 0
        },
        {
            "name": "UCB",
            "strategy_factory": lambda: UCBSelection(mode="maximize", c=2.0),
            "use_cats": True,
            "cats_label": "CATS",
            "color_idx": 1
        },
        # Epsilon Greedy
        {
            "name": "EpsilonGreedy",
            "strategy_factory": lambda: EpsilonGreedySelection(mode="maximize", epsilon=0.1, decay=0.995),
            "use_cats": False,
            "cats_label": "NoCATS",
            "color_idx": 0
        },
        {
            "name": "EpsilonGreedy",
            "strategy_factory": lambda: EpsilonGreedySelection(mode="maximize", epsilon=0.1, decay=0.995),
            "use_cats": True,
            "cats_label": "CATS",
            "color_idx": 1
        },
        # Bayes UCB
        {
            "name": "BayesUCB",
            "strategy_factory": lambda: BayesUCBSelection(mode="maximize", initial_p_high=0.90, initial_p_low=0.90),
            "use_cats": False,
            "cats_label": "NoCATS",
            "color_idx": 0
        },
        {
            "name": "BayesUCB",
            "strategy_factory": lambda: BayesUCBSelection(mode="maximize", initial_p_high=0.90, initial_p_low=0.90),
            "use_cats": True,
            "cats_label": "CATS",
            "color_idx": 1
        }
    ]

    # Storage for results
    all_dataframes = []
    recovery_stats = []

    # Run all conditions
    print(f"\n{'='*80}")
    print(f"RUNNING EXPERIMENTS")
    print(f"{'='*80}")

    for condition in conditions:
        cond_label = f"{condition['name']}_{condition['cats_label']}"
        print(f"\nCondition: {cond_label}")

        for trial in tqdm(range(NUM_TRIALS), desc=f"  {cond_label}"):
            # Create fresh strategy instance for this trial
            selection_strategy = condition['strategy_factory']()

            # Run experiment
            df = run_experiment(
                strategy_name=condition['name'],
                selection_strategy=selection_strategy,
                use_cats=condition['use_cats'],
                trial_num=trial,
                verbose=False
            )

            # Add metadata
            df["strategy"] = condition["name"]
            df["cats"] = condition['cats_label']
            df["trial"] = trial

            all_dataframes.append(df)

            # Calculate recovery
            matches, fraction = calculate_recovery(df, top_100_ref, top_n=100)
            recovery_stats.append({
                "strategy": condition["name"],
                "cats": condition['cats_label'],
                "trial": trial,
                "recovered": matches,
                "fraction": fraction,
                "best_score": df["score"].max(),
                "top_100_mean": df.nlargest(100, "score")["score"].mean()
            })

        # Print summary for this condition
        cond_recoveries = [s["recovered"] for s in recovery_stats
                         if s["strategy"] == condition["name"] and s["cats"] == condition['cats_label']]
        mean_recovery = np.mean(cond_recoveries)
        std_recovery = np.std(cond_recoveries)
        print(f"    Recovery: {mean_recovery:.1f} ± {std_recovery:.1f} out of 100")

    # Combine all results
    combined_df = pd.concat(all_dataframes, ignore_index=True)

    return combined_df, recovery_stats, top_100_ref, conditions


# ============================================================================
# Visualization
# ============================================================================

def create_strategy_comparison_figure(results_df, recovery_stats, top_100_ref, strategy_name):
    """
    Create a 2-panel comparison figure for a single strategy.

    Panel 1: Concatenated recovery (per cycle + cumulative)
    Panel 2: Mean recovery ± std with % improvement
    """

    reference_products = set(top_100_ref["Name"].values)

    # Filter data for this strategy
    strategy_data = results_df[results_df["strategy"] == strategy_name]
    strategy_recovery = [s for s in recovery_stats if s["strategy"] == strategy_name]

    # Get trial numbers
    numbered_trials = sorted(strategy_data["trial"].unique())
    all_trials = [str(t) for t in numbered_trials] + ["concat"]

    # --- Panel 1: Calculate match counts for each trial and concat ---
    match_counts = []

    for trial in all_trials:
        for cats_label in ["NoCATS", "CATS"]:
            if trial == "concat":
                # Unique recovered compounds across all trials
                trial_data = strategy_data[strategy_data["cats"] == cats_label]
                unique_names = trial_data["Name"].unique()
                match_count = len(set(unique_names).intersection(reference_products))
            else:
                # Individual trial
                trial_data = strategy_data[
                    (strategy_data["trial"] == int(trial)) &
                    (strategy_data["cats"] == cats_label)
                ]
                if len(trial_data) > 0:
                    top_100 = trial_data.nlargest(100, "score")
                    names = top_100["Name"].values
                    match_count = len(set(names).intersection(reference_products))
                else:
                    match_count = 0

            match_counts.append({
                "trial": trial,
                "cats": cats_label,
                "match_count": match_count
            })

    match_counts_df = pd.DataFrame(match_counts)

    # --- Panel 2: Calculate mean recovery ± std ---
    recovery_df = pd.DataFrame(strategy_recovery)
    mean_recovery = recovery_df.groupby("cats")["recovered"].agg(["mean", "std"]).reset_index()

    # Calculate % improvement
    no_cats_mean = mean_recovery[mean_recovery["cats"] == "NoCATS"]["mean"].values[0]
    cats_mean = mean_recovery[mean_recovery["cats"] == "CATS"]["mean"].values[0]
    improvement_pct = ((cats_mean - no_cats_mean) / no_cats_mean) * 100 if no_cats_mean > 0 else 0

    # --- Create figure with 2 panels ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
    fig.suptitle(f"{strategy_name} Selection: CATS vs No-CATS Recovery Comparison",
                 fontsize=18, fontweight='bold', y=0.98)

    # --- Panel 1: Per-trial + concatenated recovery ---
    palette = {cats: COLORS[idx] for idx, cats in enumerate(["NoCATS", "CATS"])}

    sns.barplot(data=match_counts_df, x="trial", y="match_count", hue="cats",
                palette=palette, ax=ax1, order=all_trials)

    # Add value labels
    for container in ax1.containers:
        ax1.bar_label(container, fontsize=10, padding=3)

    ax1.set_ylabel("Count of Recovered Actives (out of 100)", fontsize=13)
    ax1.set_xlabel("Trial Number", fontsize=13)
    ax1.set_title("Recovery per Trial + Cumulative (concat)", fontsize=14, fontweight='bold', pad=15)
    ax1.tick_params(axis='both', which='major', labelsize=11)
    ax1.legend(title="Condition", fontsize=11, title_fontsize=12)
    ax1.set_ylim([0, 105])
    ax1.axhline(y=100, color='red', linestyle='--', alpha=0.3, linewidth=1)
    ax1.grid(axis='y', alpha=0.3)

    # --- Panel 2: Mean recovery ± std with % improvement ---
    x_pos = np.arange(len(mean_recovery))

    bars = ax2.bar(x_pos, mean_recovery["mean"],
                   yerr=mean_recovery["std"],
                   color=[palette[cats] for cats in mean_recovery["cats"]],
                   capsize=5, alpha=0.8, error_kw={'linewidth': 2})

    # Add value labels on bars
    for i, (idx, row) in enumerate(mean_recovery.iterrows()):
        ax2.text(i, row["mean"] + row["std"] + 1.5,
                f'{row["mean"]:.1f}±{row["std"]:.1f}',
                ha='center', va='bottom', fontsize=12, fontweight='bold')

    # Add improvement annotation
    ax2.text(0.5, max(mean_recovery["mean"]) + max(mean_recovery["std"]) + 8,
            f'CATS Improvement: {improvement_pct:+.1f}%',
            ha='center', va='bottom', fontsize=13, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow' if improvement_pct > 0 else 'lightcoral', alpha=0.3))

    ax2.set_ylabel("Mean Recovered Actives (out of 100)", fontsize=13)
    ax2.set_xlabel("Condition", fontsize=13)
    ax2.set_title(f"Mean Recovery across {NUM_TRIALS} Trials", fontsize=14, fontweight='bold', pad=15)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(mean_recovery["cats"], fontsize=11)
    ax2.set_ylim([0, 105])
    ax2.axhline(y=100, color='red', linestyle='--', alpha=0.3, linewidth=1)
    ax2.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    return fig


def create_overall_summary_plot(recovery_stats):
    """Create overall summary comparing all strategies."""

    recovery_df = pd.DataFrame(recovery_stats)

    # Calculate mean recovery for each strategy and condition
    summary = recovery_df.groupby(["strategy", "cats"])["recovered"].mean().reset_index()

    # Get unique strategies
    strategies = sorted(recovery_df["strategy"].unique())

    # Create plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
    fig.suptitle("All Strategies: CATS vs No-CATS Overall Performance",
                 fontsize=18, fontweight='bold')

    # Plot 1: Grouped bar chart
    x = np.arange(len(strategies))
    width = 0.35

    no_cats_values = [summary[(summary["strategy"] == s) &
                              (summary["cats"] == "NoCATS")]["recovered"].values[0]
                     for s in strategies]
    cats_values = [summary[(summary["strategy"] == s) &
                           (summary["cats"] == "CATS")]["recovered"].values[0]
                  for s in strategies]

    bars1 = ax1.bar(x - width/2, no_cats_values, width, label='No CATS',
                    color=COLORS[0], alpha=0.8)
    bars2 = ax1.bar(x + width/2, cats_values, width, label='CATS',
                    color=COLORS[1], alpha=0.8)

    ax1.set_ylabel('Mean Recovered Actives (out of 100)', fontsize=13)
    ax1.set_xlabel('Strategy', fontsize=13)
    ax1.set_title('Mean Recovery Comparison', fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(strategies, rotation=15, ha='right')
    ax1.legend(fontsize=12)
    ax1.set_ylim([0, 105])
    ax1.axhline(y=100, color='red', linestyle='--', alpha=0.3, linewidth=1)
    ax1.grid(axis='y', alpha=0.3)

    # Add value labels
    ax1.bar_label(bars1, fmt='%.1f', padding=3, fontsize=9)
    ax1.bar_label(bars2, fmt='%.1f', padding=3, fontsize=9)

    # Plot 2: CATS improvement
    improvements = [cats_values[i] - no_cats_values[i] for i in range(len(strategies))]
    colors = [COLORS[1] if imp > 0 else COLORS[3] for imp in improvements]

    bars3 = ax2.bar(strategies, improvements, color=colors, alpha=0.8)
    ax2.set_ylabel('Improvement in Recovered Actives', fontsize=13)
    ax2.set_xlabel('Strategy', fontsize=13)
    ax2.set_title('CATS Benefit (CATS - No CATS)', fontsize=14, fontweight='bold')
    ax2.set_xticklabels(strategies, rotation=15, ha='right')
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax2.grid(axis='y', alpha=0.3)

    # Add value labels
    ax2.bar_label(bars3, fmt='%+.1f', padding=3, fontsize=9)

    plt.tight_layout()
    return fig


def create_summary_table(recovery_stats):
    """Create summary table of all strategies' performance."""

    recovery_df = pd.DataFrame(recovery_stats)

    # Overall summary
    summary = recovery_df.groupby(["strategy", "cats"]).agg({
        "recovered": ["mean", "std", "min", "max"],
        "fraction": "mean",
        "best_score": ["mean", "std"],
        "top_100_mean": ["mean", "std"]
    }).round(2)

    # Flatten columns
    summary.columns = ['_'.join(col).strip() for col in summary.columns.values]
    summary = summary.reset_index()

    # Format for display
    display_summary = pd.DataFrame({
        "Strategy": summary["strategy"],
        "CATS": summary["cats"],
        "Recovered (mean±std)": summary.apply(
            lambda x: f"{x['recovered_mean']:.1f}±{x['recovered_std']:.1f}", axis=1
        ),
        "Recovered (min-max)": summary.apply(
            lambda x: f"{int(x['recovered_min'])}-{int(x['recovered_max'])}", axis=1
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
    """Run the full strategy recovery comparison test."""

    print("="*80)
    print("SELECTION STRATEGY RECOVERY COMPARISON")
    print("="*80)
    print("\nComparing recovery of top 100 most similar compounds")
    print("across different selection strategies with and without CATS.")

    # Run all experiments
    results_df, recovery_stats, top_100_ref, conditions = run_all_experiments()

    # Create summary table
    print("\n" + "="*80)
    print("OVERALL SUMMARY")
    print("="*80)
    summary = create_summary_table(recovery_stats)
    print(summary.to_string(index=False))

    # Save summary
    summary.to_csv("strategy_comparison_summary.csv", index=False)
    print("\n✅ Summary saved to: strategy_comparison_summary.csv")

    # Save detailed results
    results_df.to_csv("strategy_comparison_detailed.csv", index=False)
    print("✅ Detailed results saved to: strategy_comparison_detailed.csv")

    # Create plots for each strategy
    print("\n📊 Generating comparison plots...")
    strategies = sorted(results_df["strategy"].unique())

    for strategy in strategies:
        fig = create_strategy_comparison_figure(results_df, recovery_stats,
                                               top_100_ref, strategy)
        filename = f"strategy_comparison_{strategy.lower()}.png"
        fig.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"✅ {strategy} comparison saved to: {filename}")

    # Create overall summary plot
    fig_summary = create_overall_summary_plot(recovery_stats)
    fig_summary.savefig("strategy_comparison_overall.png", dpi=300, bbox_inches='tight')
    print("✅ Overall comparison saved to: strategy_comparison_overall.png")

    # Analyze CATS benefit
    print("\n" + "="*80)
    print("CATS PERFORMANCE ANALYSIS")
    print("="*80)

    recovery_df = pd.DataFrame(recovery_stats)

    for strategy in strategies:
        print(f"\n{strategy}:")

        no_cats = recovery_df[(recovery_df["strategy"] == strategy) &
                              (recovery_df["cats"] == "NoCATS")]["recovered"]
        with_cats = recovery_df[(recovery_df["strategy"] == strategy) &
                                (recovery_df["cats"] == "CATS")]["recovered"]

        print(f"  No CATS:    {no_cats.mean():.1f} ± {no_cats.std():.1f}")
        print(f"  With CATS:  {with_cats.mean():.1f} ± {with_cats.std():.1f}")
        improvement = with_cats.mean() - no_cats.mean()
        improvement_pct = (improvement / no_cats.mean()) * 100 if no_cats.mean() > 0 else 0
        print(f"  Improvement: {improvement:+.1f} compounds ({improvement_pct:+.1f}%)")

    print("\n" + "="*80)
    print("✅ STRATEGY RECOVERY COMPARISON COMPLETE!")
    print("="*80)

    plt.show()


if __name__ == "__main__":
    main()
