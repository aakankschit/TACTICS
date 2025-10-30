"""
Code blocks for Jupyter notebook: Testing Different Selection Strategies on Thrombin Dataset

Copy these code blocks into your Jupyter notebook cells.
"""

# ============================================================================
# CELL 1: Setup and Imports
# ============================================================================

import sys
import json
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Add TACTICS to path
sys.path.insert(0, "../src")

from TACTICS.thompson_sampling import (
    ThompsonSampler,
    GreedySelection,
    RouletteWheelSelection,
    UCBSelection,
    EpsilonGreedySelection,
    LookupEvaluator
)

print("✅ Imports successful!")

# ============================================================================
# CELL 2: Configuration
# ============================================================================

# Dataset paths
REAGENT_FILES = [
    "input_data/acids.smi",
    "input_data/coupled_aa_sub.smi"
]

REACTION_SMARTS = "[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]>>[#6:1](=[O:2])[#7:3]"

# Reference scores (adjust to your path)
EVALUATOR_CONFIG = {
    "ref_filename": "/Users/aakankschitnandkeolyar/Desktop/TS_Chem_Space/Thrombin/Linear_amide/docking_scores/product_scores.csv"
}

# Experiment parameters
NUM_WARMUP_TRIALS = 3
NUM_TS_ITERATIONS = 1000
NUM_TO_SELECT = None  # None = use all reagents

print("✅ Configuration set!")
print(f"   - Warmup trials: {NUM_WARMUP_TRIALS}")
print(f"   - Search iterations: {NUM_TS_ITERATIONS}")

# ============================================================================
# CELL 3: Helper Function to Run a Strategy
# ============================================================================

def run_ts_strategy(strategy, strategy_name, batch_size=1, max_resamples=None):
    """
    Run Thompson Sampling with a given strategy.

    Parameters:
        strategy: Selection strategy instance
        strategy_name: Name for display
        batch_size: Number of compounds per cycle (1=single, >1=batch)
        max_resamples: Stop after N consecutive duplicates (batch mode only)

    Returns:
        pd.DataFrame: Results with scores, SMILES, and names
    """
    print(f"\n{'='*80}")
    print(f"Running: {strategy_name}")
    print(f"  Batch size: {batch_size}")
    if max_resamples:
        print(f"  Max resamples: {max_resamples}")
    print(f"{'='*80}")

    # Create sampler
    sampler = ThompsonSampler(
        selection_strategy=strategy,
        batch_size=batch_size,
        max_resamples=max_resamples
    )

    # Setup
    evaluator = LookupEvaluator(json.dumps(EVALUATOR_CONFIG))
    sampler.set_evaluator(evaluator)
    sampler.read_reagents(REAGENT_FILES, num_to_select=NUM_TO_SELECT)
    sampler.set_reaction(REACTION_SMARTS)

    print(f"Total possible products: {sampler.get_num_prods():.2e}")

    # Warmup
    print(f"\n⏳ Running warmup...")
    sampler.warm_up(num_warmup_trials=NUM_WARMUP_TRIALS)

    # Search
    print(f"⏳ Running search...")
    import time
    start = time.time()
    results = sampler.search(num_cycles=NUM_TS_ITERATIONS)
    elapsed = time.time() - start

    # Convert to DataFrame
    df = pd.DataFrame(results, columns=["score", "SMILES", "Name"])

    # Report
    print(f"\n✅ Complete!")
    print(f"   Compounds evaluated: {len(df)}")
    print(f"   Unique compounds: {df['SMILES'].nunique()}")
    print(f"   Best score: {df['score'].min():.3f}")
    print(f"   Top 10 mean: {df.nsmallest(10, 'score')['score'].mean():.3f}")
    print(f"   Time: {elapsed:.2f}s")

    return df

print("✅ Helper function defined!")

# ============================================================================
# CELL 4: Run Greedy Strategy
# ============================================================================

greedy_strategy = GreedySelection(mode="minimize")
greedy_results = run_ts_strategy(
    strategy=greedy_strategy,
    strategy_name="Greedy Selection",
    batch_size=1
)

# Show top 10
print("\n📊 Top 10 compounds:")
print(greedy_results.nsmallest(10, "score")[["score", "Name"]])

# ============================================================================
# CELL 5: Run Roulette Wheel Strategy (Single Mode)
# ============================================================================

rw_strategy = RouletteWheelSelection(
    mode="minimize",
    alpha=0.1,
    beta=0.1,
    alpha_increment=0.01,
    beta_increment=0.001,
    efficiency_threshold=0.1
)

rw_results = run_ts_strategy(
    strategy=rw_strategy,
    strategy_name="Roulette Wheel (Single)",
    batch_size=1
)

# Show top 10
print("\n📊 Top 10 compounds:")
print(rw_results.nsmallest(10, "score")[["score", "Name"]])

# ============================================================================
# CELL 6: Run Roulette Wheel Strategy (Batch Mode)
# ============================================================================

rw_batch_strategy = RouletteWheelSelection(
    mode="minimize",
    alpha=0.1,
    beta=0.1,
    alpha_increment=0.01,
    beta_increment=0.001,
    efficiency_threshold=0.1
)

rw_batch_results = run_ts_strategy(
    strategy=rw_batch_strategy,
    strategy_name="Roulette Wheel (Batch)",
    batch_size=100,
    max_resamples=1000
)

# Show top 10
print("\n📊 Top 10 compounds:")
print(rw_batch_results.nsmallest(10, "score")[["score", "Name"]])

# Check temperature adaptation
print(f"\n🌡️  Temperature adaptation:")
print(f"   Final alpha: {rw_batch_strategy.alpha:.3f}")
print(f"   Final beta: {rw_batch_strategy.beta:.3f}")

# ============================================================================
# CELL 7: Run UCB Strategy
# ============================================================================

ucb_strategy = UCBSelection(mode="minimize", c=2.0)
ucb_results = run_ts_strategy(
    strategy=ucb_strategy,
    strategy_name="UCB (c=2.0)",
    batch_size=1
)

# Show top 10
print("\n📊 Top 10 compounds:")
print(ucb_results.nsmallest(10, "score")[["score", "Name"]])

# ============================================================================
# CELL 8: Run Epsilon-Greedy Strategy
# ============================================================================

eps_strategy = EpsilonGreedySelection(
    mode="minimize",
    epsilon=0.1,
    decay=0.995
)

eps_results = run_ts_strategy(
    strategy=eps_strategy,
    strategy_name="Epsilon-Greedy (ε=0.1)",
    batch_size=1
)

# Show top 10
print("\n📊 Top 10 compounds:")
print(eps_results.nsmallest(10, "score")[["score", "Name"]])

# ============================================================================
# CELL 9: Compare All Strategies
# ============================================================================

# Combine all results
all_results = {
    "Greedy": greedy_results,
    "Roulette_Single": rw_results,
    "Roulette_Batch": rw_batch_results,
    "UCB": ucb_results,
    "Epsilon_Greedy": eps_results
}

# Create comparison table
comparison_data = []
for name, df in all_results.items():
    comparison_data.append({
        "Strategy": name,
        "N Compounds": len(df),
        "N Unique": df["SMILES"].nunique(),
        "Best Score": df["score"].min(),
        "Top 10 Mean": df.nsmallest(10, "score")["score"].mean(),
        "Top 50 Mean": df.nsmallest(50, "score")["score"].mean(),
        "Median": df["score"].median()
    })

comparison_df = pd.DataFrame(comparison_data)
print("\n" + "="*80)
print("📊 STRATEGY COMPARISON")
print("="*80)
print(comparison_df.to_string(index=False))

# Find best
best_idx = comparison_df["Best Score"].idxmin()
print(f"\n🏆 Best Strategy: {comparison_df.loc[best_idx, 'Strategy']}")
print(f"   Best Score: {comparison_df.loc[best_idx, 'Best Score']:.3f}")

# ============================================================================
# CELL 10: Visualization - Cumulative Best Score
# ============================================================================

plt.figure(figsize=(12, 6))

colors = {
    "Greedy": "#2ecc71",
    "Roulette_Single": "#e74c3c",
    "Roulette_Batch": "#e67e22",
    "UCB": "#3498db",
    "Epsilon_Greedy": "#9b59b6"
}

for name, df in all_results.items():
    cumulative_best = df["score"].cummin()
    plt.plot(range(len(cumulative_best)), cumulative_best,
             label=name, color=colors[name], linewidth=2, alpha=0.8)

plt.xlabel("Iteration", fontsize=12)
plt.ylabel("Best Score Found (lower is better)", fontsize=12)
plt.title("Cumulative Best Score Comparison - Thrombin Dataset", fontsize=14, fontweight='bold')
plt.legend(loc='best', frameon=True, shadow=True)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# ============================================================================
# CELL 11: Visualization - Distribution Comparison
# ============================================================================

fig, axes = plt.subplots(1, 2, figsize=(15, 5))

# Box plot
ax = axes[0]
data_for_box = [df["score"].values for df in all_results.values()]
bp = ax.boxplot(data_for_box, labels=all_results.keys(), patch_artist=True)
for patch, color in zip(bp['boxes'], colors.values()):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax.set_ylabel("Score", fontsize=12)
ax.set_title("Score Distributions (Box Plot)", fontsize=14, fontweight='bold')
ax.tick_params(axis='x', rotation=45)
ax.grid(True, alpha=0.3, axis='y')

# Bar chart - Top 10 performance
ax = axes[1]
strategies = list(all_results.keys())
top10_means = [df.nsmallest(10, "score")["score"].mean() for df in all_results.values()]
bars = ax.bar(strategies, top10_means, color=list(colors.values()), alpha=0.7, edgecolor='black')
ax.set_ylabel("Mean Score of Top 10", fontsize=12)
ax.set_title("Top 10 Performance Comparison", fontsize=14, fontweight='bold')
ax.tick_params(axis='x', rotation=45)
for bar, score in zip(bars, top10_means):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
            f'{score:.3f}', ha='center', va='bottom', fontsize=10)

plt.tight_layout()
plt.show()

# ============================================================================
# CELL 12: Save Results
# ============================================================================

# Save comparison table
comparison_df.to_csv("thrombin_comparison.csv", index=False)
print("✅ Saved: thrombin_comparison.csv")

# Save individual results
for name, df in all_results.items():
    filename = f"thrombin_results_{name}.csv"
    df.to_csv(filename, index=False)
    print(f"✅ Saved: {filename}")

print("\n🎉 All results saved!")

# ============================================================================
# CELL 13: Analyze Top Hits
# ============================================================================

# Get top 20 from each strategy
print("\n" + "="*80)
print("🔬 ANALYZING TOP 20 HITS FROM EACH STRATEGY")
print("="*80)

top_hits = {}
for name, df in all_results.items():
    top20 = df.nsmallest(20, "score")
    top_hits[name] = set(top20["SMILES"].values)
    print(f"\n{name}:")
    print(f"  Best score: {top20['score'].min():.3f}")
    print(f"  Top 20 mean: {top20['score'].mean():.3f}")

# Find common hits
print("\n" + "="*80)
print("🎯 OVERLAP ANALYSIS")
print("="*80)

all_top_smiles = set.union(*top_hits.values())
print(f"\nTotal unique compounds in top 20s: {len(all_top_smiles)}")

# Count how many strategies found each compound
from collections import Counter
smiles_counts = Counter()
for smiles in all_top_smiles:
    for name, hits in top_hits.items():
        if smiles in hits:
            smiles_counts[smiles] += 1

# Compounds found by all strategies
universal_hits = [s for s, c in smiles_counts.items() if c == len(all_results)]
print(f"Compounds found by ALL strategies: {len(universal_hits)}")

# Compounds found by only one strategy
unique_hits = {name: [s for s in hits if smiles_counts[s] == 1]
               for name, hits in top_hits.items()}
print(f"\nUnique hits per strategy:")
for name, hits in unique_hits.items():
    print(f"  {name}: {len(hits)} unique compounds")

print("\n✅ Analysis complete!")
