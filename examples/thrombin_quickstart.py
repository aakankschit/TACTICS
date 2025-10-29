"""
Quick Start: Testing Selection Strategies on Thrombin Dataset

Minimal example to get started quickly.
"""

import sys
sys.path.insert(0, "../src")

from TACTICS.thompson_sampling import (
    ThompsonSampler,
    RouletteWheelSelection,
    LookupEvaluator
)
import json
import pandas as pd

# ============================================================================
# Quick Configuration
# ============================================================================

# Paths
REAGENT_FILES = [
    "input_data/acids.smi",
    "input_data/coupled_aa_sub.smi"
]

REACTION_SMARTS = "[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]>>[#6:1](=[O:2])[#7:3]"

EVALUATOR_CONFIG = {
    "ref_filename": "/Users/aakankschitnandkeolyar/Desktop/TS_Chem_Space/Thrombin/Linear_amide/docking_scores/product_scores.csv"
}

# ============================================================================
# Run Thompson Sampling
# ============================================================================

print("="*80)
print("THROMBIN DATASET - QUICK START")
print("="*80)

# 1. Create selection strategy
strategy = RouletteWheelSelection(
    mode="minimize",
    alpha=0.1,
    beta=0.1
)

# 2. Create sampler
sampler = ThompsonSampler(
    selection_strategy=strategy,
    batch_size=100,          # Use batch mode for speed
    max_resamples=1000
)

# 3. Setup
evaluator = LookupEvaluator(json.dumps(EVALUATOR_CONFIG))
sampler.set_evaluator(evaluator)
sampler.read_reagents(REAGENT_FILES)
sampler.set_reaction(REACTION_SMARTS)

print(f"\nTotal possible products: {sampler.get_num_prods():.2e}")

# 4. Warmup
print("\n⏳ Running warmup...")
sampler.warm_up(num_warmup_trials=3)

# 5. Search
print("⏳ Running search (1000 iterations)...")
results = sampler.search(num_cycles=1000)

# 6. Analyze results
df = pd.DataFrame(results, columns=["score", "SMILES", "Name"])

print("\n" + "="*80)
print("RESULTS")
print("="*80)
print(f"Compounds evaluated: {len(df)}")
print(f"Unique compounds: {df['SMILES'].nunique()}")
print(f"Best score: {df['score'].min():.3f}")
print(f"Top 10 mean: {df.nsmallest(10, 'score')['score'].mean():.3f}")

print("\n📊 Top 10 compounds:")
print(df.nsmallest(10, "score")[["score", "Name"]])

# 7. Save results
df.to_csv("thrombin_quickstart_results.csv", index=False)
print("\n✅ Results saved to: thrombin_quickstart_results.csv")

# ============================================================================
# Try Different Strategy (Optional)
# ============================================================================

print("\n" + "="*80)
print("TRYING DIFFERENT STRATEGY: UCB")
print("="*80)

from TACTICS.thompson_sampling import UCBSelection

# Just swap the strategy!
ucb_strategy = UCBSelection(mode="minimize", c=2.0)
ucb_sampler = ThompsonSampler(selection_strategy=ucb_strategy, batch_size=1)

ucb_sampler.set_evaluator(LookupEvaluator(json.dumps(EVALUATOR_CONFIG)))
ucb_sampler.read_reagents(REAGENT_FILES)
ucb_sampler.set_reaction(REACTION_SMARTS)

print("\n⏳ Running warmup...")
ucb_sampler.warm_up(num_warmup_trials=3)

print("⏳ Running search...")
ucb_results = ucb_sampler.search(num_cycles=1000)

ucb_df = pd.DataFrame(ucb_results, columns=["score", "SMILES", "Name"])
print(f"\nBest score with UCB: {ucb_df['score'].min():.3f}")
print(f"Top 10 mean: {ucb_df.nsmallest(10, 'score')['score'].mean():.3f}")

print("\n🎉 Complete! Try other strategies by importing them:")
print("   - GreedySelection")
print("   - RouletteWheelSelection")
print("   - UCBSelection")
print("   - EpsilonGreedySelection")
