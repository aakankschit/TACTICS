"""
RWS and BayesUCB CATS Comparison - Clean Implementation

Compares:
1. Current RWS (with CATS)
2. Legacy RWS
3. Current BayesUCB (with CATS)

Key improvements:
- Combines warmup + search results for accurate recovery
- Simplified recovery calculation
- Clear error handling
"""

import sys
import tempfile
import os
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
    RouletteWheelSelection,
    BayesUCBSelection,
    LookupEvaluator
)
from TACTICS.thompson_sampling.warmup import StandardWarmup, EnhancedWarmup

# Import legacy RWS
from TACTICS.thompson_sampling.legacy.rws_sampling import RWSSampler as LegacyRWSSampler
from TACTICS.thompson_sampling.legacy.evaluators import LookupEvaluator as LegacyLookupEvaluator

# ============================================================================
# Configuration
# ============================================================================

# 3-component Quinazoline library
REAGENT_FILES = [
    "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/reagents/thrombin/quinazoline/niementowski_reagent_0.smi",
    "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/reagents/thrombin/quinazoline/niementowski_reagent_1.smi",
    "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/reagents/thrombin/quinazoline/niementowski_reagent_2.smi"
]

# Evaluator configs
CURRENT_EVALUATOR_CONFIG = {
    "ref_filename": "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/reagents/thrombin/quinazoline/niementowski.parquet",
    "compound_col": "Name",
    "score_col": "query_001"
}

LEGACY_EVALUATOR_CONFIG = {
    "ref_filename": "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/reagents/thrombin/quinazoline/niementowski.parquet",
    "compound_col": "Name",
    "ref_colname": "query_001"
}

REFERENCE_FILE = "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/reagents/thrombin/quinazoline/niementowski.parquet"
# Experiment parameters
NUM_WARMUP_TRIALS = 5
NUM_TS_ITERATIONS = 18500
BATCH_SIZE = 100
NUM_REPLICATES = 3

# Strategy configurations
RWS_CONFIG = {
    "alpha": 0.1,
    "beta": 0.05,
    "exploration_phase_end": 0.20,
    "transition_phase_end": 0.60,
    "min_observations": 5
}

BAYES_UCB_CONFIG = {
    "initial_p_high": 0.90,
    "initial_p_low": 0.60,
    "exploration_phase_end": 0.20,
    "transition_phase_end": 0.60,
    "min_observations": 5
}

LEGACY_RWS_CONFIG = {
    "alpha": 0.1,
    "beta": 0.05,
    "alpha_max": 0.4,
    "alpha_increment": 0.01,
    "beta_increment": 0.001,
    "efficiency_threshold": 0.1
}

# ============================================================================
# Helper Functions
# ============================================================================

def load_reference_data():
    """Load ROCS reference data and get top 100 most similar."""
    print("Loading reference (ROCS) data...")

    # Read parquet file
    ref_df = pd.read_parquet(REFERENCE_FILE)

    # Standardize column names
    score_col = CURRENT_EVALUATOR_CONFIG["score_col"]
    compound_col = CURRENT_EVALUATOR_CONFIG["compound_col"]

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


def standardize_dataframe(df):
    """Standardize dataframe to have 'Name' and 'score' columns."""
    # Convert to pandas if needed
    try:
        import polars as pl
        if isinstance(df, pl.DataFrame):
            df = df.to_pandas()
    except ImportError:
        pass

    # Ensure 'Name' column exists
    if 'Name' not in df.columns:
        if 'SMILES' in df.columns:
            df = df.rename(columns={'SMILES': 'Name'})
        elif 'Product_Code' in df.columns:
            df = df.rename(columns={'Product_Code': 'Name'})

    # Ensure 'score' column exists (lowercase)
    if 'score' not in df.columns and 'Score' in df.columns:
        df = df.rename(columns={'Score': 'score'})
    if 'score' not in df.columns and 'Scores' in df.columns:
        df = df.rename(columns={'Scores': 'score'})

    return df


def calculate_recovery_metrics(results_df, top_100_ref):
    """Calculate recovery metrics for a single run."""
    # Standardize both dataframes
    results_df = standardize_dataframe(results_df)
    top_100_ref = standardize_dataframe(top_100_ref)

    # Verify required columns exist
    if 'Name' not in results_df.columns or 'score' not in results_df.columns:
        raise ValueError(f"Results dataframe missing required columns. Has: {list(results_df.columns)}")
    if 'Name' not in top_100_ref.columns or 'score' not in top_100_ref.columns:
        raise ValueError(f"Reference dataframe missing required columns. Has: {list(top_100_ref.columns)}")

    # Get sets of compound names (ensure we have strings, not lists)
    top_100_names = set(str(name) for name in top_100_ref['Name'].tolist())
    results_names = set(str(name) for name in results_df['Name'].tolist())

    # Recovery: how many of top 100 were found
    recovered = top_100_names & results_names
    recovery_rate = len(recovered) / len(top_100_names)

    # Best score found
    best_score = results_df['score'].min()

    # Top 10 recovery
    top_10_names = set(str(name) for name in top_100_ref.head(10)['Name'].tolist())
    top_10_recovered = top_10_names & results_names
    top_10_recovery = len(top_10_recovered) / len(top_10_names)

    return {
        'recovery_rate': recovery_rate,
        'recovered_count': len(recovered),
        'best_score': best_score,
        'top_10_recovery': top_10_recovery,
        'total_evaluated': len(results_df)
    }


def combine_results(warmup_results, search_results):
    """Combine warmup and search results, handling polars/pandas."""
    try:
        import polars as pl

        # If both are polars, concat in polars
        if isinstance(warmup_results, pl.DataFrame) and isinstance(search_results, pl.DataFrame):
            combined = pl.concat([warmup_results, search_results])
            combined_pd = combined.to_pandas()
        else:
            # Convert to pandas if needed
            warmup_pd = warmup_results.to_pandas() if isinstance(warmup_results, pl.DataFrame) else warmup_results
            search_pd = search_results.to_pandas() if isinstance(search_results, pl.DataFrame) else search_results
            combined_pd = pd.concat([warmup_pd, search_pd], ignore_index=True)
    except ImportError:
        # No polars, just use pandas
        warmup_pd = warmup_results if isinstance(warmup_results, pd.DataFrame) else pd.DataFrame(warmup_results)
        search_pd = search_results if isinstance(search_results, pd.DataFrame) else pd.DataFrame(search_results)
        combined_pd = pd.concat([warmup_pd, search_pd], ignore_index=True)

    # Standardize column names
    combined_pd = standardize_dataframe(combined_pd)

    # Remove duplicates (keep best score)
    combined_pd = combined_pd.sort_values('score').drop_duplicates(subset=['Name'], keep='first')

    return combined_pd


# ============================================================================
# Strategy Run Functions
# ============================================================================

def run_current_rws(replicate_num, ref_df, top_100_ref):
    """Run current RWS implementation with CATS."""
    print(f"\n  Replicate {replicate_num + 1}/{NUM_REPLICATES}")

    # Create evaluator
    evaluator = LookupEvaluator(CURRENT_EVALUATOR_CONFIG)

    # Create strategy with CATS
    strategy = RouletteWheelSelection(
        mode="maximize",
        **RWS_CONFIG
    )

    # Create sampler
    sampler = ThompsonSampler(
        selection_strategy=strategy,
        warmup_strategy=EnhancedWarmup(),
        batch_size=BATCH_SIZE,
        use_boltzmann_weighting=True
    )

    # Setup evaluator and reagents
    sampler.set_evaluator(evaluator)
    sampler.read_reagents(REAGENT_FILES)
    sampler.set_hide_progress(True)

    # Run warmup and save results
    warmup_results = sampler.warm_up(num_warmup_trials=NUM_WARMUP_TRIALS)

    # Run search
    search_results = sampler.search(num_cycles=100000, max_evaluations=NUM_TS_ITERATIONS)

    # Combine results
    results_df = combine_results(warmup_results, search_results)

    # Calculate metrics
    metrics = calculate_recovery_metrics(results_df, top_100_ref)

    print(f"    Recovery: {metrics['recovery_rate']:.1%} ({metrics['recovered_count']}/100)")
    print(f"    Top 10 recovery: {metrics['top_10_recovery']:.1%}")
    print(f"    Best score: {metrics['best_score']:.3f}")
    print(f"    Total evaluated: {metrics['total_evaluated']:,}")

    return metrics, results_df


def run_legacy_rws(replicate_num, ref_df, top_100_ref):
    """Run legacy RWS implementation."""
    print(f"\n  Replicate {replicate_num + 1}/{NUM_REPLICATES}")

    # Create legacy sampler
    sampler = LegacyRWSSampler(processes=1, scaling=1)

    # Create and set evaluator
    evaluator = LegacyLookupEvaluator(LEGACY_EVALUATOR_CONFIG)
    sampler.set_evaluator(evaluator)

    # Read reagents
    sampler.read_reagents(REAGENT_FILES)
    sampler.set_hide_progress(True)

    # Set legacy RWS parameters
    sampler.alpha = LEGACY_RWS_CONFIG['alpha']
    sampler.beta = LEGACY_RWS_CONFIG['beta']

    # Create temp file for results
    temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv')
    temp_filename = temp_file.name
    temp_file.close()

    try:
        # Run warmup
        num_warmup_evals = sampler.warm_up(
            num_warmup_trials=NUM_WARMUP_TRIALS,
            results_filename=temp_filename
        )

        # Run search
        total_target = num_warmup_evals + NUM_TS_ITERATIONS
        percent = total_target / sampler.get_num_prods()

        sampler.search(
            percent_of_library=percent,
            stop=6000,
            min_cpds_per_core=50,
            results_filename=temp_filename
        )

        # Read results from file
        results_df = pd.read_csv(temp_filename)

    finally:
        # Clean up temp file
        if os.path.exists(temp_filename):
            os.remove(temp_filename)

    # Standardize dataframe
    results_df = standardize_dataframe(results_df)

    # Calculate metrics
    metrics = calculate_recovery_metrics(results_df, top_100_ref)

    print(f"    Recovery: {metrics['recovery_rate']:.1%} ({metrics['recovered_count']}/100)")
    print(f"    Top 10 recovery: {metrics['top_10_recovery']:.1%}")
    print(f"    Best score: {metrics['best_score']:.3f}")
    print(f"    Total evaluated: {metrics['total_evaluated']:,}")

    return metrics, results_df


def run_current_bayes_ucb(replicate_num, ref_df, top_100_ref):
    """Run current BayesUCB implementation with CATS."""
    print(f"\n  Replicate {replicate_num + 1}/{NUM_REPLICATES}")

    # Create evaluator
    evaluator = LookupEvaluator(CURRENT_EVALUATOR_CONFIG)

    # Create strategy with CATS
    strategy = BayesUCBSelection(
        mode="maximize",
        **BAYES_UCB_CONFIG
    )

    # Create sampler
    sampler = ThompsonSampler(
        selection_strategy=strategy,
        warmup_strategy=StandardWarmup(),
        batch_size=BATCH_SIZE,
        use_boltzmann_weighting=True
    )

    # Setup evaluator and reagents
    sampler.set_evaluator(evaluator)
    sampler.read_reagents(REAGENT_FILES)
    sampler.set_hide_progress(True)

    # Run warmup and save results
    warmup_results = sampler.warm_up(num_warmup_trials=NUM_WARMUP_TRIALS)

    # Run search
    search_results = sampler.search(num_cycles=100000, max_evaluations=NUM_TS_ITERATIONS)

    # Combine results
    results_df = combine_results(warmup_results, search_results)

    # Calculate metrics
    metrics = calculate_recovery_metrics(results_df, top_100_ref)

    print(f"    Recovery: {metrics['recovery_rate']:.1%} ({metrics['recovered_count']}/100)")
    print(f"    Top 10 recovery: {metrics['top_10_recovery']:.1%}")
    print(f"    Best score: {metrics['best_score']:.3f}")
    print(f"    Total evaluated: {metrics['total_evaluated']:,}")

    return metrics, results_df


# ============================================================================
# Main Comparison
# ============================================================================

def main():
    """Run comparison experiments."""
    print("="*80)
    print("RWS and BayesUCB CATS Implementation Comparison")
    print("="*80)

    # Load reference data
    ref_df, top_100_ref = load_reference_data()

    # Storage for results
    all_metrics = {
        'Current RWS (CATS)': [],
        'Legacy RWS': [],
        'Current BayesUCB (CATS)': []
    }

    # Run experiments
    for replicate in range(NUM_REPLICATES):
        print(f"\n{'='*80}")
        print(f"Replicate {replicate + 1}/{NUM_REPLICATES}")
        print(f"{'='*80}")

        # Current RWS with CATS
        print("\n[1/3] Running Current RWS (with CATS)...")
        try:
            metrics_rws, _ = run_current_rws(replicate, ref_df, top_100_ref)
            all_metrics['Current RWS (CATS)'].append(metrics_rws)
        except Exception as e:
            print(f"    ERROR: {e}")
            import traceback
            traceback.print_exc()

        # Legacy RWS
        print("\n[2/3] Running Legacy RWS...")
        try:
            metrics_legacy, _ = run_legacy_rws(replicate, ref_df, top_100_ref)
            all_metrics['Legacy RWS'].append(metrics_legacy)
        except Exception as e:
            print(f"    ERROR: {e}")
            import traceback
            traceback.print_exc()

        # Current BayesUCB with CATS
        print("\n[3/3] Running Current BayesUCB (with CATS)...")
        try:
            metrics_bayes, _ = run_current_bayes_ucb(replicate, ref_df, top_100_ref)
            all_metrics['Current BayesUCB (CATS)'].append(metrics_bayes)
        except Exception as e:
            print(f"    ERROR: {e}")
            import traceback
            traceback.print_exc()

    # ========================================================================
    # Summary Statistics
    # ========================================================================
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)

    summary_data = []
    for strategy_name, metrics_list in all_metrics.items():
        if not metrics_list:
            print(f"\nWARNING: No results for {strategy_name}")
            continue

        recovery_rates = [m['recovery_rate'] for m in metrics_list]
        top_10_rates = [m['top_10_recovery'] for m in metrics_list]
        best_scores = [m['best_score'] for m in metrics_list]

        summary_data.append({
            'Strategy': strategy_name,
            'Recovery Mean': np.mean(recovery_rates),
            'Recovery Std': np.std(recovery_rates),
            'Top 10 Mean': np.mean(top_10_rates),
            'Top 10 Std': np.std(top_10_rates),
            'Best Score Mean': np.mean(best_scores),
            'Best Score Std': np.std(best_scores)
        })

    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        print("\n" + summary_df.to_string(index=False))

        # Save detailed results
        detailed_results = []
        for strategy_name, metrics_list in all_metrics.items():
            for i, metrics in enumerate(metrics_list):
                detailed_results.append({
                    'Strategy': strategy_name,
                    'Replicate': i + 1,
                    **metrics
                })

        if detailed_results:
            detailed_df = pd.DataFrame(detailed_results)
            detailed_df.to_csv('rws_bayesucb_cats_detailed_fixed.csv', index=False)
            print(f"\nDetailed results saved to: rws_bayesucb_cats_detailed_fixed.csv")

            summary_df.to_csv('rws_bayesucb_cats_summary_fixed.csv', index=False)
            print(f"Summary results saved to: rws_bayesucb_cats_summary_fixed.csv")

            # ================================================================
            # Visualization
            # ================================================================
            print("\nGenerating comparison plots...")

            fig, axes = plt.subplots(1, 3, figsize=(15, 5))

            # Recovery rate comparison
            recovery_data = detailed_df[['Strategy', 'recovery_rate']].copy()
            recovery_data['recovery_rate'] *= 100

            sns.barplot(data=recovery_data, x='Strategy', y='recovery_rate', ax=axes[0], errorbar='sd')
            axes[0].set_title('Top 100 Recovery Rate', fontsize=12, fontweight='bold')
            axes[0].set_ylabel('Recovery Rate (%)')
            axes[0].set_xlabel('')
            axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=15, ha='right')
            axes[0].grid(axis='y', alpha=0.3)

            # Top 10 recovery comparison
            top10_data = detailed_df[['Strategy', 'top_10_recovery']].copy()
            top10_data['top_10_recovery'] *= 100

            sns.barplot(data=top10_data, x='Strategy', y='top_10_recovery', ax=axes[1], errorbar='sd')
            axes[1].set_title('Top 10 Recovery Rate', fontsize=12, fontweight='bold')
            axes[1].set_ylabel('Recovery Rate (%)')
            axes[1].set_xlabel('')
            axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=15, ha='right')
            axes[1].grid(axis='y', alpha=0.3)

            # Best score comparison
            sns.barplot(data=detailed_df, x='Strategy', y='best_score', ax=axes[2], errorbar='sd')
            axes[2].set_title('Best Score Found', fontsize=12, fontweight='bold')
            axes[2].set_ylabel('ROCS Score')
            axes[2].set_xlabel('')
            axes[2].set_xticklabels(axes[2].get_xticklabels(), rotation=15, ha='right')
            axes[2].grid(axis='y', alpha=0.3)

            plt.tight_layout()
            plt.savefig('rws_bayesucb_cats_comparison_rocs.png', dpi=300, bbox_inches='tight')
            print(f"Comparison plot saved to: rws_bayesucb_cats_comparison_fixed.png")
    else:
        print("\nERROR: No results to summarize!")

    print("\n" + "="*80)
    print("COMPARISON COMPLETE")
    print("="*80)


if __name__ == "__main__":
    main()
