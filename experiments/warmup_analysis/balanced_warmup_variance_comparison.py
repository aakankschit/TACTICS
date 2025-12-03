"""
Balanced Warmup Variance Comparison

Compares per-reagent variance vs global variance using the new BalancedWarmup strategy:

1. Current RWS + BalancedWarmup (per-reagent variance with shrinkage) - DEFAULT
2. Current RWS + BalancedWarmup (global variance)
3. Current RWS + EnhancedWarmup (legacy warmup, global variance) - BASELINE
4. Legacy RWS - REFERENCE

Key tests:
- Does per-reagent variance improve recovery?
- Does BalancedWarmup improve reproducibility (lower variance across runs)?
- How does guaranteed K observations per reagent affect performance?

Usage:
    python balanced_warmup_variance_comparison.py --output_dir /path/to/output
    python balanced_warmup_variance_comparison.py -o ./results/run_001
"""

import sys
import tempfile
import os
import argparse
from pathlib import Path
from datetime import datetime
from itertools import combinations
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from scipy import stats

# Import current implementation
from TACTICS.thompson_sampling import (
    ThompsonSampler,
    RouletteWheelSelection,
    LookupEvaluator
)
from TACTICS.thompson_sampling.warmup import BalancedWarmup, EnhancedWarmup

# Import legacy RWS
from TACTICS.thompson_sampling.legacy.rws_sampling import RWSSampler as LegacyRWSSampler
from TACTICS.thompson_sampling.legacy.evaluators import LookupEvaluator as LegacyLookupEvaluator

# ============================================================================
# Configuration
# ============================================================================

# 2-component Thrombin library
REAGENT_FILES = [
    "/Users/aakankschitnandkeolyar/Desktop/TACTICS/examples/input_files/acids.smi",
    "/Users/aakankschitnandkeolyar/Desktop/TACTICS/examples/input_files/coupled_aa_sub.smi"
]

# Evaluator configs
CURRENT_EVALUATOR_CONFIG = {
    "ref_filename": "/Users/aakankschitnandkeolyar/Desktop/TACTICS/examples/docking_scores/product_scores.csv",
    "compound_col": "Product_Code",
    "score_col": "Scores"
}

LEGACY_EVALUATOR_CONFIG = {
    "ref_filename": "/Users/aakankschitnandkeolyar/Desktop/TACTICS/examples/docking_scores/product_scores.csv",
    "ref_colname": "Scores",
    "compound_col": "Product_Code"
}

REFERENCE_FILE = "/Users/aakankschitnandkeolyar/Desktop/TACTICS/examples/docking_scores/product_scores.csv"

# Experiment parameters
OBSERVATIONS_PER_REAGENT = 5  # K parameter for BalancedWarmup
NUM_WARMUP_TRIALS = 5  # For legacy/enhanced warmup
NUM_TS_ITERATIONS = 5000
BATCH_SIZE = 100
NUM_REPLICATES = 5  # More replicates to measure variance reduction

# Strategy configurations
RWS_CONFIG = {
    "alpha": 0.1,
    "beta": 0.05,
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
    """Load brute force reference data and get top 100."""
    print("Loading reference (brute force) data...")

    ref_df = pd.read_csv(REFERENCE_FILE)

    # Standardize column names
    if "Product_Code" in ref_df.columns:
        ref_df = ref_df.rename(columns={"Product_Code": "Name", "Scores": "score"})

    # Sort and get top 100 (minimize scores)
    ref_df = ref_df.sort_values("score")
    top_100_ref = ref_df.head(100).copy()

    print(f"  Total products: {len(ref_df):,}")
    print(f"  Best score: {ref_df['score'].min():.3f}")
    print(f"  Top 100 range: {top_100_ref['score'].min():.3f} to {top_100_ref['score'].max():.3f}")

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

    # Get sets of compound names
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


def analyze_reagent_variances(sampler):
    """Analyze variance distribution across reagents after warmup."""
    variances = []
    for component_idx, reagent_list in enumerate(sampler.reagent_lists):
        for reagent in reagent_list:
            if reagent.known_var is not None:
                variances.append({
                    'component': component_idx,
                    'reagent': reagent.reagent_name,
                    'variance': reagent.known_var,
                    'std': reagent.std,
                    'mean': reagent.mean,
                    'n_samples': reagent.n_samples
                })
    return pd.DataFrame(variances)


def run_statistical_tests(detailed_df: pd.DataFrame, metric: str = 'recovery_rate') -> dict:
    """
    Run statistical tests comparing recovery rates across strategies.

    Performs:
    1. One-way ANOVA (parametric) or Kruskal-Wallis (non-parametric)
    2. Pairwise comparisons with Bonferroni correction

    Args:
        detailed_df: DataFrame with 'Strategy' and metric columns
        metric: Column name to analyze

    Returns:
        Dictionary with test results
    """
    strategies = detailed_df['Strategy'].unique()
    groups = [detailed_df[detailed_df['Strategy'] == s][metric].values for s in strategies]

    results = {
        'strategies': list(strategies),
        'n_per_group': [len(g) for g in groups],
        'means': [np.mean(g) for g in groups],
        'stds': [np.std(g) for g in groups],
    }

    # Check normality for each group (Shapiro-Wilk test)
    normality_results = []
    for i, (s, g) in enumerate(zip(strategies, groups)):
        if len(g) >= 3:  # Need at least 3 samples
            stat, p = stats.shapiro(g)
            normality_results.append({'strategy': s, 'shapiro_stat': stat, 'shapiro_p': p, 'is_normal': p > 0.05})
        else:
            normality_results.append({'strategy': s, 'shapiro_stat': np.nan, 'shapiro_p': np.nan, 'is_normal': None})

    results['normality'] = normality_results
    all_normal = all(r['is_normal'] for r in normality_results if r['is_normal'] is not None)

    # Run omnibus test
    if all_normal and len(groups) > 2:
        # One-way ANOVA
        f_stat, anova_p = stats.f_oneway(*groups)
        results['omnibus_test'] = 'One-way ANOVA'
        results['omnibus_statistic'] = f_stat
        results['omnibus_p'] = anova_p
    else:
        # Kruskal-Wallis (non-parametric)
        h_stat, kw_p = stats.kruskal(*groups)
        results['omnibus_test'] = 'Kruskal-Wallis'
        results['omnibus_statistic'] = h_stat
        results['omnibus_p'] = kw_p

    # Pairwise comparisons with Bonferroni correction
    pairwise_results = []
    n_comparisons = len(list(combinations(range(len(strategies)), 2)))

    for (i, s1), (j, s2) in combinations(enumerate(strategies), 2):
        g1, g2 = groups[i], groups[j]

        # Use t-test if both groups are normal, otherwise Mann-Whitney U
        if all_normal:
            stat, p = stats.ttest_ind(g1, g2)
            test_name = 't-test'
        else:
            stat, p = stats.mannwhitneyu(g1, g2, alternative='two-sided')
            test_name = 'Mann-Whitney U'

        # Bonferroni correction
        p_corrected = min(p * n_comparisons, 1.0)

        # Effect size (Cohen's d)
        pooled_std = np.sqrt(((len(g1)-1)*np.std(g1)**2 + (len(g2)-1)*np.std(g2)**2) / (len(g1)+len(g2)-2))
        cohens_d = (np.mean(g1) - np.mean(g2)) / pooled_std if pooled_std > 0 else 0

        pairwise_results.append({
            'group1': s1,
            'group2': s2,
            'test': test_name,
            'statistic': stat,
            'p_value': p,
            'p_corrected': p_corrected,
            'cohens_d': cohens_d,
            'significant': p_corrected < 0.05
        })

    results['pairwise'] = pairwise_results

    return results


def get_significance_symbol(p_value: float) -> str:
    """Convert p-value to significance symbol."""
    if p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return 'ns'


def add_significance_annotations(ax, detailed_df: pd.DataFrame, metric: str, stat_results: dict):
    """
    Add significance annotations to a boxplot.

    Args:
        ax: Matplotlib axis
        detailed_df: DataFrame with data
        metric: Metric being plotted
        stat_results: Results from run_statistical_tests
    """
    strategies = stat_results['strategies']
    pairwise = stat_results['pairwise']

    # Get positions for each strategy on x-axis
    strategy_positions = {s: i for i, s in enumerate(strategies)}

    # Get y-axis limits
    ymin, ymax = ax.get_ylim()
    y_range = ymax - ymin

    # Filter to only significant comparisons (or show top few)
    significant_pairs = [p for p in pairwise if p['significant']]

    # If no significant pairs, show the comparison between key strategies
    if not significant_pairs:
        # Show comparison between Balanced + Per-Reagent Var and Balanced + Global Var if they exist
        key_comparisons = [p for p in pairwise if
            ('Per-Reagent' in p['group1'] and 'Global' in p['group2']) or
            ('Global' in p['group1'] and 'Per-Reagent' in p['group2'])]
        pairs_to_show = key_comparisons[:2] if key_comparisons else pairwise[:2]
    else:
        pairs_to_show = significant_pairs[:4]  # Limit to avoid clutter

    # Add significance bars
    y_offset = 0.05 * y_range
    for idx, pair in enumerate(pairs_to_show):
        x1 = strategy_positions.get(pair['group1'], 0)
        x2 = strategy_positions.get(pair['group2'], 1)

        # Calculate bar height
        bar_height = ymax + y_offset + (idx * 0.08 * y_range)

        # Draw the bar
        ax.plot([x1, x1, x2, x2], [bar_height - 0.02*y_range, bar_height, bar_height, bar_height - 0.02*y_range],
                color='black', linewidth=1)

        # Add significance symbol
        symbol = get_significance_symbol(pair['p_corrected'])
        ax.text((x1 + x2) / 2, bar_height + 0.01*y_range, symbol,
                ha='center', va='bottom', fontsize=10, fontweight='bold')

    # Extend y-axis to accommodate annotations
    if pairs_to_show:
        new_ymax = ymax + y_offset + (len(pairs_to_show) * 0.1 * y_range)
        ax.set_ylim(ymin, new_ymax)


def print_statistical_summary(stat_results: dict, output_path: Path = None):
    """Print and optionally save statistical test results."""
    print("\n" + "="*80)
    print("STATISTICAL ANALYSIS")
    print("="*80)

    # Omnibus test
    print(f"\n{stat_results['omnibus_test']}:")
    print(f"  Statistic: {stat_results['omnibus_statistic']:.4f}")
    print(f"  p-value: {stat_results['omnibus_p']:.4e}")
    if stat_results['omnibus_p'] < 0.05:
        print("  -> Significant differences exist between groups (p < 0.05)")
    else:
        print("  -> No significant differences between groups (p >= 0.05)")

    # Pairwise comparisons
    print("\nPairwise Comparisons (Bonferroni-corrected):")
    print("-" * 80)
    print(f"{'Comparison':<50} {'p-value':>12} {'p-corrected':>12} {'Sig':>6} {'Cohen d':>10}")
    print("-" * 80)

    for pair in stat_results['pairwise']:
        comparison = f"{pair['group1']} vs {pair['group2']}"
        sig = get_significance_symbol(pair['p_corrected'])
        print(f"{comparison:<50} {pair['p_value']:>12.4e} {pair['p_corrected']:>12.4e} {sig:>6} {pair['cohens_d']:>10.3f}")

    print("-" * 80)
    print("Significance: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant")
    print(f"Effect size (Cohen's d): |d|<0.2 negligible, 0.2-0.5 small, 0.5-0.8 medium, >0.8 large")

    # Save to file if output_path provided
    if output_path is not None:
        stats_file = output_path / 'statistical_analysis.txt'
        with open(stats_file, 'w') as f:
            f.write("STATISTICAL ANALYSIS - Balanced Warmup Variance Comparison\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"{stat_results['omnibus_test']}:\n")
            f.write(f"  Statistic: {stat_results['omnibus_statistic']:.4f}\n")
            f.write(f"  p-value: {stat_results['omnibus_p']:.4e}\n\n")
            f.write("Pairwise Comparisons (Bonferroni-corrected):\n")
            f.write("-" * 80 + "\n")
            for pair in stat_results['pairwise']:
                comparison = f"{pair['group1']} vs {pair['group2']}"
                sig = get_significance_symbol(pair['p_corrected'])
                f.write(f"{comparison}: p={pair['p_corrected']:.4e} ({sig}), Cohen's d={pair['cohens_d']:.3f}\n")
        print(f"\nStatistical analysis saved to: {stats_file}")


# ============================================================================
# Strategy Run Functions
# ============================================================================

def run_balanced_per_reagent_variance(replicate_num, ref_df, top_100_ref, seed=None):
    """Run current RWS with BalancedWarmup and per-reagent variance (DEFAULT)."""
    print(f"\n  Replicate {replicate_num + 1}/{NUM_REPLICATES}")

    # Create evaluator
    evaluator = LookupEvaluator(CURRENT_EVALUATOR_CONFIG)

    # Create strategy
    strategy = RouletteWheelSelection(
        mode="minimize",
        **RWS_CONFIG
    )

    # Create warmup with per-reagent variance (DEFAULT)
    warmup = BalancedWarmup(
        observations_per_reagent=OBSERVATIONS_PER_REAGENT,
        seed=seed,
        use_per_reagent_variance=True,  # KEY: Per-reagent variance
        shrinkage_strength=3.0
    )

    # Create sampler
    sampler = ThompsonSampler(
        selection_strategy=strategy,
        warmup_strategy=warmup,
        batch_size=BATCH_SIZE,
        use_boltzmann_weighting=True
    )

    # Setup evaluator and reagents
    sampler.set_evaluator(evaluator)
    sampler.read_reagents(REAGENT_FILES)
    sampler.set_hide_progress(True)

    # Run warmup
    warmup_results = sampler.warm_up(num_warmup_trials=OBSERVATIONS_PER_REAGENT)

    # Analyze reagent variances after warmup
    variance_df = analyze_reagent_variances(sampler)

    # Run search
    search_results = sampler.search(num_cycles=100000, max_evaluations=NUM_TS_ITERATIONS)

    # Combine results
    results_df = combine_results(warmup_results, search_results)

    # Calculate metrics
    metrics = calculate_recovery_metrics(results_df, top_100_ref)

    # Add variance statistics
    if len(variance_df) > 0:
        metrics['variance_mean'] = variance_df['variance'].mean()
        metrics['variance_std'] = variance_df['variance'].std()
        metrics['variance_min'] = variance_df['variance'].min()
        metrics['variance_max'] = variance_df['variance'].max()

    print(f"    Recovery: {metrics['recovery_rate']:.1%} ({metrics['recovered_count']}/100)")
    print(f"    Top 10 recovery: {metrics['top_10_recovery']:.1%}")
    print(f"    Best score: {metrics['best_score']:.3f}")
    print(f"    Total evaluated: {metrics['total_evaluated']:,}")
    if 'variance_std' in metrics:
        print(f"    Variance spread: {metrics['variance_min']:.4f} - {metrics['variance_max']:.4f}")

    return metrics, results_df, variance_df


def run_balanced_global_variance(replicate_num, ref_df, top_100_ref, seed=None):
    """Run current RWS with BalancedWarmup and global variance."""
    print(f"\n  Replicate {replicate_num + 1}/{NUM_REPLICATES}")

    # Create evaluator
    evaluator = LookupEvaluator(CURRENT_EVALUATOR_CONFIG)

    # Create strategy
    strategy = RouletteWheelSelection(
        mode="minimize",
        **RWS_CONFIG
    )

    # Create warmup with global variance
    warmup = BalancedWarmup(
        observations_per_reagent=OBSERVATIONS_PER_REAGENT,
        seed=seed,
        use_per_reagent_variance=False,  # KEY: Global variance
        shrinkage_strength=3.0
    )

    # Create sampler
    sampler = ThompsonSampler(
        selection_strategy=strategy,
        warmup_strategy=warmup,
        batch_size=BATCH_SIZE,
        use_boltzmann_weighting=True
    )

    # Setup evaluator and reagents
    sampler.set_evaluator(evaluator)
    sampler.read_reagents(REAGENT_FILES)
    sampler.set_hide_progress(True)

    # Run warmup
    warmup_results = sampler.warm_up(num_warmup_trials=OBSERVATIONS_PER_REAGENT)

    # Analyze reagent variances after warmup
    variance_df = analyze_reagent_variances(sampler)

    # Run search
    search_results = sampler.search(num_cycles=100000, max_evaluations=NUM_TS_ITERATIONS)

    # Combine results
    results_df = combine_results(warmup_results, search_results)

    # Calculate metrics
    metrics = calculate_recovery_metrics(results_df, top_100_ref)

    # Add variance statistics
    if len(variance_df) > 0:
        metrics['variance_mean'] = variance_df['variance'].mean()
        metrics['variance_std'] = variance_df['variance'].std()
        metrics['variance_min'] = variance_df['variance'].min()
        metrics['variance_max'] = variance_df['variance'].max()

    print(f"    Recovery: {metrics['recovery_rate']:.1%} ({metrics['recovered_count']}/100)")
    print(f"    Top 10 recovery: {metrics['top_10_recovery']:.1%}")
    print(f"    Best score: {metrics['best_score']:.3f}")
    print(f"    Total evaluated: {metrics['total_evaluated']:,}")
    if 'variance_std' in metrics:
        print(f"    Variance (global): {metrics['variance_mean']:.4f} (all same)")

    return metrics, results_df, variance_df


def run_enhanced_warmup_baseline(replicate_num, ref_df, top_100_ref):
    """Run current RWS with EnhancedWarmup (legacy warmup, baseline)."""
    print(f"\n  Replicate {replicate_num + 1}/{NUM_REPLICATES}")

    # Create evaluator
    evaluator = LookupEvaluator(CURRENT_EVALUATOR_CONFIG)

    # Create strategy
    strategy = RouletteWheelSelection(
        mode="minimize",
        **RWS_CONFIG
    )

    # Create sampler with EnhancedWarmup (legacy)
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

    # Run warmup
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
    sampler = LegacyRWSSampler(processes=1, scaling=-1)

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


# ============================================================================
# Main Comparison
# ============================================================================

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Compare per-reagent vs global variance in BalancedWarmup strategy",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python balanced_warmup_variance_comparison.py --output_dir ./results/docking_run_001
    python balanced_warmup_variance_comparison.py -o ./results/$(date +%Y%m%d_%H%M%S)
        """
    )
    parser.add_argument(
        "-o", "--output_dir",
        type=str,
        default=None,
        help="Directory to save output files (CSV and PDF). Created if doesn't exist. "
             "Default: current directory with timestamp prefix."
    )
    return parser.parse_args()


def main(output_dir: str = None):
    """Run comparison experiments.

    Args:
        output_dir: Directory to save output files. If None, uses current directory.
    """
    # Setup output directory
    if output_dir is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_path = Path(f"./balanced_warmup_docking_{timestamp}")
    else:
        output_path = Path(output_dir)

    output_path.mkdir(parents=True, exist_ok=True)

    print("="*80)
    print("Balanced Warmup Variance Comparison - Docking Dataset")
    print("="*80)
    print(f"\nOutput directory: {output_path.absolute()}")
    print(f"\nKey parameters:")
    print(f"  - Observations per reagent (K): {OBSERVATIONS_PER_REAGENT}")
    print(f"  - TS iterations: {NUM_TS_ITERATIONS}")
    print(f"  - Replicates: {NUM_REPLICATES}")
    print(f"  - Shrinkage strength: 3.0")

    # Load reference data
    ref_df, top_100_ref = load_reference_data()

    # Storage for results
    all_metrics = {
        'Balanced + Per-Reagent Var': [],
        'Balanced + Global Var': [],
        'Enhanced (Baseline)': [],
        'Legacy RWS': []
    }

    all_variance_dfs = {
        'Balanced + Per-Reagent Var': [],
        'Balanced + Global Var': []
    }

    # Run experiments
    for replicate in range(NUM_REPLICATES):
        print(f"\n{'='*80}")
        print(f"Replicate {replicate + 1}/{NUM_REPLICATES}")
        print(f"{'='*80}")

        # Use same seed for balanced warmups to ensure same combinations
        seed = 42 + replicate

        # 1. Balanced + Per-reagent variance (NEW DEFAULT)
        print("\n[1/4] Running Balanced + Per-Reagent Variance (DEFAULT)...")
        try:
            metrics, _, var_df = run_balanced_per_reagent_variance(replicate, ref_df, top_100_ref, seed=seed)
            all_metrics['Balanced + Per-Reagent Var'].append(metrics)
            all_variance_dfs['Balanced + Per-Reagent Var'].append(var_df)
        except Exception as e:
            print(f"    ERROR: {e}")
            import traceback
            traceback.print_exc()

        # 2. Balanced + Global variance
        print("\n[2/4] Running Balanced + Global Variance...")
        try:
            metrics, _, var_df = run_balanced_global_variance(replicate, ref_df, top_100_ref, seed=seed)
            all_metrics['Balanced + Global Var'].append(metrics)
            all_variance_dfs['Balanced + Global Var'].append(var_df)
        except Exception as e:
            print(f"    ERROR: {e}")
            import traceback
            traceback.print_exc()

        # 3. Enhanced warmup baseline
        print("\n[3/4] Running Enhanced Warmup (Baseline)...")
        try:
            metrics, _ = run_enhanced_warmup_baseline(replicate, ref_df, top_100_ref)
            all_metrics['Enhanced (Baseline)'].append(metrics)
        except Exception as e:
            print(f"    ERROR: {e}")
            import traceback
            traceback.print_exc()

        # 4. Legacy RWS
        print("\n[4/4] Running Legacy RWS...")
        try:
            metrics, _ = run_legacy_rws(replicate, ref_df, top_100_ref)
            all_metrics['Legacy RWS'].append(metrics)
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

        # Key comparison
        print("\n" + "-"*80)
        print("KEY COMPARISON: Per-Reagent vs Global Variance")
        print("-"*80)

        per_reagent_recovery = [m['recovery_rate'] for m in all_metrics.get('Balanced + Per-Reagent Var', [])]
        global_recovery = [m['recovery_rate'] for m in all_metrics.get('Balanced + Global Var', [])]

        if per_reagent_recovery and global_recovery:
            pr_mean, pr_std = np.mean(per_reagent_recovery), np.std(per_reagent_recovery)
            gl_mean, gl_std = np.mean(global_recovery), np.std(global_recovery)

            print(f"  Per-Reagent Variance: {pr_mean:.1%} +/- {pr_std:.1%}")
            print(f"  Global Variance:      {gl_mean:.1%} +/- {gl_std:.1%}")
            print(f"  Difference:           {(pr_mean - gl_mean)*100:+.1f}%")

            # Check if per-reagent variance reduced run-to-run variability
            print(f"\n  Run-to-run variability:")
            print(f"    Per-Reagent Var std: {pr_std:.1%}")
            print(f"    Global Var std:      {gl_std:.1%}")
            if pr_std < gl_std:
                print(f"    -> Per-reagent variance REDUCED variability by {(gl_std-pr_std)*100:.1f}%")
            else:
                print(f"    -> No reduction in variability")

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
            detailed_csv = output_path / 'balanced_warmup_variance_detailed.csv'
            detailed_df.to_csv(detailed_csv, index=False)
            print(f"\nDetailed results saved to: {detailed_csv}")

            summary_csv = output_path / 'balanced_warmup_variance_summary.csv'
            summary_df.to_csv(summary_csv, index=False)
            print(f"Summary results saved to: {summary_csv}")

            # ================================================================
            # Statistical Analysis
            # ================================================================
            print("\nRunning statistical tests...")
            stat_results = run_statistical_tests(detailed_df, metric='recovery_rate')
            print_statistical_summary(stat_results, output_path)

            # ================================================================
            # Visualization
            # ================================================================
            print("\nGenerating comparison plots...")

            fig, axes = plt.subplots(2, 2, figsize=(14, 12))

            # Recovery rate comparison
            recovery_data = detailed_df[['Strategy', 'recovery_rate']].copy()
            recovery_data['recovery_rate'] *= 100

            sns.barplot(data=recovery_data, x='Strategy', y='recovery_rate', ax=axes[0, 0], errorbar='sd')
            axes[0, 0].set_title('Top 100 Recovery Rate', fontsize=12, fontweight='bold')
            axes[0, 0].set_ylabel('Recovery Rate (%)')
            axes[0, 0].set_xlabel('')
            axes[0, 0].set_xticklabels(axes[0, 0].get_xticklabels(), rotation=20, ha='right')
            axes[0, 0].grid(axis='y', alpha=0.3)

            # Top 10 recovery comparison
            top10_data = detailed_df[['Strategy', 'top_10_recovery']].copy()
            top10_data['top_10_recovery'] *= 100

            sns.barplot(data=top10_data, x='Strategy', y='top_10_recovery', ax=axes[0, 1], errorbar='sd')
            axes[0, 1].set_title('Top 10 Recovery Rate', fontsize=12, fontweight='bold')
            axes[0, 1].set_ylabel('Recovery Rate (%)')
            axes[0, 1].set_xlabel('')
            axes[0, 1].set_xticklabels(axes[0, 1].get_xticklabels(), rotation=20, ha='right')
            axes[0, 1].grid(axis='y', alpha=0.3)

            # Best score comparison
            sns.barplot(data=detailed_df, x='Strategy', y='best_score', ax=axes[1, 0], errorbar='sd')
            axes[1, 0].set_title('Best Score Found', fontsize=12, fontweight='bold')
            axes[1, 0].set_ylabel('Docking Score')
            axes[1, 0].set_xlabel('')
            axes[1, 0].set_xticklabels(axes[1, 0].get_xticklabels(), rotation=20, ha='right')
            axes[1, 0].grid(axis='y', alpha=0.3)

            # Recovery variability comparison (box plot with significance annotations)
            # Need to scale stat_results for percentage display
            stat_results_scaled = stat_results.copy()
            stat_results_scaled['means'] = [m * 100 for m in stat_results['means']]
            stat_results_scaled['stds'] = [s * 100 for s in stat_results['stds']]

            sns.boxplot(data=recovery_data, x='Strategy', y='recovery_rate', ax=axes[1, 1])
            axes[1, 1].set_title('Recovery Rate Distribution with Pairwise Comparisons', fontsize=12, fontweight='bold')
            axes[1, 1].set_ylabel('Recovery Rate (%)')
            axes[1, 1].set_xlabel('')
            axes[1, 1].set_xticklabels(axes[1, 1].get_xticklabels(), rotation=20, ha='right')
            axes[1, 1].grid(axis='y', alpha=0.3)

            # Add significance annotations to boxplot
            add_significance_annotations(axes[1, 1], recovery_data, 'recovery_rate', stat_results)

            plt.tight_layout()
            comparison_pdf = output_path / 'balanced_warmup_variance_comparison.pdf'
            plt.savefig(comparison_pdf, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"Comparison plot saved to: {comparison_pdf}")

            # ================================================================
            # Variance Distribution Analysis
            # ================================================================
            if all_variance_dfs['Balanced + Per-Reagent Var']:
                print("\nGenerating variance distribution plots...")

                fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))

                # Per-reagent variance distribution
                var_df_pr = all_variance_dfs['Balanced + Per-Reagent Var'][0]
                if len(var_df_pr) > 0:
                    sns.histplot(data=var_df_pr, x='variance', ax=axes2[0], bins=50)
                    axes2[0].set_title('Per-Reagent Variance Distribution', fontsize=12, fontweight='bold')
                    axes2[0].set_xlabel('Variance')
                    axes2[0].axvline(var_df_pr['variance'].mean(), color='r', linestyle='--', label=f'Mean: {var_df_pr["variance"].mean():.4f}')
                    axes2[0].legend()

                # Global variance (should be constant)
                var_df_gl = all_variance_dfs['Balanced + Global Var'][0]
                if len(var_df_gl) > 0:
                    sns.histplot(data=var_df_gl, x='variance', ax=axes2[1], bins=50)
                    axes2[1].set_title('Global Variance Distribution', fontsize=12, fontweight='bold')
                    axes2[1].set_xlabel('Variance')
                    axes2[1].axvline(var_df_gl['variance'].mean(), color='r', linestyle='--', label=f'Mean: {var_df_gl["variance"].mean():.4f}')
                    axes2[1].legend()

                plt.tight_layout()
                distribution_pdf = output_path / 'balanced_warmup_variance_distribution.pdf'
                plt.savefig(distribution_pdf, dpi=300, bbox_inches='tight')
                plt.close(fig2)
                print(f"Variance distribution plot saved to: {distribution_pdf}")

    else:
        print("\nERROR: No results to summarize!")

    print("\n" + "="*80)
    print("COMPARISON COMPLETE")
    print("="*80)
    print(f"\nAll outputs saved to: {output_path.absolute()}")


if __name__ == "__main__":
    args = parse_args()
    main(output_dir=args.output_dir)
