"""
Example: Using Adaptive Bayes-UCB Selection Strategy

This example demonstrates how to use the new Adaptive Bayes-UCB selection strategy
with thermal cycling for Thompson Sampling.

Bayes-UCB provides:
- Theoretically grounded confidence bounds using Student-t quantiles
- Adaptive exploration/exploitation via percentile parameters
- Thermal cycling across components
- Automatic adaptation based on sampling efficiency
"""

from TACTICS.thompson_sampling.config import ThompsonSamplingConfig
from TACTICS.thompson_sampling.strategies.config import BayesUCBConfig
from TACTICS.thompson_sampling.warmup.config import StratifiedWarmupConfig
from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig
from TACTICS.thompson_sampling.core.sampler import ThompsonSampler


def run_bayes_ucb_example():
    """Run Thompson Sampling with Adaptive Bayes-UCB strategy."""

    # Configure Bayes-UCB strategy
    bayes_ucb_config = BayesUCBConfig(
        mode="maximize",
        initial_p_high=0.90,      # Initial high percentile (heated component)
        initial_p_low=0.90,       # Initial low percentile (cooled components)
        efficiency_threshold=0.10, # 10% unique compounds triggers adaptation
        p_high_bounds=(0.85, 0.995),  # Bounds for high percentile
        p_low_bounds=(0.50, 0.90),    # Bounds for low percentile
        delta_high=0.01,          # Step size for increasing exploration
        delta_low=0.005           # Step size for decreasing exploration
    )

    # Main Thompson Sampling config (Modern approach)
    config = ThompsonSamplingConfig(
        reaction_smarts="[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]>>[#6:1](=[O:2])[#7:3]",
        reagent_file_list=[
            "data/reagents/thrombin/acids.smi",
            "data/reagents/thrombin/coupled_aa_sub.smi"
        ],
        num_ts_iterations=100,
        num_warmup_trials=3,

        # Use BayesUCB strategy
        strategy_config=bayes_ucb_config,

        # Optional: Configure warmup
        warmup_config=StratifiedWarmupConfig(),

        # Configure evaluator
        evaluator_config=LookupEvaluatorConfig(
            ref_filename="data/scores/thrombin/product_scores.csv",
            ref_colname="Score"
        ),

        # Optional: Batch sampling
        batch_size=1,  # Single compound per cycle

        # Optional: Multiprocessing
        processes=1,

        # Output
        results_filename="bayes_ucb_results.csv",
        log_filename="bayes_ucb.log"
    )

    # Create sampler from config
    print("Creating Thompson Sampler with Adaptive Bayes-UCB...")
    sampler = ThompsonSampler.from_config(config)

    # Set up reaction
    sampler.set_reaction(config.reaction_smarts)

    # Run warmup
    print("\nRunning warmup phase...")
    warmup_results = sampler.warm_up(config.num_warmup_trials)
    print(f"Warmup completed: {len(warmup_results)} evaluations")
    print(f"Best warmup score: {warmup_results['score'].max():.3f}")

    # Run search
    print("\nRunning Thompson Sampling search with Bayes-UCB...")
    search_results = sampler.search(config.num_ts_iterations)
    print(f"Search completed: {len(search_results)} evaluations")
    print(f"Best search score: {search_results['score'].max():.3f}")

    # Show final percentile values (after adaptation)
    strategy = sampler.selection_strategy
    print(f"\nFinal percentile values:")
    print(f"  p_high (heated component): {strategy.p_high:.3f}")
    print(f"  p_low (cooled components): {strategy.p_low:.3f}")

    # Clean up
    sampler.close()

    print("\nResults saved to bayes_ucb_results.csv")
    return search_results


def run_bayes_ucb_legacy_example():
    """Run Thompson Sampling with Bayes-UCB using legacy config format."""

    config = ThompsonSamplingConfig(
        reaction_smarts="[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]>>[#6:1](=[O:2])[#7:3]",
        reagent_file_list=[
            "data/reagents/thrombin/acids.smi",
            "data/reagents/thrombin/coupled_aa_sub.smi"
        ],
        num_ts_iterations=100,
        num_warmup_trials=3,

        # Legacy approach: String-based config
        selection_strategy="bayes_ucb",
        mode="maximize",
        strategy_params={
            "initial_p_high": 0.90,
            "initial_p_low": 0.90,
            "efficiency_threshold": 0.10
        },

        # Legacy evaluator config
        evaluator_class_name="LookupEvaluator",
        evaluator_arg={
            "ref_filename": "data/scores/thrombin/product_scores.csv",
            "ref_colname": "Score"
        },

        results_filename="bayes_ucb_legacy_results.csv"
    )

    print("Creating sampler with legacy config...")
    sampler = ThompsonSampler.from_config(config)
    sampler.set_reaction(config.reaction_smarts)

    print("Running warmup and search...")
    sampler.warm_up(config.num_warmup_trials)
    results = sampler.search(config.num_ts_iterations)

    print(f"Completed: {len(results)} evaluations")
    print(f"Best score: {results['score'].max():.3f}")

    sampler.close()
    return results


def direct_instantiation_example():
    """Example: Direct instantiation of BayesUCBSelection (without config)."""

    from TACTICS.thompson_sampling.strategies import BayesUCBSelection
    from TACTICS.thompson_sampling.core.evaluators import LookupEvaluator

    # Create strategy directly
    strategy = BayesUCBSelection(
        mode="maximize",
        initial_p_high=0.95,
        initial_p_low=0.80,
        efficiency_threshold=0.10
    )

    # Create sampler
    sampler = ThompsonSampler(selection_strategy=strategy)

    # Configure sampler
    sampler.read_reagents([
        "data/reagents/thrombin/acids.smi",
        "data/reagents/thrombin/coupled_aa_sub.smi"
    ])
    sampler.set_reaction(
        "[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]>>[#6:1](=[O:2])[#7:3]"
    )
    sampler.set_evaluator(LookupEvaluator({
        "ref_filename": "data/scores/thrombin/product_scores.csv"
    }))

    # Run sampling
    print("Running with direct instantiation...")
    sampler.warm_up(num_warmup_trials=3)
    results = sampler.search(num_cycles=50)

    print(f"Best score: {results['score'].max():.3f}")

    sampler.close()
    return results


if __name__ == "__main__":
    print("=" * 60)
    print("Adaptive Bayes-UCB Selection Strategy Example")
    print("=" * 60)

    # Run modern config example
    print("\n1. Modern Config Approach")
    print("-" * 60)
    results = run_bayes_ucb_example()

    print("\n\n2. Legacy Config Approach")
    print("-" * 60)
    results_legacy = run_bayes_ucb_legacy_example()

    print("\n\n3. Direct Instantiation")
    print("-" * 60)
    results_direct = direct_instantiation_example()

    print("\n" + "=" * 60)
    print("All examples completed successfully!")
    print("=" * 60)
