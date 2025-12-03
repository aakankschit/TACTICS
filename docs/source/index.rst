.. TACTICS documentation master file, created by
   sphinx-quickstart on Thu Apr 24 15:14:03 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to TACTICS documentation!
=================================

TACTICS (Thompson Sampling-Assisted Chemical Targeting and Iterative Compound Selection for Drug Discovery) is a comprehensive Python package for Thompson Sampling-based optimization of chemical combinatorial libraries with a unified, flexible architecture.

Key Features
------------

* **Unified Thompson Sampler**: Single sampler class that accepts different selection strategies
* **Multiple Selection Strategies**: Greedy, Roulette Wheel (thermal cycling), UCB, Epsilon-Greedy, Bayes-UCB
* **Flexible Warmup Strategies**: Balanced (recommended), Standard, Enhanced (legacy)
* **Multiple Evaluators**: Lookup, Database, Fingerprint, ML models, ROCS, FRED docking
* **Parallel Evaluation**: Built-in multiprocessing support for expensive evaluators
* **Pydantic Configuration**: Type-safe configuration with validation and presets
* **Configuration Presets**: Pre-configured setups for common use cases
* **Extensible Design**: Easy integration of custom strategies, warmup methods, and evaluators
* **Library Enumeration**: Efficient generation of combinatorial reaction products
* **SMARTS Toolkit**: Advanced pattern validation, multi-SMARTS routing, and multi-step synthesis
* **Library Analysis**: Comprehensive analysis and visualization tools

Quick Start
-----------

Simple Out-of-the-Box Usage (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The easiest way to get started with TACTICS is using the ``run_ts()`` convenience function with presets:

.. code-block:: python

    from TACTICS.thompson_sampling import run_ts, get_preset
    from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

    # 1. Create evaluator config
    evaluator = LookupEvaluatorConfig(ref_filename="scores.csv")

    # 2. Get a preset configuration
    config = get_preset(
        "fast_exploration",  # Quick screening with epsilon-greedy
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        reagent_file_list=["acids.smi", "amines.smi"],
        evaluator_config=evaluator,
        mode="minimize",  # Use "minimize" for docking scores
        num_iterations=1000
    )

    # 3. Run and get results (returns polars DataFrame)
    results_df = run_ts(config)

**Available Presets:**

*Modern Presets* (use BalancedWarmup):

* ``fast_exploration`` - Epsilon-greedy strategy, quick screening
* ``parallel_batch`` - Batch processing with multiprocessing (for slow evaluators)
* ``conservative_exploit`` - Greedy strategy, focus on best reagents
* ``balanced_sampling`` - UCB strategy with theoretical guarantees
* ``diverse_coverage`` - Maximum diversity exploration

*Legacy Presets* (for replicating published RWS results):

* ``legacy_rws_maximize`` - Original RWS algorithm with Boltzmann weighting (maximize mode)
* ``legacy_rws_minimize`` - Original RWS algorithm with Boltzmann weighting (minimize mode, e.g., docking)

Parallel Batch Processing
~~~~~~~~~~~~~~~~~~~~~~~~~~

For expensive evaluators (docking, ML models), use parallel batch mode:

.. code-block:: python

    from TACTICS.thompson_sampling import run_ts, get_preset
    from TACTICS.thompson_sampling.core.evaluator_config import FredEvaluatorConfig

    # Configure slow evaluator
    evaluator = FredEvaluatorConfig(design_unit_file="receptor.oedu")

    # Get parallel batch preset
    config = get_preset(
        "parallel_batch",
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        reagent_file_list=["acids.smi", "amines.smi"],
        evaluator_config=evaluator,
        mode="minimize",  # Docking scores
        batch_size=100,   # Sample 100 compounds per cycle
        processes=8       # Use 8 CPU cores
    )

    results_df = run_ts(config)

Advanced Usage: Direct Sampler Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For maximum flexibility, use ``ThompsonSampler`` directly:

.. code-block:: python

    from TACTICS.thompson_sampling import ThompsonSampler
    from TACTICS.thompson_sampling.presets import get_preset
    from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

    # Use a preset configuration
    config = get_preset(
        "balanced_sampling",
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        reagent_file_list=["reagents1.smi", "reagents2.smi"],
        evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv")
    )

    # Create sampler from config
    sampler = ThompsonSampler.from_config(config)
    sampler.set_reaction(config.reaction_smarts)

    # Run optimization with full control
    warmup_results = sampler.warm_up(num_warmup_trials=config.num_warmup_trials)
    search_results = sampler.search(num_cycles=config.num_ts_iterations)

    # Cleanup
    sampler.close()

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   thompson_sampling
   api_reference
   configuration
   library_enumeration
   library_analysis

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

