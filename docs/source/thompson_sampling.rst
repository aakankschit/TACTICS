Thompson Sampling Module
========================

The Thompson Sampling module implements a unified, flexible Thompson sampling framework for chemical library exploration with support for multiple selection strategies, evaluators, and warmup approaches.

Overview
--------

The package provides a **unified ThompsonSampler** that accepts different strategies:

* **Selection Strategies**: Control how reagents are selected (greedy, roulette wheel, UCB, epsilon-greedy)
* **Warmup Strategies**: Determine how reagent priors are initialized (standard, stratified, enhanced, Latin hypercube)
* **Evaluators**: Score compounds using different methods (lookup, database, fingerprint, ML models, docking)

This modular design allows easy experimentation with different optimization approaches and simple integration of custom strategies.

Core Sampler
------------

.. autoclass:: TACTICS.thompson_sampling.core.sampler.ThompsonSampler
   :members:
   :undoc-members:
   :show-inheritance:

   .. automethod:: __init__
   .. automethod:: from_config
   .. automethod:: warm_up
   .. automethod:: search
   .. automethod:: evaluate
   .. automethod:: evaluate_batch
   .. automethod:: close

   The ThompsonSampler is the main class for Thompson Sampling optimization:

   * **Flexible Strategy System**: Use any selection strategy by passing a SelectionStrategy instance
   * **Parallel Evaluation**: Built-in support for multiprocessing (processes > 1)
   * **Batch Sampling**: Sample multiple compounds per cycle (batch_size > 1)
   * **Warmup Strategies**: Initialize priors using different sampling approaches
   * **Progress Tracking**: Optional progress bars and logging
   * **Config-based Setup**: Create sampler from Pydantic configs using ``from_config()``

Selection Strategies
-------------------

Selection strategies determine how reagents are chosen during the search phase.

.. autoclass:: TACTICS.thompson_sampling.strategies.base_strategy.SelectionStrategy
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: TACTICS.thompson_sampling.strategies.greedy_selection.GreedySelection
   :members:
   :undoc-members:
   :show-inheritance:

   Simple greedy selection (argmax/argmin of sampled scores):

   * Direct optimization toward best-performing reagents
   * Fastest convergence but may get stuck in local optima
   * Best for: Simple optimization landscapes, limited computational budgets

.. autoclass:: TACTICS.thompson_sampling.strategies.roulette_wheel.RouletteWheelSelection
   :members:
   :undoc-members:
   :show-inheritance:

   Roulette wheel selection with adaptive thermal cycling:

   * Better exploration/exploitation balance via thermal cycling
   * Adaptive temperature control based on sampling efficiency
   * Component rotation for systematic exploration
   * Best for: Complex multi-modal landscapes, large libraries

.. autoclass:: TACTICS.thompson_sampling.strategies.ucb_selection.UCBSelection
   :members:
   :undoc-members:
   :show-inheritance:

   Upper Confidence Bound (UCB) selection:

   * Balances exploitation and exploration via confidence bounds
   * Deterministic selection based on UCB values
   * Best for: Situations where deterministic behavior is preferred

.. autoclass:: TACTICS.thompson_sampling.strategies.epsilon_greedy.EpsilonGreedy
   :members:
   :undoc-members:
   :show-inheritance:

   Epsilon-greedy selection:

   * Simple exploration strategy: random selection with probability ε
   * Greedy selection with probability 1-ε
   * Best for: Baseline comparisons, simple exploration needs

Warmup Strategies
-----------------

Warmup strategies determine how reagent combinations are sampled to initialize posteriors.

.. autoclass:: TACTICS.thompson_sampling.warmup.base.WarmupStrategy
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: TACTICS.thompson_sampling.warmup.standard.StandardWarmup
   :members:
   :undoc-members:
   :show-inheritance:

   Standard warmup: Each reagent tested num_trials times with random partners

   * Simple and straightforward
   * Ensures all reagents are evaluated
   * Expected evaluations: (sum of reagents) × num_trials

.. autoclass:: TACTICS.thompson_sampling.warmup.stratified.StratifiedWarmup
   :members:
   :undoc-members:
   :show-inheritance:

   Stratified warmup: Ensures balanced coverage across all components

   * More uniform coverage than standard warmup
   * Reduces bias from component imbalances

.. autoclass:: TACTICS.thompson_sampling.warmup.enhanced.EnhancedWarmup
   :members:
   :undoc-members:
   :show-inheritance:

   Enhanced warmup: Anchor-based approach for better diversity

   * Uses anchor compounds to explore diverse regions
   * Multiple anchor strategies (random, max_variance, etc.)
   * Better for complex libraries

.. autoclass:: TACTICS.thompson_sampling.warmup.latin_hypercube.LatinHypercubeWarmup
   :members:
   :undoc-members:
   :show-inheritance:

   Latin hypercube sampling for space-filling warmup

   * Ensures even coverage of the search space
   * Particularly useful for large combinatorial libraries

Evaluator Classes
----------------

.. autoclass:: TACTICS.thompson_sampling.core.evaluators.LookupEvaluator
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: TACTICS.thompson_sampling.core.evaluators.ROCSEvaluator
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: TACTICS.thompson_sampling.core.evaluators.DBEvaluator
   :members:
   :undoc-members:
   :show-inheritance:

Utility Functions
----------------

.. autofunction:: TACTICS.thompson_sampling.utils.ts_utils.read_reagents
.. autofunction:: TACTICS.thompson_sampling.utils.ts_utils.create_reagents
.. autofunction:: TACTICS.thompson_sampling.utils.ts_logger.get_logger

Legacy Interface
----------------

For backward compatibility, the legacy interface is still available:

.. autofunction:: TACTICS.thompson_sampling.legacy.ts_main.run_ts
   :noindex:

.. autoclass:: TACTICS.thompson_sampling.legacy.standard_thompson_sampling.StandardThompsonSampler
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: TACTICS.thompson_sampling.legacy.enhanced_thompson_sampling.EnhancedThompsonSampler
   :members:
   :undoc-members:
   :show-inheritance:

Configuration-Based Approach
-----------------------------

The modern way to use Thompson Sampling is through Pydantic configuration:

.. code-block:: python

    from TACTICS.thompson_sampling import ThompsonSampler
    from TACTICS.thompson_sampling.config import ThompsonSamplingConfig
    from TACTICS.thompson_sampling.strategies.config import RouletteWheelConfig
    from TACTICS.thompson_sampling.warmup.config import StratifiedWarmupConfig
    from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

    # Create configuration
    config = ThompsonSamplingConfig(
        reaction_smarts="[#6:1](=[O:2])[OH].[#7:3]>>[#6:1](=[O:2])[#7:3]",
        reagent_file_list=["acids.smi", "amines.smi"],
        num_ts_iterations=1000,
        num_warmup_trials=3,
        strategy_config=RouletteWheelConfig(mode="maximize", alpha=0.1, beta=0.1),
        warmup_config=StratifiedWarmupConfig(),
        evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv")
    )

    # Create sampler from config (handles all setup automatically)
    sampler = ThompsonSampler.from_config(config)
    sampler.set_reaction(config.reaction_smarts)

    # Run optimization
    warmup_results = sampler.warm_up(num_warmup_trials=config.num_warmup_trials)
    search_results = sampler.search(num_cycles=config.num_ts_iterations)
    sampler.close()

Using Configuration Presets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For common use cases, use pre-configured presets:

.. code-block:: python

    from TACTICS.thompson_sampling.presets import get_preset
    from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

    # Use a preset configuration
    config = get_preset(
        "parallel_batch",  # Fast exploration, parallel batch, etc.
        reaction_smarts="[#6:1](=[O:2])[OH].[#7:3]>>[#6:1](=[O:2)][#7:3]",
        reagent_file_list=["acids.smi", "amines.smi"],
        evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv"),
        batch_size=100,
        processes=4
    )

    sampler = ThompsonSampler.from_config(config)
    sampler.set_reaction(config.reaction_smarts)
    warmup_results = sampler.warm_up(num_warmup_trials=config.num_warmup_trials)
    search_results = sampler.search(num_cycles=config.num_ts_iterations)
    sampler.close()

Convenience Function: run_ts()
-------------------------------

For the simplest out-of-the-box usage, use the ``run_ts()`` convenience wrapper:

.. autofunction:: TACTICS.thompson_sampling.main.run_ts

The ``run_ts()`` function provides automatic setup, logging, file saving, and cleanup:

Simple Usage with Presets
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from TACTICS.thompson_sampling import run_ts, get_preset
    from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

    # 1. Create evaluator config
    evaluator = LookupEvaluatorConfig(ref_filename="scores.csv")

    # 2. Get preset
    config = get_preset(
        "fast_exploration",
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        reagent_file_list=["acids.smi", "amines.smi"],
        evaluator_config=evaluator,
        mode="minimize"  # For docking scores
    )

    # 3. Run and get results (returns polars DataFrame)
    results_df = run_ts(config)

With Warmup Results
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from TACTICS.thompson_sampling import run_ts, get_preset
    from TACTICS.thompson_sampling.core.evaluator_config import DBEvaluatorConfig

    evaluator = DBEvaluatorConfig(db_filename="scores.db")
    config = get_preset(
        "balanced_sampling",
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        reagent_file_list=["acids.smi", "amines.smi"],
        evaluator_config=evaluator,
        num_iterations=1000
    )

    # Get both search and warmup results
    search_df, warmup_df = run_ts(config, return_warmup=True)

    # Analyze warmup phase
    print(warmup_df.head())

    # Analyze search results
    print(search_df.sort("score").head(10))

Parallel Batch Processing
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from TACTICS.thompson_sampling import run_ts, get_preset
    from TACTICS.thompson_sampling.core.evaluator_config import FredEvaluatorConfig

    # For slow evaluators (docking, ML models)
    evaluator = FredEvaluatorConfig(design_unit_file="receptor.oedu")

    config = get_preset(
        "parallel_batch",
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        reagent_file_list=["acids.smi", "amines.smi"],
        evaluator_config=evaluator,
        mode="minimize",  # Docking scores
        batch_size=100,   # Sample 100 per cycle
        processes=8       # Use 8 CPU cores
    )

    results_df = run_ts(config)

Custom Configuration
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from TACTICS.thompson_sampling import ThompsonSamplingConfig, run_ts
    from TACTICS.thompson_sampling.strategies.config import EpsilonGreedyConfig
    from TACTICS.thompson_sampling.warmup.config import StratifiedWarmupConfig
    from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

    config = ThompsonSamplingConfig(
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        reagent_file_list=["acids.smi", "amines.smi"],
        num_ts_iterations=1000,
        num_warmup_trials=5,
        strategy_config=EpsilonGreedyConfig(mode="maximize", epsilon=0.2, decay=0.995),
        warmup_config=StratifiedWarmupConfig(),
        evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv"),
        results_filename="results.csv",
        log_filename="run.log"
    )

    results_df = run_ts(config)

Benefits of run_ts()
~~~~~~~~~~~~~~~~~~~~

The ``run_ts()`` function provides several benefits:

* **Automatic logging**: Sets up logging from config
* **Automatic file saving**: Saves results if filename provided
* **Progress display**: Shows top results unless hidden
* **Resource cleanup**: Automatically closes multiprocessing pools
* **Simple API**: 3-line workflow for common cases
* **Flexible**: Works with presets or custom configs
* **Polars DataFrames**: Returns efficient polars DataFrames

For maximum flexibility and control, use ``ThompsonSampler`` directly instead of ``run_ts()``.

See the :doc:`configuration` page for more details on configuration options and presets.

Usage Examples
--------------

Basic Thompson Sampling with Greedy Selection:

.. code-block:: python

    from TACTICS.thompson_sampling.core.sampler import ThompsonSampler
    from TACTICS.thompson_sampling.strategies import GreedySelection
    from TACTICS.thompson_sampling.core.evaluators import LookupEvaluator

    # Create greedy selection strategy
    strategy = GreedySelection(mode="maximize")

    # Create sampler
    sampler = ThompsonSampler(selection_strategy=strategy)

    # Setup
    sampler.read_reagents(["reagents1.smi", "reagents2.smi"])
    sampler.set_reaction("[C:1]=[O:2]>>[C:1][O:2]")
    sampler.set_evaluator(LookupEvaluator({"ref_filename": "scores.csv"}))

    # Run optimization
    warmup_results = sampler.warm_up(num_warmup_trials=3)
    search_results = sampler.search(num_cycles=100)

    # Cleanup
    sampler.close()

Roulette Wheel Selection with Thermal Cycling:

.. code-block:: python

    from TACTICS.thompson_sampling.strategies import RouletteWheelSelection

    # Create roulette wheel strategy with thermal cycling
    strategy = RouletteWheelSelection(
        mode="maximize",
        alpha=0.1,          # Initial heating temperature
        beta=0.1,           # Initial cooling temperature
        scaling=1.0
    )

    sampler = ThompsonSampler(
        selection_strategy=strategy,
        batch_size=10,      # Sample 10 compounds per cycle
        processes=4         # Use 4 CPU cores for parallel evaluation
    )

    # Setup and run as before
    sampler.read_reagents(["reagents1.smi", "reagents2.smi"])
    sampler.set_reaction("[C:1]=[O:2]>>[C:1][O:2]")
    sampler.set_evaluator(LookupEvaluator({"ref_filename": "scores.csv"}))

    warmup_results = sampler.warm_up(num_warmup_trials=3)
    search_results = sampler.search(num_cycles=100)

    sampler.close()

Custom Warmup Strategy:

.. code-block:: python

    from TACTICS.thompson_sampling.warmup import EnhancedWarmup

    # Use enhanced warmup with anchor compounds
    warmup_strategy = EnhancedWarmup(
        anchor_strategy="max_variance",
        num_anchors=5
    )

    strategy = GreedySelection(mode="maximize")

    sampler = ThompsonSampler(
        selection_strategy=strategy,
        warmup_strategy=warmup_strategy  # Use custom warmup
    )

    # Run as usual
    sampler.read_reagents(["reagents1.smi", "reagents2.smi"])
    sampler.set_reaction("[C:1]=[O:2]>>[C:1][O:2]")
    sampler.set_evaluator(LookupEvaluator({"ref_filename": "scores.csv"}))

    warmup_results = sampler.warm_up(num_warmup_trials=3)
    search_results = sampler.search(num_cycles=100)

    sampler.close()

Creating Custom Strategies:

.. code-block:: python

    from TACTICS.thompson_sampling.strategies.base_strategy import SelectionStrategy
    import numpy as np

    class MyCustomStrategy(SelectionStrategy):
        """Custom selection strategy"""

        def select_reagent(self, reagent_list, disallow_mask=None, **kwargs):
            rng = kwargs.get('rng', np.random.default_rng())

            # Sample from posterior distributions
            scores = self.prepare_scores(reagent_list, rng)

            # Apply custom selection logic
            # ... your implementation ...

            # Apply disallow mask
            if disallow_mask:
                scores[np.array(list(disallow_mask))] = np.nan

            # Return selected index
            if self.mode in ["maximize", "maximize_boltzmann"]:
                return np.nanargmax(scores)
            else:
                return np.nanargmin(scores)

    # Use custom strategy
    strategy = MyCustomStrategy(mode="maximize")
    sampler = ThompsonSampler(selection_strategy=strategy)

Strategy Selection Guide
-----------------------

**GreedySelection:**
* Best for: Simple optimization landscapes, limited budgets
* Pros: Fast convergence, simple
* Cons: Can get stuck in local optima

**RouletteWheelSelection:**
* Best for: Complex multi-modal landscapes, large libraries
* Pros: Better exploration/exploitation balance, adaptive
* Cons: More complex, requires parameter tuning

**UCBSelection:**
* Best for: Deterministic optimization needs
* Pros: Theoretically grounded, good exploration
* Cons: Less stochastic exploration than Thompson Sampling

**EpsilonGreedy:**
* Best for: Baseline comparisons, simple needs
* Pros: Very simple, easy to understand
* Cons: Less sophisticated than other methods 