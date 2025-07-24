Thompson Sampling Module
========================

The Thompson Sampling module implements Thompson sampling algorithms for chemical library exploration, featuring modern Pydantic configuration and improved package organization.

The package provides two main Thompson Sampling strategies:

* **Standard Thompson Sampler**: Uses greedy selection mode for straightforward optimization
* **Enhanced Thompson Sampler**: Uses thermal cycling for better exploration/exploitation balance

Configuration Models
-------------------

.. autoclass:: PRISMS.thompson_sampling.config.StandardSamplerConfig
   :members:
   :undoc-members:
   :show-inheritance:

The StandardSamplerConfig is used for basic Thompson Sampling with greedy selection:

* **sampler_type**: Must be "standard"
* **ts_mode**: One of "maximize", "minimize", "maximize_boltzmann", "minimize_boltzmann"
* **evaluator_class_name**: Name of the evaluator class to use
* **evaluator_arg**: Argument passed to the evaluator constructor
* **reaction_smarts**: SMARTS string defining the reaction
* **num_ts_iterations**: Number of Thompson Sampling iterations
* **reagent_file_list**: List of reagent file paths
* **num_warmup_trials**: Number of warmup trials
* **results_filename**: Optional output file path
* **log_filename**: Optional log file path

.. autoclass:: PRISMS.thompson_sampling.config.EnhancedSamplerConfig
   :members:
   :undoc-members:
   :show-inheritance:

The EnhancedSamplerConfig uses thermal cycling for improved exploration:

* **sampler_type**: Must be "enhanced"
* **processes**: Number of parallel processes (must be > 0)
* **scaling**: Scaling factor for scores
* **percent_of_library**: Fraction of library to search (0 < x ≤ 1)
* **minimum_no_of_compounds_per_core**: Minimum compounds per core (must be > 0)
* **stopping_criteria**: Stopping criteria threshold (must be > 0)

Main Interface
-------------

.. autofunction:: PRISMS.thompson_sampling.main.run_ts
   :noindex:

Core Samplers
-------------

.. autoclass:: PRISMS.thompson_sampling.core.standard_sampler.StandardThompsonSampler
   :members:
   :undoc-members:
   :show-inheritance:

   .. automethod:: __init__
   .. automethod:: warm_up
   .. automethod:: search
   .. automethod:: evaluate

   The StandardThompsonSampler implements basic Thompson Sampling with greedy selection:
   
   * Simple greedy selection based on sampled scores
   * Direct optimization toward best-performing reagents
   * Suitable for straightforward optimization problems
   * Lower computational overhead

.. autoclass:: PRISMS.thompson_sampling.core.enhanced_sampler.EnhancedThompsonSampler
   :members:
   :undoc-members:
   :show-inheritance:

   .. automethod:: __init__
   .. automethod:: warm_up
   .. automethod:: search
   .. automethod:: evaluate

   The EnhancedThompsonSampler implements thermal cycling for better exploration:
   
   * Thermal cycling for exploration/exploitation balance
   * Parallel processing support
   * Adaptive temperature control
   * Batch processing optimization
   * Enhanced progress tracking
   * Better for complex, multi-modal search landscapes

Evaluator Classes
----------------

.. autoclass:: PRISMS.thompson_sampling.core.evaluators.LookupEvaluator
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: PRISMS.thompson_sampling.core.evaluators.ROCSEvaluator
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: PRISMS.thompson_sampling.core.evaluators.DBEvaluator
   :members:
   :undoc-members:
   :show-inheritance:

Utility Functions
----------------

.. autofunction:: PRISMS.thompson_sampling.utils.ts_utils.read_reagents
.. autofunction:: PRISMS.thompson_sampling.utils.ts_utils.create_reagents
.. autofunction:: PRISMS.thompson_sampling.utils.ts_logger.get_logger

Legacy Interface
---------------

For backward compatibility, the legacy interface is still available:

.. autofunction:: PRISMS.thompson_sampling.legacy.ts_main.run_ts
   :noindex:

.. autoclass:: PRISMS.thompson_sampling.legacy.standard_thompson_sampling.StandardThompsonSampler
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: PRISMS.thompson_sampling.legacy.enhanced_thompson_sampling.EnhancedThompsonSampler
   :members:
   :undoc-members:
   :show-inheritance:

Usage Examples
-------------

Standard Thompson Sampling (Greedy Selection):

.. code-block:: python

    from PRISMS.thompson_sampling import StandardSamplerConfig, run_ts

    # Create configuration for greedy selection
    config = StandardSamplerConfig(
        sampler_type="standard",
        ts_mode="maximize",
        evaluator_class_name="DBEvaluator",
        evaluator_arg="scores.csv",
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        num_ts_iterations=100,
        reagent_file_list=["aldehydes.smi", "amines.smi"],
        num_warmup_trials=3,
        results_filename="results.csv"
    )

    # Run Thompson Sampling with greedy selection
    results_df = run_ts(config)

Enhanced Thompson Sampling (Thermal Cycling):

.. code-block:: python

    from PRISMS.thompson_sampling import EnhancedSamplerConfig

    # Create configuration for thermal cycling
    config = EnhancedSamplerConfig(
        sampler_type="enhanced",
        processes=4,
        scaling=1.0,
        percent_of_library=0.1,
        minimum_no_of_compounds_per_core=10,
        stopping_criteria=1000,
        evaluator_class_name="DBEvaluator",
        evaluator_arg="scores.csv",
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        num_ts_iterations=100,
        reagent_file_list=["aldehydes.smi", "amines.smi"],
        num_warmup_trials=3
    )

    # Run Thompson Sampling with thermal cycling
    results_df = run_ts(config)

Configuration Validation:

.. code-block:: python

    from pydantic import ValidationError

    try:
        config = StandardSamplerConfig(
            sampler_type="invalid",  # ❌ ValidationError
            ts_mode="maximize",
            # ...
        )
    except ValidationError as e:
        print(f"Configuration error: {e}")

Strategy Selection Guide
-----------------------

**Use Standard Thompson Sampling (Greedy Selection) when:**
* You have a straightforward optimization problem
* Computational resources are limited
* You want direct, greedy optimization
* The search landscape is relatively simple

**Use Enhanced Thompson Sampling (Thermal Cycling) when:**
* You have complex, multi-modal search landscapes
* You want better exploration/exploitation balance
* You have access to parallel processing
* The search space has many local optima 