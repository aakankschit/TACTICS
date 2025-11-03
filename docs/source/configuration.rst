Configuration System
====================

The TACTICS package uses Pydantic v2 for robust configuration management with automatic validation and type checking.

Overview
--------

The configuration system provides:

* **Type Safety**: Automatic type checking and conversion
* **Validation**: Constraint validation for all parameters with detailed error messages
* **Nested Configs**: Hierarchical configuration with strategy, warmup, and evaluator configs
* **Documentation**: Self-documenting configuration models
* **IDE Support**: Full autocompletion and error detection
* **Presets**: Pre-configured setups for common use cases
* **Backward Compatibility**: Legacy string-based configs still supported

Modern Configuration Approach
------------------------------

The modern configuration system uses nested Pydantic models for each component:

.. code-block:: python

    from TACTICS.thompson_sampling import ThompsonSampler
    from TACTICS.thompson_sampling.config import ThompsonSamplingConfig
    from TACTICS.thompson_sampling.strategies.config import RouletteWheelConfig
    from TACTICS.thompson_sampling.warmup.config import StratifiedWarmupConfig
    from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

    # Create nested configuration
    config = ThompsonSamplingConfig(
        reaction_smarts="[#6:1](=[O:2])[OH].[#7:3]>>[#6:1](=[O:2])[#7:3]",
        reagent_file_list=["acids.smi", "amines.smi"],
        num_ts_iterations=1000,
        num_warmup_trials=3,
        strategy_config=RouletteWheelConfig(
            mode="maximize",
            alpha=0.1,
            beta=0.1
        ),
        warmup_config=StratifiedWarmupConfig(),
        evaluator_config=LookupEvaluatorConfig(
            ref_filename="scores.csv",
            ref_colname="Scores"
        ),
        batch_size=1,
        processes=1
    )

    # Create sampler from config
    sampler = ThompsonSampler.from_config(config)
    sampler.set_reaction(config.reaction_smarts)
    sampler.warm_up(num_warmup_trials=config.num_warmup_trials)
    results = sampler.search(num_cycles=config.num_ts_iterations)

Configuration Models
--------------------

Main Configuration
~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.config.ThompsonSamplingConfig
   :members:
   :undoc-members:
   :show-inheritance:

The main configuration class supports both modern (nested configs) and legacy (string-based) approaches:

* **Modern**: Use strategy_config, warmup_config, evaluator_config
* **Legacy**: Use selection_strategy, evaluator_class_name, evaluator_arg

Strategy Configurations
~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.strategies.config.GreedyConfig
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: TACTICS.thompson_sampling.strategies.config.RouletteWheelConfig
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: TACTICS.thompson_sampling.strategies.config.UCBConfig
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: TACTICS.thompson_sampling.strategies.config.EpsilonGreedyConfig
   :members:
   :undoc-members:
   :show-inheritance:

Warmup Configurations
~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.warmup.config.StandardWarmupConfig
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: TACTICS.thompson_sampling.warmup.config.StratifiedWarmupConfig
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: TACTICS.thompson_sampling.warmup.config.EnhancedWarmupConfig
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: TACTICS.thompson_sampling.warmup.config.LatinHypercubeWarmupConfig
   :members:
   :undoc-members:
   :show-inheritance:

Evaluator Configurations
~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.core.evaluator_config.LookupEvaluatorConfig
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: TACTICS.thompson_sampling.core.evaluator_config.DBEvaluatorConfig
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: TACTICS.thompson_sampling.core.evaluator_config.FPEvaluatorConfig
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: TACTICS.thompson_sampling.core.evaluator_config.ROCSEvaluatorConfig
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: TACTICS.thompson_sampling.core.evaluator_config.FredEvaluatorConfig
   :members:
   :undoc-members:
   :show-inheritance:

Configuration Presets
---------------------

TACTICS provides pre-configured setups for common use cases:

.. autoclass:: TACTICS.thompson_sampling.presets.ConfigPresets
   :members:
   :undoc-members:
   :show-inheritance:

.. autofunction:: TACTICS.thompson_sampling.presets.get_preset

Preset Examples
~~~~~~~~~~~~~~~

Fast Exploration
^^^^^^^^^^^^^^^^

.. code-block:: python

    from TACTICS.thompson_sampling.presets import ConfigPresets
    from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

    config = ConfigPresets.fast_exploration(
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        reagent_file_list=["acids.smi", "amines.smi"],
        evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv"),
        num_iterations=1000
    )

    sampler = ThompsonSampler.from_config(config)

Parallel Batch Processing
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    config = ConfigPresets.parallel_batch(
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        reagent_file_list=["acids.smi", "amines.smi"],
        evaluator_config=ROCSEvaluatorConfig(query_molfile="ref.sdf"),
        batch_size=100,
        processes=8
    )

Conservative Exploitation
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    config = ConfigPresets.conservative_exploit(
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        reagent_file_list=["acids.smi", "amines.smi"],
        evaluator_config=DBEvaluatorConfig(
            db_filename="scores.db",
            db_prefix="compound_"
        )
    )

Diverse Coverage
^^^^^^^^^^^^^^^^

.. code-block:: python

    config = ConfigPresets.diverse_coverage(
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        reagent_file_list=["acids.smi", "amines.smi"],
        evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv")
    )

Minimize Mode (for Docking)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    config = ConfigPresets.minimize_mode(
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        reagent_file_list=["acids.smi", "amines.smi"],
        evaluator_config=FredEvaluatorConfig(design_unit_file="receptor.oedu"),
        strategy="roulette_wheel"
    )

Using get_preset
^^^^^^^^^^^^^^^^

.. code-block:: python

    from TACTICS.thompson_sampling.presets import get_preset

    config = get_preset(
        "balanced_sampling",
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        reagent_file_list=["acids.smi", "amines.smi"],
        evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv")
    )

Factory Functions
-----------------

.. autofunction:: TACTICS.thompson_sampling.factories.create_strategy
.. autofunction:: TACTICS.thompson_sampling.factories.create_warmup
.. autofunction:: TACTICS.thompson_sampling.factories.create_evaluator

Factory Usage
~~~~~~~~~~~~~

.. code-block:: python

    from TACTICS.thompson_sampling.factories import (
        create_strategy,
        create_warmup,
        create_evaluator
    )
    from TACTICS.thompson_sampling.strategies.config import RouletteWheelConfig
    from TACTICS.thompson_sampling.warmup.config import LatinHypercubeWarmupConfig
    from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

    # Create components from configs
    strategy = create_strategy(RouletteWheelConfig(mode="maximize", alpha=0.1))
    warmup = create_warmup(LatinHypercubeWarmupConfig())
    evaluator = create_evaluator(LookupEvaluatorConfig(ref_filename="scores.csv"))

    # Use with sampler
    sampler = ThompsonSampler(
        selection_strategy=strategy,
        warmup_strategy=warmup
    )
    sampler.set_evaluator(evaluator)

Validation Examples
-------------------

Type Validation
~~~~~~~~~~~~~~~

.. code-block:: python

    from TACTICS.thompson_sampling.config import ThompsonSamplingConfig
    from TACTICS.thompson_sampling.strategies.config import GreedyConfig
    from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig
    from pydantic import ValidationError

    # Valid configuration
    config = ThompsonSamplingConfig(
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        reagent_file_list=["acids.smi", "amines.smi"],
        num_ts_iterations=100,
        strategy_config=GreedyConfig(mode="maximize"),
        evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv")
    )

    # Invalid type raises ValidationError
    try:
        config = ThompsonSamplingConfig(
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            reagent_file_list=["acids.smi"],
            num_ts_iterations="not_an_integer",  # ❌ Invalid type
            strategy_config=GreedyConfig(),
            evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv")
        )
    except ValidationError as e:
        print(f"Type error: {e}")

Constraint Validation
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from TACTICS.thompson_sampling.strategies.config import RouletteWheelConfig

    # Invalid constraints raise ValidationError
    try:
        config = RouletteWheelConfig(
            mode="maximize",
            alpha=-0.1,  # ❌ Must be > 0
            efficiency_threshold=1.5  # ❌ Must be ≤ 1
        )
    except ValidationError as e:
        print(f"Constraint error: {e}")

Nested Config Validation
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    try:
        # Cannot mix modern and legacy configs
        config = ThompsonSamplingConfig(
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            reagent_file_list=["acids.smi"],
            num_ts_iterations=100,
            # Modern
            strategy_config=GreedyConfig(mode="maximize"),
            # Legacy (conflicts!)
            evaluator_class_name="LookupEvaluator",
            evaluator_arg="{}"
        )
    except ValueError as e:
        print(f"Config conflict: {e}")

JSON/YAML Configuration
-----------------------

Save Configuration to JSON
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import json
    from TACTICS.thompson_sampling.config import ThompsonSamplingConfig
    from TACTICS.thompson_sampling.strategies.config import RouletteWheelConfig
    from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

    config = ThompsonSamplingConfig(
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        reagent_file_list=["acids.smi", "amines.smi"],
        num_ts_iterations=1000,
        strategy_config=RouletteWheelConfig(mode="maximize", alpha=0.1),
        evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv")
    )

    # Save to JSON
    with open("config.json", "w") as f:
        json.dump(config.model_dump(), f, indent=2)

Load Configuration from JSON
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import json
    from TACTICS.thompson_sampling.config import ThompsonSamplingConfig

    # Load from JSON file
    with open("config.json", "r") as f:
        data = json.load(f)

    # Parse with validation
    config = ThompsonSamplingConfig.model_validate(data)

    # Use configuration
    sampler = ThompsonSampler.from_config(config)

Example JSON Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: json

    {
        "reaction_smarts": "[#6:1](=[O:2])[OH].[#7:3]>>[#6:1](=[O:2])[#7:3]",
        "reagent_file_list": ["acids.smi", "amines.smi"],
        "num_ts_iterations": 1000,
        "num_warmup_trials": 3,
        "strategy_config": {
            "strategy_type": "roulette_wheel",
            "mode": "maximize",
            "alpha": 0.1,
            "beta": 0.1,
            "scaling": 1.0,
            "alpha_increment": 0.01,
            "beta_increment": 0.001,
            "efficiency_threshold": 0.1
        },
        "warmup_config": {
            "warmup_type": "stratified"
        },
        "evaluator_config": {
            "evaluator_type": "lookup",
            "ref_filename": "scores.csv",
            "ref_colname": "Scores"
        },
        "batch_size": 100,
        "processes": 4,
        "results_filename": "results.csv"
    }

Legacy Configuration (Backward Compatibility)
----------------------------------------------

The legacy string-based configuration is still supported for backward compatibility:

.. code-block:: python

    from TACTICS.thompson_sampling.config import ThompsonSamplingConfig

    # Legacy approach (still works)
    config = ThompsonSamplingConfig(
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        reagent_file_list=["acids.smi", "amines.smi"],
        num_ts_iterations=100,
        evaluator_class_name="LookupEvaluator",
        evaluator_arg='{"ref_filename": "scores.csv"}',
        selection_strategy="greedy",
        mode="maximize"
    )

    sampler = ThompsonSampler.from_config(config)

Configuration Best Practices
-----------------------------

1. **Use Modern Configs**: Prefer nested configs for better validation and type safety
2. **Use Presets**: Start with presets and customize as needed
3. **Validate Early**: Catch configuration errors before running expensive computations
4. **Save Configs**: Save configurations to JSON for reproducibility
5. **Document Parameters**: Use clear parameter names and add comments
6. **Test Configurations**: Include configuration tests in your test suite

Preset Selection Guide
----------------------

Choose the right preset for your use case:

**fast_exploration**
    * Quick initial screening
    * Fast evaluators (Lookup, Database)
    * Epsilon-greedy strategy
    * Single-compound mode

**parallel_batch**
    * Slow evaluators (ROCS, Fred, ML)
    * Large-scale screening
    * Multi-core systems
    * Batch processing

**conservative_exploit**
    * Hit optimization
    * Focus on best reagents
    * Well-characterized space
    * Greedy strategy

**balanced_sampling**
    * General-purpose screening
    * Theoretical guarantees
    * UCB strategy
    * Diverse exploration

**diverse_coverage**
    * Maximum diversity
    * Reagents ordered by properties
    * Latin Hypercube warmup
    * Exploration-heavy

**minimize_mode**
    * Docking scores (lower is better)
    * Binding energies
    * Any minimization problem

Error Handling
--------------

The configuration system provides detailed error messages:

.. code-block:: python

    from pydantic import ValidationError

    try:
        config = ThompsonSamplingConfig(
            reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
            reagent_file_list=["acids.smi"],
            num_ts_iterations=-100,  # Invalid: must be positive
            batch_size=0,  # Invalid: must be > 0
            processes=0,  # Invalid: must be > 0
            strategy_config=GreedyConfig(mode="invalid_mode"),  # Invalid mode
            evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv")
        )
    except ValidationError as e:
        print("Configuration errors:")
        for error in e.errors():
            print(f"  Field: {error['loc']}")
            print(f"  Error: {error['msg']}")
            print(f"  Input: {error['input']}")
            print(f"  Type: {error['type']}")

This will output detailed error information including:

* Field location in the configuration hierarchy
* Clear error message
* Invalid input value
* Error type for programmatic handling
