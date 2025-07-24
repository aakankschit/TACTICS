Configuration System
====================

The TACTICS package uses Pydantic for robust configuration management with automatic validation and type checking.

Overview
--------

The configuration system provides:

* **Type Safety**: Automatic type checking and conversion
* **Validation**: Constraint validation for all parameters
* **Documentation**: Self-documenting configuration models
* **IDE Support**: Full autocompletion and error detection
* **Error Handling**: Clear error messages for invalid configurations

Configuration Models
-------------------

Standard Sampler Configuration (Greedy Selection)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.config.StandardSamplerConfig
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

Enhanced Sampler Configuration (Thermal Cycling)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.config.EnhancedSamplerConfig
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

Validation Examples
------------------

Type Validation
~~~~~~~~~~~~~~~

.. code-block:: python

    from TACTICS.thompson_sampling import StandardSamplerConfig
    from pydantic import ValidationError

    # Valid configuration
    config = StandardSamplerConfig(
        sampler_type="standard",
        ts_mode="maximize",
        evaluator_class_name="DBEvaluator",
        evaluator_arg="scores.csv",
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        num_ts_iterations=100,
        reagent_file_list=["aldehydes.smi", "amines.smi"],
        num_warmup_trials=3
    )

    # Invalid type raises ValidationError
    try:
        config = StandardSamplerConfig(
            sampler_type="standard",
            ts_mode="maximize",
            num_ts_iterations="not_an_integer",  # ❌ Invalid type
            # ...
        )
    except ValidationError as e:
        print(f"Type error: {e}")

Constraint Validation
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from TACTICS.thompson_sampling import EnhancedSamplerConfig

    # Invalid constraints raise ValidationError
    try:
        config = EnhancedSamplerConfig(
            sampler_type="enhanced",
            processes=-1,  # ❌ Must be > 0
            percent_of_library=1.5,  # ❌ Must be ≤ 1
            minimum_no_of_compounds_per_core=0,  # ❌ Must be > 0
            # ...
        )
    except ValidationError as e:
        print(f"Constraint error: {e}")

Enumeration Validation
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # Invalid enumeration values raise ValidationError
    try:
        config = StandardSamplerConfig(
            sampler_type="invalid",  # ❌ Must be "standard"
            ts_mode="invalid_mode",  # ❌ Must be valid mode
            # ...
        )
    except ValidationError as e:
        print(f"Enumeration error: {e}")

Usage Patterns
--------------

Direct Configuration
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from TACTICS.thompson_sampling import StandardSamplerConfig, run_ts

    # Create configuration directly
    config = StandardSamplerConfig(
        sampler_type="standard",
        ts_mode="maximize",
        evaluator_class_name="DBEvaluator",
        evaluator_arg="scores.csv",
        reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]",
        num_ts_iterations=100,
        reagent_file_list=["aldehydes.smi", "amines.smi"],
        num_warmup_trials=3
    )

    # Run with configuration
    results = run_ts(config)

JSON Configuration
~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import json
    from TACTICS.thompson_sampling import StandardSamplerConfig

    # Load from JSON file
    with open("config.json", "r") as f:
        data = json.load(f)

    # Parse with validation
    config = StandardSamplerConfig.parse_obj(data)

    # Run with configuration
    results = run_ts(config)

Example JSON configuration:

.. code-block:: json

    {
        "sampler_type": "standard",
        "ts_mode": "maximize",
        "evaluator_class_name": "DBEvaluator",
        "evaluator_arg": "scores.csv",
        "reaction_smarts": "[C:1]=[O:2]>>[C:1][O:2]",
        "num_ts_iterations": 100,
        "reagent_file_list": ["aldehydes.smi", "amines.smi"],
        "num_warmup_trials": 3,
        "results_filename": "results.csv"
    }

Configuration Best Practices
---------------------------

1. **Use Type Hints**: Always specify types for better IDE support
2. **Validate Early**: Catch configuration errors before running expensive computations
3. **Document Parameters**: Use docstrings to explain parameter purposes
4. **Provide Defaults**: Use sensible defaults for optional parameters
5. **Test Configurations**: Include configuration tests in your test suite

Error Handling
--------------

The configuration system provides detailed error messages:

.. code-block:: python

    try:
        config = StandardSamplerConfig(
            sampler_type="invalid",
            ts_mode="maximize",
            # ... other required fields
        )
    except ValidationError as e:
        print("Configuration errors:")
        for error in e.errors():
            print(f"  {error['loc']}: {error['msg']}")
            print(f"    Input: {error['input']}")
            print(f"    Type: {error['type']}")

This will output detailed error information including:
* Field location in the configuration
* Error message
* Invalid input value
* Expected type 