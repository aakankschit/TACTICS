.. TACTICS documentation master file, created by
   sphinx-quickstart on Thu Apr 24 15:14:03 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to TACTICS documentation!
=================================

TACTICS (Thompson Sampling-Assisted Chemical Targeting and Iterative Compound Selection for Drug Discovery) is a comprehensive Python package for Thompson Sampling-based optimization of chemical combinatorial libraries, featuring modern Pydantic configuration and improved package organization.

Key Features
------------

* **Modern Configuration**: Pydantic-based configuration with type safety and validation
* **Two Thompson Sampling Strategies**: Standard (greedy selection) and Enhanced (thermal cycling)
* **Library Enumeration**: Efficient generation of combinatorial reaction products
* **Library Analysis**: Comprehensive analysis and visualization tools
* **Advanced Search Methods**: Including DisallowTracker-enhanced strategies for guaranteed uniqueness
* **Clean Package Structure**: Well-organized modules with clear separation of concerns

Quick Start
-----------

.. code-block:: python

    from TACTICS.thompson_sampling import StandardSamplerConfig, run_ts

    # Create configuration using Pydantic models
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

    # Run Thompson Sampling
    results_df = run_ts(config)

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   configuration
   thompson_sampling
   library_enumeration
   library_analysis

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

