API Reference
=============

This page provides detailed API documentation for all major components of TACTICS.

Core Components
---------------

Thompson Sampler
~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.core.sampler.ThompsonSampler
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __del__

Reagent
~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.core.reagent.Reagent
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

Selection Strategies
--------------------

Base Strategy
~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.strategies.base_strategy.SelectionStrategy
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

Greedy Selection
~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.strategies.greedy_selection.GreedySelection
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

Roulette Wheel Selection
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.strategies.roulette_wheel.RouletteWheelSelection
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

   **Key Methods:**

   .. automethod:: select_reagent
   .. automethod:: select_batch
   .. automethod:: update_temperature
   .. automethod:: rotate_component
   .. automethod:: reset_temperature

UCB Selection
~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.strategies.ucb_selection.UCBSelection
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

Epsilon-Greedy Selection
~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.strategies.epsilon_greedy.EpsilonGreedySelection
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

Bayes-UCB Selection
~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.strategies.bayes_ucb_selection.BayesUCBSelection
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

   Adaptive Bayes-UCB selection strategy with thermal cycling:

   * Uses Student-t quantiles for Bayesian confidence bounds
   * Adaptive percentile parameters (p_high, p_low) for exploration
   * Thermal cycling for escaping local optima
   * Requires scipy for quantile calculations

   **Key Methods:**

   .. automethod:: select_reagent
   .. automethod:: select_batch
   .. automethod:: update_percentiles
   .. automethod:: rotate_component

Warmup Strategies
-----------------

Base Warmup Strategy
~~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.warmup.base.WarmupStrategy
   :members:
   :undoc-members:
   :show-inheritance:

   **Key Methods:**

   .. automethod:: generate_warmup_combinations
   .. automethod:: get_expected_evaluations
   .. automethod:: get_name
   .. automethod:: get_description

Balanced Warmup (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.warmup.balanced.BalancedWarmup
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

   Recommended warmup strategy with exactly K observations per reagent:

   * Guarantees uniform coverage across all reagents
   * Uses James-Stein shrinkage for per-reagent variance estimation
   * Stratified partner selection for better prior initialization

Standard Warmup
~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.warmup.standard.StandardWarmup
   :members:
   :undoc-members:
   :show-inheritance:

Enhanced Warmup (Legacy)
~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.warmup.enhanced.EnhancedWarmup
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

   Legacy warmup with stochastic parallel pairing:

   * Used in original RWS algorithm (Zhao et al. 2025)
   * Parallel pairing of reagents across components with shuffling
   * Required for legacy_rws_maximize and legacy_rws_minimize presets

Evaluators
----------

Base Evaluator
~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.core.evaluators.Evaluator
   :members:
   :undoc-members:
   :show-inheritance:

   **Abstract Methods:**

   .. automethod:: evaluate
   .. autoproperty:: counter

Lookup Evaluator
~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.core.evaluators.LookupEvaluator
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

   Simple evaluation class that looks up values from a CSV file.
   This is primarily used for testing and benchmarking.

   **Parameters:**

   * **input_dictionary** (dict): Dictionary with 'ref_filename' key pointing to CSV file

   **CSV Format:**

   * Must contain 'Product_Code' and 'Scores' columns
   * Product_Code should match product names from reaction

Database Evaluator
~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.core.evaluators.DBEvaluator
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

   Evaluator that looks up values from a SQLite database.
   More efficient than LookupEvaluator for large datasets.

   **Parameters:**

   * **input_dictionary** (dict): Dictionary with keys:
      - 'db_filename': Path to SQLite database
      - 'db_prefix': Prefix for database keys

Fingerprint Evaluator
~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.core.evaluators.FPEvaluator
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

   Evaluator that calculates Morgan fingerprint Tanimoto similarity to a reference molecule.

   **Parameters:**

   * **input_dict** (dict): Dictionary with 'query_smiles' key

Molecular Weight Evaluator
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.core.evaluators.MWEvaluator
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

   Simple evaluator that returns molecular weight.
   Primarily used for development and testing.

ROCS Evaluator
~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.core.evaluators.ROCSEvaluator
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

   Evaluator that performs 3D shape-based alignment using OpenEye ROCS.

   **Requires:** OpenEye Toolkit license

   **Parameters:**

   * **input_dict** (dict): Dictionary with 'query_molfile' key pointing to reference structure

FRED Evaluator
~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.core.evaluators.FredEvaluator
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

   Evaluator that performs molecular docking using OpenEye FRED.

   **Requires:** OpenEye Toolkit license

   **Parameters:**

   * **input_dict** (dict): Dictionary with 'design_unit_file' key pointing to .oedu file

ML Classifier Evaluator
~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.core.evaluators.MLClassifierEvaluator
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

   Evaluator that uses a trained scikit-learn classifier to score compounds.

   **Parameters:**

   * **input_dict** (dict): Dictionary with 'model_filename' key pointing to pickled model

Utility Functions
-----------------

Reagent Utilities
~~~~~~~~~~~~~~~~~

.. autofunction:: TACTICS.thompson_sampling.utils.ts_utils.read_reagents

.. autofunction:: TACTICS.thompson_sampling.utils.ts_utils.create_reagents

Logging Utilities
~~~~~~~~~~~~~~~~~

.. autofunction:: TACTICS.thompson_sampling.utils.ts_logger.get_logger

Parallel Evaluation
-------------------

Parallel Evaluator
~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.core.parallel_evaluator.ParallelEvaluator
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

   Handles parallel evaluation of reagent combinations using multiprocessing.

   **Key Methods:**

   .. automethod:: evaluate_batch
   .. automethod:: close

Legacy Components
-----------------

Disallow Tracker
~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.thompson_sampling.legacy.disallow_tracker.DisallowTracker
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

   Tracks and prevents selection of disallowed reagent combinations.

   **Key Methods:**

   .. automethod:: get_disallowed_selection_mask
   .. automethod:: update
   .. automethod:: retire_one_synthon
