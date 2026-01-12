Configuration System
====================

The TACTICS package uses Pydantic v2 for robust configuration management with automatic 
validation and type checking. The configuration system uses ``SynthesisPipeline`` as the
single source of truth for reactions and reagents.

Configuration Hierarchy
-----------------------

.. graphviz::

    digraph Configuration {
        rankdir=TB;
        node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=10];
        edge [fontname="Helvetica", fontsize=9];
        nodesep=0.2;
        ranksep=0.4;

        // Main config at top
        ThompsonSamplingConfig [label="ThompsonSamplingConfig", fillcolor="#FFD700"];

        // Strategy section header
        strategy_label [label="strategy_config", shape=plaintext, fontname="Helvetica Bold"];
        
        // Strategy configs - stacked vertically
        GreedyConfig [label="GreedyConfig", fillcolor="#FFB6C1"];
        RouletteWheelConfig [label="RouletteWheelConfig", fillcolor="#FFB6C1"];
        UCBConfig [label="UCBConfig", fillcolor="#FFB6C1"];
        EpsilonGreedyConfig [label="EpsilonGreedyConfig", fillcolor="#FFB6C1"];
        BayesUCBConfig [label="BayesUCBConfig", fillcolor="#FFB6C1"];

        // Warmup section header
        warmup_label [label="warmup_config", shape=plaintext, fontname="Helvetica Bold"];
        
        // Warmup configs - stacked vertically
        BalancedWarmupConfig [label="BalancedWarmupConfig", fillcolor="#00CED1"];
        StandardWarmupConfig [label="StandardWarmupConfig", fillcolor="#AFEEEE"];
        EnhancedWarmupConfig [label="EnhancedWarmupConfig", fillcolor="#AFEEEE"];

        // Evaluator section header
        evaluator_label [label="evaluator_config", shape=plaintext, fontname="Helvetica Bold"];
        
        // Evaluator configs - stacked vertically
        LookupEvaluatorConfig [label="LookupEvaluatorConfig", fillcolor="#DDA0DD"];
        DBEvaluatorConfig [label="DBEvaluatorConfig", fillcolor="#DDA0DD"];
        FPEvaluatorConfig [label="FPEvaluatorConfig", fillcolor="#DDA0DD"];
        MWEvaluatorConfig [label="MWEvaluatorConfig", fillcolor="#DDA0DD"];
        ROCSEvaluatorConfig [label="ROCSEvaluatorConfig", fillcolor="#D8BFD8"];
        FredEvaluatorConfig [label="FredEvaluatorConfig", fillcolor="#D8BFD8"];
        MLClassifierConfig [label="MLClassifierEvaluatorConfig", fillcolor="#D8BFD8"];

        // Force vertical layout with invisible edges
        ThompsonSamplingConfig -> strategy_label [style=invis];
        strategy_label -> GreedyConfig [style=invis];
        GreedyConfig -> RouletteWheelConfig [style=invis];
        RouletteWheelConfig -> UCBConfig [style=invis];
        UCBConfig -> EpsilonGreedyConfig [style=invis];
        EpsilonGreedyConfig -> BayesUCBConfig [style=invis];
        BayesUCBConfig -> warmup_label [style=invis];
        warmup_label -> BalancedWarmupConfig [style=invis];
        BalancedWarmupConfig -> StandardWarmupConfig [style=invis];
        StandardWarmupConfig -> EnhancedWarmupConfig [style=invis];
        EnhancedWarmupConfig -> evaluator_label [style=invis];
        evaluator_label -> LookupEvaluatorConfig [style=invis];
        LookupEvaluatorConfig -> DBEvaluatorConfig [style=invis];
        DBEvaluatorConfig -> FPEvaluatorConfig [style=invis];
        FPEvaluatorConfig -> MWEvaluatorConfig [style=invis];
        MWEvaluatorConfig -> ROCSEvaluatorConfig [style=invis];
        ROCSEvaluatorConfig -> FredEvaluatorConfig [style=invis];
        FredEvaluatorConfig -> MLClassifierConfig [style=invis];

        // Visible connections from main config to section headers
        ThompsonSamplingConfig -> strategy_label [style=bold];
        ThompsonSamplingConfig -> warmup_label [style=bold];
        ThompsonSamplingConfig -> evaluator_label [style=bold];

        // Visible connections from headers to first option (one-of relationship)
        strategy_label -> GreedyConfig [style=dashed, label="one of"];
        strategy_label -> RouletteWheelConfig [style=dashed];
        strategy_label -> UCBConfig [style=dashed];
        strategy_label -> EpsilonGreedyConfig [style=dashed];
        strategy_label -> BayesUCBConfig [style=dashed];

        warmup_label -> BalancedWarmupConfig [style=dashed, label="one of"];
        warmup_label -> StandardWarmupConfig [style=dashed];
        warmup_label -> EnhancedWarmupConfig [style=dashed];

        evaluator_label -> LookupEvaluatorConfig [style=dashed, label="one of"];
        evaluator_label -> DBEvaluatorConfig [style=dashed];
        evaluator_label -> FPEvaluatorConfig [style=dashed];
        evaluator_label -> MWEvaluatorConfig [style=dashed];
        evaluator_label -> ROCSEvaluatorConfig [style=dashed];
        evaluator_label -> FredEvaluatorConfig [style=dashed];
        evaluator_label -> MLClassifierConfig [style=dashed];
    }

Quick Start
-----------

**Modern approach (recommended):**

.. code-block:: python
   :caption: Using nested Pydantic configs with SynthesisPipeline

   from TACTICS.library_enumeration import SynthesisPipeline
   from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef
   from TACTICS.thompson_sampling.config import ThompsonSamplingConfig
   from TACTICS.thompson_sampling.strategies.config import RouletteWheelConfig
   from TACTICS.thompson_sampling.warmup.config import BalancedWarmupConfig
   from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

   # Create synthesis pipeline (single source of truth)
   rxn_config = ReactionConfig(
       reactions=[ReactionDef(
           reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
           step_index=0
       )],
       reagent_file_list=["acids.smi", "amines.smi"]
   )
   pipeline = SynthesisPipeline(rxn_config)

   # Create Thompson Sampling config
   config = ThompsonSamplingConfig(
       synthesis_pipeline=pipeline,
       num_ts_iterations=1000,
       strategy_config=RouletteWheelConfig(mode="maximize", alpha=0.1),
       warmup_config=BalancedWarmupConfig(observations_per_reagent=3),
       evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv"),
       batch_size=10
   )

**Using presets:**

.. code-block:: python
   :caption: Simplified preset-based configuration

   from TACTICS.library_enumeration import SynthesisPipeline
   from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef
   from TACTICS.thompson_sampling import get_preset
   from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

   # Create synthesis pipeline
   rxn_config = ReactionConfig(
       reactions=[ReactionDef(
           reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
           step_index=0
       )],
       reagent_file_list=["acids.smi", "amines.smi"]
   )
   pipeline = SynthesisPipeline(rxn_config)

   # Get preset configuration
   config = get_preset(
       "fast_exploration",
       synthesis_pipeline=pipeline,
       evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv"),
       num_iterations=1000
   )


.. _thompson-sampling-config:

ThompsonSamplingConfig
----------------------

.. rst-class:: class-core

The main configuration class for Thompson Sampling optimization.

.. admonition:: Dependencies
   :class: dependencies

   Accepts these nested config objects:

   - :ref:`SynthesisPipeline <synthesis-pipeline>` - via ``synthesis_pipeline`` (required)
   - :ref:`Strategy Configs <strategy-configs>` - via ``strategy_config``
   - :ref:`Warmup Configs <warmup-configs>` - via ``warmup_config``
   - :ref:`Evaluator Configs <evaluator-configs>` - via ``evaluator_config``

.. list-table:: Core Parameters
   :header-rows: 1
   :widths: 24 18 8 50

   * - Parameter
     - Type
     - Required
     - Description
   * - ``synthesis_pipeline``
     - ``SynthesisPipeline``
     - Yes
     - Single source of truth for reactions and reagents.
   * - ``num_ts_iterations``
     - ``int``
     - Yes
     - Maximum sampling cycles.
   * - ``num_warmup_trials``
     - ``int``
     - No
     - Warmup trials per reagent. Default: 3.

.. list-table:: Component Configs (Modern)
   :header-rows: 1
   :widths: 24 22 8 46

   * - Parameter
     - Type
     - Required
     - Description
   * - ``strategy_config``
     - ``StrategyConfig``
     - Yes*
     - Selection strategy configuration.
   * - ``warmup_config``
     - ``WarmupConfig``
     - No
     - Warmup strategy. Default: BalancedWarmupConfig.
   * - ``evaluator_config``
     - ``EvaluatorConfig``
     - Yes*
     - Evaluator configuration.

.. list-table:: Batch & Performance
   :header-rows: 1
   :widths: 22 15 8 55

   * - Parameter
     - Type
     - Required
     - Description
   * - ``batch_size``
     - ``int``
     - No
     - Compounds to sample per cycle. Default: 1.
   * - ``processes``
     - ``int``
     - No
     - CPU cores for parallel evaluation. Default: 1.
   * - ``min_cpds_per_core``
     - ``int``
     - No
     - Min compounds per core before batch evaluation. Default: 10.
   * - ``max_resamples``
     - ``int``
     - No
     - Stop after this many consecutive duplicates.

.. list-table:: Output & Advanced
   :header-rows: 1
   :widths: 28 15 8 49

   * - Parameter
     - Type
     - Required
     - Description
   * - ``results_filename``
     - ``str``
     - No
     - Output CSV path. Default: ``"results.csv"``.
   * - ``log_filename``
     - ``str``
     - No
     - Log file path.
   * - ``hide_progress``
     - ``bool``
     - No
     - Hide progress bars. Default: False.
   * - ``use_boltzmann_weighting``
     - ``bool``
     - No
     - Legacy RWS Boltzmann updates. Default: False.
   * - ``auto_detect_smarts_compatibility``
     - ``bool``
     - No
     - Auto-detect reagent SMARTS compatibility. Default: False.
   * - ``deprotect_for_compatibility``
     - ``bool``
     - No
     - Apply deprotection during detection. Default: False.
   * - ``desalt_for_compatibility``
     - ``bool``
     - No
     - Apply desalting during detection. Default: False.


.. _strategy-configs:

Strategy Configurations
-----------------------

All strategy configs are used with ``ThompsonSamplingConfig.strategy_config``.

.. _greedy-config:

GreedyConfig
~~~~~~~~~~~~

.. rst-class:: class-config

Configuration for greedy (argmax/argmin) selection.

**Creates:** :ref:`GreedySelection <greedy-selection>`

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 18 10 52

   * - Parameter
     - Type
     - Required
     - Description
   * - ``strategy_type``
     - ``Literal["greedy"]``
     - Auto
     - Set automatically.
   * - ``mode``
     - ``str``
     - No
     - ``"maximize"`` or ``"minimize"``. Default: ``"maximize"``.

.. _roulette-wheel-config:

RouletteWheelConfig
~~~~~~~~~~~~~~~~~~~

.. rst-class:: class-config

Configuration for roulette wheel selection with thermal cycling and CATS.

**Creates:** :ref:`RouletteWheelSelection <roulette-wheel-selection>`

.. list-table:: Parameters
   :header-rows: 1
   :widths: 24 15 8 53

   * - Parameter
     - Type
     - Required
     - Description
   * - ``strategy_type``
     - ``Literal``
     - Auto
     - Set automatically.
   * - ``mode``
     - ``str``
     - No
     - ``"maximize"``, ``"minimize"``, or Boltzmann variants.
   * - ``alpha``
     - ``float``
     - No
     - Base temperature for heated component. Default: 0.1.
   * - ``beta``
     - ``float``
     - No
     - Base temperature for cooled components. Default: 0.05.
   * - ``exploration_phase_end``
     - ``float``
     - No
     - Fraction before CATS starts [0, 1]. Default: 0.20.
   * - ``transition_phase_end``
     - ``float``
     - No
     - Fraction when CATS fully applied [0, 1]. Default: 0.60.
   * - ``min_observations``
     - ``int``
     - No
     - Min observations for criticality trust. Default: 5.

.. _ucb-config:

UCBConfig
~~~~~~~~~

.. rst-class:: class-config

Configuration for Upper Confidence Bound selection.

**Creates:** :ref:`UCBSelection <ucb-selection>`

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 18 10 52

   * - Parameter
     - Type
     - Required
     - Description
   * - ``strategy_type``
     - ``Literal["ucb"]``
     - Auto
     - Set automatically.
   * - ``mode``
     - ``str``
     - No
     - ``"maximize"`` or ``"minimize"``. Default: ``"maximize"``.
   * - ``c``
     - ``float``
     - No
     - Exploration parameter. Higher = more exploration. Default: 2.0.

.. _epsilon-greedy-config:

EpsilonGreedyConfig
~~~~~~~~~~~~~~~~~~~

.. rst-class:: class-config

Configuration for epsilon-greedy selection with decay.

**Creates:** :ref:`EpsilonGreedySelection <epsilon-greedy-selection>`

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
   * - ``strategy_type``
     - ``Literal``
     - Auto
     - Set automatically.
   * - ``mode``
     - ``str``
     - No
     - ``"maximize"`` or ``"minimize"``. Default: ``"maximize"``.
   * - ``epsilon``
     - ``float``
     - No
     - Initial exploration probability [0, 1]. Default: 0.1.
   * - ``decay``
     - ``float``
     - No
     - Decay rate per iteration (0, 1]. Default: 0.995.

.. _bayes-ucb-config:

BayesUCBConfig
~~~~~~~~~~~~~~

.. rst-class:: class-config

Configuration for Bayesian UCB with CATS integration.

**Creates:** :ref:`BayesUCBSelection <bayes-ucb-selection>`

.. list-table:: Parameters
   :header-rows: 1
   :widths: 24 15 8 53

   * - Parameter
     - Type
     - Required
     - Description
   * - ``strategy_type``
     - ``Literal``
     - Auto
     - Set automatically.
   * - ``mode``
     - ``str``
     - No
     - ``"maximize"`` or ``"minimize"``. Default: ``"maximize"``.
   * - ``initial_p_high``
     - ``float``
     - No
     - Percentile for heated component [0.5, 0.999]. Default: 0.90.
   * - ``initial_p_low``
     - ``float``
     - No
     - Percentile for cooled components [0.5, 0.999]. Default: 0.60.
   * - ``exploration_phase_end``
     - ``float``
     - No
     - Fraction before CATS starts. Default: 0.20.
   * - ``transition_phase_end``
     - ``float``
     - No
     - Fraction when CATS fully applied. Default: 0.60.
   * - ``min_observations``
     - ``int``
     - No
     - Min observations for criticality. Default: 5.


.. _warmup-configs:

Warmup Configurations
---------------------

All warmup configs are used with ``ThompsonSamplingConfig.warmup_config``.

.. _balanced-warmup-config:

BalancedWarmupConfig (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. rst-class:: class-config

Configuration for balanced warmup with per-reagent variance estimation.

**Creates:** :ref:`BalancedWarmup <balanced-warmup>`

.. list-table:: Parameters
   :header-rows: 1
   :widths: 26 15 8 51

   * - Parameter
     - Type
     - Required
     - Description
   * - ``warmup_type``
     - ``Literal``
     - Auto
     - Set automatically.
   * - ``observations_per_reagent``
     - ``int``
     - No
     - Observations per reagent. Default: 3.
   * - ``use_per_reagent_variance``
     - ``bool``
     - No
     - Use per-reagent variance estimation. Default: True.
   * - ``shrinkage_strength``
     - ``float``
     - No
     - James-Stein shrinkage strength. Default: 3.0.
   * - ``seed``
     - ``int``
     - No
     - Random seed for reproducibility.

.. _standard-warmup-config:

StandardWarmupConfig
~~~~~~~~~~~~~~~~~~~~

.. rst-class:: class-config

Configuration for standard warmup with random partners.

**Creates:** :ref:`StandardWarmup <standard-warmup>`

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 18 10 52

   * - Parameter
     - Type
     - Required
     - Description
   * - ``warmup_type``
     - ``Literal``
     - Auto
     - Set automatically.
   * - ``seed``
     - ``int``
     - No
     - Random seed for reproducibility.

.. _enhanced-warmup-config:

EnhancedWarmupConfig
~~~~~~~~~~~~~~~~~~~~

.. rst-class:: class-config

Configuration for legacy enhanced warmup (stochastic parallel pairing).

**Creates:** :ref:`EnhancedWarmup <enhanced-warmup>`

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 18 10 52

   * - Parameter
     - Type
     - Required
     - Description
   * - ``warmup_type``
     - ``Literal``
     - Auto
     - Set automatically.
   * - ``seed``
     - ``int``
     - No
     - Random seed for reproducibility.


.. _evaluator-configs:

Evaluator Configurations
------------------------

All evaluator configs are used with ``ThompsonSamplingConfig.evaluator_config``.

.. _lookup-evaluator-config:

LookupEvaluatorConfig
~~~~~~~~~~~~~~~~~~~~~

.. rst-class:: class-config

Configuration for CSV-based score lookup.

**Creates:** :ref:`LookupEvaluator <lookup-evaluator>`

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
   * - ``evaluator_type``
     - ``Literal``
     - Auto
     - Set automatically.
   * - ``ref_filename``
     - ``str``
     - Yes
     - Path to CSV file with scores.
   * - ``score_col``
     - ``str``
     - No
     - Score column name. Default: ``"Scores"``.
   * - ``compound_col``
     - ``str``
     - No
     - Compound ID column. Default: ``"Product_Code"``.

.. _db-evaluator-config:

DBEvaluatorConfig
~~~~~~~~~~~~~~~~~

.. rst-class:: class-config

Configuration for SQLite database lookup.

**Creates:** :ref:`DBEvaluator <db-evaluator>`

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
   * - ``evaluator_type``
     - ``Literal``
     - Auto
     - Set automatically.
   * - ``db_filename``
     - ``str``
     - Yes
     - Path to SQLite database.
   * - ``db_prefix``
     - ``str``
     - No
     - Key prefix for lookups. Default: ``""``.

.. _fp-evaluator-config:

FPEvaluatorConfig
~~~~~~~~~~~~~~~~~

.. rst-class:: class-config

Configuration for fingerprint similarity evaluation.

**Creates:** :ref:`FPEvaluator <fp-evaluator>`

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
   * - ``evaluator_type``
     - ``Literal``
     - Auto
     - Set automatically.
   * - ``query_smiles``
     - ``str``
     - Yes
     - Reference molecule SMILES.
   * - ``radius``
     - ``int``
     - No
     - Morgan fingerprint radius. Default: 2.
   * - ``n_bits``
     - ``int``
     - No
     - Fingerprint bit length. Default: 2048.

.. _rocs-evaluator-config:

ROCSEvaluatorConfig
~~~~~~~~~~~~~~~~~~~

.. rst-class:: class-config

Configuration for 3D shape similarity (requires OpenEye).

**Creates:** :ref:`ROCSEvaluator <rocs-evaluator>`

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
   * - ``evaluator_type``
     - ``Literal``
     - Auto
     - Set automatically.
   * - ``query_molfile``
     - ``str``
     - Yes
     - Path to reference structure (.sdf).
   * - ``max_confs``
     - ``int``
     - No
     - Max conformers to generate. Default: 50.

.. _fred-evaluator-config:

FredEvaluatorConfig
~~~~~~~~~~~~~~~~~~~

.. rst-class:: class-config

Configuration for molecular docking (requires OpenEye).

**Creates:** :ref:`FredEvaluator <fred-evaluator>`

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
   * - ``evaluator_type``
     - ``Literal``
     - Auto
     - Set automatically.
   * - ``design_unit_file``
     - ``str``
     - Yes
     - Path to receptor file (.oedu).
   * - ``max_confs``
     - ``int``
     - No
     - Max conformers to generate. Default: 100.

.. _ml-classifier-evaluator-config:

MLClassifierEvaluatorConfig
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. rst-class:: class-config

Configuration for ML model-based evaluation.

**Creates:** :ref:`MLClassifierEvaluator <ml-classifier-evaluator>`

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
   * - ``evaluator_type``
     - ``Literal``
     - Auto
     - Set automatically.
   * - ``model_filename``
     - ``str``
     - Yes
     - Path to pickled sklearn model.


.. _configuration-presets:

Configuration Presets
---------------------

.. rst-class:: class-core

TACTICS provides pre-configured setups for common use cases via ``get_preset()``.

.. admonition:: Dependencies
   :class: dependencies

   Returns fully configured :ref:`ThompsonSamplingConfig <thompson-sampling-config>` objects with preset
   strategy, warmup, and evaluator configurations.

**Available Presets:**

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Preset Name
     - Description
   * - ``fast_exploration``
     - Quick screening with epsilon-greedy (epsilon=0.2, decay=0.995).
   * - ``parallel_batch``
     - Batch processing for slow evaluators with RouletteWheel.
   * - ``conservative_exploit``
     - Hit optimization with greedy strategy.
   * - ``balanced_sampling``
     - General-purpose with UCB (c=2.0).
   * - ``diverse_coverage``
     - Maximum diversity with high-temperature RouletteWheel.
   * - ``legacy_rws_maximize``
     - Original RWS algorithm for maximize mode.
   * - ``legacy_rws_minimize``
     - Original RWS algorithm for minimize mode (docking).

**Example usage:**

.. code-block:: python

   from TACTICS.library_enumeration import SynthesisPipeline
   from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef
   from TACTICS.thompson_sampling import get_preset
   from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

   # Create synthesis pipeline
   rxn_config = ReactionConfig(
       reactions=[ReactionDef(
           reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
           step_index=0
       )],
       reagent_file_list=["acids.smi", "amines.smi"]
   )
   pipeline = SynthesisPipeline(rxn_config)

   # Fast exploration for initial screening
   config = get_preset(
       "fast_exploration",
       synthesis_pipeline=pipeline,
       evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv"),
       num_iterations=1000
   )

   # Parallel batch for docking
   from TACTICS.thompson_sampling.core.evaluator_config import FredEvaluatorConfig

   config = get_preset(
       "parallel_batch",
       synthesis_pipeline=pipeline,
       evaluator_config=FredEvaluatorConfig(design_unit_file="receptor.oedu"),
       mode="minimize",
       batch_size=100
   )


.. _factory-functions:

Factory Functions
-----------------

.. rst-class:: class-core

Factory functions create component instances from configurations.

.. admonition:: Dependencies
   :class: dependencies

   - ``create_strategy()`` - takes :ref:`Strategy Configs <strategy-configs>`
   - ``create_warmup()`` - takes :ref:`Warmup Configs <warmup-configs>`
   - ``create_evaluator()`` - takes :ref:`Evaluator Configs <evaluator-configs>`

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Function
     - Description
   * - ``create_strategy(config)``
     - Create SelectionStrategy from config.
   * - ``create_warmup(config)``
     - Create WarmupStrategy from config.
   * - ``create_evaluator(config)``
     - Create Evaluator from config.

**Example:**

.. code-block:: python

   from TACTICS.thompson_sampling.factories import (
       create_strategy,
       create_warmup,
       create_evaluator
   )
   from TACTICS.thompson_sampling.strategies.config import RouletteWheelConfig
   from TACTICS.thompson_sampling.warmup.config import BalancedWarmupConfig
   from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

   strategy = create_strategy(RouletteWheelConfig(mode="maximize", alpha=0.1))
   warmup = create_warmup(BalancedWarmupConfig(observations_per_reagent=5))
   evaluator = create_evaluator(LookupEvaluatorConfig(ref_filename="scores.csv"))


JSON/YAML Configuration
-----------------------

.. note::

   ``ThompsonSamplingConfig`` requires a ``SynthesisPipeline`` object which cannot be
   serialized directly to JSON. For reproducibility, serialize the ``ReactionConfig`` 
   and evaluator/strategy configs separately.

**Save ReactionConfig to JSON:**

.. code-block:: python

   import json
   from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef

   rxn_config = ReactionConfig(
       reactions=[ReactionDef(
           reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
           step_index=0
       )],
       reagent_file_list=["acids.smi", "amines.smi"]
   )

   with open("reaction_config.json", "w") as f:
       json.dump(rxn_config.model_dump(), f, indent=2)

**Load and create pipeline:**

.. code-block:: python

   import json
   from TACTICS.library_enumeration import SynthesisPipeline
   from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig

   with open("reaction_config.json", "r") as f:
       data = json.load(f)

   rxn_config = ReactionConfig.model_validate(data)
   pipeline = SynthesisPipeline(rxn_config)

**Example ReactionConfig JSON:**

.. code-block:: json

   {
       "reactions": [
           {
               "reaction_smarts": "[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
               "step_index": 0,
               "pattern_id": null,
               "description": null,
               "deprotections": []
           }
       ],
       "reagent_file_list": ["acids.smi", "amines.smi"],
       "step_inputs": null,
       "step_modes": null,
       "protecting_groups": null
   }


Validation Examples
-------------------

**Type validation:**

.. code-block:: python

   from TACTICS.library_enumeration import SynthesisPipeline
   from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef
   from TACTICS.thompson_sampling.config import ThompsonSamplingConfig
   from TACTICS.thompson_sampling.strategies.config import GreedyConfig
   from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig
   from pydantic import ValidationError

   # Create pipeline
   rxn_config = ReactionConfig(
       reactions=[ReactionDef(reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]", step_index=0)],
       reagent_file_list=["acids.smi", "amines.smi"]
   )
   pipeline = SynthesisPipeline(rxn_config)

   try:
       config = ThompsonSamplingConfig(
           synthesis_pipeline=pipeline,
           num_ts_iterations="not_an_integer",  # Invalid!
           strategy_config=GreedyConfig(),
           evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv")
       )
   except ValidationError as e:
       print(f"Validation error: {e}")

**Constraint validation:**

.. code-block:: python

   from TACTICS.thompson_sampling.strategies.config import RouletteWheelConfig
   from pydantic import ValidationError

   try:
       config = RouletteWheelConfig(
           mode="maximize",
           alpha=-0.1,  # Invalid: must be > 0
       )
   except ValidationError as e:
       print(f"Constraint error: {e}")


Best Practices
--------------

1. **Use SynthesisPipeline** - Single source of truth for reactions and reagents
2. **Use presets** - Start with presets and customize as needed
3. **Validate early** - Pydantic catches errors before expensive computations
4. **Save ReactionConfig** - JSON export for reproducibility
5. **Use BalancedWarmup** - Recommended for most use cases
6. **Choose evaluator wisely** - Use ``processes=1`` for fast evaluators (Lookup, DB)
