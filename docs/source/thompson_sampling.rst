Thompson Sampling
=================

The Thompson Sampling module implements a unified, flexible Thompson Sampling framework 
for chemical library exploration. It provides pluggable selection strategies, warmup 
approaches, and evaluators for efficiently screening ultra-large combinatorial libraries.

The module follows a composition-based architecture where the core ``ThompsonSampler`` 
class accepts pluggable components:

- **Selection Strategies** - How to choose reagents during search
- **Warmup Strategies** - How to initialize priors before search
- **Evaluators** - How to score generated compounds

Module Architecture
-------------------

.. graphviz::

    digraph ThompsonSamplingModule {
        rankdir=TB;
        node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=10];
        edge [fontname="Helvetica", fontsize=9];
        nodesep=0.3;
        ranksep=0.5;

        // Core at top
        subgraph cluster_core {
            label="Core";
            style=filled;
            color="#FFF8DC";
            
            ThompsonSampler [label="ThompsonSampler", fillcolor="#FFD700"];
            
            subgraph cluster_core_helpers {
                label="";
                style=invis;
                rank=same;
                Reagent [label="Reagent", fillcolor="#FFFACD"];
                ParallelEvaluator [label="ParallelEvaluator", fillcolor="#FFFACD"];
            }
        }

        // Strategies - vertical list
        subgraph cluster_strategies {
            label="Selection Strategies";
            style=filled;
            color="#FFE4E1";

            BaseStrategy [label="SelectionStrategy (ABC)", fillcolor="white"];
            Greedy [label="GreedySelection", fillcolor="#FFB6C1"];
            RouletteWheel [label="RouletteWheelSelection", fillcolor="#FFB6C1"];
            UCB [label="UCBSelection", fillcolor="#FFB6C1"];
            EpsilonGreedy [label="EpsilonGreedySelection", fillcolor="#FFB6C1"];
            BayesUCB [label="BayesUCBSelection", fillcolor="#FFB6C1"];
            
            // Force vertical
            BaseStrategy -> Greedy [style=invis];
            Greedy -> RouletteWheel [style=invis];
            RouletteWheel -> UCB [style=invis];
            UCB -> EpsilonGreedy [style=invis];
            EpsilonGreedy -> BayesUCB [style=invis];
        }

        // Warmup - vertical list
        subgraph cluster_warmup {
            label="Warmup Strategies";
            style=filled;
            color="#E0FFFF";

            BaseWarmup [label="WarmupStrategy (ABC)", fillcolor="white"];
            Balanced [label="BalancedWarmup", fillcolor="#00CED1"];
            Standard [label="StandardWarmup", fillcolor="#AFEEEE"];
            Enhanced [label="EnhancedWarmup", fillcolor="#AFEEEE"];
            
            // Force vertical
            BaseWarmup -> Balanced [style=invis];
            Balanced -> Standard [style=invis];
            Standard -> Enhanced [style=invis];
        }

        // Evaluators - vertical list
        subgraph cluster_evaluators {
            label="Evaluators";
            style=filled;
            color="#E6E6FA";

            BaseEvaluator [label="Evaluator (ABC)", fillcolor="white"];
            Lookup [label="LookupEvaluator", fillcolor="#DDA0DD"];
            DB [label="DBEvaluator", fillcolor="#DDA0DD"];
            FP [label="FPEvaluator", fillcolor="#DDA0DD"];
            MW [label="MWEvaluator", fillcolor="#DDA0DD"];
            ROCS [label="ROCSEvaluator", fillcolor="#D8BFD8"];
            FRED [label="FredEvaluator", fillcolor="#D8BFD8"];
            ML [label="MLClassifierEvaluator", fillcolor="#D8BFD8"];
            
            // Force vertical
            BaseEvaluator -> Lookup [style=invis];
            Lookup -> DB [style=invis];
            DB -> FP [style=invis];
            FP -> MW [style=invis];
            MW -> ROCS [style=invis];
            ROCS -> FRED [style=invis];
            FRED -> ML [style=invis];
        }

        // Inheritance arrows (visible)
        BaseStrategy -> Greedy [style=dashed, constraint=false];
        BaseStrategy -> RouletteWheel [style=dashed, constraint=false];
        BaseStrategy -> UCB [style=dashed, constraint=false];
        BaseStrategy -> EpsilonGreedy [style=dashed, constraint=false];
        BaseStrategy -> BayesUCB [style=dashed, constraint=false];

        BaseWarmup -> Balanced [style=dashed, constraint=false];
        BaseWarmup -> Standard [style=dashed, constraint=false];
        BaseWarmup -> Enhanced [style=dashed, constraint=false];

        BaseEvaluator -> Lookup [style=dashed, constraint=false];
        BaseEvaluator -> DB [style=dashed, constraint=false];
        BaseEvaluator -> FP [style=dashed, constraint=false];
        BaseEvaluator -> MW [style=dashed, constraint=false];
        BaseEvaluator -> ROCS [style=dashed, constraint=false];
        BaseEvaluator -> FRED [style=dashed, constraint=false];
        BaseEvaluator -> ML [style=dashed, constraint=false];

        // Core connections
        ThompsonSampler -> BaseStrategy [style=bold];
        ThompsonSampler -> BaseWarmup [style=bold];
        ThompsonSampler -> BaseEvaluator [style=bold];
    }

Quick Start
-----------

**Using presets (recommended):**

.. code-block:: python
   :caption: Simplest usage with presets

   from TACTICS.library_enumeration import SynthesisPipeline
   from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef
   from TACTICS.thompson_sampling import run_ts, get_preset
   from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

   # 1. Create synthesis pipeline (single source of truth)
   rxn_config = ReactionConfig(
       reactions=[ReactionDef(
           reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
           step_index=0
       )],
       reagent_file_list=["acids.smi", "amines.smi"]
   )
   pipeline = SynthesisPipeline(rxn_config)

   # 2. Get preset configuration
   config = get_preset(
       "fast_exploration",
       synthesis_pipeline=pipeline,
       evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv"),
       mode="minimize",
       num_iterations=1000
   )

   # 3. Run and get results
   results_df = run_ts(config)
   print(results_df.sort("score").head(10))

**Direct sampler control:**

.. code-block:: python
   :caption: Manual sampler setup

   from TACTICS.library_enumeration import SynthesisPipeline
   from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef
   from TACTICS.thompson_sampling.core.sampler import ThompsonSampler
   from TACTICS.thompson_sampling.strategies import RouletteWheelSelection
   from TACTICS.thompson_sampling.warmup import BalancedWarmup
   from TACTICS.thompson_sampling.factories import create_evaluator
   from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

   # 1. Create synthesis pipeline
   rxn_config = ReactionConfig(
       reactions=[ReactionDef(
           reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
           step_index=0
       )],
       reagent_file_list=["acids.smi", "amines.smi"]
   )
   pipeline = SynthesisPipeline(rxn_config)

   # 2. Create components
   strategy = RouletteWheelSelection(mode="maximize", alpha=0.1, beta=0.05)
   warmup = BalancedWarmup(observations_per_reagent=3)
   evaluator = create_evaluator(LookupEvaluatorConfig(ref_filename="scores.csv"))

   # 3. Create sampler with pipeline
   sampler = ThompsonSampler(
       synthesis_pipeline=pipeline,
       selection_strategy=strategy,
       warmup_strategy=warmup,
       batch_size=10
   )

   # 4. Set evaluator and run
   sampler.set_evaluator(evaluator)
   warmup_df = sampler.warm_up(num_warmup_trials=3)
   results_df = sampler.search(num_cycles=1000)
   sampler.close()


.. _thompson-sampler:

ThompsonSampler
---------------

.. rst-class:: class-core

**The main class for Thompson Sampling optimization.**

The ``ThompsonSampler`` is the central orchestrator that coordinates selection strategies,
warmup strategies, and evaluators to efficiently explore combinatorial chemical libraries.

.. admonition:: Dependencies
   :class: dependencies

   Requires these components:

   - :ref:`SynthesisPipeline <synthesis-pipeline>` - single source of truth for reactions and reagents
   - :ref:`SelectionStrategy <selection-strategy>` - for reagent selection during search
   - :ref:`WarmupStrategy <warmup-strategy>` - for initializing priors (optional, defaults to StandardWarmup)
   - :ref:`Evaluator <evaluator-base>` - for scoring compounds (set via ``set_evaluator()``)

**Depends on:** :ref:`SynthesisPipeline <synthesis-pipeline>`, :ref:`SelectionStrategy <selection-strategy>`, :ref:`WarmupStrategy <warmup-strategy>`, :ref:`Evaluator <evaluator-base>`

Constructor
~~~~~~~~~~~

.. list-table:: Parameters
   :header-rows: 1
   :widths: 22 18 10 50

   * - Parameter
     - Type
     - Required
     - Description
   * - ``synthesis_pipeline``
     - ``SynthesisPipeline``
     - Yes
     - Pipeline containing reaction config and reagent files (single source of truth).
   * - ``selection_strategy``
     - ``SelectionStrategy``
     - Yes
     - Selection strategy instance (Greedy, RouletteWheel, UCB, etc.).
   * - ``warmup_strategy``
     - ``WarmupStrategy``
     - No
     - Warmup strategy. Default: StandardWarmup().
   * - ``batch_size``
     - ``int``
     - No
     - Compounds to sample per cycle. Default: 1.
   * - ``processes``
     - ``int``
     - No
     - CPU cores for parallel evaluation. Default: 1 (sequential).
   * - ``min_cpds_per_core``
     - ``int``
     - No
     - Min compounds per core before batch evaluation. Default: 10.
   * - ``max_resamples``
     - ``int``
     - No
     - Stop after this many consecutive duplicates. Default: None.
   * - ``log_filename``
     - ``str``
     - No
     - Path for log file output.
   * - ``product_library_file``
     - ``str``
     - No
     - Pre-enumerated product CSV for testing mode.
   * - ``use_boltzmann_weighting``
     - ``bool``
     - No
     - Use Boltzmann-weighted updates (legacy RWS). Default: False.

Factory Method: from_config
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a sampler from a Pydantic configuration.

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 25 10 45

   * - Parameter
     - Type
     - Required
     - Description
   * - ``config``
     - ``ThompsonSamplingConfig``
     - Yes
     - Configuration with strategy, warmup, and evaluator settings.

**Returns**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Type
     - Description
   * - ``ThompsonSampler``
     - Configured sampler ready for warmup and search.

**Example**

.. code-block:: python

   from TACTICS.library_enumeration import SynthesisPipeline
   from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef
   from TACTICS.thompson_sampling.core.sampler import ThompsonSampler
   from TACTICS.thompson_sampling.config import ThompsonSamplingConfig
   from TACTICS.thompson_sampling.strategies.config import RouletteWheelConfig
   from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

   # Create synthesis pipeline
   rxn_config = ReactionConfig(
       reactions=[ReactionDef(reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]", step_index=0)],
       reagent_file_list=["acids.smi", "amines.smi"]
   )
   pipeline = SynthesisPipeline(rxn_config)

   # Create Thompson Sampling config
   config = ThompsonSamplingConfig(
       synthesis_pipeline=pipeline,
       num_ts_iterations=1000,
       strategy_config=RouletteWheelConfig(mode="maximize"),
       evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv")
   )

   sampler = ThompsonSampler.from_config(config)

Core Methods
~~~~~~~~~~~~

warm_up
^^^^^^^

Initialize reagent posteriors with warmup evaluations.

.. list-table:: Parameters
   :header-rows: 1
   :widths: 22 15 10 53

   * - Parameter
     - Type
     - Required
     - Description
   * - ``num_warmup_trials``
     - ``int``
     - No
     - Trials per reagent. Default: 3.

**Returns**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Type
     - Description
   * - ``polars.DataFrame``
     - Warmup results with columns: ``score``, ``SMILES``, ``Name``.

search
^^^^^^

Run the main Thompson Sampling search loop.

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
   * - ``num_cycles``
     - ``int``
     - No
     - Maximum sampling cycles. Default: 100.
   * - ``max_evaluations``
     - ``int``
     - No
     - Stop after this many unique evaluations.

**Returns**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Type
     - Description
   * - ``polars.DataFrame``
     - Search results with columns: ``score``, ``SMILES``, ``Name``.

evaluate
^^^^^^^^

Evaluate a single reagent combination.

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 18 10 52

   * - Parameter
     - Type
     - Required
     - Description
   * - ``choice_list``
     - ``list[int]``
     - Yes
     - Reagent indices for each component.

**Returns**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Type
     - Description
   * - ``tuple[str, str, float]``
     - (product_smiles, product_name, score).

Setup Methods
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Method
     - Description
   * - ``set_evaluator(evaluator)``
     - Set the scoring evaluator.
   * - ``load_product_library(library_file)``
     - Load pre-enumerated products for testing.
   * - ``close()``
     - Cleanup multiprocessing resources.

.. note::

   The ``synthesis_pipeline`` is now passed to the constructor and is the single source
   of truth for reactions and reagents. The old ``read_reagents()`` and ``set_reaction()``
   methods have been removed.


Selection Strategies
--------------------

Selection strategies determine how reagents are chosen during the search phase.
All strategies implement the ``SelectionStrategy`` abstract base class.

.. _selection-strategy:

.. graphviz::

    digraph StrategySelection {
        rankdir=TB;
        node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=10];
        nodesep=0.2;
        ranksep=0.3;

        SelectReagent [label="Select Reagent", fillcolor="#FFFACD"];

        Greedy [label="Greedy (argmax/argmin)", fillcolor="#90EE90"];
        RouletteWheel [label="RouletteWheel (Thermal Cycling)", fillcolor="#ADD8E6"];
        UCB [label="UCB (Confidence Bounds)", fillcolor="#FFB6C1"];
        BayesUCB [label="BayesUCB (Student-t + CATS)", fillcolor="#E6E6FA"];
        EpsilonGreedy [label="EpsilonGreedy (Random/Greedy)", fillcolor="#FFFACD"];

        // Force vertical layout
        SelectReagent -> Greedy;
        Greedy -> RouletteWheel [style=invis];
        RouletteWheel -> UCB [style=invis];
        UCB -> BayesUCB [style=invis];
        BayesUCB -> EpsilonGreedy [style=invis];
        
        SelectReagent -> RouletteWheel;
        SelectReagent -> UCB;
        SelectReagent -> BayesUCB;
        SelectReagent -> EpsilonGreedy;
    }

SelectionStrategy (Base Class)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. rst-class:: class-fundamental

Abstract base class for all selection strategies. Extend this to create custom strategies.

.. list-table:: Required Methods
   :header-rows: 1
   :widths: 40 60

   * - Method
     - Description
   * - ``select_reagent(reagent_list, disallow_mask, rng, ...)``
     - Select one reagent from the list.
   * - ``select_batch(reagent_list, batch_size, ...)``
     - Select multiple reagents (optional override).

.. _greedy-selection:

GreedySelection
~~~~~~~~~~~~~~~

.. rst-class:: class-config

Simple greedy selection using argmax/argmin of sampled scores.

**Extends:** :ref:`SelectionStrategy <selection-strategy>`

- Fast convergence but may get stuck in local optima
- Best for: Simple optimization landscapes, limited budgets

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 18 10 52

   * - Parameter
     - Type
     - Required
     - Description
   * - ``mode``
     - ``str``
     - No
     - ``"maximize"`` or ``"minimize"``. Default: ``"maximize"``.

**Example**

.. code-block:: python

   from TACTICS.thompson_sampling.strategies import GreedySelection

   strategy = GreedySelection(mode="maximize")
   # For docking scores (lower is better)
   strategy = GreedySelection(mode="minimize")

.. _roulette-wheel-selection:

RouletteWheelSelection
~~~~~~~~~~~~~~~~~~~~~~

.. rst-class:: class-config

Roulette wheel selection with thermal cycling and Component-Aware Thompson Sampling (CATS).

**Extends:** :ref:`SelectionStrategy <selection-strategy>`

- Boltzmann-weighted selection with adaptive temperature control
- Component rotation for systematic exploration
- CATS: Shannon entropy-based criticality analysis
- Best for: Complex multi-modal landscapes, large libraries

.. list-table:: Parameters
   :header-rows: 1
   :widths: 22 15 10 53

   * - Parameter
     - Type
     - Required
     - Description
   * - ``mode``
     - ``str``
     - No
     - ``"maximize"``, ``"minimize"``, ``"maximize_boltzmann"``, or ``"minimize_boltzmann"``.
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
     - Fraction before CATS starts. Default: 0.20.
   * - ``transition_phase_end``
     - ``float``
     - No
     - Fraction when CATS fully applied. Default: 0.60.
   * - ``min_observations``
     - ``int``
     - No
     - Min observations before trusting criticality. Default: 5.

**Example**

.. code-block:: python

   from TACTICS.thompson_sampling.strategies import RouletteWheelSelection

   # Standard thermal cycling
   strategy = RouletteWheelSelection(
       mode="maximize",
       alpha=0.1,
       beta=0.05
   )

   # Higher exploration
   strategy = RouletteWheelSelection(
       mode="maximize",
       alpha=0.2,
       beta=0.1
   )

.. _ucb-selection:

UCBSelection
~~~~~~~~~~~~

.. rst-class:: class-config

Upper Confidence Bound selection with deterministic behavior.

**Extends:** :ref:`SelectionStrategy <selection-strategy>`

- Balances exploitation and exploration via confidence bounds
- Best for: Situations requiring deterministic, reproducible behavior

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
   * - ``mode``
     - ``str``
     - No
     - ``"maximize"`` or ``"minimize"``. Default: ``"maximize"``.
   * - ``c``
     - ``float``
     - No
     - Exploration parameter. Higher = more exploration. Default: 2.0.

**Example**

.. code-block:: python

   from TACTICS.thompson_sampling.strategies import UCBSelection

   strategy = UCBSelection(mode="maximize", c=2.0)
   # Higher exploration
   strategy = UCBSelection(mode="maximize", c=4.0)

.. _epsilon-greedy-selection:

EpsilonGreedySelection
~~~~~~~~~~~~~~~~~~~~~~

.. rst-class:: class-config

Simple exploration strategy with decaying epsilon.

**Extends:** :ref:`SelectionStrategy <selection-strategy>`

- Random selection with probability epsilon, greedy otherwise
- Best for: Baseline comparisons, simple exploration needs

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
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
     - Decay rate per iteration. Default: 0.995.

**Example**

.. code-block:: python

   from TACTICS.thompson_sampling.strategies import EpsilonGreedySelection

   # 20% exploration with decay
   strategy = EpsilonGreedySelection(
       mode="maximize",
       epsilon=0.2,
       decay=0.995
   )

.. _bayes-ucb-selection:

BayesUCBSelection
~~~~~~~~~~~~~~~~~

.. rst-class:: class-config

Bayesian UCB with Student-t quantiles and CATS integration.

**Extends:** :ref:`SelectionStrategy <selection-strategy>`

- Theoretically grounded Bayesian confidence bounds
- Percentile-based thermal cycling (analog to temperature)
- Component-aware exploration based on Shannon entropy
- Best for: Complex landscapes, escaping local optima
- Requires: scipy

.. list-table:: Parameters
   :header-rows: 1
   :widths: 24 15 8 53

   * - Parameter
     - Type
     - Required
     - Description
   * - ``mode``
     - ``str``
     - No
     - ``"maximize"`` or ``"minimize"``. Default: ``"maximize"``.
   * - ``initial_p_high``
     - ``float``
     - No
     - Base percentile for heated component [0.5, 0.999]. Default: 0.90.
   * - ``initial_p_low``
     - ``float``
     - No
     - Base percentile for cooled components [0.5, 0.999]. Default: 0.60.
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
     - Min observations before trusting criticality. Default: 5.

**Example**

.. code-block:: python

   from TACTICS.thompson_sampling.strategies import BayesUCBSelection

   strategy = BayesUCBSelection(mode="maximize")

   # More aggressive exploration
   strategy = BayesUCBSelection(
       mode="maximize",
       initial_p_high=0.95,
       initial_p_low=0.70,
       exploration_phase_end=0.25
   )


Warmup Strategies
-----------------

Warmup strategies determine how reagent combinations are sampled to initialize posteriors
before the main search begins.

.. graphviz::

    digraph WarmupStrategies {
        rankdir=TB;
        node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=10];
        nodesep=0.3;
        ranksep=0.4;

        Init [label="Initialize Priors", fillcolor="#FFFACD"];

        Balanced [label="BalancedWarmup (Recommended)\nK obs per reagent", fillcolor="#00CED1"];
        Standard [label="StandardWarmup\nRandom partners", fillcolor="#AFEEEE"];
        Enhanced [label="EnhancedWarmup\nParallel pairing (Legacy RWS)", fillcolor="#AFEEEE"];

        Init -> Balanced;
        Init -> Standard;
        Init -> Enhanced;
    }

.. _warmup-strategy:

WarmupStrategy (Base Class)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. rst-class:: class-fundamental

Abstract base class for warmup strategies.

.. list-table:: Required Methods
   :header-rows: 1
   :widths: 45 55

   * - Method
     - Description
   * - ``generate_warmup_combinations(reagent_lists, num_trials, disallow_tracker)``
     - Generate list of combinations to evaluate.
   * - ``get_expected_evaluations(reagent_lists, num_trials)``
     - Estimate number of evaluations.
   * - ``get_name()``
     - Return strategy name.

.. _balanced-warmup:

BalancedWarmup (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. rst-class:: class-config

Balanced warmup guaranteeing exactly K observations per reagent with stratified partners.

**Extends:** :ref:`WarmupStrategy <warmup-strategy>`

- Uniform coverage across all reagents
- Per-reagent variance estimation with James-Stein shrinkage
- Reduces bias from random sampling
- Best for: Most use cases, especially asymmetric component sizes

.. list-table:: Parameters
   :header-rows: 1
   :widths: 26 15 8 51

   * - Parameter
     - Type
     - Required
     - Description
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

**Example**

.. code-block:: python

   from TACTICS.thompson_sampling.warmup import BalancedWarmup

   warmup = BalancedWarmup(observations_per_reagent=5)

   # With per-reagent variance
   warmup = BalancedWarmup(
       observations_per_reagent=5,
       use_per_reagent_variance=True,
       shrinkage_strength=3.0
   )

.. _standard-warmup:

StandardWarmup
~~~~~~~~~~~~~~

.. rst-class:: class-config

Standard warmup testing each reagent with random partners.

**Extends:** :ref:`WarmupStrategy <warmup-strategy>`

- Simple and straightforward
- Ensures all reagents evaluated
- Expected evaluations: sum(reagent_counts) * num_trials

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
   * - ``seed``
     - ``int``
     - No
     - Random seed for reproducibility.

.. _enhanced-warmup:

EnhancedWarmup (Legacy)
~~~~~~~~~~~~~~~~~~~~~~~

.. rst-class:: class-config

Stochastic parallel pairing with shuffling from the original RWS algorithm.

**Extends:** :ref:`WarmupStrategy <warmup-strategy>`

- Parallel pairing of reagents across components
- Required for replicating legacy RWS results
- Best for: legacy_rws_maximize and legacy_rws_minimize presets

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
   * - ``seed``
     - ``int``
     - No
     - Random seed for reproducibility.


Evaluators
----------

Evaluators score compounds based on various criteria. Choose based on your data source
and computational requirements.

.. graphviz::

    digraph Evaluators {
        rankdir=TB;
        node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=10];
        nodesep=0.2;
        ranksep=0.3;

        Evaluate [label="Evaluate Compound", fillcolor="#FFFACD"];

        // Fast evaluators
        fast_label [label="Fast", shape=plaintext, fontsize=9];
        Lookup [label="LookupEvaluator (CSV)", fillcolor="#90EE90"];
        DB [label="DBEvaluator (SQLite)", fillcolor="#90EE90"];

        // Computational evaluators
        compute_label [label="Computational", shape=plaintext, fontsize=9];
        FP [label="FPEvaluator (Fingerprints)", fillcolor="#ADD8E6"];
        MW [label="MWEvaluator (Mol Weight)", fillcolor="#ADD8E6"];

        // Slow evaluators
        slow_label [label="Slow (use processes>1)", shape=plaintext, fontsize=9];
        ROCS [label="ROCSEvaluator (3D Shape)", fillcolor="#E6E6FA"];
        FRED [label="FredEvaluator (Docking)", fillcolor="#E6E6FA"];
        ML [label="MLClassifierEvaluator", fillcolor="#E6E6FA"];

        // Force vertical layout
        Evaluate -> fast_label [style=invis];
        fast_label -> Lookup [style=invis];
        Lookup -> DB [style=invis];
        DB -> compute_label [style=invis];
        compute_label -> FP [style=invis];
        FP -> MW [style=invis];
        MW -> slow_label [style=invis];
        slow_label -> ROCS [style=invis];
        ROCS -> FRED [style=invis];
        FRED -> ML [style=invis];

        // Visible connections
        Evaluate -> Lookup;
        Evaluate -> DB;
        Evaluate -> FP;
        Evaluate -> MW;
        Evaluate -> ROCS;
        Evaluate -> FRED;
        Evaluate -> ML;
    }

.. _evaluator-base:

Evaluator (Base Class)
~~~~~~~~~~~~~~~~~~~~~~

.. rst-class:: class-fundamental

Abstract base class for all evaluators.

.. list-table:: Required Methods
   :header-rows: 1
   :widths: 35 65

   * - Method
     - Description
   * - ``evaluate(input)``
     - Score a compound (accepts Mol or product_name depending on evaluator).
   * - ``counter`` (property)
     - Number of evaluations performed.

.. _lookup-evaluator:

LookupEvaluator
~~~~~~~~~~~~~~~

.. rst-class:: class-config

Fast evaluator that looks up pre-computed scores from a CSV file.

**Extends:** :ref:`Evaluator <evaluator-base>`

- Use for: Pre-computed scores, benchmarking
- Recommendation: Use ``processes=1`` (parallel overhead exceeds lookup time)

.. list-table:: Config Parameters (LookupEvaluatorConfig)
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
   * - ``ref_filename``
     - ``str``
     - Yes
     - Path to CSV file with scores.
   * - ``score_col``
     - ``str``
     - No
     - Column name for scores. Default: ``"Scores"``.
   * - ``compound_col``
     - ``str``
     - No
     - Column name for compound IDs. Default: ``"Product_Code"``.

**Example**

.. code-block:: python

   from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig
   from TACTICS.thompson_sampling.factories import create_evaluator

   config = LookupEvaluatorConfig(
       ref_filename="scores.csv",
       score_col="binding_affinity"
   )
   evaluator = create_evaluator(config)

.. _db-evaluator:

DBEvaluator
~~~~~~~~~~~

.. rst-class:: class-config

Fast evaluator using SQLite database for large datasets.

**Extends:** :ref:`Evaluator <evaluator-base>`

- Use for: Large pre-computed datasets (millions of compounds)
- Recommendation: Use ``processes=1``

.. list-table:: Config Parameters (DBEvaluatorConfig)
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
   * - ``db_filename``
     - ``str``
     - Yes
     - Path to SQLite database.
   * - ``db_prefix``
     - ``str``
     - No
     - Key prefix for lookups. Default: ``""``.

.. _fp-evaluator:

FPEvaluator
~~~~~~~~~~~

.. rst-class:: class-config

Evaluator using Morgan fingerprint Tanimoto similarity.

**Extends:** :ref:`Evaluator <evaluator-base>`

- Use for: Similarity-based virtual screening
- Returns: Tanimoto similarity [0, 1]

.. list-table:: Config Parameters (FPEvaluatorConfig)
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
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

.. _mw-evaluator:

MWEvaluator
~~~~~~~~~~~

.. rst-class:: class-config

Simple evaluator returning molecular weight. Primarily for testing.

**Extends:** :ref:`Evaluator <evaluator-base>`

.. _rocs-evaluator:

ROCSEvaluator
~~~~~~~~~~~~~

.. rst-class:: class-config

3D shape-based evaluator using OpenEye ROCS.

**Extends:** :ref:`Evaluator <evaluator-base>`

- Use for: Shape-based virtual screening
- Requires: OpenEye Toolkit license
- Recommendation: Use ``processes>1`` for parallel evaluation

.. list-table:: Config Parameters (ROCSEvaluatorConfig)
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
   * - ``query_molfile``
     - ``str``
     - Yes
     - Path to reference structure (.sdf).
   * - ``max_confs``
     - ``int``
     - No
     - Max conformers to generate. Default: 50.

.. _fred-evaluator:

FredEvaluator
~~~~~~~~~~~~~

.. rst-class:: class-config

Molecular docking evaluator using OpenEye FRED.

**Extends:** :ref:`Evaluator <evaluator-base>`

- Use for: Structure-based virtual screening
- Requires: OpenEye Toolkit license
- Recommendation: Use ``processes>1`` for parallel evaluation
- Mode: ``minimize`` (lower docking scores = better)

.. list-table:: Config Parameters (FredEvaluatorConfig)
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
   * - ``design_unit_file``
     - ``str``
     - Yes
     - Path to receptor file (.oedu).
   * - ``max_confs``
     - ``int``
     - No
     - Max conformers to generate. Default: 100.

.. _ml-classifier-evaluator:

MLClassifierEvaluator
~~~~~~~~~~~~~~~~~~~~~

.. rst-class:: class-config

Evaluator using a trained scikit-learn classifier.

**Extends:** :ref:`Evaluator <evaluator-base>`

- Use for: ML-based scoring with trained models
- Requires: scikit-learn, trained model pickle file

.. list-table:: Config Parameters (MLClassifierEvaluatorConfig)
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
   * - ``model_filename``
     - ``str``
     - Yes
     - Path to pickled sklearn model.


.. _run-ts:

Convenience Function: run_ts()
------------------------------

.. rst-class:: class-core

The simplest way to run Thompson Sampling optimization.

.. admonition:: Dependencies
   :class: dependencies

   - :ref:`ThompsonSamplingConfig <thompson-sampling-config>` - full configuration object
   - Internally creates :ref:`ThompsonSampler <thompson-sampler>` and runs warmup + search

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 22 10 48

   * - Parameter
     - Type
     - Required
     - Description
   * - ``config``
     - ``ThompsonSamplingConfig``
     - Yes
     - Full configuration object.
   * - ``return_warmup``
     - ``bool``
     - No
     - Also return warmup results. Default: False.

**Returns**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Type
     - Description
   * - ``DataFrame`` or ``tuple``
     - Search results, or (search_df, warmup_df) if return_warmup=True.

**Example**

.. code-block:: python

   from TACTICS.library_enumeration import SynthesisPipeline
   from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef
   from TACTICS.thompson_sampling import run_ts, get_preset
   from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

   # Create synthesis pipeline
   rxn_config = ReactionConfig(
       reactions=[ReactionDef(reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]", step_index=0)],
       reagent_file_list=["acids.smi", "amines.smi"]
   )
   pipeline = SynthesisPipeline(rxn_config)

   # Get preset configuration
   config = get_preset(
       "balanced_sampling",
       synthesis_pipeline=pipeline,
       evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv"),
       num_iterations=1000
   )

   # Get both results
   search_df, warmup_df = run_ts(config, return_warmup=True)
   print(f"Warmup: {len(warmup_df)}, Search: {len(search_df)}")


Strategy Selection Guide
------------------------

Choose the right strategy based on your use case:

.. list-table::
   :header-rows: 1
   :widths: 22 25 25 28

   * - Strategy
     - Best For
     - Pros
     - Cons
   * - **Greedy**
     - Simple landscapes, limited budgets
     - Fast convergence
     - Can get stuck in local optima
   * - **RouletteWheel**
     - Complex multi-modal landscapes
     - Thermal cycling, CATS, adaptive
     - More parameters to tune
   * - **UCB**
     - Deterministic optimization needs
     - Theoretically grounded
     - Less stochastic
   * - **BayesUCB**
     - Complex landscapes, escaping optima
     - Bayesian bounds, CATS
     - Requires scipy
   * - **EpsilonGreedy**
     - Baseline comparisons
     - Very simple
     - Less sophisticated


Evaluator Selection Guide
-------------------------

Choose based on your data source and computational requirements:

**Fast Evaluators** (use ``processes=1``):

- ``LookupEvaluator``: Pre-computed scores in CSV
- ``DBEvaluator``: Pre-computed scores in SQLite

**Computational Evaluators**:

- ``FPEvaluator``: Fingerprint similarity (fast)
- ``MWEvaluator``: Molecular weight (testing only)

**Slow Evaluators** (use ``processes>1``):

- ``ROCSEvaluator``: 3D shape similarity (requires OpenEye)
- ``FredEvaluator``: Molecular docking (requires OpenEye)
- ``MLClassifierEvaluator``: ML model predictions

See the :doc:`configuration` page for preset configurations and detailed examples.
