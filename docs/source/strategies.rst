.. _strategies-guide:

Strategies Guide
================

This guide helps you choose the right TACTICS configuration for your use case.
If you just want the best method, use ``get_preset()`` with no arguments — it
defaults to the best-performing configuration from large-scale benchmarking.

.. contents:: On this page
   :local:
   :depth: 2

Which Preset Should I Use?
--------------------------

TACTICS provides five presets. Each corresponds to a method configuration
validated across 21 combinatorial libraries (114,450+ total trials).

.. list-table::
   :header-rows: 1
   :widths: 20 22 18 18 22

   * - Preset
     - Selection Strategy
     - Warmup
     - Recovery
     - When to use
   * - **recommended** (default)
     - Top-Two TS
     - Enhanced
     - **86.1%**
     - New work. Best overall.
   * - **recommended_rws**
     - RWS / CATS
     - Enhanced
     - 85.5%
     - When you want the original TACTICS method.
   * - **baseline**
     - Greedy
     - Balanced
     - ~81%
     - Reference point to measure strategy gains.
   * - **legacy_rws**
     - RWS (round-robin)
     - Enhanced
     - ~76%
     - Reproducing Zhao et al. 2025 results. Pass ``mode="minimize"`` for docking.

Recovery is mean top-100 recovery across all benchmark libraries and queries.
Both ``recommended`` and ``recommended_rws`` default to ``batch_size=100``.

Basic Usage
~~~~~~~~~~~

.. code-block:: python

    from TACTICS.thompson_sampling import ThompsonSampler, get_preset
    from TACTICS.library_enumeration import SynthesisPipeline
    from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef
    from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

    pipeline = SynthesisPipeline(ReactionConfig(
        reactions=[ReactionDef(reaction_smarts="...", step_index=0)],
        reagent_file_list=["acids.smi", "amines.smi"]
    ))

    # Best method — no need to specify a preset name
    config = get_preset(
        synthesis_pipeline=pipeline,
        evaluator_config=LookupEvaluatorConfig(ref_filename="scores.csv"),
    )

    sampler = ThompsonSampler.from_config(config)
    warmup_df = sampler.warm_up(num_warmup_trials=config.num_warmup_trials)
    results_df = sampler.search(num_cycles=config.num_ts_iterations)
    sampler.close()

Docking (Minimize Mode)
~~~~~~~~~~~~~~~~~~~~~~~

For scoring functions where lower is better (docking scores, binding energies):

.. code-block:: python

    config = get_preset(
        synthesis_pipeline=pipeline,
        evaluator_config=FredEvaluatorConfig(design_unit_file="receptor.oedu"),
        mode="minimize",
    )

Custom Batch Size
~~~~~~~~~~~~~~~~~

Both ``recommended`` and ``recommended_rws`` accept a ``batch_size`` parameter
(default 100). Adjust based on your evaluator speed:

.. code-block:: python

    # Slow evaluator — larger batches for better throughput
    config = get_preset(
        synthesis_pipeline=pipeline,
        evaluator_config=evaluator,
        batch_size=500,
    )

    # Fast evaluator — single compound per cycle for tightest feedback loop
    config = get_preset(
        synthesis_pipeline=pipeline,
        evaluator_config=evaluator,
        batch_size=1,
    )


Understanding the Methods
-------------------------

This section explains *why* the recommended methods work. Skip this if you
just want to run TACTICS — the presets handle everything.

How Thompson Sampling Works
~~~~~~~~~~~~~~~~~~~~~~~~~~~

TACTICS maintains a Bayesian posterior distribution (Normal-Normal conjugate)
for each reagent in the combinatorial library. At each iteration:

1. **Sample** from each reagent's posterior
2. **Select** reagents to form a compound (using the selection strategy)
3. **Evaluate** the compound (using the evaluator)
4. **Update** the posteriors with the observed score

The selection strategy determines *how* reagents are chosen from the sampled
posteriors. This is where the methods differ.

Top-Two Thompson Sampling (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Top-Two TS targets **best-arm identification** rather than cumulative regret
minimization. At each step, it samples posteriors twice:

- The **leader** is the reagent with the best first sample
- The **challenger** is the reagent with the best second sample (excluding the leader)
- With probability :math:`\beta`, the challenger is selected instead of the leader

This mechanism allocates exploration budget toward reagents where posterior
uncertainty matters — exactly the reagents whose top-K membership is ambiguous.

TACTICS adds **adaptive per-component thermal cycling** on top of TT-TS:
posteriors are scaled up (heated) or down (cooled) to control exploration
intensity. The adaptive mechanism tracks per-component disagreement rates
via exponential moving averages and self-tunes the heated scale, which is
critical for imbalanced libraries (e.g., 130 acids x 3,844 amines).

Roulette Wheel Selection / CATS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RWS uses **Boltzmann-weighted selection** with thermal cycling. Instead of
always picking the best-sampled reagent (greedy), it converts sampled
posteriors into selection probabilities via a Boltzmann distribution and
samples from that.

**CATS** (Component-Aware Thompson Sampling) adds **GMIC-weighted component
rotation**: instead of cycling through components in round-robin order, it
measures each component's GMIC criticality score and preferentially heats
the most "flexible" component (the one where the top-K set is least certain).

This is the mechanism responsible for the largest gains over legacy methods
(+4.3 pts on 2-component libraries).

GMIC Criticality
~~~~~~~~~~~~~~~~

GMIC (Generalized Mean Information Coefficient) measures how concentrated
or spread the posterior means are for a component's reagents. A high GMIC
indicates a component where a few reagents dominate — this component is
"solved" and needs less exploration. A low GMIC indicates a flat landscape
where more exploration would help.

TACTICS uses GMIC to allocate thermal cycling budget: components with lower
GMIC get heated more often, directing search toward unsolved parts of the
chemical space.

Enhanced vs Balanced Warmup
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Enhanced warmup** (recommended) uses stochastic parallel pairing: all
reagents are shuffled and paired exhaustively in each warmup trial. On
imbalanced libraries, this naturally over-samples the small component,
pre-solving its ranking during warmup. GMIC then detects this differential
posterior quality and redirects search to the unsolved large component.

**Balanced warmup** guarantees exactly K observations per reagent via
stratified partner selection. It is useful for isolating the framework's
warmup contribution (Balanced-Greedy provides +1.5 pts over Legacy-Greedy
on 2-component libraries), but Enhanced warmup is universally better when
paired with GMIC rotation.

Boltzmann Weighting
~~~~~~~~~~~~~~~~~~~

When ``use_boltzmann_weighting=True`` (enabled in all recommended presets),
the posterior update weights observations exponentially: good scores have
more influence on the posterior than bad scores. This accelerates
convergence toward high-scoring reagents.


Custom Configuration
--------------------

For full control, construct a ``ThompsonSamplingConfig`` directly:

.. code-block:: python

    from TACTICS.thompson_sampling import ThompsonSamplingConfig, ThompsonSampler
    from TACTICS.thompson_sampling.strategies.config import TopTwoConfig
    from TACTICS.thompson_sampling.warmup.config import EnhancedWarmupConfig
    from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

    config = ThompsonSamplingConfig(
        synthesis_pipeline=pipeline,
        num_ts_iterations=5000,
        num_warmup_trials=5,
        strategy_config=TopTwoConfig(
            mode="maximize",
            beta=0.5,              # Challenger selection probability
            heated_scale=1.5,      # Posterior inflation for heated component
            cooled_scale=0.75,     # Posterior deflation for cooled component
            adaptive_disagreement=True,  # Per-component adaptive scaling
        ),
        warmup_config=EnhancedWarmupConfig(),
        evaluator_config=LookupEvaluatorConfig(
            ref_filename="scores.csv",
            score_col="binding_affinity",
        ),
        batch_size=100,
        use_boltzmann_weighting=True,
    )

    sampler = ThompsonSampler.from_config(config)
    warmup_df = sampler.warm_up(num_warmup_trials=config.num_warmup_trials)
    results_df = sampler.search(num_cycles=config.num_ts_iterations)
    sampler.close()

See the :doc:`configuration` reference for all available parameters.

Selection Strategies Reference
------------------------------

All strategies inherit from ``SelectionStrategy`` and are interchangeable
via the ``strategy_config`` parameter.

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Strategy
     - Tier
     - Description
   * - :ref:`TopTwoSelection <top-two-selection>`
     - Recommended
     - Top-Two TS with adaptive thermal cycling and GMIC rotation.
       Best overall (86.1% recovery).
   * - :ref:`RouletteWheelSelection <roulette-wheel-selection>`
     - Recommended
     - CATS with Boltzmann thermal cycling and GMIC rotation.
       Close second (85.5% recovery).
   * - :ref:`GreedySelection <greedy-selection>`
     - Baseline
     - Pure argmax. Simplest strategy — sample posteriors, pick the best.
   * - :ref:`UCBSelection <ucb-selection>`
     - Baseline
     - Upper Confidence Bound. Consistently outperformed by TopTwo/RWS.
   * - :ref:`EpsilonGreedySelection <epsilon-greedy-selection>`
     - Baseline
     - Epsilon-greedy with decay. Consistently outperformed by TopTwo/RWS.
   * - :ref:`BayesUCBSelection <bayes-ucb-selection>`
     - Baseline
     - Bayesian UCB with Student-t quantiles. Consistently outperformed by TopTwo/RWS.

See the :doc:`thompson_sampling` reference for detailed parameter documentation
for each strategy.

Warmup Strategies Reference
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Strategy
     - Tier
     - Description
   * - :ref:`EnhancedWarmup <enhanced-warmup>`
     - Recommended
     - Stochastic parallel pairing. Universally optimal — used by all
       recommended presets.
   * - :ref:`BalancedWarmup <balanced-warmup>`
     - Alternative
     - Exactly K observations per reagent. Useful for controlled experiments
       isolating warmup contribution.
   * - :ref:`StandardWarmup <standard-warmup>`
     - Baseline
     - Random partner selection. For comparison only.
