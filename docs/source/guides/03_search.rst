.. _guide-search:

3. Search
=========

**What this block is.** The search is a loop: sample each reagent's
posterior, pick one reagent per component, make the product, score it,
update the posteriors. A *preset* fixes every choice in that loop to the
combination that won the benchmarks; a *config* lets you change any of
them. Either way the entry point is the same object,
:class:`~TACTICS.thompson_sampling.config.ThompsonSamplingConfig`, and the
same five calls.

Core
----

Library (Block 1) plus scorer (Block 2), through the recommended preset:

.. literalinclude:: /snippets/search_minimal.py
   :language: python
   :start-after: # [start:search]
   :end-before: # [end:search]
   :dedent:

The five calls:

.. list-table::
   :widths: 30 70

   * - ``get_preset(...)``
     - Returns a fully populated
       :class:`~TACTICS.thompson_sampling.config.ThompsonSamplingConfig`.
       Everything you pass overrides the preset's default for that field.
   * - ``ThompsonSampler.from_config(config)``
     - Builds the strategy, warmup and evaluator from the config, reads the
       reagent files, and checks which SMARTS pattern each reagent matches.
   * - ``sampler.warm_up(num_warmup_trials=…)``
     - Scores enough products to give every reagent a starting posterior.
       Returns those products as a DataFrame.
   * - ``sampler.search(num_cycles=…)``
     - The loop. Returns every product scored during the search.
   * - ``sampler.close()``
     - Shuts down the worker pool (a no-op with ``processes=1``, but always
       call it).

**The parameters that matter first**

- ``num_iterations`` (→ ``num_ts_iterations``): cycles. Each cycle samples
  ``batch_size`` products, so the evaluation budget is roughly
  ``num_iterations × batch_size`` plus warmup. 1,000 × 100 is the benchmark
  setting for a ~500 k library.
- ``batch_size``: products per cycle. Larger batches update the posteriors
  less often but suit parallel scoring; 50–100 is the usual range.
- ``num_warmup_trials``: how many times each reagent is tried before the
  search starts (preset default 5). See *Warmup* below for what one trial
  costs.
- ``seed``: makes reagent selection and component rotation reproducible.
  Set it on the config after ``get_preset`` (``config.seed = 42``) or in a
  hand-built config.

**What it produces:** a Polars DataFrame with one row per scored product —
``score``, ``SMILES``, ``Name`` — in evaluation order. Sort it, join it to
your own tables, or ``results.write_parquet("run.parquet")``. Nothing is
written to disk unless you do.

.. admonition:: Imports
   :class: note

   ``TACTICS`` and ``TACTICS.thompson_sampling`` export the same names, and
   both resolve them lazily — ``import TACTICS`` takes ~40 ms, and RDKit,
   SciPy and the rest load on first use. ``from TACTICS import
   ThompsonSampler, get_preset, TopTwoConfig, LookupEvaluatorConfig`` is the
   whole import surface for most scripts. Deep paths
   (``TACTICS.thompson_sampling.core.sampler.ThompsonSampler``) also work and
   are what the reference pages use.

Build on it
-----------

Which preset
~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 25 15 40

   * - Preset
     - Strategy + warmup
     - Top-100 recovery
     - When
   * - ``recommended`` *(default)*
     - Top-Two TS + Enhanced
     - **86.1 %**
     - New work. Best overall across 21 libraries / 114 k trials.
   * - ``recommended_rws``
     - Roulette wheel (CATS) + Enhanced
     - 85.5 %
     - The original TACTICS method; sometimes wins on particular libraries.
       If you are benchmarking, run both.
   * - ``baseline``
     - Greedy + Balanced (K = 5)
     - —
     - What the warmup alone buys you (+1.5 pts over random warmup on
       2-component libraries). Batch size 1.

Both ``recommended`` presets run with ``use_boltzmann_weighting=True`` — the
posterior update that weights good observations more heavily — and five
warmup trials. Presets take ``num_iterations``, ``batch_size``, ``mode`` and
``output_dir`` (where the run log goes); anything else you set on the
returned config.

Hand-built config
~~~~~~~~~~~~~~~~~

The config is a plain Pydantic model. Build one when you want a specific
strategy, warmup or parameter the preset does not expose:

.. literalinclude:: /snippets/search_handbuilt_config.py
   :language: python
   :start-after: # [start:config]
   :end-before: # [end:config]
   :dedent:

Every field is documented under :class:`~TACTICS.thompson_sampling.config.ThompsonSamplingConfig`.
Strategy and warmup configs reject unknown fields, so a typo raises
``ValidationError`` rather than being ignored.

Choosing a strategy
~~~~~~~~~~~~~~~~~~~

TACTICS keeps a Normal posterior per reagent. Every cycle it draws from
each posterior and picks one reagent per component; the *strategy* is the
rule for turning draws into a pick. All strategies share one idea,
**thermal cycling**: one component at a time is "heated" (its selection
made more exploratory) while the rest are "cooled" (more exploitative), and
the heated component rotates.

**Top-Two TS** — :class:`~TACTICS.thompson_sampling.strategies.config.TopTwoConfig`
   Draws twice. If the two draws disagree about the best reagent, takes the
   challenger with probability ``beta`` (0.5). Heating scales the posterior
   *standard deviation* (``heated_scale`` 1.5, ``cooled_scale`` 0.75), and
   by default the heated scale adapts per component from how often the two
   draws disagree (``adaptive_disagreement``). Targets *finding the top set*
   rather than average reward, which is what recovery measures.

**Roulette wheel / CATS** — :class:`~TACTICS.thompson_sampling.strategies.config.RouletteWheelConfig`
   Turns draws into selection probabilities with a Boltzmann softmax at
   temperature ``alpha`` (heated, 0.1) or ``beta`` (cooled, 0.05). CATS
   modulates the heated temperature by how "solved" the component looks
   (GMIC, below), gated on the posteriors having stabilised
   (``divergence_threshold``).

**Baselines** — ``GreedyConfig``, ``UCBConfig``, ``EpsilonGreedyConfig``, ``BayesUCBConfig``
   Argmax of the draws; UCB1; ε-greedy with decay; Bayes-UCB on Student-t
   quantiles with CATS. Kept for comparison. Greedy, UCB and ε-greedy do no
   thermal cycling and record no diagnostics. All are in
   :doc:`/reference/search`.

**Which component to heat.** Both recommended strategies rotate the heated
component by *GMIC* (Gaussian Mutual Information Criticality:
``0.5·log(1 + var(means) / mean(variances))``). A high-GMIC component has
clear winners and is left cool; a low-GMIC one is heated more often. This
is what lets the search spend its budget on the component that is still
undecided — the largest single gain over round-robin rotation.

Choosing a warmup
~~~~~~~~~~~~~~~~~

Before the first posterior exists, every reagent needs a few observations.

- :class:`~TACTICS.thompson_sampling.warmup.config.EnhancedWarmupConfig`
  *(default)* — each trial shuffles every component and pairs reagents
  exhaustively. One trial costs ``max(component sizes)`` products, so the
  small component is over-sampled: on 130 acids × 3,844 amines, five trials
  give every amine 5 observations and every acid ~150. That pre-solves the
  small component, and GMIC rotation then spends the search on the large
  one.
- :class:`~TACTICS.thompson_sampling.warmup.config.BalancedWarmupConfig`
  — exactly ``observations_per_reagent`` (K, default 5) per reagent with
  stratified partners; ``seed`` and James–Stein-shrunk per-reagent variance.
  Costs ``sum(component sizes) × K``. Use it when you want the warmup
  contribution held constant across an experiment.

.. _guide-search-direct:

Driving the sampler directly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The config is convenience; the sampler works without it. Useful in tests
and in pipelines that already hold strategy/evaluator objects:

.. literalinclude:: /snippets/search_direct.py
   :language: python
   :start-after: # [start:direct]
   :end-before: # [end:direct]
   :dedent:

``max_evaluations`` caps the number of scored products regardless of
``num_cycles``. Without a config the sampler cannot rebuild the evaluator
in workers, so ``processes > 1`` needs ``set_evaluator(evaluator,
evaluator_config=…)`` — see :doc:`04_scale`.

.. admonition:: Gotchas
   :class: warning

   - ``mode`` defaults to ``"maximize"``. Docking scores need
     ``mode="minimize"`` or you optimise for the worst binders.
   - ``get_preset("fast_exploration")`` and other names from 1.x do not
     exist; the three above are the full list.
   - ``seed`` reproduces selection and rotation. Enhanced warmup pairs with
     the standard-library ``random`` module and is not seeded by it — seed
     ``random`` yourself for a bit-identical warmup, or use Balanced warmup
     with ``seed=``.
   - ``search()`` returns only products scored *during the search*; warmup
     products come back from ``warm_up()``. No product is ever scored twice
     (a disallow tracker prevents resampling), so concatenating the two
     frames is every evaluation the run made.
   - Always ``close()``. With ``processes > 1`` an unclosed pool keeps
     worker processes alive.

Reference
---------

:doc:`/reference/search` — the config, presets, sampler, every strategy and
warmup with their fields.

:doc:`/theory/index` — derivations of Thompson Sampling, CATS, TT-TS and
the warmup analysis.

Interactive: ``marimo edit tutorials/thompson_sampling_tutorial.py`` runs a
strategy × warmup comparison on the thrombin data with recovery charts.

**Next:** :doc:`04_scale` — when scoring is slow.
