.. _guide-extending:

Extending TACTICS
=================

The blocks are built on three abstract classes. Each is small; a working
implementation is a page of code. This page shows the contract the sampler
relies on for each, and what to reuse.

A scoring function (no subclass needed)
---------------------------------------

Most custom scoring needs no new class:
:class:`~TACTICS.thompson_sampling.core.evaluator_config.CustomEvaluatorConfig`
wraps any ``Callable[[Mol], float]`` (see :doc:`guides/02_scoring`). Write
an evaluator class when the scorer holds state that must be built once per
worker — a loaded model, a receptor, a database handle.

A new evaluator class
---------------------

Subclass :class:`~TACTICS.thompson_sampling.core.evaluators.Evaluator`;
implement ``evaluate(mol) -> float`` and the ``counter`` property. Return
``NaN`` for "could not score" — the sampler skips it.

.. code-block:: python

   import numpy as np
   from TACTICS import Evaluator

   class ModelEvaluator(Evaluator):
       def __init__(self, input_dict):
           import joblib
           self.model = joblib.load(input_dict["model_path"])   # once per process
           self.num_evaluations = 0

       @property
       def counter(self):
           return self.num_evaluations

       def evaluate(self, mol):
           self.num_evaluations += 1
           try:
               return float(self.model.predict_one(mol))
           except Exception:
               return np.nan

To use it with ``processes > 1`` it needs a picklable config the workers
can rebuild from. Add a Pydantic model next to the others in
``core/evaluator_config.py`` with an ``evaluator_type`` discriminator, add
it to ``EvaluatorConfigType`` in ``thompson_sampling/config.py``, and add
an ``isinstance`` branch in
:func:`~TACTICS.thompson_sampling.factories.create_evaluator`. For a
single-process run, ``sampler.set_evaluator(ModelEvaluator({...}))`` on a
directly-constructed sampler is enough.

A selection strategy
--------------------

Subclass :class:`~TACTICS.thompson_sampling.strategies.base_strategy.SelectionStrategy`.
The one abstract method is

.. code-block:: python

   def select_reagent(self, reagent_list, disallow_mask=None, **kwargs) -> int

The sampler calls it once per component per cycle with, in ``kwargs``:
``rng`` (a ``numpy.random.Generator`` — use it, never the global RNG, or
``seed`` stops working), ``component_idx``, ``iteration``,
``current_cycle`` and ``total_cycles``. ``reagent_list`` holds
:class:`~TACTICS.thompson_sampling.core.reagent.Reagent` objects with
``mean``, ``std``, ``n_samples`` and ``sample(rng)``. ``disallow_mask`` is
a set of indices that must not be returned. ``self.mode`` is ``"maximize"``
or ``"minimize"``.

That is a complete strategy. Greedy is nine lines:

.. literalinclude:: ../../src/TACTICS/thompson_sampling/strategies/greedy_selection.py
   :language: python
   :pyobject: GreedySelection

**Optional hooks.** The sampler probes for these with ``hasattr`` after each
cycle, so a strategy opts in by defining them:

.. list-table::
   :widths: 40 60

   * - ``rotate_component_weighted(n_components, reagent_lists, rng)``
     - Choose the next heated component. Preferred over
       ``rotate_component(n_components)`` (round-robin) when both exist.
   * - ``adapt_temperatures(n_unique, n_attempted)``
     - React to sampling efficiency (duplicates rising → heat up).
   * - ``adapt_heated_scale()``
     - Per-cycle self-tuning (TT-TS uses it for the disagreement EMA).
   * - ``get_component_state(reagent_list, component_idx, current_cycle, total_cycles) -> dict | None``
     - What ``get_diagnostics()`` records each cycle. Return a flat dict;
       include ``current_cycle``, ``component_idx`` and ``criticality`` so
       the analysis functions work. Return ``None`` to record nothing.
   * - ``get_component_criticality(reagent_list) -> float | None``
     - A scalar "how solved is this component"; used by the diagnostics.

**Reuse the thermal-cycling machinery.** Two mixins in
``strategies/_thermal.py`` implement the heated-component index, both
rotation methods, and GMIC:

.. literalinclude:: /snippets/extending_strategy.py
   :language: python
   :start-after: # [start:strategy]
   :end-before: # [end:strategy]

With the mixin, ``rotate_component_weighted`` (GMIC-weighted),
``rotate_component`` and ``get_component_criticality`` come for free.
:class:`~TACTICS.thompson_sampling.strategies._thermal.ThermalCyclingMixin`
alone gives the rotation scaffold with a ``_rotation_flexibility`` hook for
a different weighting.

To make it configurable, add a Pydantic model in ``strategies/config.py``
(subclass ``_StrictModel`` so typos raise), add it to ``StrategyConfigType``,
and add a branch to
:func:`~TACTICS.thompson_sampling.factories.create_strategy`.

A warmup strategy
-----------------

Subclass :class:`~TACTICS.thompson_sampling.warmup.base.WarmupStrategy`
and implement

.. code-block:: python

   def generate_warmup_combinations(self, reagent_lists, num_warmup_trials, disallow_tracker) -> list[list[int]]
   def get_name(self) -> str

returning a list of reagent-index tuples (one index per component) to
evaluate before the search. Override ``get_expected_evaluations`` if the
count is not ``sum(len(rl)) * num_warmup_trials``, so the progress bar is
right. The sampler also reads two optional attributes:
``use_per_reagent_variance`` (default off) and ``shrinkage_strength``
(default 3.0) — set them to opt into the James–Stein per-reagent variance
that :class:`~TACTICS.thompson_sampling.warmup.balanced.BalancedWarmup`
uses. Draw randomness from a generator you own (``seed`` attribute) so runs
are reproducible.

Wire it into ``warmup/config.py`` → ``WarmupConfigType`` →
:func:`~TACTICS.thompson_sampling.factories.create_warmup` the same way.

How the package is laid out
---------------------------

- ``library_enumeration/`` — Block 1. ``SynthesisPipeline`` wraps the
  Pydantic ``ReactionConfig``; ``smarts_toolkit/_validator.py`` is the
  chemistry checker.
- ``thompson_sampling/core/`` — the sampler, reagent posteriors, evaluators,
  the disallow tracker and the parallel evaluator.
- ``thompson_sampling/strategies/``, ``warmup/`` — one file per class plus a
  ``config.py`` of Pydantic models; ``factories.py`` maps config → object.
- ``library_analysis/`` — Block 6, behind the ``viz`` extra.

The package ``__init__`` files re-export names **lazily** (PEP 562, via
:func:`TACTICS._lazy.install`): ``import TACTICS`` costs ~40 ms and nothing
heavy loads until a name is used. Every re-export in a hub is a
``{name: submodule}`` entry; add yours there and to the ``TYPE_CHECKING``
block beside it. ``tests/test_import_time.py`` fails the suite if a
config-only import starts pulling in RDKit, SciPy or the plotting stack —
so import heavy things inside the function that needs them, as
``evaluators.py`` does for OpenEye.

Tests live in ``tests/``; ``pytest tests/`` runs in under a minute on the
bundled data. The docs' code examples are files under
``docs/source/snippets/`` executed by ``tests/test_doc_snippets.py``, so a
new feature's example is also its test.

Reference
---------

:doc:`/reference/extending` — the three base classes and the mixins.
