.. _guide-inspect:

5. Inspect
==========

**What this block is.** The results DataFrame tells you *which* products
scored well. The posteriors and the per-cycle trajectory tell you *why* —
which component the search decided was solved, when, and which reagents
carry the signal. That is the SAR the search has learned, and it is
available from the same sampler object once the search ends.

Core
----

Turn on ``track_diagnostics`` before building the sampler, then ask for
three things after ``search()``:

.. literalinclude:: /snippets/inspect_diagnostics.py
   :language: python
   :start-after: # [start:track]
   :end-before: # [end:track]
   :dedent:

.. list-table::
   :widths: 30 70

   * - ``get_diagnostics()``
     - One row per (cycle, component): what the strategy knew and did at
       that moment. Columns depend on the strategy (below); every schema
       shares ``current_cycle``, ``component_idx`` and ``criticality``.
       Empty for Greedy/UCB/ε-greedy, which keep no component state.
   * - ``get_posterior_landscape()``
     - One row per reagent — ``component_idx``, ``reagent_name``, ``mean``,
       ``std``, ``n_samples``. The learned SAR, as a table. Works without
       ``track_diagnostics``.
   * - ``get_sar_summary()``
     - A dict: per-component entropy, concentration, Gini, top-1/top-5
       share, SNR, the dominant reagent and a convergence verdict;
       ``landscape_type`` (``structured_SAR`` / ``diffuse_SAR`` / ``mixed`` /
       ``insufficient_data``); ``summary_text`` in prose; and, when
       diagnostics were tracked, ``convergence_dynamics`` per component.

**What it produces:** two Polars DataFrames and a dict. Write the frames to
Parquet — everything in *Build on it* and in :doc:`06_visualise` runs on
the saved files without a sampler.

Build on it
-----------

Reading the trajectory
~~~~~~~~~~~~~~~~~~~~~~

Each recommended strategy records what drives its decisions. The columns
worth watching:

**Top-Two TS** (12 columns)
   ``gmic`` and ``criticality`` (the component's GMIC), ``is_heated``,
   ``heated_scale`` / ``cooled_scale`` / ``effective_scale`` (what the
   posterior std was multiplied by), ``disagreement_ema`` (how often the two
   draws disagreed — the signal the adaptive scale follows),
   ``n_active_reagents``.

**Roulette wheel / CATS** (18 columns)
   ``gmic`` with its parts ``signal_var`` / ``mean_noise_var``;
   ``divergence`` against ``divergence_threshold`` and the resulting
   ``is_stable`` / ``cats_mode`` (``"diversity"`` until posteriors settle,
   then GMIC-driven); ``base_temp`` → ``cats_multiplier`` →
   ``final_temperature``; ``ema_relative_gmic``.

**Bayes-UCB** (18 columns)
   ``criticality`` from the z-score softmax with ``snr``,
   ``participation_ratio`` / ``effective_n``, ``sharpening_factor``, the
   observation-gated ``criticality_weight`` and its ``decay``, and the
   percentile pipeline ``base_temp`` → ``cats_multiplier`` →
   ``final_temperature``.

A component whose ``gmic`` climbs and stays high early is one the search
considered solved; one that stays flat is where the budget went.

Analysis functions
~~~~~~~~~~~~~~~~~~

:mod:`TACTICS.thompson_sampling.diagnostics` is Polars in, Polars out, no
sampler needed:

.. literalinclude:: /snippets/inspect_diagnostics.py
   :language: python
   :start-after: # [start:analyse]
   :end-before: # [end:analyse]
   :dedent:

- :func:`~TACTICS.thompson_sampling.diagnostics.compute_posterior_entropy`
  — from the landscape: per-component entropy, concentration and SNR.
  Pass the same ``mode`` as the search.
- :func:`~TACTICS.thompson_sampling.diagnostics.compute_convergence_point`
  — from the diagnostics: the first cycle each component's criticality
  crossed ``threshold`` and the cycle it stayed there.
- :func:`~TACTICS.thompson_sampling.diagnostics.compare_trajectory_vs_snapshot`
  — joins the two: does the trajectory add information the final snapshot
  does not?
- :func:`~TACTICS.thompson_sampling.diagnostics.compute_disagreement_convergence`
  and :func:`~TACTICS.thompson_sampling.diagnostics.compute_scale_adaptation`
  — Top-Two-specific: when disagreement settled, how far the heated scale
  moved.
- :func:`~TACTICS.thompson_sampling.diagnostics.format_sar_report` — the
  summary dict as text.

Many runs
~~~~~~~~~

To compare methods or replicates, save each run's diagnostics with a
``replicate`` (and ``method``) column added, concatenate, and hand the
result to the plots in :doc:`06_visualise` — they average across
replicates and draw confidence bands.

.. admonition:: Gotchas
   :class: warning

   - ``track_diagnostics`` must be set **before** ``from_config`` /
     construction. Setting it afterwards records nothing.
   - The diagnostics schema is strategy-specific; code that reads
     ``disagreement_ema`` will fail on an RWS run. Branch on
     ``"disagreement_ema" in df.columns`` or on the strategy you ran.
   - ``get_sar_summary()`` needs observations: with very few cycles the
     verdict is ``insufficient_data`` and the per-component numbers are
     ``NaN``.
   - ``compute_posterior_entropy`` is direction-aware — pass ``mode`` or a
     minimize run reads as inverted.

Reference
---------

:doc:`/reference/inspect` — the three accessors, the six analysis
functions, and the full column lists.

:doc:`/theory/search_performance_metrics` — the metrics behind these
numbers.

**Next:** :doc:`06_visualise` — pictures of all of this.
