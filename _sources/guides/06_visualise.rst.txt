.. _guide-visualise:

6. Visualise
============

**What this block is.** Two families of charts over the frames Blocks 3
and 5 produce. *Benchmarks* answer "which method found more of the known
hits, and how fast" — Altair, interactive HTML. *Diagnostic plots* show the
mechanism — what GMIC, temperature and disagreement did over the run —
matplotlib figures for a paper. Both live in ``TACTICS.library_analysis``
and need the ``viz`` extra:

.. code-block:: bash

   pip install "chem-tactics[viz]"

Importing ``TACTICS.library_analysis`` does not load matplotlib or Altair;
they are imported when you build a chart, and a missing extra raises an
``ImportError`` that says so.

Core
----

Compare two presets on the thrombin data against its exhaustive scores:

.. literalinclude:: /snippets/visualise_benchmarks.py
   :language: python
   :start-after: # [start:benchmarks]
   :end-before: # [end:benchmarks]
   :dedent:

:class:`~TACTICS.library_analysis.visualization.TS_Benchmarks` takes, per
method, a list of ``search()`` DataFrames — one per replicate — and a
reference frame of the ground truth. It counts how many of the ``top_n``
best reference compounds each replicate found, and:

- ``plot_barplot_TS_results`` — hits per replicate and per method, with
  the reference total as the ceiling;
- ``plot_line_performance_with_error_bars`` — fraction of the top-*N*
  recovered, for each *N* in ``top_ns``, with ±1 sd bands;
- ``stripplot_TS_results`` — every scored product's score, jittered by
  replicate and coloured by method.

Each returns an Altair chart (display it in a notebook) and writes HTML
when ``save_path`` is given.

**What it produces:** interactive HTML files, or chart objects for a
notebook.

Build on it
-----------

Diagnostic plots
~~~~~~~~~~~~~~~~

:mod:`TACTICS.library_analysis.diagnostic_plots` draws the trajectories
from :doc:`05_inspect`. Every function takes Polars frames and returns a
``matplotlib.figure.Figure``:

.. literalinclude:: /snippets/visualise_diagnostic_plot.py
   :language: python
   :start-after: # [start:plot]
   :end-before: # [end:plot]
   :dedent:

**Single-run trajectories** — one ``get_diagnostics()`` frame, any strategy
that records ``criticality``:

- :func:`~TACTICS.library_analysis.diagnostic_plots.plot_criticality_trajectory`
  — criticality per component over cycles.
- :func:`~TACTICS.library_analysis.diagnostic_plots.plot_snr_trajectory`
  — SNR with the threshold line (Bayes-UCB frames).
- :func:`~TACTICS.library_analysis.diagnostic_plots.plot_temperature_decomposition`
  — base temperature → CATS multiplier → final temperature for one
  component (RWS / Bayes-UCB frames).

**Replicate-averaged mechanism plots** — require a ``replicate`` column:

- :func:`~TACTICS.library_analysis.diagnostic_plots.plot_rws_diagnostic`
  — GMIC with 95 % bands on the left axis, CATS multiplier on the right,
  diversity-mode shading.
- :func:`~TACTICS.library_analysis.diagnostic_plots.plot_ttts_diagnostic`
  — disagreement EMA and adapted ``heated_scale``.

**Layered panels** — require the *search* frame too, with a ``phase``
column (``"search"``), ``replicate``, and ``method`` for the last two, plus
the reagent files so product names map back to components:

- :func:`~TACTICS.library_analysis.diagnostic_plots.plot_reagent_usage_action_panel`
  — heating timeline over the fraction of each component's reagents tried
  per batch, new vs revisited.
- :func:`~TACTICS.library_analysis.diagnostic_plots.plot_gmic_directed_exploration`
  — *Layer 1*: which component GMIC chose to explore, and what came of it,
  one row per method.
- :func:`~TACTICS.library_analysis.diagnostic_plots.plot_adaptive_intensity`
  — *Layer 2*: how each method tuned exploration *within* the component
  (``cats_multiplier`` for RWS, ``heated_scale`` for TT-TS).

The column contract, since nothing adds these for you:

.. code-block:: python

   import polars as pl

   search = results.with_columns(
       pl.lit("search").alias("phase"),      # warm_up() rows would be "warmup"
       pl.lit(0).alias("replicate"),
       pl.lit("TT-TS").alias("method"),
   )
   diag = diagnostics.with_columns(pl.lit(0).alias("replicate"), pl.lit("TT-TS").alias("method"))

Stack several runs with ``pl.concat`` and pass ``replicate=None`` to
average, or ``replicate=2`` to show one.

Saving
~~~~~~

Altair: ``save_path="chart.html"``. matplotlib: ``fig.savefig("fig.png",
dpi=150, bbox_inches="tight")`` or ``fig.savefig("fig.pdf")`` for a
manuscript.

.. admonition:: Gotchas
   :class: warning

   - ``TS_Benchmarks(no_of_cycles=...)`` is the number of **replicate
     runs** per method, not search cycles. The lists in ``TS_runs_data``
     must each have exactly that many DataFrames.
   - The reference frame needs columns ``score`` and ``Name`` **in that
     order** — the same shape as ``search()`` output minus ``SMILES``; a
     different order raises a Polars ``ShapeError``.
   - ``sort_type`` on ``TS_Benchmarks`` must match the search ``mode``
     (``"minimize"`` for docking) or the "top" reference compounds are the
     worst ones.
   - The replicate-averaged and layered plots read from benchmark output
     with ``replicate`` / ``phase`` / ``method`` columns that a single run
     does not carry — add them as above.
   - ``plt.show()`` is not called for you; in a script, save the figure.

Reference
---------

:doc:`/reference/visualise` — ``TS_Benchmarks`` and the eight plot
functions with every argument.

Interactive: ``marimo edit tutorials/diagnostic_benchmark_plots.py`` renders
the mechanism plots over the published benchmark output;
``tutorials/thompson_sampling_tutorial.py`` builds a ``TS_Benchmarks``
comparison live.

**Next:** :doc:`/extending` — write your own strategy, warmup or scorer
class.
