Library Analysis
================

Post-run analysis and visualisation for Thompson Sampling benchmark output.
Two things live here:

- :ref:`TS_Benchmarks <ts-benchmarks>` — Altair charts comparing recovery
  across cycles and strategies.
- :ref:`diagnostic_plots <diagnostic-plots>` — matplotlib figures over the
  per-cycle diagnostics a run records when ``track_diagnostics=True``.

Both need the optional ``viz`` extra:

.. code-block:: bash

   pip install "chem-tactics[viz]"

Importing ``TACTICS.library_analysis`` (or ``diagnostic_plots``) does not load
matplotlib or altair; they are imported when a chart is first built, and a
missing extra raises an ``ImportError`` that says so.

.. _ts-benchmarks:

TS_Benchmarks
-------------

.. rst-class:: class-core

Comprehensive benchmarking and visualization for Thompson Sampling results across
multiple cycles and search strategies.

**Key Features:**

- Automatic data generation during initialization
- Consistent color schemes across all plots
- Multi-cycle analysis with statistics
- Reference comparison with performance metrics
- Multiple visualization types (strip, bar, line plots)

Constructor
~~~~~~~~~~~

All required data is automatically generated during initialization.

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 22 10 48

   * - Parameter
     - Type
     - Required
     - Description
   * - ``no_of_cycles``
     - ``int``
     - Yes
     - Number of cycles to analyze.
   * - ``methods_list``
     - ``list[str]``
     - Yes
     - List of method names (search strategies).
   * - ``TS_runs_data``
     - ``dict``
     - Yes
     - Maps method names to lists of DataFrames (one per cycle).
   * - ``reference_data``
     - ``DataFrame``
     - No
     - Ground truth reference data for comparison.
   * - ``top_n``
     - ``int``
     - No
     - Top products for bar plot. Default: 100.
   * - ``sort_type``
     - ``str``
     - No
     - ``"minimize"`` or ``"maximize"``. Default: ``"minimize"``.
   * - ``top_ns``
     - ``list[int]``
     - No
     - Top-N values for line plot. Default: [50, 100, 200, 300, 400, 500].

**Automatic Data Storage:**

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Attribute
     - Description
   * - ``combined_df_top_n``
     - Top N compounds from each method/cycle (for stripplot).
   * - ``combined_df_all``
     - All compounds from each method/cycle.
   * - ``bar_plot_df``
     - Hit recovery data for bar plots.
   * - ``line_plot_df``
     - Raw performance data across cycles.
   * - ``grouped_stats``
     - Statistical summaries with mean, std, error bounds.
   * - ``actual_methods``
     - Methods found in data (for validation).

Visualization Methods
~~~~~~~~~~~~~~~~~~~~~

stripplot_TS_results
^^^^^^^^^^^^^^^^^^^^

Generate strip plot showing score distributions across cycles and methods.

.. list-table:: Parameters
   :header-rows: 1
   :widths: 18 15 10 57

   * - Parameter
     - Type
     - Required
     - Description
   * - ``width``
     - ``int``
     - No
     - Plot width in pixels (auto-calculated if None).
   * - ``height``
     - ``int``
     - No
     - Plot height in pixels (auto-calculated if None).
   * - ``save_path``
     - ``str``
     - No
     - Path to save (.html, .png, .svg).
   * - ``show_plot``
     - ``bool``
     - No
     - Display in Jupyter. Default: True.
   * - ``legend_position``
     - ``str``
     - No
     - Position of legend: ``"right"`` (default) or ``"bottom"`` for horizontal legend below plot.

**Returns**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Type
     - Description
   * - ``altair.Chart``
     - Altair chart object (or None if saved).

plot_barplot_TS_results
^^^^^^^^^^^^^^^^^^^^^^^

Create grouped bar plot showing reference hit recovery by method and cycle.

.. list-table:: Parameters
   :header-rows: 1
   :widths: 18 15 10 57

   * - Parameter
     - Type
     - Required
     - Description
   * - ``width``
     - ``int``
     - No
     - Plot width in pixels (auto-calculated if None).
   * - ``height``
     - ``int``
     - No
     - Plot height in pixels. Default: 400.
   * - ``save_path``
     - ``str``
     - No
     - Path to save (.html, .png, .svg).
   * - ``show_plot``
     - ``bool``
     - No
     - Display in Jupyter. Default: True.
   * - ``legend_position``
     - ``str``
     - No
     - Position of legend: ``"right"`` (default) or ``"bottom"`` for horizontal legend below plot.
   * - ``dark_mode``
     - ``bool``
     - No
     - Use white text for bar labels (for dark backgrounds). Default: False.

plot_line_performance_with_error_bars
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create line plot with error bars showing mean performance across top-N cutoffs.

.. list-table:: Parameters
   :header-rows: 1
   :widths: 18 15 10 57

   * - Parameter
     - Type
     - Required
     - Description
   * - ``width``
     - ``int``
     - No
     - Plot width in pixels. Default: 800.
   * - ``height``
     - ``int``
     - No
     - Plot height in pixels. Default: 500.
   * - ``save_path``
     - ``str``
     - No
     - Path to save (.html, .png, .svg).
   * - ``show_plot``
     - ``bool``
     - No
     - Display in Jupyter. Default: True.
   * - ``legend_position``
     - ``str``
     - No
     - Position of legend: ``"right"`` (default) or ``"bottom"`` for horizontal legend below plot.

Other Methods
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Method
     - Description
   * - ``gen_TS_runs_data(top_n, sort_type)``
     - Internal: Generate combined datasets.
   * - ``get_barplot_TS_results_data(top_n)``
     - Internal: Generate bar plot data.
   * - ``gen_line_plot_performance_data(top_ns)``
     - Internal: Generate line plot data.

Complete Example
~~~~~~~~~~~~~~~~

.. code-block:: python
   :caption: Full benchmarking workflow

   import polars as pl
   from TACTICS.library_analysis.visualization import TS_Benchmarks

   # Load results from multiple runs
   rw_runs = [pl.read_csv(f"rw_cycle_{i}.csv") for i in range(10)]
   ucb_runs = [pl.read_csv(f"ucb_cycle_{i}.csv") for i in range(10)]
   greedy_runs = [pl.read_csv(f"greedy_cycle_{i}.csv") for i in range(10)]

   # Load reference data
   reference = pl.read_csv("reference_scores.csv")

   # Create benchmarks (all data generated automatically)
   benchmarks = TS_Benchmarks(
       no_of_cycles=10,
       methods_list=["RouletteWheel", "BayesUCB", "Greedy"],
       TS_runs_data={
           "RouletteWheel": rw_runs,
           "BayesUCB": ucb_runs,
           "Greedy": greedy_runs,
       },
       reference_data=reference,
       top_n=100,
       sort_type="minimize",
       top_ns=[25, 50, 100, 200, 300]
   )

   # Generate all visualizations
   strip_chart = benchmarks.stripplot_TS_results(
       width=800, height=500, save_path="strip_plot.html",
       legend_position="right"  # or "bottom" for horizontal legend
   )

   bar_chart = benchmarks.plot_barplot_TS_results(
       width=700, height=400, save_path="bar_plot.html",
       legend_position="bottom",  # horizontal legend below plot
       dark_mode=False  # set True for white text on dark backgrounds
   )

   line_chart = benchmarks.plot_line_performance_with_error_bars(
       width=900, height=600, save_path="line_plot.html",
       legend_position="right"
   )

   # Access computed statistics
   print(f"Methods analyzed: {benchmarks.actual_methods}")


.. _diagnostic-plots:

Diagnostic Plots
----------------

.. rst-class:: class-core

``TACTICS.library_analysis.diagnostic_plots`` turns the diagnostics DataFrame
from ``sampler.get_diagnostics()`` (and, for the reagent-usage panels, the
search results DataFrame) into matplotlib figures. Every function accepts
Polars DataFrames and returns a ``matplotlib.figure.Figure``.

**Requires:** ``track_diagnostics=True`` in ``ThompsonSamplingConfig``.

Single-mechanism trajectories (RWS / CATS)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Function
     - Description
   * - ``plot_criticality_trajectory(diagnostics_df, figsize=(10, 4))``
     - Per-component criticality over cycles; the diffuse band (≤ 0.3) is shaded.
   * - ``plot_snr_trajectory(diagnostics_df, figsize=(10, 4))``
     - Signal-to-noise evolution with a threshold line at SNR = 1.
   * - ``plot_temperature_decomposition(diagnostics_df, component_idx, figsize=(10, 6))``
     - Full temperature pipeline (base temperature, CATS multiplier, final
       temperature) for one component.

Combined method diagnostics
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Function
     - Description
   * - ``plot_rws_diagnostic(diagnostics_df, title="", figsize=(12, 5))``
     - GMIC with confidence bands on the left axis, CATS multiplier on the
       right, and the divergence gate.
   * - ``plot_ttts_diagnostic(diagnostics_df, title="", figsize=(12, 5))``
     - Disagreement EMA on the left axis, adaptive ``heated_scale`` on the right.

Reagent-usage and layered mechanism panels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These take the search results DataFrame, the diagnostics DataFrame, and the
reagent files so they can map product names back to components.

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Function
     - Description
   * - ``plot_reagent_usage_action_panel(search_df, diag_df, reagent_files, batch_size=100, ...)``
     - Two panels: heating timeline, and per-batch fraction of each component's
       reagents sampled (revisits solid, new hatched), with heated-component
       opacity encoding.
   * - ``plot_gmic_directed_exploration(search_df, diag_df, reagent_files, methods, *, ...)``
     - Layer 1 — how GMIC directs *which* component to explore
       (signal → decision → result), one row per method.
   * - ``plot_adaptive_intensity(search_df, diag_df, reagent_files, methods, *, ...)``
     - Layer 2 — how each method tunes exploration *within* a component
       (RWS ``cats_multiplier`` vs TT-TS ``heated_scale``).

**Example**

.. code-block:: python

   from TACTICS.library_analysis.diagnostic_plots import plot_rws_diagnostic

   # After running search with track_diagnostics=True
   diag_df = sampler.get_diagnostics()
   fig = plot_rws_diagnostic(diag_df, title="Thrombin — recommended_rws")
   fig.savefig("diagnostics.png", dpi=150, bbox_inches="tight")
