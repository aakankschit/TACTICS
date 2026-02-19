Library Analysis
================

The Library Analysis module provides tools for analyzing Thompson Sampling results,
identifying top building blocks, and creating visualizations for benchmarking studies.

Module Architecture
-------------------

.. graphviz::

    digraph LibraryAnalysis {
        rankdir=TB;
        node [shape=box, style="rounded,filled", fontname="Helvetica"];

        subgraph cluster_input {
            label="Input Data";
            style=filled;
            color=lightblue;

            Results [label="Results DataFrame\n(Polars)", shape=cylinder, fillcolor=lightblue];
            Reference [label="Reference Data\n(Ground Truth)", shape=cylinder, fillcolor=lightblue];
        }

        subgraph cluster_analysis {
            label="Analysis";
            style=filled;
            color=lightyellow;

            LibraryAnalysis [label="LibraryAnalysis\n- Top building blocks\n- Overlap analysis", fillcolor=gold];
        }

        subgraph cluster_viz {
            label="Visualization";
            style=filled;
            color=lightgreen;

            LibraryVisualization [label="LibraryVisualization\n- Structure grids\n- Comparison plots", fillcolor=lightgreen];
            TS_Benchmarks [label="TS_Benchmarks\n- Strip plots\n- Bar plots\n- Line plots", fillcolor=lightgreen];
        }

        subgraph cluster_output {
            label="Outputs";
            style=filled;
            color=lavender;

            Plots [label="Interactive Plots\n(Altair)", fillcolor=lavender];
            MolGrids [label="Molecule Grids\n(RDKit)", fillcolor=lavender];
            Stats [label="Statistics\n& Summaries", fillcolor=lavender];
        }

        Results -> LibraryAnalysis;
        Reference -> TS_Benchmarks;
        LibraryAnalysis -> LibraryVisualization;
        LibraryVisualization -> Plots;
        LibraryVisualization -> MolGrids;
        TS_Benchmarks -> Plots;
        TS_Benchmarks -> Stats;
    }

Quick Start
-----------

**Analyze results and find top building blocks:**

.. code-block:: python
   :caption: Basic analysis workflow

   import polars as pl
   from TACTICS.library_analysis import LibraryAnalysis, LibraryVisualization

   # Load results
   results_df = pl.read_csv("results.csv")

   # Create analysis with SMILES files for building block lookup
   analysis = LibraryAnalysis(
       df=results_df,
       smiles_files=["acids.smi", "amines.smi"],
       product_code_column="Product_Code",
       score_column="Scores"
   )
   
   # Find top building blocks in top 100 compounds
   counters, total = analysis.find_top_building_blocks(cutoff=100)

   # Visualize top building blocks
   viz = LibraryVisualization(analysis)
   viz.visualize_top_building_blocks(top_n=20)

**Benchmark Thompson Sampling methods:**

.. code-block:: python
   :caption: Benchmarking multiple methods

   from TACTICS.library_analysis.visualization import TS_Benchmarks

   # TS_Benchmarks auto-generates all data during initialization
   benchmarks = TS_Benchmarks(
       no_of_cycles=10,
       methods_list=["roulette_wheel", "bayes_ucb", "greedy"],
       TS_runs_data={
           "roulette_wheel": [cycle1_df, cycle2_df, ...],
           "bayes_ucb": [cycle1_df, cycle2_df, ...],
           "greedy": [cycle1_df, cycle2_df, ...],
       },
       reference_data=reference_df,
       top_n=100,
       sort_type="minimize"
   )

   # All data pre-calculated - ready to plot immediately
   benchmarks.stripplot_TS_results()
   benchmarks.plot_barplot_TS_results()
   benchmarks.plot_line_performance_with_error_bars()


.. _library-analysis:

LibraryAnalysis
---------------

.. rst-class:: class-core

Analyzes chemical libraries to identify top building blocks and their overlap.

.. admonition:: Dependencies
   :class: dependencies

   - Polars or Pandas DataFrame with product codes and scores
   - SMILES files (.smi) containing building block structures
   - Used by :ref:`LibraryVisualization <library-visualization>` for visualizations

Constructor
~~~~~~~~~~~

.. list-table:: Parameters
   :header-rows: 1
   :widths: 22 20 10 48

   * - Parameter
     - Type
     - Required
     - Description
   * - ``df``
     - ``DataFrame``
     - Yes
     - Polars or Pandas DataFrame with product codes and scores.
   * - ``smiles_files``
     - ``str | list[str]``
     - Yes
     - Path(s) to .smi file(s) containing SMILES and building block codes.
   * - ``product_code_column``
     - ``str``
     - No
     - Name of product code column. Default: ``"Product_Code"``.
   * - ``score_column``
     - ``str``
     - No
     - Name of score column. Default: ``"Scores"``.

Methods
~~~~~~~

find_top_building_blocks
^^^^^^^^^^^^^^^^^^^^^^^^

Identify and count building blocks in top products.

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
   * - ``cutoff``
     - ``int``
     - Yes
     - Number of top products to analyze.
   * - ``sort_scores_by``
     - ``str``
     - No
     - ``"ascending"`` or ``"descending"``. Default: ``"ascending"``.
   * - ``top_n``
     - ``int``
     - No
     - Number of top building blocks per position. Default: 100.

**Returns**

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Type
     - Description
   * - ``Tuple[List[Counter], int]``
     - (position_counters, total_molecules)

**Example**

.. code-block:: python

   from TACTICS.library_analysis import LibraryAnalysis

   analysis = LibraryAnalysis(
       df=results_df,
       smiles_files=["acids.smi", "amines.smi"],
       product_code_column="Product_Code",
       score_column="Scores"
   )
   counters, total = analysis.find_top_building_blocks(cutoff=100)

check_overlap
^^^^^^^^^^^^^

Check overlap between current cutoff and a new cutoff value.

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 15 10 55

   * - Parameter
     - Type
     - Required
     - Description
   * - ``new_cutoff``
     - ``int``
     - Yes
     - New cutoff value to compare against.

**Returns**

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Type
     - Description
   * - ``list[tuple]``
     - Overlap information for each position.

compare_analysis_overlap
^^^^^^^^^^^^^^^^^^^^^^^^

Compare building block overlap with another LibraryAnalysis instance.

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 20 10 50

   * - Parameter
     - Type
     - Required
     - Description
   * - ``other_analysis``
     - ``LibraryAnalysis``
     - Yes
     - Another analysis to compare with.
   * - ``top_n``
     - ``int``
     - No
     - Number of top building blocks. Default: 100.

**Returns**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Type
     - Description
   * - ``list``
     - Overlap results for each position.


.. _library-visualization:

LibraryVisualization
--------------------

.. rst-class:: class-core

Creates visualizations from LibraryAnalysis instances.

.. admonition:: Dependencies
   :class: dependencies

   - Requires :ref:`LibraryAnalysis <library-analysis>` instance
   - RDKit for molecular structure rendering

**Depends on:** :ref:`LibraryAnalysis <library-analysis>`

Constructor
~~~~~~~~~~~

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 20 10 50

   * - Parameter
     - Type
     - Required
     - Description
   * - ``analysis``
     - ``LibraryAnalysis``
     - Yes
     - Analysis instance containing data to visualize.

Methods
~~~~~~~

visualize_top_building_blocks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Display top building blocks using RDKit molecular drawings.

.. list-table:: Parameters
   :header-rows: 1
   :widths: 24 18 10 48

   * - Parameter
     - Type
     - Required
     - Description
   * - ``show_overlap``
     - ``bool``
     - No
     - Show overlapping building blocks. Default: False.
   * - ``mols_per_row``
     - ``int``
     - No
     - Molecules per row in grid. Default: 5.
   * - ``sub_img_size``
     - ``tuple``
     - No
     - Size of each molecule image. Default: (300, 300).
   * - ``comparison_analysis``
     - ``LibraryAnalysis``
     - No
     - Another analysis for overlap comparison.
   * - ``top_n``
     - ``int``
     - No
     - Number of top building blocks. Default: 20.

plot_top_products_comparison
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Generate subplots comparing overlap between analysis instances and a reference.

.. list-table:: Parameters
   :header-rows: 1
   :widths: 22 20 10 48

   * - Parameter
     - Type
     - Required
     - Description
   * - ``analysis_instances``
     - ``list``
     - Yes
     - List of LibraryAnalysis instances to compare.
   * - ``reference_instance``
     - ``LibraryAnalysis``
     - Yes
     - Reference analysis instance.
   * - ``top_n``
     - ``int``
     - No
     - Top products to consider. Default: 100.
   * - ``title``
     - ``str``
     - No
     - Plot title.
   * - ``figsize``
     - ``tuple``
     - No
     - Figure size (width, height). Default: (15, 10).
   * - ``analysis_labels``
     - ``list``
     - No
     - Labels for each analysis group.
   * - ``save_path``
     - ``str``
     - No
     - Path to save the plot.


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
   * - ``get_performance_summary()``
     - Get dict containing all stored data and charts.
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
   summary = benchmarks.get_performance_summary()


Workflow Overview
-----------------

.. graphviz::

    digraph Workflow {
        rankdir=TB;
        node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=10];
        nodesep=0.3;
        ranksep=0.4;

        TSRun [label="Thompson Sampling Results", fillcolor="#ADD8E6"];
        Analysis [label="LibraryAnalysis\n(find top BBs)", fillcolor="#FFD700"];
        Viz [label="LibraryVisualization\n(create plots)", fillcolor="#90EE90"];
        Benchmark [label="TS_Benchmarks\n(compare methods)", fillcolor="#E6E6FA"];

        TSRun -> Analysis;
        Analysis -> Viz;
        TSRun -> Benchmark;
    }

Typical analysis workflow:

1. **Run Thompson Sampling** - Generate results DataFrames
2. **LibraryAnalysis** - Identify top building blocks and their frequencies
3. **LibraryVisualization** - Create molecular structure grids
4. **TS_Benchmarks** - Compare multiple methods with statistical analysis
