Library Analysis Module
========================

The Library Analysis module provides tools for analyzing chemical libraries and their building blocks, 
as well as visualization capabilities for chemical analysis and Thompson Sampling benchmark results.

LibraryAnalysis Class
---------------------

The LibraryAnalysis class provides tools for analyzing chemical libraries and their building blocks,
including identification of top building blocks and overlap analysis between different cutoffs.

.. autoclass:: PRISMS.library_analysis.library_analysis.LibraryAnalysis
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:

**Class Methods:**

.. class:: LibraryAnalysis(polars_dataframe, top_n_cutoff=100)
   
   Initialize the LibraryAnalysis class.

   :param polars.DataFrame polars_dataframe: DataFrame containing library data with 'product_code' and 'score' columns
   :param int top_n_cutoff: Number of top products to consider for analysis

   .. method:: find_top_building_blocks()
      
      Identify and count the most frequently occurring building blocks in top products.
      
      :returns: None (results stored in class attributes)

   .. method:: check_overlap(new_cutoff)
      
      Check overlap between current cutoff and a new cutoff value.
      
      :param int new_cutoff: New cutoff value to compare against
      :returns: List of tuples containing overlap information for each position
      :rtype: list

   .. method:: visualize_top_building_blocks(show_overlap=False, mols_per_row=5, sub_img_size=(300, 300), comparison_analysis=None, top_n=20)
      
      Visualize top building blocks using RDKit molecular drawings.
      
      :param bool show_overlap: Whether to show overlapping building blocks
      :param int mols_per_row: Number of molecules per row in the grid
      :param tuple sub_img_size: Size of each molecule image
      :param LibraryAnalysis comparison_analysis: Another analysis instance for comparison
      :param int top_n: Number of top building blocks to visualize

   .. method:: compare_analysis_overlap(other_analysis, top_n=100)
      
      Compare building block overlap with another LibraryAnalysis instance.
      
      :param LibraryAnalysis other_analysis: Another LibraryAnalysis instance to compare with
      :param int top_n: Number of top building blocks to consider
      :returns: List of overlap results for each position
      :rtype: list

LibraryVisualization Class
---------------------------

The LibraryVisualization class provides tools for creating visualizations from LibraryAnalysis instances,
including chemical structure visualization and comparative analysis plots.

.. autoclass:: PRISMS.library_analysis.visualization.LibraryVisualization
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:

**Class Methods:**

.. class:: LibraryVisualization(analysis)
   
   Initialize the LibraryVisualization class.

   :param LibraryAnalysis analysis: LibraryAnalysis instance containing data to visualize

   .. method:: plot_top_products_comparison(analysis_instances, reference_instance, top_n=100, title="Top Products Comparison", figsize=(15, 10), analysis_labels=None, save_path=None)
      
      Generate subplots comparing overlap of top products between analysis instances and a reference.
      
      :param list analysis_instances: List of LibraryAnalysis instances to compare
      :param LibraryAnalysis reference_instance: Reference LibraryAnalysis instance
      :param int top_n: Number of top products to consider
      :param str title: Plot title
      :param tuple figsize: Figure size (width, height)
      :param list analysis_labels: Labels for each analysis group
      :param str save_path: Path to save the plot

   .. method:: visualize_top_building_blocks(show_overlap=False, mols_per_row=5, sub_img_size=(300, 300), comparison_analysis=None, top_n=20)
      
      Visualize top building blocks for each position using RDKit.
      
      :param bool show_overlap: Whether to show overlapping building blocks
      :param int mols_per_row: Number of molecules to display per row
      :param tuple sub_img_size: Size of each molecule image
      :param LibraryAnalysis comparison_analysis: Another LibraryAnalysis instance for overlap comparison
      :param int top_n: Number of top building blocks to consider

TS_Benchmarks Class
-------------------

The TS_Benchmarks class provides comprehensive tools for benchmarking and visualizing Thompson Sampling (TS) 
results across multiple cycles and search strategies. All required data is automatically generated during 
initialization, making it extremely user-friendly with a simple one-step setup process.

**Key Features:**

* **🤖 Automatic Data Generation**: All datasets generated during initialization - no manual steps required
* **🎨 Consistent Color Schemes**: Unified colors across all plot types for professional appearance
* **📊 Multi-cycle Analysis**: Compare different TS methods across multiple cycles with comprehensive statistics
* **📈 Reference Comparison**: Benchmark against ground truth reference data with detailed performance metrics
* **🎯 Multiple Visualization Types**: Strip plots, bar plots, and line plots with error bars
* **⚙️ Configurable Parameters**: Customizable top_n, sort_type, and analysis points via constructor
* **🧹 Clean Output**: Organized summaries with emojis and comprehensive statistics

**New Simplified Workflow:**

1. **One-step initialization**: ``TS_Benchmarks(...)`` - Everything generated automatically!
2. **Direct plotting**: Call any plotting method immediately
   - Distribution analysis: ``stripplot_TS_results()``
   - Hit recovery: ``plot_barplot_TS_results()``
   - Performance trends: ``plot_line_performance_with_error_bars()``

**Major Improvements:**

* **🚀 Streamlined Usage**: No more manual data generation calls - everything happens in ``__init__``
* **🔧 Enhanced Constructor**: Support for custom ``top_n``, ``sort_type``, and ``top_ns`` parameters
* **📊 Pre-calculated Statistics**: Grouped statistics and error bar data generated during initialization
* **🎨 Color Consistency**: Same method gets same color across stripplot, bar plot, and line plot
* **✨ Better User Experience**: Clear progress indicators and comprehensive summaries

.. autoclass:: PRISMS.library_analysis.visualization.TS_Benchmarks
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:

**Class Methods:**

.. class:: TS_Benchmarks(no_of_cycles, methods_list, TS_runs_data, reference_data=None, top_n=100, sort_type="minimize", top_ns=None)
   
   Initialize the TS_Benchmarks class and automatically generate all required data for Thompson Sampling analysis.
   
   **🔄 Automatic Data Generation**: This constructor automatically generates all datasets including:
   TS runs data, bar plot data, line plot data, and grouped statistics for error bars.

   :param int no_of_cycles: Number of cycles to analyze
   :param list methods_list: List of method names (search strategies)
   :param dict TS_runs_data: Dictionary mapping method names to lists of DataFrames (one per cycle)
   :param DataFrame reference_data: Optional reference/ground truth data for comparison
   :param int top_n: Number of top products to consider for bar plot analysis (default: 100)
   :param str sort_type: Type of sorting ("minimize" or "maximize", default: "minimize")
   :param list top_ns: List of top N values for line plot analysis (default: [50, 100, 200, 300, 400, 500])

   .. method:: stripplot_TS_results(width=None, height=None, save_path=None, show_plot=True)
      
      🎨 Generate a strip plot showing score distributions across cycles and methods.
      **Data is pre-generated during initialization** - ready to plot immediately!
      
      :param int width: Plot width in pixels (auto-calculated if None)
      :param int height: Plot height in pixels (auto-calculated if None)
      :param str save_path: Path to save the plot (supports .html, .png, .svg)
      :param bool show_plot: Whether to display the plot in Jupyter
      :returns: Altair chart object or None
      :rtype: altair.Chart or None

   .. method:: plot_barplot_TS_results(width=None, height=None, save_path=None, show_plot=True)
      
      📊 Create a grouped bar plot showing reference hit recovery by method and cycle.
      **Data is pre-generated during initialization** - ready to plot immediately!
      
      :param int width: Plot width in pixels (auto-calculated if None)
      :param int height: Plot height in pixels (default: 400)
      :param str save_path: Path to save the plot (supports .html, .png, .svg)
      :param bool show_plot: Whether to display the plot in Jupyter
      :returns: Altair chart object or None
      :rtype: altair.Chart or None

   .. method:: plot_line_performance_with_error_bars(width=None, height=None, save_path=None, show_plot=True)
      
      📈 Create line plot with error bars showing mean performance trends across different top-N cutoffs.
      **Data and grouped statistics are pre-generated during initialization** - ready to plot immediately!
      
      :param int width: Plot width in pixels (default: 800)
      :param int height: Plot height in pixels (default: 500)
      :param str save_path: Path to save the plot (supports .html, .png, .svg)
      :param bool show_plot: Whether to display the plot in Jupyter
      :returns: Altair chart object or None
      :rtype: altair.Chart or None

**Internal/Advanced Methods:**

These methods are automatically called during initialization but can be accessed for advanced usage:

   .. method:: gen_TS_runs_data(top_n=100, sort_type="minimize")
      
      🔧 **Internal method** - Generate combined datasets from all TS runs.
      **Automatically called during initialization.**
      
      :param int top_n: Number of top products to consider per method/cycle
      :param str sort_type: "minimize" (ascending) or "maximize" (descending) sorting
      :returns: Tuple of (combined_df_top_n, combined_df_all)
      :rtype: tuple

   .. method:: get_barplot_TS_results_data(top_n=100)
      
      🔧 **Internal method** - Generate data for bar plot showing fraction of reference compounds found.
      **Automatically called during initialization.**
      
      :param int top_n: Number of top reference compounds to use as benchmark
      :returns: DataFrame with hit counts for bar plotting
      :rtype: polars.DataFrame

   .. method:: gen_line_plot_performance_data(top_ns=None)
      
      🔧 **Internal method** - Generate performance data for line plots across different top_n cutoffs.
      **Automatically called during initialization.**
      
      :param list top_ns: List of top_n values to test (default: [50,100,200,300,400,500])
      :returns: DataFrame with performance fractions
      :rtype: polars.DataFrame

   .. method:: get_performance_summary()
      
      Get summary of all processed performance data and chart components.
      
      :returns: Dictionary containing all stored data and charts
      :rtype: dict

**Automatic Data Storage:**

🤖 **All data is automatically generated and stored during initialization**:

* ``combined_df_top_n``: Top N compounds from each method/cycle (ready for stripplot)
* ``combined_df_all``: All compounds from each method/cycle (complete dataset)
* ``bar_plot_df``: Hit recovery data for bar plots (pre-calculated)
* ``line_plot_df``: Raw performance data across cycles (individual cycle data)
* ``grouped_stats``: Statistical summaries with mean, std, and error bounds (for line plot)
* ``grouped_stats_caps``: Error bar cap positions (for error bar styling)
* ``unique_top_ns``: Sorted list of top-N values analyzed
* ``actual_methods``: Methods found in data (for validation)

**Color Consistency:**

🎨 **Unified color scheme across all plots** - each method gets the same color in stripplot, bar plot, and line plot via the internal ``_get_color_scheme()`` method.

**New Simplified Usage Pattern:**

.. code-block:: python

   # ✨ NEW: One-step initialization with automatic data generation
   ts_benchmarks = TS_Benchmarks(
       no_of_cycles=10,
       methods_list=['thompson_standard', 'thompson_boltzmann', 'random_baseline'],
       TS_runs_data={
           'thompson_standard': [cycle1_df, cycle2_df, ...], 
           'thompson_boltzmann': [cycle1_df, cycle2_df, ...],
           'random_baseline': [cycle1_df, cycle2_df, ...]
       },
       reference_data=reference_df,
       top_n=100,  # Custom analysis size
       sort_type="minimize",  # Lower scores = better
       top_ns=[25, 50, 100, 200]  # Custom line plot points
   )
   # 🎉 All data automatically generated above!

   # 🚀 Ready to plot immediately - no additional data generation needed!
   strip_plot = ts_benchmarks.stripplot_TS_results(width=800, height=500)
   bar_plot = ts_benchmarks.plot_barplot_TS_results(width=700, height=400)
   line_plot = ts_benchmarks.plot_line_performance_with_error_bars(width=900, height=600)

   # 📊 Access all pre-calculated data
   print(f"Methods analyzed: {ts_benchmarks.actual_methods}")
   print(f"Performance range: {ts_benchmarks.grouped_stats['mean'].min():.3f} to {ts_benchmarks.grouped_stats['mean'].max():.3f}")
   
   # 🎯 Get comprehensive summary
   summary = ts_benchmarks.get_performance_summary()

**Migration from Old Usage:**

.. code-block:: python

   # ❌ OLD: Multi-step process (no longer needed)
   # ts_benchmarks = TS_Benchmarks(...)
   # ts_benchmarks.gen_TS_runs_data(top_n=100)
   # ts_benchmarks.get_barplot_TS_results_data(top_n=100)
   # ts_benchmarks.gen_line_plot_performance_data(top_ns=[50,100,200])
   # ts_benchmarks.plot_line_performance_with_error_bars()

   # ✅ NEW: Single-step process
   ts_benchmarks = TS_Benchmarks(..., top_n=100, top_ns=[50,100,200])
   ts_benchmarks.plot_line_performance_with_error_bars()  # Ready immediately! 