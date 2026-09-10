Visualise
=========

Block 6 — charts. See the :doc:`guide </guides/06_visualise>`. Requires
``pip install chem-tactics[viz]``.

Benchmark comparison (Altair)
-----------------------------

.. autoclass:: TACTICS.library_analysis.visualization.TS_Benchmarks
   :members: stripplot_TS_results, plot_barplot_TS_results, plot_line_performance_with_error_bars, gen_TS_runs_data, get_barplot_TS_results_data, gen_line_plot_performance_data

Diagnostic plots (matplotlib)
-----------------------------

.. automodule:: TACTICS.library_analysis.diagnostic_plots
   :members: plot_criticality_trajectory, plot_snr_trajectory, plot_temperature_decomposition, plot_rws_diagnostic, plot_ttts_diagnostic, plot_reagent_usage_action_panel, plot_gmic_directed_exploration, plot_adaptive_intensity
