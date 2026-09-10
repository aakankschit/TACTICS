Tutorials
=========

Interactive `marimo <https://marimo.io>`_ notebooks under ``tutorials/`` in
the repository. They are not rendered here; run one with

.. code-block:: bash

   pip install "chem-tactics[tutorials]"
   marimo edit tutorials/<name>.py      # or: marimo run ... for app mode

Working with the package
------------------------

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Thompson Sampling benchmark
      :link: https://github.com/aakankschit/TACTICS/blob/main/tutorials/thompson_sampling_tutorial.py

      ``thompson_sampling_tutorial.py`` — pick strategies and warmups, run
      them on the bundled thrombin library, and compare recovery with
      ``TS_Benchmarks``. Blocks 1–3, 6.

   .. grid-item-card:: Library enumeration
      :link: https://github.com/aakankschit/TACTICS/blob/main/tutorials/library_enumeration_tutorial.py

      ``library_enumeration_tutorial.py`` — single-step, multi-step and
      alternative-SMARTS pipelines, full enumeration, then a search.
      Block 1, with a taste of 3.

   .. grid-item-card:: Reaction config builder
      :link: https://github.com/aakankschit/TACTICS/blob/main/tutorials/reaction_config_builder.py

      ``reaction_config_builder.py`` — build a ``ReactionConfig`` in a form,
      validate it against your reagent files, see the failures. Block 1.

   .. grid-item-card:: Custom evaluator tester
      :link: https://github.com/aakankschit/TACTICS/blob/main/tutorials/custom_evaluator_tester.py

      ``custom_evaluator_tester.py`` — paste a scoring function, check it
      compiles and behaves, run a short search with it. Block 2.

   .. grid-item-card:: Diagnostic benchmark plots
      :link: https://github.com/aakankschit/TACTICS/blob/main/tutorials/diagnostic_benchmark_plots.py

      ``diagnostic_benchmark_plots.py`` — the mechanism plots
      (``plot_rws_diagnostic``, ``plot_ttts_diagnostic``, the layered
      panels) over the published diagnostic benchmark. Blocks 5–6.
      Needs the benchmark parquet under ``outputs/``.

   .. grid-item-card:: Interactive SAR explorer
      :link: https://github.com/aakankschit/TACTICS/blob/main/tutorials/interactive_sar_explorer.py

      ``interactive_sar_explorer.py`` — hover a point, see the structure;
      per-component oracle GMIC over the ground-truth scores. Needs
      ``data/scores``.

Manuscript figures
------------------

Reproduce the paper's figures from the benchmark output (not the package
API). Each expects the benchmark parquet files under ``outputs/``.

- ``manuscript_plots_ROCS.py`` — ROCS libraries: aggregate top-*N*
  recovery, per-library breakdowns, Tukey HSD, budget sensitivity.
- ``manuscript_plots_docking.py`` — the eight docking libraries: method ×
  library heatmap, significance, budget sensitivity.
- ``manuscript_sar_plots.py`` — reagent score landscapes from the
  brute-force ground truth.
