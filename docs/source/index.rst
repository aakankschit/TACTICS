.. TACTICS documentation master file

.. image:: _static/images/TACTICS_logo.png
   :alt: TACTICS
   :align: center
   :width: 420px

TACTICS
=======

**Thompson Sampling over combinatorial libraries: define, score, search.**

You have reagent files, a reaction, and a way to score a product. TACTICS
finds the best products in a library of millions by scoring a few thousand
— maintaining a posterior per reagent and choosing what to make next from
what it has learned. Six building blocks, composed in order:

.. grid:: 1 2 3 3
   :gutter: 3

   .. grid-item-card:: 1 · Library
      :link: guides/01_library
      :link-type: doc

      A reaction SMARTS and one reagent file per component →
      ``SynthesisPipeline``. Validate coverage; multi-step; deprotection.

   .. grid-item-card:: 2 · Scoring
      :link: guides/02_scoring
      :link-type: doc

      A product in, a number out. A lookup table, your own Python
      function, fingerprints, docking, a model.

   .. grid-item-card:: 3 · Search
      :link: guides/03_search
      :link-type: doc

      ``get_preset → from_config → warm_up → search``. Presets, hand-built
      configs, choosing a strategy and warmup.

   .. grid-item-card:: 4 · Scale
      :link: guides/04_scale
      :link-type: doc

      Worker processes for slow scorers; pre-enumerated products;
      budgeting a run on a cluster.

   .. grid-item-card:: 5 · Inspect
      :link: guides/05_inspect
      :link-type: doc

      What the search learned: per-cycle diagnostics, the posterior
      landscape, a SAR summary, analysis functions.

   .. grid-item-card:: 6 · Visualise
      :link: guides/06_visualise
      :link-type: doc

      Recovery benchmarks (Altair) and mechanism plots (matplotlib).
      ``pip install chem-tactics[viz]``.

.. grid:: 1 2 4 4
   :gutter: 3

   .. grid-item-card:: Extend
      :link: extending
      :link-type: doc

      Your own strategy, warmup or evaluator class.

   .. grid-item-card:: Reference
      :link: reference/index
      :link-type: doc

      Every class and function, generated from the code.

   .. grid-item-card:: Theory
      :link: theory/index
      :link-type: doc

      The mathematics: Thompson Sampling, CATS, TT-TS, Bayes-UCB.

   .. grid-item-card:: Tutorials
      :link: tutorials
      :link-type: doc

      Interactive marimo notebooks.

Ten minutes
-----------

.. code-block:: bash

   pip install chem-tactics

Then run the bundled thrombin library (130 acids × 3,844 amines) against
its precomputed docking scores — this is also the README quickstart, and
the test suite executes it:

.. literalinclude:: /snippets/search_minimal.py
   :language: python
   :start-after: # [start:search]
   :end-before: # [end:search]
   :dedent:

``results`` is a Polars DataFrame of every product scored — ``score``,
``SMILES``, ``Name``. From here, :doc:`guides/01_library` is where you swap
in your own reagents and reaction, and :doc:`guides/02_scoring` is where
you swap in your own scorer.

Install extras
--------------

.. code-block:: bash

   pip install "chem-tactics[viz]"        # matplotlib + altair: Block 6
   pip install "chem-tactics[tutorials]"  # marimo + viz
   pip install "chem-tactics[openeye]"    # ROCS / FRED evaluators (licence required)
   pip install -e ".[test,docs]"          # development

Requires Python 3.11+. Core dependencies are RDKit, NumPy, SciPy, Polars,
Pydantic, tqdm, sqlitedict and joblib; everything else is optional and
imported lazily.

.. toctree::
   :hidden:
   :caption: Guides

   guides/01_library
   guides/02_scoring
   guides/03_search
   guides/04_scale
   guides/05_inspect
   guides/06_visualise
   extending
   tutorials

.. toctree::
   :hidden:
   :caption: Reference

   reference/index

.. toctree::
   :hidden:
   :caption: Theory

   theory/index

.. toctree::
   :hidden:
   :caption: Project

   changelog
   GitHub <https://github.com/aakankschit/TACTICS>
