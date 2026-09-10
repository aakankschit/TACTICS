.. _guide-scoring:

2. Scoring
==========

**What this block is.** A product goes in, a number comes out. The sampler
only ever sees that number; it never knows whether it came from a docking
run, a lookup table or a line of Python. You describe the scorer with a
small config object, and the sampler builds (and, in parallel runs,
rebuilds) the evaluator from it.

Core
----

The quickest scorer is a table of precomputed values keyed by product name.
The bundled dataset ships one — docking scores for all 500 k thrombin
amides:

.. code-block:: python

   from TACTICS.thompson_sampling import LookupEvaluatorConfig

   evaluator = LookupEvaluatorConfig(
       ref_filename=str(data / "product_scores.parquet"),   # .csv or .parquet
       compound_col="Product_Code",                          # default
       score_col="Scores",                                   # default
   )

Two decisions travel with every scorer:

- **Direction.** ``mode="minimize"`` for docking scores and free energies
  (lower is better); ``mode="maximize"`` for probabilities, similarities,
  read counts. It is set on the *strategy* — via ``get_preset(mode=...)`` or
  ``TopTwoConfig(mode=...)`` — because it is the search that needs to know.
- **What "missing" means.** A product the table does not contain scores
  ``NaN`` and is skipped. For sparse DEL-style tables where absence means
  *not a binder*, set ``default_score=0.0`` so absence is a real observation.

**What it produces:** an ``*EvaluatorConfig`` that you pass as
``evaluator_config=`` in Block 3.

Build on it
-----------

Your own Python function
~~~~~~~~~~~~~~~~~~~~~~~~

Most real screens score with code you already have. Wrap any callable that
takes an RDKit ``Mol`` and returns a ``float``:

.. literalinclude:: /snippets/scoring_custom_function.py
   :language: python
   :start-after: # [start:function]
   :end-before: # [end:function]

.. literalinclude:: /snippets/scoring_custom_function.py
   :language: python
   :start-after: # [start:run]
   :end-before: # [end:run]
   :dedent:

Scores are cached by canonical SMILES; an exception inside your function
yields ``NaN`` for that product and the run continues. Two rules:

- Define the function at **module level** (not a lambda, not inside another
  function). With ``processes > 1`` each worker rebuilds the evaluator from
  the config, which means pickling a reference to your function.
- Expensive setup (loading a model, reading a receptor) belongs in a
  module-level cache, not in the function body — it runs once per product.

The other built-in scorers
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 40 10 20

   * - Config
     - Scores
     - Cost
     - Needs
   * - :class:`~TACTICS.thompson_sampling.core.evaluator_config.LookupEvaluatorConfig`
     - value from a CSV/Parquet table, by product name
     - instant
     - the table
   * - :class:`~TACTICS.thompson_sampling.core.evaluator_config.DBEvaluatorConfig`
     - value from a ``sqlitedict`` database, by product name
     - instant
     - the ``.db`` file
   * - :class:`~TACTICS.thompson_sampling.core.evaluator_config.CustomEvaluatorConfig`
     - your function
     - yours
     - a picklable callable
   * - :class:`~TACTICS.thompson_sampling.core.evaluator_config.FPEvaluatorConfig`
     - Morgan-fingerprint Tanimoto to a query SMILES
     - fast
     - a query SMILES
   * - :class:`~TACTICS.thompson_sampling.core.evaluator_config.MWEvaluatorConfig`
     - molecular weight (a smoke test)
     - fast
     - nothing
   * - :class:`~TACTICS.thompson_sampling.core.evaluator_config.MLClassifierEvaluatorConfig`
     - positive-class probability from a pickled scikit-learn model on a Morgan FP
     - fast
     - the ``.pkl``
   * - :class:`~TACTICS.thompson_sampling.core.evaluator_config.ROCSEvaluatorConfig`
     - shape + colour overlay to a 3D query
     - slow
     - ``[openeye]``, licence
   * - :class:`~TACTICS.thompson_sampling.core.evaluator_config.FredEvaluatorConfig`
     - FRED docking score into an ``.oedu`` receptor
     - slow
     - ``[openeye]``, licence

The two lookup scorers key on the product **name** and skip synthesis
altogether, so ``SMILES`` in the results reads ``FAIL`` for them; the others
need the structure and populate it. For the slow ones, see :doc:`04_scale`.

Using an evaluator without a config
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every config has a class behind it (``LookupEvaluatorConfig`` →
:class:`~TACTICS.thompson_sampling.core.evaluators.LookupEvaluator`) that
takes the same keys as a dict. Construct one directly when driving the
sampler by hand (see :doc:`03_search`) or when
you just want to score a molecule:

.. code-block:: python

   from TACTICS import FPEvaluator
   from rdkit import Chem

   scorer = FPEvaluator({"query_smiles": "CC(=O)Nc1ccc(O)cc1"})
   scorer.evaluate(Chem.MolFromSmiles("CC(=O)Nc1ccc(OC)cc1"))   # 0.593

.. admonition:: Gotchas
   :class: warning

   - **Direction lives on the strategy, not the evaluator.** A docking table
     with the default ``mode="maximize"`` will happily find you the *worst*
     binders.
   - **Lookup tables are keyed by exact product name** — ``<acid>_<amine>``
     built from the reagent-file names. Rename a reagent and every score for
     it goes missing (→ ``NaN``, skipped; or ``default_score``).
   - **A lambda or nested function fails with** ``processes > 1``. Module
     level only.
   - **``FPEvaluatorConfig`` has no radius/bits knobs** — it is fixed at
     radius 2, 2048 bits. Write a ``CustomEvaluatorConfig`` for anything else.
   - Passing a key the config does not have is silently ignored for evaluator
     configs (they are permissive); check the field list in the reference.

Reference
---------

:doc:`/reference/scoring` — every config and evaluator class with its keys.

Interactive: ``marimo edit tutorials/custom_evaluator_tester.py`` lets you
paste a scoring function and run it against a live search.

**Next:** :doc:`03_search` — put the library and the scorer together and run.
