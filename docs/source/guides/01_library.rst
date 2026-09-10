.. _guide-library:

1. Library
==========

**What this block is.** A combinatorial library is a reaction plus one reagent
file per reactant position. TACTICS never enumerates the whole library up
front; it makes each product on demand from a tuple of reagents. This block
is how you tell it what "make a product" means.

Core
----

One reaction SMARTS and two ``.smi`` files describe the bundled thrombin
library (130 acids × 3,844 amines ≈ 500 k amides):

.. literalinclude:: /snippets/library_define.py
   :language: python
   :start-after: # [start:define]
   :end-before: # [end:define]
   :dedent:

**Reagent files** are whitespace-separated ``SMILES name`` lines, one reagent
per line. The *name* becomes part of every product's name
(``<acid>_<amine>``), which is how scores are keyed later — keep names unique
and free of underscores.

**Reaction SMARTS** map atoms across the arrow; the *order of reactants* in
the SMARTS is the order of files in ``reagent_file_list``.

Make one product to check the chemistry does what you think:

.. literalinclude:: /snippets/library_define.py
   :language: python
   :start-after: # [start:single]
   :end-before: # [end:single]
   :dedent:

And check every reagent against its template before you screen — an acid
file with 10 % non-acids silently wastes 10 % of your budget:

.. literalinclude:: /snippets/library_define.py
   :language: python
   :start-after: # [start:validate]
   :end-before: # [end:validate]
   :dedent:

:meth:`~TACTICS.library_enumeration.smarts_toolkit.config.ReactionDef.validate_reaction`
returns a :class:`~TACTICS.library_enumeration.smarts_toolkit._validator.ValidationResult`
with the incompatible reagents by position, unparseable SMILES, duplicates,
and reagents carrying protecting groups or salt fragments. Pass
``deprotect=True`` / ``desalt=True`` to re-check after cleaning them up.

**What it produces:** a :class:`~TACTICS.library_enumeration.synthesis_pipeline.SynthesisPipeline`.
Everything downstream takes it as ``synthesis_pipeline=``.

Build on it
-----------

Multi-step synthesis
~~~~~~~~~~~~~~~~~~~~

A step can consume the product of an earlier step. ``step_inputs`` says, for
each step, where each reactant comes from — a reagent file or a previous
step's product:

.. literalinclude:: /snippets/library_multistep.py
   :language: python
   :start-after: # [start:multistep]
   :end-before: # [end:multistep]
   :dedent:

The sampler treats every reagent *file* as a component, so this config has
three components regardless of how many steps use them.

Alternative SMARTS at one step
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Give a step several :class:`~TACTICS.library_enumeration.smarts_toolkit.config.ReactionDef`
entries (same ``step_index``) and set ``step_modes={0: "alternative"}``. The
pipeline tries each pattern in order until one applies — useful when one
reagent file mixes, say, primary and secondary amines that need different
templates. :meth:`~TACTICS.library_enumeration.synthesis_pipeline.SynthesisPipeline.auto_detect_compatibility`
(run for you by ``ThompsonSampler.from_config``) works out which pattern
each reagent matches so the search does not try the wrong one.

Protecting groups and salts
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Reagent files from vendors carry Boc/Fmoc/Cbz groups and counter-ions.
Two tools:

- **Detection** — ``validate_reaction(...)`` lists ``protected_reagents`` and
  ``multi_fragment_reagents``. Ten common groups are built in
  (:data:`~TACTICS.library_enumeration.smarts_toolkit.constants.DEFAULT_PROTECTING_GROUPS`);
  add your own with :class:`~TACTICS.library_enumeration.smarts_toolkit.config.ProtectingGroupInfo`
  and ``ReactionConfig(protecting_groups=[...])``.
- **Removal during synthesis** — a
  :class:`~TACTICS.library_enumeration.smarts_toolkit.config.DeprotectionSpec`
  on a ``ReactionDef`` strips a group from a reactant (``target=1``) or from
  the step's product (``target="product"``) before the next step, as in the
  multi-step example above.

Enumerate everything (no Thompson Sampling)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a library small enough to score exhaustively, or to produce a file for
another tool, the pipeline enumerates on its own:

.. literalinclude:: /snippets/library_enumerate_all.py
   :language: python
   :start-after: # [start:enumerate]
   :end-before: # [end:enumerate]
   :dedent:

``enumerate_library(n_jobs=8)`` parallelises across reagent combinations.
:func:`~TACTICS.library_enumeration.file_writer.write_enumerated_library`
writes ``csv``, ``smi`` or ``sdf``;
:func:`~TACTICS.library_enumeration.file_writer.write_products_chunked`
splits large outputs into numbered files.

.. admonition:: Gotchas
   :class: warning

   - **Reactant order = file order.** The first reactant in the SMARTS is
     ``reagent_file_list[0]``. Swapping them gives 0 % coverage, not an error.
   - **Product names come from reagent names.** ``LookupEvaluator`` and
     ``DBEvaluator`` key on ``<name1>_<name2>``; a reagent name containing an
     underscore breaks the join.
   - **SMARTS that parse are not SMARTS that match.** ``ReactionDef`` only
     checks that the SMARTS parses; run ``validate_reaction`` against the real
     files to see coverage.
   - **Salts and protecting groups are not removed unless you ask.** Use
     ``deprotect=True`` / ``desalt=True`` in validation and ``DeprotectionSpec``
     in synthesis.

Reference
---------

:doc:`/reference/library` — :class:`~TACTICS.library_enumeration.synthesis_pipeline.SynthesisPipeline`,
:class:`~TACTICS.library_enumeration.smarts_toolkit.config.ReactionDef`,
:class:`~TACTICS.library_enumeration.smarts_toolkit.config.ReactionConfig`,
:class:`~TACTICS.library_enumeration.smarts_toolkit.config.StepInput`,
:class:`~TACTICS.library_enumeration.smarts_toolkit.config.DeprotectionSpec`,
:class:`~TACTICS.library_enumeration.smarts_toolkit._validator.ValidationResult`,
:class:`~TACTICS.library_enumeration.enumeration_utils.EnumerationResult`.

Interactive: ``marimo edit tutorials/reaction_config_builder.py`` builds and
validates a config step by step; ``tutorials/library_enumeration_tutorial.py``
walks single-step, multi-step and alternative-SMARTS cases.

**Next:** :doc:`02_scoring` — how a product gets a number.
