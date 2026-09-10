Library
=======

Block 1 — describe the combinatorial library. See the :doc:`guide </guides/01_library>`.

Pipeline
--------

.. autoclass:: TACTICS.library_enumeration.synthesis_pipeline.SynthesisPipeline
   :members:

Reaction definition
-------------------

.. automodule:: TACTICS.library_enumeration.smarts_toolkit.config
   :members: ReactionDef, ReactionConfig, StepInput, DeprotectionSpec, InputSource, ProtectingGroupInfo

Validation result
-----------------

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit._validator.ValidationResult
   :members:

Enumeration results and helpers
-------------------------------

.. automodule:: TACTICS.library_enumeration.enumeration_utils
   :members: EnumerationResult, EnumerationError, AutoDetectionResult, read_reagent_file, results_to_dataframe, failures_to_dataframe, summarize_failures

.. automodule:: TACTICS.library_enumeration.file_writer
   :members: write_enumerated_library, write_products_chunked

.. automodule:: TACTICS.library_enumeration.generate_products
   :members: enumerate_products, generate_all_combinations

Protecting groups and salts
---------------------------

.. automodule:: TACTICS.library_enumeration.smarts_toolkit.constants
   :members: DEFAULT_PROTECTING_GROUPS, DEFAULT_SALT_FRAGMENTS, PROTECTING_GROUP_MAP, get_protecting_group, get_all_protecting_group_names
