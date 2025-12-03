Library Enumeration Module
==========================

The Library Enumeration module provides tools for generating and managing chemical libraries,
including advanced SMARTS pattern validation, multi-SMARTS routing, and multi-step synthesis support.

LibraryEnumerator Class
-----------------------

.. autoclass:: TACTICS.library_enumeration.generate_products.LibraryEnumerator
   :members:
   :undoc-members:
   :show-inheritance:

   .. automethod:: __init__
   .. automethod:: enumerate_library
   .. automethod:: get_product_smiles

Utility Functions
-----------------

.. autofunction:: TACTICS.library_enumeration.enumeration_utils.find_reactants_from_product_code
.. autofunction:: TACTICS.library_enumeration.enumeration_utils.write_products_to_files

SMARTS Toolkit
--------------

The SMARTS Toolkit provides comprehensive tools for troubleshooting, validation, and advanced
reaction pattern handling including multi-SMARTS routing and multi-step synthesis.

SMARTSValidator
~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.SMARTSValidator
   :members:
   :undoc-members:
   :show-inheritance:

   Validates SMARTS patterns and checks for common issues:

   * Syntax validation
   * Atom mapping verification
   * Reactant/product template checking

ReagentCompatibilityAnalyzer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.ReagentCompatibilityAnalyzer
   :members:
   :undoc-members:
   :show-inheritance:

   Analyzes reagent compatibility with SMARTS patterns:

   * Tests which reagents match a given pattern
   * Identifies incompatible reagent combinations
   * Suggests alternative SMARTS patterns

ExceptionFinder
~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.ExceptionFinder
   :members:
   :undoc-members:
   :show-inheritance:

   Identifies reagents that fail to react with a given SMARTS pattern.

SMARTSVisualizer
~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.SMARTSVisualizer
   :members:
   :undoc-members:
   :show-inheritance:

   Visualization tools for SMARTS patterns and reaction results.

Multi-SMARTS Routing
--------------------

For reactions where a single SMARTS pattern doesn't cover all reagent types,
the SMARTS Router provides automatic pattern selection.

SMARTSRouter
~~~~~~~~~~~~

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.SMARTSRouter
   :members:
   :undoc-members:
   :show-inheritance:

   Routes reagents to appropriate SMARTS patterns based on compatibility.

SMARTSPatternConfig
~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.SMARTSPatternConfig
   :members:
   :undoc-members:
   :show-inheritance:

   Configuration for individual SMARTS patterns in a router.

AlternativeSMARTSConfig
~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.AlternativeSMARTSConfig
   :members:
   :undoc-members:
   :show-inheritance:

   Configuration for alternative SMARTS patterns (primary + fallbacks).

Multi-Step Synthesis
--------------------

Support for multi-step synthesis routes with protecting groups and intermediate tracking.

ReactionSequence
~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.ReactionSequence
   :members:
   :undoc-members:
   :show-inheritance:

   Orchestrates multi-step synthesis with intermediate handling.

MultiStepSynthesisConfig
~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.MultiStepSynthesisConfig
   :members:
   :undoc-members:
   :show-inheritance:

   Top-level configuration for multi-step synthesis routes.

ReactionStepConfig
~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.ReactionStepConfig
   :members:
   :undoc-members:
   :show-inheritance:

   Configuration for a single step in a multi-step synthesis.

ReactionInputConfig
~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.ReactionInputConfig
   :members:
   :undoc-members:
   :show-inheritance:

   Configuration for reaction inputs (reagents or intermediates).

ProtectingGroupConfig
~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.ProtectingGroupConfig
   :members:
   :undoc-members:
   :show-inheritance:

   Configuration for protecting group handling.

Usage Examples
--------------

Basic Multi-SMARTS Setup
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from TACTICS.library_enumeration.smarts_toolkit import (
        SMARTSRouter,
        AlternativeSMARTSConfig,
        SMARTSPatternConfig
    )

    # Configure alternative SMARTS patterns
    config = AlternativeSMARTSConfig(
        primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
        alternatives=[
            SMARTSPatternConfig(
                pattern_id="secondary_amine",
                reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]"
            )
        ]
    )

    router = SMARTSRouter(config)

Multi-Step Synthesis
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from TACTICS.library_enumeration.smarts_toolkit import (
        ReactionSequence,
        MultiStepSynthesisConfig,
        ReactionStepConfig,
        AlternativeSMARTSConfig
    )

    # Configure multi-step synthesis
    config = MultiStepSynthesisConfig(
        alternative_smarts=AlternativeSMARTSConfig(
            primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"
        ),
        reagent_file_list=["acids.smi", "amines.smi"]
    )

    sequence = ReactionSequence(config) 