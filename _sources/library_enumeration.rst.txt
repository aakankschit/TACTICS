Library Enumeration
===================

The Library Enumeration module provides tools for generating combinatorial chemical products
from reagent libraries using reaction SMARTS patterns. It supports single-step reactions,
alternative SMARTS routing, multi-step synthesis pipelines, and protecting group deprotection.


Module Architecture
-------------------

The following diagram shows the class hierarchy and dependencies:

.. graphviz::

    digraph LibraryEnumeration {
        rankdir=TB;
        node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=10];
        edge [fontname="Helvetica", fontsize=9];
        nodesep=0.2;
        ranksep=0.4;

        // Fundamental types
        fundamental_label [label="Fundamental Types", shape=plaintext, fontname="Helvetica Bold"];
        InputSource [label="InputSource (enum)", fillcolor="#DCDCDC"];
        ProtectingGroupInfo [label="ProtectingGroupInfo", fillcolor="#DCDCDC"];
        DeprotectionSpec [label="DeprotectionSpec", fillcolor="#DCDCDC"];
        StepInput [label="StepInput", fillcolor="#DCDCDC"];

        // Core reaction definition
        core_label [label="Reaction Definition", shape=plaintext, fontname="Helvetica Bold"];
        ReactionDef [label="ReactionDef (core class)", fillcolor="#007bff", fontcolor="white"];

        // Configuration
        config_label [label="Configuration", shape=plaintext, fontname="Helvetica Bold"];
        ReactionConfig [label="ReactionConfig", fillcolor="#17a2b8", fontcolor="white"];

        // Pipeline
        pipeline_label [label="Pipeline", shape=plaintext, fontname="Helvetica Bold"];
        SynthesisPipeline [label="SynthesisPipeline", fillcolor="#007bff", fontcolor="white"];

        // Results
        result_label [label="Result Types", shape=plaintext, fontname="Helvetica Bold"];
        ValidationResult [label="ValidationResult", fillcolor="#28a745", fontcolor="white"];
        EnumerationResult [label="EnumerationResult", fillcolor="#28a745", fontcolor="white"];
        EnumerationError [label="EnumerationError", fillcolor="#28a745", fontcolor="white"];
        AutoDetectionResult [label="AutoDetectionResult", fillcolor="#28a745", fontcolor="white"];

        // File Utilities
        utility_label [label="File Utilities", shape=plaintext, fontname="Helvetica Bold"];
        WriteLibrary [label="write_enumerated_library()", fillcolor="#ffc107"];
        ResultsDF [label="results_to_dataframe()", fillcolor="#ffc107"];
        FailuresDF [label="failures_to_dataframe()", fillcolor="#ffc107"];

        // Force vertical layout
        fundamental_label -> InputSource [style=invis];
        InputSource -> ProtectingGroupInfo [style=invis];
        ProtectingGroupInfo -> DeprotectionSpec [style=invis];
        DeprotectionSpec -> StepInput [style=invis];
        StepInput -> core_label [style=invis];
        core_label -> ReactionDef [style=invis];
        ReactionDef -> config_label [style=invis];
        config_label -> ReactionConfig [style=invis];
        ReactionConfig -> pipeline_label [style=invis];
        pipeline_label -> SynthesisPipeline [style=invis];
        SynthesisPipeline -> result_label [style=invis];
        result_label -> ValidationResult [style=invis];
        ValidationResult -> EnumerationResult [style=invis];
        EnumerationResult -> EnumerationError [style=invis];
        EnumerationError -> AutoDetectionResult [style=invis];
        AutoDetectionResult -> utility_label [style=invis];
        utility_label -> WriteLibrary [style=invis];
        WriteLibrary -> ResultsDF [style=invis];
        ResultsDF -> FailuresDF [style=invis];

        // Visible dependencies
        StepInput -> InputSource [label="uses", style=dashed, constraint=false];
        ReactionDef -> DeprotectionSpec [label="contains", style=dashed, constraint=false];
        ReactionDef -> ValidationResult [label="returns", constraint=false];
        ReactionConfig -> ReactionDef [label="contains 1+", style=bold, constraint=false];
        ReactionConfig -> StepInput [label="uses", style=dashed, constraint=false];
        ReactionConfig -> ProtectingGroupInfo [label="uses (opt)", style=dashed, constraint=false];
        SynthesisPipeline -> ReactionConfig [label="constructed from", style=bold, constraint=false];
        SynthesisPipeline -> EnumerationResult [label="returns", constraint=false];
        SynthesisPipeline -> AutoDetectionResult [label="returns (opt)", constraint=false];
        EnumerationResult -> WriteLibrary [label="input", style=dashed, constraint=false];
        EnumerationResult -> ResultsDF [label="input", style=dashed, constraint=false];
    }

Quick Start
-----------

**Single reaction with validation:**

.. code-block:: python
   :caption: Define and validate a reaction

   from TACTICS.library_enumeration import ReactionDef, ReactionConfig, SynthesisPipeline

   # Define a reaction
   rxn = ReactionDef(
       reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
       step_index=0,
       description="Amide coupling"
   )

   # Create configuration and pipeline
   config = ReactionConfig(
       reactions=[rxn],
       reagent_file_list=["acids.smi", "amines.smi"]
   )
   pipeline = SynthesisPipeline(config)

   # Validate against reagent files
   result = rxn.validate_reaction(reagent_files=["acids.smi", "amines.smi"])
   print(f"Coverage: {result.coverage_stats}")


**Enumerate products with SynthesisPipeline:**

.. code-block:: python
   :caption: Simple enumeration

   from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef

   # Create configuration
   config = ReactionConfig(
       reactions=[
           ReactionDef(
               reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
               step_index=0,
               description="Amide coupling",
           )
       ],
       reagent_file_list=["acids.smi", "amines.smi"]
   )

   # Create pipeline directly from config
   pipeline = SynthesisPipeline(config)

   # Enumerate a single product from SMILES
   result = pipeline.enumerate_single_from_smiles(["CC(=O)O", "CCN"])

   if result.success:
       print(f"Product: {result.product_smiles}")  # CCNC(C)=O


**With alternative SMARTS patterns:**

.. code-block:: python
   :caption: Auto-routing to compatible patterns

   from TACTICS.library_enumeration import (
       SynthesisPipeline, ReactionConfig, ReactionDef,
       StepInput, InputSource
   )

   # Define primary and alternative patterns
   config = ReactionConfig(
       reactions=[
           ReactionDef(
               reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
               step_index=0,
               pattern_id="primary",
               description="Primary amines"
           ),
           ReactionDef(
               reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]",
               step_index=0,
               pattern_id="secondary",
               description="Secondary amines"
           ),
       ],
       reagent_file_list=["acids.smi", "amines.smi"],
       step_inputs={
           0: [
               StepInput(source=InputSource.REAGENT_FILE, file_index=0),
               StepInput(source=InputSource.REAGENT_FILE, file_index=1)
           ]
       },
       step_modes={0: "alternative"}  # Mark step 0 as having alternatives
   )

   pipeline = SynthesisPipeline(config)

   # Pipeline automatically tries alternative patterns at runtime
   result = pipeline.enumerate_single_from_smiles(["CC(=O)O", "CCN(C)C"])
   print(f"Pattern used: {result.patterns_used}")  # {0: "secondary"}


Fundamental Types
-----------------

These are the basic building blocks used by higher-level classes.

.. _input-source:

InputSource
~~~~~~~~~~~

.. rst-class:: class-fundamental

Enum specifying the source of an input for a reaction step.

.. list-table:: Values
   :header-rows: 1
   :widths: 30 70

   * - Value
     - Description
   * - ``REAGENT_FILE``
     - Input comes from a reagent file (use ``file_index`` in :ref:`StepInput <step-input>`)
   * - ``PREVIOUS_STEP``
     - Input comes from output of a previous step (use ``step_index`` in :ref:`StepInput <step-input>`)

**Example**

.. code-block:: python

   from TACTICS.library_enumeration import InputSource

   source = InputSource.REAGENT_FILE
   source = InputSource.PREVIOUS_STEP


.. _protecting-group-info:

ProtectingGroupInfo
~~~~~~~~~~~~~~~~~~~

.. rst-class:: class-fundamental

Dataclass defining a protecting group for detection and optional removal.

.. admonition:: Dependencies
   :class: dependencies

   None - this is a standalone dataclass.

.. list-table:: Parameters
   :header-rows: 1
   :widths: 22 18 10 50

   * - Parameter
     - Type
     - Required
     - Description
   * - ``name``
     - ``str``
     - Yes
     - Human-readable name (e.g., "Boc", "Fmoc")
   * - ``smarts``
     - ``str``
     - Yes
     - SMARTS pattern to detect the group
   * - ``deprotection_smarts``
     - ``str``
     - No
     - Reaction SMARTS for removal (optional)

**Example**

.. code-block:: python

   from TACTICS.library_enumeration import ProtectingGroupInfo

   boc = ProtectingGroupInfo(
       name="Boc",
       smarts="[NX3][C](=O)OC(C)(C)C",
       deprotection_smarts="[N:1][C](=O)OC(C)(C)C>>[N:1]"
   )


.. _deprotection-spec:

DeprotectionSpec
~~~~~~~~~~~~~~~~

.. rst-class:: class-fundamental

Pydantic model specifying a deprotection to apply during synthesis.
Deprotections can target either a reactant (before the reaction) or the product (after the reaction).

.. admonition:: Dependencies
   :class: dependencies

   None - uses protecting group names as strings.

.. list-table:: Parameters
   :header-rows: 1
   :widths: 22 25 10 43

   * - Parameter
     - Type
     - Required
     - Description
   * - ``group``
     - ``str``
     - Yes
     - Name of protecting group (e.g., "Boc", "Fmoc")
   * - ``target``
     - ``int`` or ``"product"``
     - Yes
     - Reactant index (int >= 0) for pre-reaction, or ``"product"`` for post-reaction

.. list-table:: Properties
   :header-rows: 1
   :widths: 30 20 50

   * - Property
     - Type
     - Description
   * - ``is_product_deprotection``
     - ``bool``
     - True if this deprotects the product (after reaction)
   * - ``reactant_index``
     - ``int | None``
     - The reactant index if targeting a reactant, None if targeting product

**Example: Reactant Deprotection**

.. code-block:: python

   from TACTICS.library_enumeration import DeprotectionSpec

   # Remove Boc from the second reactant (index 1) BEFORE reaction
   deprot = DeprotectionSpec(group="Boc", target=1)

   print(deprot.is_product_deprotection)  # False
   print(deprot.reactant_index)           # 1

**Example: Product Deprotection**

.. code-block:: python

   from TACTICS.library_enumeration import DeprotectionSpec

   # Remove Fmoc from the product AFTER reaction
   deprot = DeprotectionSpec(group="Fmoc", target="product")

   print(deprot.is_product_deprotection)  # True
   print(deprot.reactant_index)           # None


.. _step-input:

StepInput
~~~~~~~~~

.. rst-class:: class-fundamental

Pydantic model specifying where an input comes from for a reaction step.

.. admonition:: Dependencies
   :class: dependencies

   - :ref:`InputSource <input-source>` - enum specifying the source type

**Depends on:** :ref:`InputSource <input-source>`

.. list-table:: Parameters
   :header-rows: 1
   :widths: 20 18 12 50

   * - Parameter
     - Type
     - Required
     - Description
   * - ``source``
     - ``InputSource``
     - Yes
     - Either ``REAGENT_FILE`` or ``PREVIOUS_STEP``
   * - ``file_index``
     - ``int``
     - Conditional
     - Index into ``reagent_file_list``. Required if source is ``REAGENT_FILE``
   * - ``step_index``
     - ``int``
     - Conditional
     - Index of previous step. Required if source is ``PREVIOUS_STEP``

**Example**

.. code-block:: python

   from TACTICS.library_enumeration import StepInput, InputSource

   # Input from first reagent file
   input1 = StepInput(source=InputSource.REAGENT_FILE, file_index=0)

   # Input from output of step 0
   input2 = StepInput(source=InputSource.PREVIOUS_STEP, step_index=0)


Reaction Definition
-------------------

.. _reaction-def:

ReactionDef
~~~~~~~~~~~

.. rst-class:: class-core

**The core class for defining a single reaction with built-in validation and visualization.**

This is the fundamental building block of the SMARTS toolkit. Each ``ReactionDef`` represents
a single chemical reaction that can validate reagent compatibility, visualize template matches,
and be combined into multi-step syntheses via :ref:`ReactionConfig <reaction-config>`.

.. admonition:: Dependencies
   :class: dependencies

   - :ref:`DeprotectionSpec <deprotection-spec>` - for deprotections (optional)
   - :ref:`ProtectingGroupInfo <protecting-group-info>` - used during validation (optional)
   - Returns :ref:`ValidationResult <validation-result>` from validation methods

**Depends on:** :ref:`DeprotectionSpec <deprotection-spec>` (optional), :ref:`ProtectingGroupInfo <protecting-group-info>` (optional)

Constructor Parameters
^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 18 10 52

   * - Parameter
     - Type
     - Required
     - Description
   * - ``reaction_smarts``
     - ``str``
     - Yes
     - Reaction SMARTS string (validated on creation)
   * - ``step_index``
     - ``int``
     - No
     - Step in the sequence (0 = first step). Default: 0
   * - ``pattern_id``
     - ``str``
     - No
     - Identifier for alternatives (auto-generated if not provided)
   * - ``description``
     - ``str``
     - No
     - Human-readable description
   * - ``deprotections``
     - ``list[DeprotectionSpec]``
     - No
     - Deprotections to apply for this reaction. Default: []

Properties
^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Property
     - Type
     - Description
   * - ``num_reactants``
     - ``int``
     - Number of reactants in the reaction
   * - ``is_validated``
     - ``bool``
     - True if ``validate_reaction()`` has been called
   * - ``coverage_stats``
     - ``dict[int, float]``
     - Coverage percentage per position (0-100)
   * - ``validation_result``
     - ``ValidationResult``
     - Cached validation result (None if not validated)

Validation Methods
^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Method
     - Description
   * - ``validate_reaction(reagent_files=None, reagent_smiles=None, ...)``
     - Validate reaction against reagent files or SMILES lists. Returns :ref:`ValidationResult <validation-result>`
   * - ``get_compatible_reagents(position)``
     - Get list of ``(smiles, name)`` tuples for compatible reagents at position
   * - ``get_incompatible_reagents(position)``
     - Get list of ``(smiles, name)`` tuples for incompatible reagents at position
   * - ``get_reactant_template(position)``
     - Get RDKit Mol template for a specific position
   * - ``summary()``
     - Human-readable validation summary string

**validate_reaction parameters:**

.. list-table::
   :header-rows: 1
   :widths: 22 18 10 50

   * - Parameter
     - Type
     - Required
     - Description
   * - ``reagent_files``
     - ``list[str]``
     - No*
     - Paths to reagent files (.smi format)
   * - ``reagent_smiles``
     - ``list[list[tuple]]``
     - No*
     - Lists of ``(smiles, name)`` tuples per position
   * - ``protecting_groups``
     - ``list[ProtectingGroupInfo]``
     - No
     - Custom protecting groups for detection
   * - ``deprotect``
     - ``bool``
     - No
     - Apply deprotection during validation. Default: False
   * - ``desalt``
     - ``bool``
     - No
     - Remove salt fragments during validation. Default: False
   * - ``test_reactions``
     - ``bool``
     - No
     - Run test reactions to verify products form. Default: False

\* Must provide either ``reagent_files`` or ``reagent_smiles``

Visualization Methods
^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Method
     - Description
   * - ``visualize_template_match(smiles, position, ...)``
     - Visualize which atoms in a molecule match the reaction template
   * - ``visualize_reaction(size=(800, 200))``
     - Visualize the reaction scheme

Example
^^^^^^^

.. code-block:: python
   :caption: Complete ReactionDef workflow

   from TACTICS.library_enumeration import ReactionDef, DeprotectionSpec

   # Define reaction with deprotection on reactant
   rxn = ReactionDef(
       reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
       step_index=0,
       pattern_id="amide_coupling",
       description="Amide bond formation",
       deprotections=[DeprotectionSpec(group="Boc", target=1)]
   )

   # Validate against reagent files
   result = rxn.validate_reaction(
       reagent_files=["acids.smi", "amines.smi"],
       deprotect=True,
       test_reactions=True
   )

   # Check coverage
   print(f"Position 0 coverage: {rxn.coverage_stats[0]:.1f}%")
   print(f"Position 1 coverage: {rxn.coverage_stats[1]:.1f}%")
   print(f"Reaction success rate: {result.reaction_success_rate:.1f}%")

   # Get compatible reagents for position 1
   compatible = rxn.get_compatible_reagents(1)
   print(f"Found {len(compatible)} compatible amines")

   # Troubleshoot a problematic reagent
   img = rxn.visualize_template_match("CC(C)(C)OC(=O)NCCn", position=1)


Configuration
-------------

.. _reaction-config:

ReactionConfig
~~~~~~~~~~~~~~

.. rst-class:: class-config

**Container for synthesis configuration holding one or more ReactionDef objects.**

Use ``ReactionConfig`` to define complex syntheses including:

- Single reactions (one ``ReactionDef``)
- Alternative SMARTS patterns (multiple ``ReactionDef`` objects with same ``step_index``)
- Multi-step syntheses (multiple ``ReactionDef`` objects with different ``step_index`` values)
- Protecting group deprotections (via ``ReactionDef.deprotections``)

.. admonition:: Dependencies
   :class: dependencies

   Composes multiple lower-level classes:

   - :ref:`ReactionDef <reaction-def>` - one or more reaction definitions (required)
   - :ref:`StepInput <step-input>` - input sources for multi-step synthesis (required for multi-step)
   - :ref:`ProtectingGroupInfo <protecting-group-info>` - custom protecting groups (optional)

**Depends on:** :ref:`ReactionDef <reaction-def>`, :ref:`StepInput <step-input>`, :ref:`ProtectingGroupInfo <protecting-group-info>`

Constructor Parameters
^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 22 22 10 46

   * - Parameter
     - Type
     - Required
     - Description
   * - ``reactions``
     - ``list[ReactionDef]``
     - Yes
     - List of reaction definitions (minimum 1)
   * - ``reagent_file_list``
     - ``list[str]``
     - No
     - Paths to reagent files. Default: []
   * - ``step_inputs``
     - ``dict[int, list[StepInput]]``
     - Conditional
     - Mapping of step_index to input sources. Required if multiple reactions
   * - ``step_modes``
     - ``dict[int, str]``
     - No
     - Mark steps with alternative patterns: ``{step_index: "alternative"}``
   * - ``protecting_groups``
     - ``list[ProtectingGroupInfo]``
     - No
     - Custom protecting group definitions

Properties
^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Property
     - Type
     - Description
   * - ``num_steps``
     - ``int``
     - Number of unique steps
   * - ``is_multi_step``
     - ``bool``
     - True if more than one step
   * - ``steps_with_alternatives``
     - ``list[int]``
     - Step indices that have alternative patterns
   * - ``step_indices``
     - ``list[int]``
     - Sorted list of all step indices

Methods
^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Method
     - Description
   * - ``get_reactions_for_step(step_index)``
     - Get all ReactionDef objects for a given step
   * - ``get_primary_reaction(step_index)``
     - Get the primary reaction (pattern_id='primary' or first)
   * - ``get_inputs_for_step(step_index)``
     - Get StepInput configuration for a step
   * - ``has_alternatives_at_step(step_index)``
     - Check if step has alternative SMARTS patterns
   * - ``validate_all(deprotect=False, desalt=False)``
     - Validate all reactions, returns ``{step: {pattern_id: ValidationResult}}``

Example: Single Reaction
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python
   :caption: Single reaction configuration

   from TACTICS.library_enumeration import ReactionDef, ReactionConfig, SynthesisPipeline

   rxn = ReactionDef(
       reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
       description="Amide coupling"
   )

   config = ReactionConfig(
       reactions=[rxn],
       reagent_file_list=["acids.smi", "amines.smi"]
   )

   # Create pipeline directly from config
   pipeline = SynthesisPipeline(config)

Example: Multi-Step Synthesis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python
   :caption: Two-step synthesis with intermediate

   from TACTICS.library_enumeration import (
       ReactionDef, ReactionConfig, StepInput, InputSource, DeprotectionSpec,
       SynthesisPipeline
   )

   # Step 0: Amide coupling
   step0 = ReactionDef(
       reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
       step_index=0,
       description="Amide coupling"
   )

   # Step 1: Reductive amination (with Boc deprotection on reactant)
   step1 = ReactionDef(
       reaction_smarts="[NH2:1].[CH:2]=O>>[NH:1][CH2:2]",
       step_index=1,
       description="Reductive amination",
       deprotections=[DeprotectionSpec(group="Boc", target=0)]  # Deprotect first input
   )

   config = ReactionConfig(
       reactions=[step0, step1],
       reagent_file_list=["acids.smi", "amines.smi", "aldehydes.smi"],
       step_inputs={
           0: [
               StepInput(source=InputSource.REAGENT_FILE, file_index=0),
               StepInput(source=InputSource.REAGENT_FILE, file_index=1)
           ],
           1: [
               StepInput(source=InputSource.PREVIOUS_STEP, step_index=0),
               StepInput(source=InputSource.REAGENT_FILE, file_index=2)
           ],
       }
   )

   pipeline = SynthesisPipeline(config)

Example: Alternative SMARTS
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python
   :caption: Multiple patterns for the same step with runtime fallback

   from TACTICS.library_enumeration import (
       ReactionDef, ReactionConfig, StepInput, InputSource, SynthesisPipeline
   )

   # Primary pattern for primary amines
   primary = ReactionDef(
       reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
       step_index=0,
       pattern_id="primary",  # Required: first alternative should be "primary"
       description="Primary amines"
   )

   # Alternative for secondary amines
   secondary = ReactionDef(
       reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]",
       step_index=0,
       pattern_id="secondary",
       description="Secondary amines"
   )

   config = ReactionConfig(
       reactions=[primary, secondary],
       reagent_file_list=["acids.smi", "amines.smi"],
       step_inputs={
           0: [
               StepInput(source=InputSource.REAGENT_FILE, file_index=0),
               StepInput(source=InputSource.REAGENT_FILE, file_index=1)
           ]
       },
       step_modes={0: "alternative"}  # Mark step 0 as having alternatives
   )

   pipeline = SynthesisPipeline(config)

   # Runtime pattern fallback: pipeline automatically tries patterns until one succeeds
   # No need to call auto_detect_compatibility()


Pipeline
--------

.. _synthesis-pipeline:

SynthesisPipeline
~~~~~~~~~~~~~~~~~

.. rst-class:: class-core

**Main entry point for executing syntheses and enumerating products.**

``SynthesisPipeline`` orchestrates the synthesis process, handling:

- Single reactions from SMARTS strings
- Automatic routing to compatible alternative patterns (runtime fallback)
- Multi-step syntheses with intermediate tracking
- Batch enumeration with optional parallelization
- Integration with Thompson Sampling via multiprocessing support

.. admonition:: Dependencies
   :class: dependencies

   - :ref:`ReactionConfig <reaction-config>` - passed to constructor
   - :ref:`ReactionDef <reaction-def>` - accessed via config
   - Returns :ref:`EnumerationResult <enumeration-result>`, :ref:`EnumerationError <enumeration-error>`, :ref:`AutoDetectionResult <auto-detection-result>`

**Depends on:** :ref:`ReactionConfig <reaction-config>`

Constructor
^^^^^^^^^^^

.. code-block:: python

   from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef

   config = ReactionConfig(
       reactions=[ReactionDef(reaction_smarts="...", step_index=0)],
       reagent_file_list=["reagents.smi"]
   )

   pipeline = SynthesisPipeline(config)

.. list-table:: Constructor Parameters
   :header-rows: 1
   :widths: 20 20 10 50

   * - Parameter
     - Type
     - Required
     - Description
   * - ``config``
     - ``ReactionConfig``
     - Yes
     - Reaction configuration with reactions and reagent files

Enumeration Methods
^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Method
     - Description
   * - ``enumerate_single(reagent_mols, reagent_keys=None, ...)``
     - Enumerate single product from RDKit Mol objects
   * - ``enumerate_single_from_smiles(smiles_list, reagent_keys=None, ...)``
     - Enumerate single product from SMILES strings
   * - ``enumerate_library(n_jobs=1, ...)``
     - Enumerate all combinations from reagent files
   * - ``enumerate_batch(combinations, n_jobs=1, ...)``
     - Enumerate specific combinations, optionally in parallel

**enumerate_single_from_smiles parameters:**

.. list-table::
   :header-rows: 1
   :widths: 22 22 10 46

   * - Parameter
     - Type
     - Required
     - Description
   * - ``smiles_list``
     - ``list[str]``
     - Yes
     - SMILES strings for each reagent position
   * - ``reagent_keys``
     - ``list[str]``
     - No
     - Keys for compatibility lookup (names or identifiers)
   * - ``store_intermediates``
     - ``bool``
     - No
     - Store intermediate products for multi-step. Default: False

**Returns:** :ref:`EnumerationResult <enumeration-result>`

Utility Methods
^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Method
     - Description
   * - ``get_validator(step_index=0, pattern_id="primary")``
     - Get internal validator for troubleshooting
   * - ``validate_all()``
     - Validate all steps and patterns
   * - ``auto_detect_compatibility(reagent_lists=None, ...)``
     - Pre-detect which reagents work with which patterns (optional optimization)
   * - ``register_compatibility(position, reagent_key, patterns, ...)``
     - Manually register a reagent's compatible patterns
   * - ``get_compatible_patterns(reagent_keys, step_index=0)``
     - Find pattern compatible with all reagents
   * - ``get_compatibility_map()``
     - Get full compatibility mapping
   * - ``prepare_worker_data()``
     - Serialize for multiprocessing workers
   * - ``from_worker_data(data)``
     - Reconstruct pipeline in worker process (class method)

Properties
^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Property
     - Type
     - Description
   * - ``num_steps``
     - ``int``
     - Number of reaction steps
   * - ``num_components``
     - ``int``
     - Number of reagent files/positions required
   * - ``is_multi_step``
     - ``bool``
     - True if multi-step synthesis
   * - ``has_alternatives``
     - ``bool``
     - True if any step has alternative patterns
   * - ``pattern_ids``
     - ``dict[int, list[str]]``
     - Available pattern IDs at each step
   * - ``reactions``
     - ``list[ReactionDef]``
     - All reactions from config
   * - ``reagent_file_list``
     - ``list[str]``
     - Reagent file paths from config

Example: Basic Enumeration
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python
   :caption: Enumerate products from SMILES

   from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef

   config = ReactionConfig(
       reactions=[
           ReactionDef(
               reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
               step_index=0
           )
       ],
       reagent_file_list=["acids.smi", "amines.smi"]
   )

   pipeline = SynthesisPipeline(config)

   # Single enumeration
   result = pipeline.enumerate_single_from_smiles(["CC(=O)O", "CCN"])

   if result.success:
       print(f"Product: {result.product_smiles}")
       print(f"Pattern used: {result.patterns_used}")
   else:
       print(f"Failed: {result.error}")

Example: Library Enumeration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python
   :caption: Enumerate all products with parallelization

   from TACTICS.library_enumeration import SynthesisPipeline, ReactionConfig, ReactionDef, results_to_dataframe

   config = ReactionConfig(
       reactions=[
           ReactionDef(
               reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
               step_index=0
           )
       ],
       reagent_file_list=["acids.smi", "amines.smi"]
   )

   pipeline = SynthesisPipeline(config)

   # Enumerate all combinations with 4 parallel workers
   all_results = pipeline.enumerate_library(n_jobs=4, show_progress=True)

   # Analyze results
   successes = [r for r in all_results if r.success]
   failures = [r for r in all_results if not r.success]
   print(f"Success rate: {len(successes)/len(all_results)*100:.1f}%")

   # Convert to Polars DataFrame
   df = results_to_dataframe(all_results)


Protecting Groups and Deprotection
----------------------------------

The SMARTS toolkit includes built-in support for common protecting groups and allows
defining custom groups for specialized chemistry.

.. _built-in-protecting-groups:

Built-in Protecting Groups
~~~~~~~~~~~~~~~~~~~~~~~~~~

The following protecting groups are available by default:

.. list-table::
   :header-rows: 1
   :widths: 15 18 35 32

   * - Name
     - Protects
     - Detection SMARTS
     - Result
   * - Boc
     - Amine (N)
     - ``[NX3][C](=O)OC(C)(C)C``
     - Free amine
   * - Fmoc
     - Amine (N)
     - ``[NX3]C(=O)OCC1c2ccccc2-c2ccccc12``
     - Free amine
   * - Cbz
     - Amine (N)
     - ``[NX3]C(=O)OCc1ccccc1``
     - Free amine
   * - Acetamide
     - Amine (N)
     - ``[NX3][C](=O)[CH3]``
     - Free amine
   * - TBS
     - Alcohol (O)
     - ``[OX2][Si](C)(C)C(C)(C)C``
     - Free alcohol
   * - O-Benzyl
     - Alcohol (O)
     - ``[OX2]Cc1ccccc1``
     - Free alcohol
   * - Trityl
     - Amine/Alcohol
     - ``[NX3,OX2]C(c1ccccc1)(c1ccccc1)c1ccccc1``
     - Free N/O
   * - tBu-ester
     - Carboxylic acid
     - ``[CX3](=O)OC(C)(C)C``
     - Free acid
   * - Me-ester
     - Carboxylic acid
     - ``[CX3](=O)O[CH3]``
     - Free acid
   * - Et-ester
     - Carboxylic acid
     - ``[CX3](=O)OCC``
     - Free acid

**Accessing protecting group information:**

.. code-block:: python

   from TACTICS.library_enumeration.smarts_toolkit.constants import (
       get_all_protecting_group_names,
       get_protecting_group,
   )

   # List all available protecting groups
   print(get_all_protecting_group_names())
   # ['Boc', 'Fmoc', 'Cbz', 'Acetamide', 'TBS', 'O-Benzyl', 'Trityl', 'tBu-ester', 'Me-ester', 'Et-ester']

   # Get details for a specific group
   boc = get_protecting_group("Boc")
   print(f"Name: {boc.name}")
   print(f"Detection SMARTS: {boc.smarts}")
   print(f"Deprotection SMARTS: {boc.deprotection_smarts}")

Reactant Deprotection
~~~~~~~~~~~~~~~~~~~~~

Deprotect a reactant **before** the reaction runs using ``target=<reactant_index>``:

.. code-block:: python

   from TACTICS.library_enumeration import (
       ReactionConfig, ReactionDef, DeprotectionSpec, SynthesisPipeline
   )

   config = ReactionConfig(
       reactions=[
           ReactionDef(
               reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
               step_index=0,
               deprotections=[
                   # Remove Boc from second reactant (index 1) before reaction
                   DeprotectionSpec(group="Boc", target=1),
               ],
           )
       ],
       reagent_file_list=["acids.smi", "boc_amines.smi"],
   )

   pipeline = SynthesisPipeline(config)

   # Boc-protected amine will be deprotected before amide coupling
   result = pipeline.enumerate_single_from_smiles([
       "CC(=O)O",                    # Acetic acid
       "CC(C)(C)OC(=O)NCCN"          # Boc-ethylenediamine
   ])

Product Deprotection
~~~~~~~~~~~~~~~~~~~~

Deprotect the product **after** the reaction runs using ``target="product"``:

.. code-block:: python

   from TACTICS.library_enumeration import (
       ReactionConfig, ReactionDef, DeprotectionSpec, SynthesisPipeline
   )

   config = ReactionConfig(
       reactions=[
           ReactionDef(
               reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
               step_index=0,
               deprotections=[
                   # Remove Fmoc from product after reaction
                   DeprotectionSpec(group="Fmoc", target="product"),
               ],
           )
       ],
       reagent_file_list=["acids.smi", "amines.smi"],
   )

   pipeline = SynthesisPipeline(config)

Custom Protecting Groups
~~~~~~~~~~~~~~~~~~~~~~~~

Define custom protecting groups using :ref:`ProtectingGroupInfo <protecting-group-info>`:

.. code-block:: python

   from TACTICS.library_enumeration import (
       ReactionConfig, ReactionDef, DeprotectionSpec, ProtectingGroupInfo,
       SynthesisPipeline
   )

   # Define a custom protecting group
   alloc = ProtectingGroupInfo(
       name="Alloc",
       smarts="[NX3]C(=O)OCC=C",
       deprotection_smarts="[N:1]C(=O)OCC=C>>[N:1]"
   )

   config = ReactionConfig(
       reactions=[
           ReactionDef(
               reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
               step_index=0,
               deprotections=[
                   DeprotectionSpec(group="Alloc", target="product"),
               ],
           ),
       ],
       reagent_file_list=["acids.smi", "amines.smi"],
       protecting_groups=[alloc],  # Register custom group
   )

   pipeline = SynthesisPipeline(config)


Result Types
------------

.. _validation-result:

ValidationResult
~~~~~~~~~~~~~~~~

.. rst-class:: class-result

**Comprehensive results from SMARTS validation.**

Returned by :ref:`ReactionDef.validate_reaction() <reaction-def>`.

.. list-table:: Attributes
   :header-rows: 1
   :widths: 28 25 47

   * - Attribute
     - Type
     - Description
   * - ``compatible_reagents``
     - ``dict[int, list[tuple]]``
     - ``{position: [(smiles, name), ...]}`` - reagents that match
   * - ``incompatible_reagents``
     - ``dict[int, list[tuple]]``
     - ``{position: [(smiles, name), ...]}`` - reagents that don't match
   * - ``invalid_smiles``
     - ``dict[int, list[tuple]]``
     - Unparseable SMILES entries
   * - ``duplicate_smiles``
     - ``dict[int, list[tuple]]``
     - Duplicate entries
   * - ``protected_reagents``
     - ``dict[int, list[tuple]]``
     - ``{position: [(smiles, name, [groups]), ...]}`` - with protecting groups
   * - ``multi_fragment_reagents``
     - ``dict[int, list[tuple]]``
     - ``{position: [(smiles, name, [fragments]), ...]}`` - with salts
   * - ``coverage_stats``
     - ``dict[int, float]``
     - Percent compatible per position (0-100)
   * - ``reaction_success_rate``
     - ``float``
     - Percent of test reactions that succeeded
   * - ``error_messages``
     - ``list[str]``
     - Critical errors
   * - ``warnings``
     - ``list[str]``
     - Non-critical warnings

**Methods and Properties:**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Member
     - Description
   * - ``is_valid()``
     - True if all positions have >0% coverage and no critical errors
   * - ``total_compatible``
     - Total compatible reagents across all positions
   * - ``total_incompatible``
     - Total incompatible reagents


.. _enumeration-result:

EnumerationResult
~~~~~~~~~~~~~~~~~

.. rst-class:: class-result

**Complete result from pipeline enumeration.**

Returned by :ref:`SynthesisPipeline <synthesis-pipeline>` enumeration methods.

.. list-table:: Attributes
   :header-rows: 1
   :widths: 22 25 53

   * - Attribute
     - Type
     - Description
   * - ``product``
     - ``Mol | None``
     - Final product as RDKit Mol object
   * - ``product_smiles``
     - ``str | None``
     - SMILES of final product
   * - ``product_name``
     - ``str | None``
     - Name derived from reagent keys
   * - ``patterns_used``
     - ``dict[int, str]``
     - ``{step_index: pattern_id}`` for each step
   * - ``intermediates``
     - ``dict[int, Mol] | None``
     - ``{step_index: mol}`` if ``store_intermediates=True``
   * - ``error``
     - ``EnumerationError | None``
     - Error details if failed

**Properties:**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Property
     - Description
   * - ``success``
     - True if enumeration succeeded (product is not None and no error)


.. _enumeration-error:

EnumerationError
~~~~~~~~~~~~~~~~

.. rst-class:: class-result

**Details about a failed enumeration.**

Contained in :ref:`EnumerationResult.error <enumeration-result>` when enumeration fails.

.. list-table:: Attributes
   :header-rows: 1
   :widths: 22 25 53

   * - Attribute
     - Type
     - Description
   * - ``step_index``
     - ``int``
     - Which step failed
   * - ``pattern_id``
     - ``str | None``
     - Which pattern was attempted
   * - ``error_type``
     - ``str``
     - One of: ``no_compatible_pattern``, ``reaction_failed``, ``invalid_input``, ``deprotection_failed``
   * - ``message``
     - ``str``
     - Human-readable description
   * - ``reagent_smiles``
     - ``list[str]``
     - Input SMILES that caused the failure

**Example: Handling Errors**

.. code-block:: python

   result = pipeline.enumerate_single_from_smiles(["CC(=O)O", "CCC"])  # Propane has no amine

   if not result.success:
       err = result.error
       print(f"Failed at step {err.step_index}")
       print(f"Error type: {err.error_type}")
       print(f"Message: {err.message}")
       print(f"Reagents: {err.reagent_smiles}")


.. _auto-detection-result:

AutoDetectionResult
~~~~~~~~~~~~~~~~~~~

.. rst-class:: class-result

**Results from automatic SMARTS compatibility detection.**

Returned by :ref:`SynthesisPipeline.auto_detect_compatibility() <synthesis-pipeline>`.

.. list-table:: Attributes
   :header-rows: 1
   :widths: 25 30 45

   * - Attribute
     - Type
     - Description
   * - ``pattern_results``
     - ``dict[int, dict[str, ValidationResult]]``
     - ``{step_index: {pattern_id: ValidationResult}}``
   * - ``compatibility_map``
     - ``dict[tuple, set[str]]``
     - ``{(step_index, position, reagent_key): {pattern_ids}}``
   * - ``coverage_by_pattern``
     - ``dict[int, dict[str, dict[int, float]]]``
     - ``{step_index: {pattern_id: {position: coverage%}}}``
   * - ``unmatched_reagents``
     - ``dict[int, list[str]]``
     - ``{position: [reagent_names]}`` - reagents that match no pattern
   * - ``warnings``
     - ``list[str]``
     - Warning messages


Constants
---------

The module provides default protecting groups and salt fragments for common use cases.

DEFAULT_PROTECTING_GROUPS
~~~~~~~~~~~~~~~~~~~~~~~~~

List of 10 common protecting groups (see :ref:`built-in-protecting-groups` for details).

DEFAULT_SALT_FRAGMENTS
~~~~~~~~~~~~~~~~~~~~~~

List of ~25 common salt/counterion fragments including:

- Halides (Cl⁻, Br⁻, I⁻)
- Metal cations (Na⁺, K⁺, Li⁺, Ca²⁺, Mg²⁺)
- Organic acids (TFA, acetate, formate)
- Ammonium salts

Utility Functions
~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Function
     - Description
   * - ``get_protecting_group(name)``
     - Get :ref:`ProtectingGroupInfo <protecting-group-info>` by name. Raises ``KeyError`` if not found.
   * - ``get_all_protecting_group_names()``
     - Get list of all default protecting group names.
   * - ``results_to_dataframe(results)``
     - Convert list of EnumerationResult to Polars DataFrame.
   * - ``failures_to_dataframe(results)``
     - Convert failed EnumerationResult objects to DataFrame.
   * - ``summarize_failures(results)``
     - Get summary statistics of enumeration failures.


