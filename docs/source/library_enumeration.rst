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

The SMARTS Toolkit (version 0.2.0) provides comprehensive tools for troubleshooting, validation, and advanced
reaction pattern handling including multi-SMARTS routing and multi-step synthesis.

The toolkit is organized into several key components:

- **Validation**: :class:`SMARTSValidator` for pattern validation and reagent compatibility analysis
- **Visualization**: :class:`SMARTSVisualizer` for Jupyter-compatible debugging visualizations
- **Multi-SMARTS Routing**: :class:`SMARTSRouter` for handling alternative reaction patterns
- **Multi-Step Synthesis**: :class:`ReactionSequence` for orchestrating complex synthesis routes

SMARTSValidator
~~~~~~~~~~~~~~~

The core validation class for analyzing SMARTS patterns against reagent libraries.

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.SMARTSValidator
   :members:
   :undoc-members:
   :show-inheritance:

**Key Features:**

* Validates reaction SMARTS syntax and structure
* Loads reagents from ``.smi`` and ``.csv`` files
* Detects protecting groups and salt fragments (with customizable lists)
* Optional import-time deprotection and desalting
* Template matching analysis (which reagents are compatible)
* Reaction testing with sampling
* Coverage statistics per reagent position
* Exports compatible reagents to files

**Basic Usage Example:**

.. code-block:: python

    from TACTICS.library_enumeration.smarts_toolkit import SMARTSValidator

    # Define reaction SMARTS (amide coupling)
    reaction_smarts = "[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"

    # Initialize validator with reagent files
    validator = SMARTSValidator(
        reaction_smarts=reaction_smarts,
        reagent_files=["acids.smi", "amines.smi"]
    )

    # Run validation
    result = validator.validate(test_reactions=True, sample_size=100)

    # Check coverage statistics
    for position, coverage in result.coverage_stats.items():
        print(f"Position {position}: {coverage:.1f}% compatible")

    # Check for errors and warnings
    if result.error_messages:
        print("Errors:", result.error_messages)
    if result.warnings:
        print("Warnings:", result.warnings)

**Using Deprotection and Desalting:**

.. code-block:: python

    from TACTICS.library_enumeration.smarts_toolkit import SMARTSValidator

    # Enable automatic deprotection and desalting during import
    validator = SMARTSValidator(
        reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
        reagent_files=["acids.smi", "amines.smi"],
        deprotect_on_import=True,              # Remove protecting groups
        deprotect_groups=["Boc", "Fmoc"],      # Only these groups (None = all)
        desalt_on_import=True                  # Remove salt counterions
    )

    result = validator.validate()

    # Or deprotect/desalt individual SMILES manually
    deprotected, removed_groups = validator.deprotect_smiles(
        "CC(C)(C)OC(=O)NCc1ccccc1",
        groups_to_remove=["Boc"]
    )
    print(f"Deprotected: {deprotected}, removed: {removed_groups}")

    desalted, removed_salts = validator.desalt_smiles("CCN.[Cl-]")
    print(f"Desalted: {desalted}, removed: {removed_salts}")

**Exporting Compatible Reagents:**

.. code-block:: python

    from TACTICS.library_enumeration.smarts_toolkit import SMARTSValidator

    validator = SMARTSValidator(
        reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
        reagent_files=["acids.smi", "amines.smi"]
    )

    # Export compatible reagents to files
    created_files = validator.export_compatible_reagents(
        output_dir="./compatible_reagents",
        prefix=["acid", "amine"],  # Naming prefix per position
        deprotect=True,            # Apply deprotection to exported SMILES
        desalt=True                # Apply desalting to exported SMILES
    )
    print(f"Created files: {created_files}")

**Generating Compatibility Reports:**

.. code-block:: python

    from TACTICS.library_enumeration.smarts_toolkit import SMARTSValidator

    validator = SMARTSValidator(
        reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
        reagent_files=["acids.smi", "amines.smi"]
    )

    # Get compatibility report as Polars DataFrame
    report = validator.generate_compatibility_report()
    print(report)

    # Or use the property shortcut
    print(validator.compatibility_report)

ValidationResult
~~~~~~~~~~~~~~~~

Container dataclass for validation results.

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.ValidationResult
   :members:
   :undoc-members:
   :show-inheritance:

**Attributes:**

* ``compatible_reagents``: Dict mapping position -> list of (SMILES, name) tuples
* ``incompatible_reagents``: Dict mapping position -> list of (SMILES, name) tuples
* ``invalid_smiles``: Dict mapping position -> list of invalid (SMILES, name) tuples
* ``duplicate_smiles``: Dict mapping position -> list of duplicate (SMILES, name) tuples
* ``protected_reagents``: Dict mapping position -> list of (SMILES, name, [protecting_group_names])
* ``multi_fragment_reagents``: Dict mapping position -> list of (SMILES, name, [detected_salt_names])
* ``coverage_stats``: Dict mapping position -> coverage percentage (0-100)
* ``reaction_success_rate``: Float percentage of successful test reactions
* ``error_messages``: List of error strings
* ``warnings``: List of warning strings

ProtectingGroupInfo
~~~~~~~~~~~~~~~~~~~

Definition for a protecting group with detection and deprotection patterns.

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.ProtectingGroupInfo
   :members:
   :undoc-members:
   :show-inheritance:

**Example - Adding Custom Protecting Groups:**

.. code-block:: python

    from TACTICS.library_enumeration.smarts_toolkit import (
        SMARTSValidator,
        ProtectingGroupInfo
    )

    # Define a custom protecting group
    custom_pg = ProtectingGroupInfo(
        name="MyProtectingGroup",
        smarts="[NX3]C(=O)CC",  # Detection pattern
        deprotection_smarts="[N:1]C(=O)CC>>[N:1]"  # Deprotection reaction
    )

    validator = SMARTSValidator(
        reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
        reagent_files=["acids.smi", "amines.smi"],
        additional_protecting_groups=[custom_pg]
    )

Default Protecting Groups and Salt Fragments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The toolkit provides default lists of common protecting groups and salt fragments:

.. autodata:: TACTICS.library_enumeration.smarts_toolkit.DEFAULT_PROTECTING_GROUPS

**Default protecting groups include:** Boc, Fmoc, Cbz, Acetamide, TBS, O-Bn, Trityl, tBu-ester, Me-ester, Et-ester

.. autodata:: TACTICS.library_enumeration.smarts_toolkit.DEFAULT_SALT_FRAGMENTS

**Default salt fragments include:** halides (F⁻, Cl⁻, Br⁻, I⁻), common counterions (triflate, mesylate, tosylate),
cations (Na⁺, K⁺, Li⁺), and common solvents (water, methanol, ethanol).

**Example - Adding Custom Salt Fragments:**

.. code-block:: python

    from TACTICS.library_enumeration.smarts_toolkit import SMARTSValidator

    custom_salts = [
        ("[BF4-]", "tetrafluoroborate"),
        ("[PF6-]", "hexafluorophosphate"),
    ]

    validator = SMARTSValidator(
        reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
        reagent_files=["acids.smi", "amines.smi"],
        additional_salt_fragments=custom_salts
    )

SMARTSVisualizer
~~~~~~~~~~~~~~~~

Jupyter-compatible visualization tools for SMARTS pattern debugging.

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.SMARTSVisualizer
   :members:
   :undoc-members:
   :show-inheritance:

**Key Methods:**

* ``visualize_reaction()``: Display the reaction SMARTS pattern
* ``visualize_compatible_reagents(position, max_molecules)``: Grid of compatible molecules with template highlighting
* ``visualize_incompatible_reagents(position, max_molecules)``: Grid of incompatible molecules
* ``generate_summary_plot()``: Matplotlib bar/pie charts with statistics

**Usage Example (Jupyter Notebook):**

.. code-block:: python

    from TACTICS.library_enumeration.smarts_toolkit import (
        SMARTSValidator,
        SMARTSVisualizer
    )

    # Initialize validator
    validator = SMARTSValidator(
        reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
        reagent_files=["acids.smi", "amines.smi"]
    )

    # Create visualizer
    viz = SMARTSVisualizer(
        validator,
        img_size=(300, 300),
        mols_per_row=4
    )

    # Visualize the reaction SMARTS
    viz.visualize_reaction()

    # Show compatible reagents at position 0 with template highlighting
    viz.visualize_compatible_reagents(position=0, max_molecules=20)

    # Show incompatible reagents at position 1
    viz.visualize_incompatible_reagents(position=1, max_molecules=10)

    # Generate summary statistics plot
    viz.generate_summary_plot()

Multi-SMARTS Routing
--------------------

For reactions where a single SMARTS pattern doesn't cover all reagent types,
the SMARTS Router provides automatic pattern selection based on reagent compatibility.

SMARTSRouter
~~~~~~~~~~~~

Routes reagent combinations to appropriate SMARTS patterns based on compatibility.

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.SMARTSRouter
   :members:
   :undoc-members:
   :show-inheritance:

**Key Features:**

* Manages multiple SMARTS patterns for a single reaction step (primary + alternatives)
* Determines which SMARTS pattern to use based on reagent compatibility
* Executes reactions with the selected SMARTS
* Handles enumeration failures gracefully
* Pre-compiles all reactions for efficiency
* Caches compatibility information

**Basic Usage Example:**

.. code-block:: python

    from rdkit import Chem
    from TACTICS.library_enumeration.smarts_toolkit import (
        SMARTSRouter,
        AlternativeSMARTSConfig,
        SMARTSPatternConfig
    )

    # Configure alternative SMARTS patterns
    # Primary handles primary amines, alternative handles secondary amines
    config = AlternativeSMARTSConfig(
        primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
        alternatives=[
            SMARTSPatternConfig(
                pattern_id="secondary_amine",
                reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]",
                description="For secondary amines"
            )
        ]
    )

    # Create router
    router = SMARTSRouter(config)

    # Register which patterns each reagent is compatible with
    router.register_reagent_compatibility("acetic_acid", {"primary", "secondary_amine"})
    router.register_reagent_compatibility("methylamine", {"primary"})
    router.register_reagent_compatibility("dimethylamine", {"secondary_amine"})

    # Enumerate products - router automatically selects the right pattern
    acid_mol = Chem.MolFromSmiles("CC(=O)O")
    amine_mol = Chem.MolFromSmiles("CNC")  # Secondary amine

    product, pattern_used = router.enumerate(
        reagent_mols=[acid_mol, amine_mol],
        reagent_keys=["acetic_acid", "dimethylamine"]
    )

    if product:
        print(f"Product: {Chem.MolToSmiles(product)}")
        print(f"Pattern used: {pattern_used}")

**Finding Compatible Patterns:**

.. code-block:: python

    from TACTICS.library_enumeration.smarts_toolkit import (
        SMARTSRouter,
        AlternativeSMARTSConfig,
        SMARTSPatternConfig
    )

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

    # Register compatibility
    router.register_reagent_compatibility("reagent_A", {"primary", "secondary_amine"})
    router.register_reagent_compatibility("reagent_B", {"primary"})
    router.register_reagent_compatibility("reagent_C", {"secondary_amine"})

    # Find which pattern works for a combination
    pattern = router.find_compatible_smarts(["reagent_A", "reagent_B"])
    print(f"Compatible pattern: {pattern}")  # "primary"

    pattern = router.find_compatible_smarts(["reagent_A", "reagent_C"])
    print(f"Compatible pattern: {pattern}")  # "secondary_amine"

    # Get compatibility summary
    summary = router.get_compatibility_summary()
    print(f"Reagents per pattern: {summary}")

SMARTSPatternConfig
~~~~~~~~~~~~~~~~~~~

Configuration for a single SMARTS pattern in a router.

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.SMARTSPatternConfig
   :members:
   :undoc-members:
   :show-inheritance:

**Attributes:**

* ``pattern_id``: Unique identifier (e.g., "primary", "secondary_amine")
* ``reaction_smarts``: The reaction SMARTS string
* ``description``: Optional human-readable description
* ``component_compatibility``: Optional dict mapping component_idx to list of substructure SMARTS

**Example:**

.. code-block:: python

    from TACTICS.library_enumeration.smarts_toolkit import SMARTSPatternConfig

    pattern = SMARTSPatternConfig(
        pattern_id="secondary_amine",
        reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]",
        description="For secondary amines that lack a free NH2",
        component_compatibility={
            1: ["[NH]"]  # Component 1 must match this SMARTS
        }
    )

AlternativeSMARTSConfig
~~~~~~~~~~~~~~~~~~~~~~~

Configuration for alternative SMARTS patterns (primary + fallbacks).

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.AlternativeSMARTSConfig
   :members:
   :undoc-members:
   :show-inheritance:

**Attributes:**

* ``primary_smarts``: The default/primary reaction SMARTS
* ``alternatives``: Optional list of :class:`SMARTSPatternConfig` alternative patterns

**Single SMARTS (backwards compatible):**

.. code-block:: python

    from TACTICS.library_enumeration.smarts_toolkit import AlternativeSMARTSConfig

    config = AlternativeSMARTSConfig(
        primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"
    )

    print(config.has_alternatives())  # False

**Multiple SMARTS patterns:**

.. code-block:: python

    from TACTICS.library_enumeration.smarts_toolkit import (
        AlternativeSMARTSConfig,
        SMARTSPatternConfig
    )

    config = AlternativeSMARTSConfig(
        primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
        alternatives=[
            SMARTSPatternConfig(
                pattern_id="secondary_amine",
                reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]"
            ),
            SMARTSPatternConfig(
                pattern_id="tertiary_amine",
                reaction_smarts="[C:1](=O)[OH].[N:2]>>[C:1](=O)[N:2]"
            )
        ]
    )

    print(config.has_alternatives())  # True

    # Get all patterns in priority order (primary first)
    all_patterns = config.get_all_patterns()
    for p in all_patterns:
        print(f"{p.pattern_id}: {p.reaction_smarts}")

Multi-Step Synthesis
--------------------

Support for multi-step synthesis routes with intermediate tracking and protecting group handling.

ReactionSequence
~~~~~~~~~~~~~~~~

Unified handler for single-step and multi-step synthesis.

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.ReactionSequence
   :members:
   :undoc-members:
   :show-inheritance:

**Key Features:**

* Supports three modes: simple SMARTS, SMARTS with alternatives, multi-step
* Uses :class:`SMARTSRouter` at each step for alternative pattern support
* Tracks intermediate products between steps
* Supports protecting group tracking
* Validates required/expected substructures

**Single-Step Usage:**

.. code-block:: python

    from rdkit import Chem
    from TACTICS.library_enumeration.smarts_toolkit import (
        ReactionSequence,
        MultiStepSynthesisConfig
    )

    # Single-step synthesis (mode 1)
    config = MultiStepSynthesisConfig(
        reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
        reagent_file_list=["acids.smi", "amines.smi"]
    )

    sequence = ReactionSequence(config)

    # Enumerate product
    acid_mol = Chem.MolFromSmiles("CC(=O)O")
    amine_mol = Chem.MolFromSmiles("CCN")

    product, patterns_used = sequence.enumerate(
        reagent_mols=[acid_mol, amine_mol],
        reagent_keys=["acetic_acid", "ethylamine"]
    )

    if product:
        print(f"Product: {Chem.MolToSmiles(product)}")
        print(f"Patterns used: {patterns_used}")

**Multi-Step Synthesis:**

.. code-block:: python

    from rdkit import Chem
    from TACTICS.library_enumeration.smarts_toolkit import (
        ReactionSequence,
        MultiStepSynthesisConfig,
        ReactionStepConfig,
        ReactionInputConfig,
        ReactionInputSource,
        AlternativeSMARTSConfig
    )

    # Define a 2-step synthesis
    step1 = ReactionStepConfig(
        step_idx=0,
        smarts_config=AlternativeSMARTSConfig(
            primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"
        ),
        inputs=[
            ReactionInputConfig(source=ReactionInputSource.REAGENT_FILE, component_idx=0),
            ReactionInputConfig(source=ReactionInputSource.REAGENT_FILE, component_idx=1)
        ],
        description="Amide coupling step 1"
    )

    step2 = ReactionStepConfig(
        step_idx=1,
        smarts_config=AlternativeSMARTSConfig(
            primary_smarts="[N:1].[C:2]=O>>[N:1][C:2]"  # Reductive amination
        ),
        inputs=[
            ReactionInputConfig(source=ReactionInputSource.PREVIOUS_STEP, step_idx=0),
            ReactionInputConfig(source=ReactionInputSource.REAGENT_FILE, component_idx=2)
        ],
        description="Reductive amination step 2"
    )

    config = MultiStepSynthesisConfig(
        reagent_file_list=["acids.smi", "amines.smi", "aldehydes.smi"],
        steps=[step1, step2]
    )

    sequence = ReactionSequence(config)

    # Check properties
    print(f"Number of components: {sequence.num_components}")
    print(f"Number of steps: {sequence.num_steps}")

    # Enumerate
    reagent_mols = [
        Chem.MolFromSmiles("CC(=O)O"),    # acid
        Chem.MolFromSmiles("NCCN"),       # diamine
        Chem.MolFromSmiles("CC=O")        # aldehyde
    ]

    product, patterns_used = sequence.enumerate(reagent_mols)

    if product:
        print(f"Final product: {Chem.MolToSmiles(product)}")
        for step, pattern in patterns_used.items():
            print(f"  Step {step}: {pattern}")

MultiStepSynthesisConfig
~~~~~~~~~~~~~~~~~~~~~~~~

Top-level configuration for single or multi-step synthesis routes.

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.MultiStepSynthesisConfig
   :members:
   :undoc-members:
   :show-inheritance:

**Supports three modes:**

1. **Single SMARTS** (backwards compatible): Provide ``reaction_smarts``
2. **Single with alternatives**: Provide ``alternative_smarts``
3. **Multi-step**: Provide ``steps`` list

**Mode 1 - Single SMARTS:**

.. code-block:: python

    from TACTICS.library_enumeration.smarts_toolkit import MultiStepSynthesisConfig

    config = MultiStepSynthesisConfig(
        reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
        reagent_file_list=["acids.smi", "amines.smi"]
    )

    print(config.is_multi_step())  # False

**Mode 2 - Single with alternatives:**

.. code-block:: python

    from TACTICS.library_enumeration.smarts_toolkit import (
        MultiStepSynthesisConfig,
        AlternativeSMARTSConfig,
        SMARTSPatternConfig
    )

    config = MultiStepSynthesisConfig(
        alternative_smarts=AlternativeSMARTSConfig(
            primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
            alternatives=[
                SMARTSPatternConfig(
                    pattern_id="secondary",
                    reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]"
                )
            ]
        ),
        reagent_file_list=["acids.smi", "amines.smi"]
    )

**Mode 3 - Multi-step:**

.. code-block:: python

    from TACTICS.library_enumeration.smarts_toolkit import (
        MultiStepSynthesisConfig,
        ReactionStepConfig,
        ReactionInputConfig,
        ReactionInputSource,
        AlternativeSMARTSConfig
    )

    config = MultiStepSynthesisConfig(
        reagent_file_list=["component_a.smi", "component_b.smi", "component_c.smi"],
        steps=[
            ReactionStepConfig(
                step_idx=0,
                smarts_config=AlternativeSMARTSConfig(primary_smarts="..."),
                inputs=[
                    ReactionInputConfig(source=ReactionInputSource.REAGENT_FILE, component_idx=0),
                    ReactionInputConfig(source=ReactionInputSource.REAGENT_FILE, component_idx=1)
                ]
            ),
            ReactionStepConfig(
                step_idx=1,
                smarts_config=AlternativeSMARTSConfig(primary_smarts="..."),
                inputs=[
                    ReactionInputConfig(source=ReactionInputSource.PREVIOUS_STEP, step_idx=0),
                    ReactionInputConfig(source=ReactionInputSource.REAGENT_FILE, component_idx=2)
                ]
            )
        ]
    )

ReactionStepConfig
~~~~~~~~~~~~~~~~~~

Configuration for a single step in a multi-step synthesis.

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.ReactionStepConfig
   :members:
   :undoc-members:
   :show-inheritance:

**Attributes:**

* ``step_idx``: 0-indexed step number
* ``smarts_config``: :class:`AlternativeSMARTSConfig` for this step
* ``inputs``: List of :class:`ReactionInputConfig` specifying inputs
* ``description``: Optional human-readable description
* ``required_substructure``: Optional SMARTS that must be present in inputs
* ``expected_substructure``: Optional SMARTS that should appear in output

**Example with substructure validation:**

.. code-block:: python

    from TACTICS.library_enumeration.smarts_toolkit import (
        ReactionStepConfig,
        ReactionInputConfig,
        ReactionInputSource,
        AlternativeSMARTSConfig
    )

    step = ReactionStepConfig(
        step_idx=0,
        smarts_config=AlternativeSMARTSConfig(
            primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]"
        ),
        inputs=[
            ReactionInputConfig(source=ReactionInputSource.REAGENT_FILE, component_idx=0),
            ReactionInputConfig(source=ReactionInputSource.REAGENT_FILE, component_idx=1)
        ],
        description="Amide coupling",
        required_substructure="[NH2]",      # Input must have primary amine
        expected_substructure="C(=O)N"       # Product should have amide bond
    )

ReactionInputConfig
~~~~~~~~~~~~~~~~~~~

Configuration for reaction inputs (reagents or intermediates).

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.ReactionInputConfig
   :members:
   :undoc-members:
   :show-inheritance:

**Attributes:**

* ``source``: :class:`ReactionInputSource` (``REAGENT_FILE`` or ``PREVIOUS_STEP``)
* ``component_idx``: Required when source is ``REAGENT_FILE`` (which reagent file)
* ``step_idx``: Required when source is ``PREVIOUS_STEP`` (which step's output)

**Examples:**

.. code-block:: python

    from TACTICS.library_enumeration.smarts_toolkit import (
        ReactionInputConfig,
        ReactionInputSource
    )

    # Input from reagent file (component 0)
    input_from_file = ReactionInputConfig(
        source=ReactionInputSource.REAGENT_FILE,
        component_idx=0
    )

    # Input from previous step's output
    input_from_step = ReactionInputConfig(
        source=ReactionInputSource.PREVIOUS_STEP,
        step_idx=0
    )

ReactionInputSource
~~~~~~~~~~~~~~~~~~~

Enum defining the source of input for a reaction step.

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.ReactionInputSource
   :members:
   :undoc-members:
   :show-inheritance:

**Values:**

* ``REAGENT_FILE``: Input comes from a reagent ``.smi`` file
* ``PREVIOUS_STEP``: Input comes from the output of a previous synthesis step

ProtectingGroupConfig
~~~~~~~~~~~~~~~~~~~~~

Configuration for tracking a protecting group through multi-step synthesis.

.. autoclass:: TACTICS.library_enumeration.smarts_toolkit.ProtectingGroupConfig
   :members:
   :undoc-members:
   :show-inheritance:

**Attributes:**

* ``name``: Human-readable name (e.g., "Boc", "Fmoc")
* ``protected_smarts``: SMARTS pattern to detect protected form
* ``deprotected_smarts``: SMARTS pattern to detect deprotected form
* ``component_idx``: Which reagent component has this protecting group
* ``protection_removed_at_step``: Which step removes the protection

**Example:**

.. code-block:: python

    from TACTICS.library_enumeration.smarts_toolkit import ProtectingGroupConfig

    boc_config = ProtectingGroupConfig(
        name="Boc",
        protected_smarts="[NH]C(=O)OC(C)(C)C",
        deprotected_smarts="[NH2]",
        component_idx=0,
        protection_removed_at_step=1  # Deprotection happens at step 1
    )

Command-Line Interface
----------------------

The SMARTS toolkit includes a command-line interface for validation and testing.

**Basic Usage:**

.. code-block:: bash

    python -m TACTICS.library_enumeration.smarts_toolkit.cli \
        "[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]" \
        acids.smi amines.smi \
        --validate

**With Reaction Testing:**

.. code-block:: bash

    python -m TACTICS.library_enumeration.smarts_toolkit.cli \
        "[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]" \
        acids.smi amines.smi \
        --validate \
        --test-reactions \
        --sample-size 1000

**Export Compatible Reagents:**

.. code-block:: bash

    python -m TACTICS.library_enumeration.smarts_toolkit.cli \
        "[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]" \
        acids.smi amines.smi \
        --export-compatible ./compatible_reagents

**CLI Arguments:**

* ``reaction_smarts``: Required. The reaction SMARTS pattern to validate.
* ``reagent_files``: Required. One or more reagent files (``.smi`` or ``.csv``).
* ``--validate``: Run validation analysis.
* ``--test-reactions``: Test actual reaction combinations (can be slow).
* ``--sample-size``: Number of combinations to test (default: 100).
* ``--export-compatible``: Directory to export compatible reagents.
* ``--output``: Output JSON file for results (default: ``smarts_analysis.json``).

Complete Workflow Example
-------------------------

This example shows a complete workflow for validating a reaction, visualizing results,
and exporting compatible reagents:

.. code-block:: python

    from rdkit import Chem
    from TACTICS.library_enumeration.smarts_toolkit import (
        # Validation
        SMARTSValidator,
        ValidationResult,
        # Visualization
        SMARTSVisualizer,
        # Multi-SMARTS routing
        SMARTSRouter,
        AlternativeSMARTSConfig,
        SMARTSPatternConfig,
        # Multi-step synthesis
        ReactionSequence,
        MultiStepSynthesisConfig
    )

    # 1. Validate the SMARTS pattern against reagent libraries
    validator = SMARTSValidator(
        reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
        reagent_files=["acids.smi", "amines.smi"],
        desalt_on_import=True  # Clean up salt forms
    )

    result = validator.validate(test_reactions=True, sample_size=500)

    print("=== Validation Results ===")
    for pos, coverage in result.coverage_stats.items():
        print(f"Position {pos}: {coverage:.1f}% coverage")
    print(f"Reaction success rate: {result.reaction_success_rate:.1f}%")

    # 2. Visualize results (in Jupyter)
    viz = SMARTSVisualizer(validator)
    viz.visualize_reaction()
    viz.visualize_compatible_reagents(position=0, max_molecules=10)
    viz.generate_summary_plot()

    # 3. Export compatible reagents for use in TACTICS
    validator.export_compatible_reagents(
        output_dir="./cleaned_reagents",
        prefix=["acid", "amine"],
        desalt=True,
        deprotect=True
    )

    # 4. Set up multi-SMARTS routing for production
    config = AlternativeSMARTSConfig(
        primary_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
        alternatives=[
            SMARTSPatternConfig(
                pattern_id="secondary",
                reaction_smarts="[C:1](=O)[OH].[NH:2]>>[C:1](=O)[N:2]"
            )
        ]
    )

    router = SMARTSRouter(config)

    # Register compatibility based on validation
    for smiles, name in result.compatible_reagents.get(1, []):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Check if primary or secondary amine
            if mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]")):
                router.register_reagent_compatibility(name, {"primary"})
            elif mol.HasSubstructMatch(Chem.MolFromSmarts("[NH]")):
                router.register_reagent_compatibility(name, {"secondary"})

    # 5. Use in enumeration
    acid = Chem.MolFromSmiles("CC(=O)O")
    amine = Chem.MolFromSmiles("CNC")

    product, pattern_used = router.enumerate(
        [acid, amine],
        ["acetic_acid", "methylamine"]
    )

    if product:
        print(f"Product: {Chem.MolToSmiles(product)}")
        print(f"Used pattern: {pattern_used}")
