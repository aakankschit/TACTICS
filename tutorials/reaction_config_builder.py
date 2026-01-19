# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "polars",
#     "rdkit",
# ]
# ///
"""
Reaction Configuration Builder

An interactive notebook for building ReactionConfig objects through a UI interface.
Supports single-step, multi-step, and alternative SMARTS configurations.

Run as app: marimo run notebooks/reaction_config_builder.py
Edit mode:  marimo edit notebooks/reaction_config_builder.py
"""

import marimo

__generated_with = "0.18.2"
app = marimo.App(width="full", app_title="TACTICS ReactionConfig Builder")


@app.cell
def _():
    """Imports and project setup."""
    import marimo as mo
    import sys
    from pathlib import Path

    # Add TACTICS project paths
    try:
        project_root = Path(__file__).parent.parent.resolve()
    except NameError:
        project_root = Path("/Users/aakankschitnandkeolyar/Desktop/TACTICS")

    sys.path.insert(0, str(project_root))
    sys.path.insert(0, str(project_root / "src"))

    # TACTICS Library Enumeration imports
    from TACTICS.library_enumeration import SynthesisPipeline, results_to_dataframe
    from TACTICS.library_enumeration.enumeration_utils import (
        failures_to_dataframe,
        summarize_failures,
        read_reagent_file,
    )
    from TACTICS.library_enumeration.smarts_toolkit import (
        ReactionConfig,
        ReactionDef,
        StepInput,
        InputSource,
        DeprotectionSpec,
        DEFAULT_PROTECTING_GROUPS,
    )
    import polars as pl
    from functools import reduce
    import operator

    # Extract protecting group names for dropdown
    PROTECTING_GROUP_NAMES = [pg.name for pg in DEFAULT_PROTECTING_GROUPS]

    # RDKit for visualization
    from rdkit import Chem
    from rdkit.Chem import Draw

    return (
        Chem,
        DeprotectionSpec,
        Draw,
        failures_to_dataframe,
        InputSource,
        mo,
        operator,
        Path,
        pl,
        PROTECTING_GROUP_NAMES,
        read_reagent_file,
        reduce,
        ReactionConfig,
        ReactionDef,
        results_to_dataframe,
        StepInput,
        summarize_failures,
        SynthesisPipeline,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Reaction Configuration Builder

    Build a `ReactionConfig` interactively for use with TACTICS library enumeration.

    **Workflow:**
    1. Select the number of synthesis steps
    2. Configure SMARTS patterns for each step
    3. Specify reagent file paths
    4. Configure step inputs (for multi-step reactions)
    5. Click **Build ReactionConfig** to validate and preview
    6. Click **Enumerate Full Library** to generate all products as a DataFrame
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Step 1: Number of Steps

    Select how many reaction steps in your synthesis.
    """)
    return


@app.cell
def _(mo):
    """Number of steps dropdown."""
    num_steps_dropdown = mo.ui.dropdown(
        options=["1", "2", "3", "4", "5"],
        value="1",
        label="Number of Steps"
    )
    return (num_steps_dropdown,)


@app.cell
def _(num_steps_dropdown):
    """Display the dropdown."""
    num_steps_dropdown
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Step 2: Reaction SMARTS Configuration

    Define the SMARTS pattern for each reaction step.
    Optionally add alternative SMARTS for substrate-specific routing.
    """)
    return


@app.cell
def _(mo, num_steps_dropdown, PROTECTING_GROUP_NAMES):
    """Create dynamic widgets for each step."""
    _num_steps = int(num_steps_dropdown.value)

    # Create widgets for each step
    step_smarts_inputs = []
    step_description_inputs = []
    step_alt_checkboxes = []
    step_alt_smarts_inputs = []
    step_deprot_group_dropdowns = []
    step_deprot_target_dropdowns = []

    # Deprotection options: "None" + all protecting group names
    _deprot_options = ["None"] + PROTECTING_GROUP_NAMES

    for _i in range(_num_steps):
        step_smarts_inputs.append(
            mo.ui.text(
                value="",
                placeholder="e.g., [C:1](=O)O.[N:2]>>[C:1](=O)[N:2]",
                label=f"Step {_i} SMARTS",
                full_width=True
            )
        )
        step_description_inputs.append(
            mo.ui.text(
                value="",
                placeholder="e.g., Amide coupling",
                label=f"Description"
            )
        )
        step_alt_checkboxes.append(
            mo.ui.checkbox(value=False, label="Add alternative SMARTS")
        )
        step_alt_smarts_inputs.append(
            mo.ui.text(
                value="",
                placeholder="Alternative SMARTS pattern",
                label=f"Step {_i} Alt SMARTS",
                full_width=True
            )
        )
        # Deprotection widgets - dropdown with "None" option to disable
        step_deprot_group_dropdowns.append(
            mo.ui.dropdown(
                options=_deprot_options,
                value="None",
                label="Deprotection"
            )
        )
        step_deprot_target_dropdowns.append(
            mo.ui.dropdown(
                options=["Reactant 0", "Reactant 1", "Product"],
                value="Reactant 0",
                label="Deprotect Target"
            )
        )

    return (
        step_smarts_inputs,
        step_description_inputs,
        step_alt_checkboxes,
        step_alt_smarts_inputs,
        step_deprot_group_dropdowns,
        step_deprot_target_dropdowns,
    )


@app.cell
def _(mo, num_steps_dropdown, step_smarts_inputs, step_description_inputs, step_alt_checkboxes, step_alt_smarts_inputs, step_deprot_group_dropdowns, step_deprot_target_dropdowns):
    """Display step configuration widgets."""
    _num_steps = int(num_steps_dropdown.value)

    _step_sections = []
    for _i in range(_num_steps):
        # Alternative SMARTS row (only if checkbox is checked)
        _alt_row = mo.hstack([step_alt_smarts_inputs[_i]], justify="start") if step_alt_checkboxes[_i].value else None

        # Deprotection row - always show both dropdowns
        # User selects "None" in the first dropdown to skip deprotection
        _deprot_row = mo.hstack([
            step_deprot_group_dropdowns[_i],
            step_deprot_target_dropdowns[_i],
        ], gap=2)

        # Build section
        _section_items = [
            mo.md(f"### Step {_i}"),
            step_smarts_inputs[_i],
            mo.hstack([step_description_inputs[_i], step_alt_checkboxes[_i]], justify="start", gap=2),
        ]
        if _alt_row:
            _section_items.append(_alt_row)
        _section_items.append(_deprot_row)

        _section = mo.vstack(_section_items, gap=1)
        _step_sections.append(_section)

    mo.vstack(_step_sections, gap=2)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Step 3: Reagent Files

    Specify the paths to reagent SMILES files. By default, paths to the bundled
    Thrombin dataset are provided. Replace these with your own file paths as needed.
    """)
    return


@app.cell
def _(mo):
    """Number of reagent files dropdown."""
    num_files_dropdown = mo.ui.dropdown(
        options=["2", "3", "4", "5", "6"],
        value="2",
        label="Number of Reagent Files"
    )
    return (num_files_dropdown,)


@app.cell
def _(num_files_dropdown):
    """Display the dropdown."""
    num_files_dropdown
    return


@app.cell
def _(mo, num_files_dropdown):
    """Create reagent file path inputs."""
    import importlib.resources

    _num_files = int(num_files_dropdown.value)

    # Get bundled data paths for default examples
    _data_files = importlib.resources.files("TACTICS.data.thrombin")
    _example_paths = [
        str(_data_files / "acids.smi"),
        str(_data_files / "coupled_aa_sub.smi"),
        str(_data_files / "amino_acids_no_fmoc.smi"),
        "",
        "",
        "",
    ]

    reagent_file_inputs = []
    for _i in range(_num_files):
        reagent_file_inputs.append(
            mo.ui.text(
                value=_example_paths[_i] if _i < len(_example_paths) else "",
                placeholder=f"Path to reagent file {_i}",
                label=f"File {_i}",
                full_width=True
            )
        )

    return (reagent_file_inputs,)


@app.cell
def _(mo, reagent_file_inputs):
    """Display reagent file inputs."""
    mo.vstack([
        mo.hstack([inp], justify="start") for inp in reagent_file_inputs
    ], gap=1)
    return


@app.cell(hide_code=True)
def _(mo, num_steps_dropdown):
    """Step inputs section header - only shown for multi-step."""
    _num_steps = int(num_steps_dropdown.value)
    if _num_steps > 1:
        _output = mo.md(r"""
        ---
        ## Step 4: Step Inputs Configuration

        For multi-step reactions, configure where each step gets its inputs from.
        Each input can come from a **Reagent File** or from a **Previous Step**.
        """)
    else:
        _output = mo.md("")
    _output
    return


@app.cell
def _(mo, num_steps_dropdown, num_files_dropdown):
    """Create step inputs configuration widgets."""
    _num_steps = int(num_steps_dropdown.value)
    _num_files = int(num_files_dropdown.value)

    # For each step, create input source configuration
    # Each step typically has 2 inputs (for a bimolecular reaction)
    step_input_widgets = {}

    if _num_steps > 1:
        for _step_idx in range(_num_steps):
            step_input_widgets[_step_idx] = {
                "input0_source": mo.ui.dropdown(
                    options=["REAGENT_FILE", "PREVIOUS_STEP"],
                    value="REAGENT_FILE" if _step_idx == 0 else "PREVIOUS_STEP",
                    label=f"Step {_step_idx} Input 0 Source"
                ),
                "input0_index": mo.ui.dropdown(
                    options=[str(_j) for _j in range(max(_num_files, _num_steps))],
                    value="0" if _step_idx == 0 else str(_step_idx - 1),
                    label="Index"
                ),
                "input1_source": mo.ui.dropdown(
                    options=["REAGENT_FILE", "PREVIOUS_STEP"],
                    value="REAGENT_FILE",
                    label=f"Step {_step_idx} Input 1 Source"
                ),
                "input1_index": mo.ui.dropdown(
                    options=[str(_j) for _j in range(max(_num_files, _num_steps))],
                    value=str(min(_step_idx + 1, _num_files - 1)),
                    label="Index"
                ),
            }

    return (step_input_widgets,)


@app.cell
def _(mo, num_steps_dropdown, step_input_widgets):
    """Display step inputs configuration."""
    _num_steps = int(num_steps_dropdown.value)

    if _num_steps > 1:
        _sections = []
        for _step_idx in range(_num_steps):
            _widgets = step_input_widgets[_step_idx]
            _section = mo.vstack([
                mo.md(f"**Step {_step_idx} Inputs:**"),
                mo.hstack([
                    mo.vstack([
                        mo.md("Input 0:"),
                        _widgets["input0_source"],
                        _widgets["input0_index"],
                    ]),
                    mo.vstack([
                        mo.md("Input 1:"),
                        _widgets["input1_source"],
                        _widgets["input1_index"],
                    ]),
                ], gap=4),
            ], gap=1)
            _sections.append(_section)

        _output = mo.vstack(_sections, gap=2)
    else:
        _output = mo.md("")
    _output
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Step 5: Build Configuration

    Click the button below to build the `ReactionConfig` and preview sample products.
    """)
    return


@app.cell
def _(mo):
    """Build button."""
    build_button = mo.ui.run_button(label="Build ReactionConfig")
    return (build_button,)


@app.cell
def _(build_button):
    """Display build button."""
    build_button
    return


@app.cell
def _(
    mo,
    build_button,
    num_steps_dropdown,
    num_files_dropdown,
    step_smarts_inputs,
    step_description_inputs,
    step_alt_checkboxes,
    step_alt_smarts_inputs,
    step_deprot_group_dropdowns,
    step_deprot_target_dropdowns,
    reagent_file_inputs,
    step_input_widgets,
    ReactionConfig,
    ReactionDef,
    StepInput,
    InputSource,
    DeprotectionSpec,
    Path,
):
    """Build the ReactionConfig from widget values."""
    # Initialize outputs
    config = None
    config_error = None

    # Only build if button clicked
    if build_button.value:
        _num_steps = int(num_steps_dropdown.value)
        _num_files = int(num_files_dropdown.value)

        # Collect reactions
        _reactions = []
        _step_modes = {}
        _validation_error = None

        for _i in range(_num_steps):
            if _validation_error:
                break

            _smarts = step_smarts_inputs[_i].value.strip()
            if not _smarts:
                _validation_error = f"Step {_i} SMARTS is required"
                break

            _desc = step_description_inputs[_i].value.strip() or None

            # Build deprotections list if a protecting group is selected (not "None")
            _deprotections = []
            _group = step_deprot_group_dropdowns[_i].value
            if _group != "None":
                _target_str = step_deprot_target_dropdowns[_i].value
                # Convert target string to appropriate value
                if _target_str == "Product":
                    _target = "product"
                elif _target_str == "Reactant 0":
                    _target = 0
                elif _target_str == "Reactant 1":
                    _target = 1
                else:
                    _target = 0
                _deprotections.append(DeprotectionSpec(group=_group, target=_target))

            # Primary reaction
            _reactions.append(
                ReactionDef(
                    reaction_smarts=_smarts,
                    step_index=_i,
                    pattern_id="primary" if step_alt_checkboxes[_i].value else None,
                    description=_desc,
                    deprotections=_deprotections,
                )
            )

            # Alternative reaction if enabled
            if step_alt_checkboxes[_i].value:
                _alt_smarts = step_alt_smarts_inputs[_i].value.strip()
                if _alt_smarts:
                    _reactions.append(
                        ReactionDef(
                            reaction_smarts=_alt_smarts,
                            step_index=_i,
                            pattern_id="alternative",
                            description=f"{_desc} (alternative)" if _desc else "Alternative",
                            deprotections=_deprotections,
                        )
                    )
                    _step_modes[_i] = "alternative"

        if _validation_error:
            config_error = _validation_error
        else:
            # Collect reagent files
            _reagent_files = []
            for _i in range(_num_files):
                _path = reagent_file_inputs[_i].value.strip()
                if _path:
                    _reagent_files.append(_path)

            if not _reagent_files:
                config_error = "At least one reagent file is required"
            else:
                # Build step_inputs
                _step_inputs = None
                if _num_steps > 1:
                    _step_inputs = {}
                    for _step_idx in range(_num_steps):
                        _widgets = step_input_widgets[_step_idx]
                        _inputs = []

                        # Input 0
                        _src0 = _widgets["input0_source"].value
                        _idx0 = int(_widgets["input0_index"].value)
                        if _src0 == "REAGENT_FILE":
                            _inputs.append(StepInput(source=InputSource.REAGENT_FILE, file_index=_idx0))
                        else:
                            _inputs.append(StepInput(source=InputSource.PREVIOUS_STEP, step_index=_idx0))

                        # Input 1
                        _src1 = _widgets["input1_source"].value
                        _idx1 = int(_widgets["input1_index"].value)
                        if _src1 == "REAGENT_FILE":
                            _inputs.append(StepInput(source=InputSource.REAGENT_FILE, file_index=_idx1))
                        else:
                            _inputs.append(StepInput(source=InputSource.PREVIOUS_STEP, step_index=_idx1))

                        _step_inputs[_step_idx] = _inputs

                # Build the config
                try:
                    config = ReactionConfig(
                        reactions=_reactions,
                        reagent_file_list=_reagent_files,
                        step_inputs=_step_inputs,
                        step_modes=_step_modes if _step_modes else None,
                    )
                except Exception as e:
                    config_error = str(e)

    config, config_error


@app.cell
def _(mo, config, config_error):
    """Display the generated config."""
    mo.stop(config is None and config_error is None, mo.md("*Click 'Build ReactionConfig' to generate the configuration*"))

    if config_error:
        _output = mo.md(f"""
        ## Configuration Error

        ```
        {config_error}
        ```
        """)
    else:
        _output = mo.md(f"""
        ## Generated ReactionConfig

        **Properties:**
        - Number of steps: {config.num_steps}
        - Is multi-step: {config.is_multi_step}
        - Steps with alternatives: {config.steps_with_alternatives}

        **Reactions:**
        """)
    _output
    return


@app.cell
def _(mo, config, config_error):
    """Display reaction details."""
    mo.stop(config is None or config_error is not None)

    _lines = []
    for _i, _rxn in enumerate(config.reactions):
        _lines.append(f"- **Reaction {_i}** (step {_rxn.step_index}): `{_rxn.reaction_smarts[:50]}...`")
        if _rxn.description:
            _lines.append(f"  - Description: {_rxn.description}")
        if _rxn.pattern_id:
            _lines.append(f"  - Pattern ID: {_rxn.pattern_id}")
        if _rxn.deprotections:
            for _deprot in _rxn.deprotections:
                _target_str = f"Reactant {_deprot.target}" if isinstance(_deprot.target, int) else "Product"
                _lines.append(f"  - Deprotection: Remove **{_deprot.group}** from {_target_str}")

    mo.md("\n".join(_lines))
    return


@app.cell
def _(mo, config, config_error, read_reagent_file, reduce, operator, Path):
    """Display reagent files and library statistics."""
    mo.stop(config is None or config_error is not None)

    # Resolve file paths and count reagents
    try:
        _project_root = Path(__file__).parent.parent.resolve()
    except NameError:
        _project_root = Path("/Users/aakankschitnandkeolyar/Desktop/TACTICS")

    _lines = ["**Reagent Files:**"]
    _reagent_counts = []

    for _i, _f in enumerate(config.reagent_file_list):
        _p = Path(_f)
        if not _p.is_absolute():
            _p = _project_root / _f

        # Try to count reagents if file exists
        _count = 0
        if _p.exists():
            try:
                _reagents = read_reagent_file(str(_p))
                _count = len(_reagents)
            except Exception:
                _count = 0

        _reagent_counts.append(_count)
        _lines.append(f"- File {_i}: `{_f}` — **{_count:,}** reagents")

    # Calculate expected total products
    if _reagent_counts and all(_c > 0 for _c in _reagent_counts):
        _expected_products = reduce(operator.mul, _reagent_counts, 1)
        _lines.append("")
        _lines.append(f"**Expected Library Size:** {_expected_products:,} products")
    else:
        _expected_products = 0

    mo.md("\n".join(_lines))
    return


@app.cell(hide_code=True)
def _(mo, config, config_error):
    """Sample enumeration header."""
    mo.stop(config is None or config_error is not None)

    mo.md(r"""
    ---
    ## Sample Products

    Enumerating a few sample products from the configured library.
    """)
    return


@app.cell
def _(mo, config, config_error, SynthesisPipeline, results_to_dataframe, summarize_failures, failures_to_dataframe, Chem, Draw, Path):
    """Generate and display sample products."""
    mo.stop(config is None or config_error is not None, mo.md("*Build the ReactionConfig first to see sample products*"))

    # Resolve file paths
    try:
        _project_root = Path(__file__).parent.parent.resolve()
    except NameError:
        _project_root = Path("/Users/aakankschitnandkeolyar/Desktop/TACTICS")

    # Update reagent file paths to absolute
    _abs_files = []
    for _f in config.reagent_file_list:
        _p = Path(_f)
        if not _p.is_absolute():
            _p = _project_root / _f
        _abs_files.append(str(_p))

    # Check if files exist
    _missing_files = [_f for _f in _abs_files if not Path(_f).exists()]
    if _missing_files:
        mo.stop(True, mo.md(f"**Error:** Files not found: {_missing_files}"))

    # Create config with absolute paths
    from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig as _RC

    _config_with_abs_paths = _RC(
        reactions=config.reactions,
        reagent_file_list=_abs_files,
        step_inputs=config.step_inputs,
        step_modes=config.step_modes,
    )

    # Create pipeline and enumerate
    try:
        _pipeline = SynthesisPipeline(_config_with_abs_paths)
        _results = _pipeline.enumerate_library()

        # Get success/failure summary
        _summary = summarize_failures(_results)
        sample_total = _summary["total"]
        sample_successes = _summary["successes"]
        sample_failures_count = _summary["failures"]
        sample_failure_rate = _summary["failure_rate"]

        # Convert successful results to DataFrame
        _products_df = results_to_dataframe(_results)

        # Get failures DataFrame
        sample_failures_df = failures_to_dataframe(_results)

        # Detect duplicate SMILES
        if len(_products_df) > 0:
            # Find SMILES that appear more than once
            _smiles_counts = _products_df.group_by("product_smiles").agg(
                pl.count().alias("count"),
                pl.col("product_name").alias("reagent_combinations")
            ).filter(pl.col("count") > 1)

            if len(_smiles_counts) > 0:
                # Get all rows with duplicate SMILES
                _duplicate_smiles = _smiles_counts["product_smiles"].to_list()
                sample_duplicates_df = _products_df.filter(
                    pl.col("product_smiles").is_in(_duplicate_smiles)
                ).sort("product_smiles")
                sample_unique_smiles_count = len(_products_df["product_smiles"].unique())
            else:
                sample_duplicates_df = None
                sample_unique_smiles_count = len(_products_df)
        else:
            sample_duplicates_df = None
            sample_unique_smiles_count = 0

        # Get first 6 successful products for display
        if len(_products_df) > 0:
            _sample_df = _products_df.head(6)
            _sample_smiles = _sample_df["product_smiles"].to_list()
            _sample_names = _sample_df["product_name"].to_list()

            # Convert to molecules
            _mols = []
            _legends = []
            for _i, (_smiles, _name) in enumerate(zip(_sample_smiles, _sample_names)):
                _mol = Chem.MolFromSmiles(_smiles)
                if _mol:
                    _mols.append(_mol)
                    # Use product name (reagent combination) as legend
                    _legends.append(_name if _name else f"Product {_i+1}")

            if _mols:
                _img = Draw.MolsToGridImage(_mols, molsPerRow=3, subImgSize=(300, 200), legends=_legends)
                sample_image = _img
            else:
                sample_image = None
        else:
            sample_image = None

        sample_error = None

    except Exception as e:
        sample_image = None
        sample_total = 0
        sample_successes = 0
        sample_failures_count = 0
        sample_failure_rate = 0.0
        sample_failures_df = None
        sample_duplicates_df = None
        sample_unique_smiles_count = 0
        sample_error = str(e)

    return sample_image, sample_total, sample_successes, sample_failures_count, sample_failure_rate, sample_failures_df, sample_duplicates_df, sample_unique_smiles_count, sample_error


@app.cell
def _(mo, config, config_error, sample_image, sample_total, sample_successes, sample_failures_count, sample_failure_rate, sample_failures_df, sample_duplicates_df, sample_unique_smiles_count, sample_error):
    """Display sample products, failure summary, and duplicates."""
    mo.stop(config is None or config_error is not None)

    if sample_error:
        _output = mo.callout(
            mo.md(f"""
            **Enumeration Error:**
            ```
            {sample_error}
            ```
            """),
            kind="danger"
        )
    else:
        # Calculate duplicate info
        _num_duplicates = sample_successes - sample_unique_smiles_count if sample_successes > 0 else 0

        # Build summary section
        _summary_items = [
            mo.md(f"""
            **Enumeration Summary:**
            - Total attempted: **{sample_total:,}**
            - Successful: **{sample_successes:,}** ✓
            - Failed: **{sample_failures_count:,}** ✗
            - Failure rate: **{sample_failure_rate:.1f}%**
            - Unique SMILES: **{sample_unique_smiles_count:,}** ({_num_duplicates:,} duplicates)
            """),
        ]

        # Add sample products image if available
        if sample_image is not None:
            _summary_items.append(mo.md("**Sample Products:**"))
            _summary_items.append(sample_image)

        # Add duplicates section if there are any
        if sample_duplicates_df is not None and len(sample_duplicates_df) > 0:
            _summary_items.append(mo.md("---"))
            _summary_items.append(
                mo.callout(
                    mo.md(f"**{_num_duplicates:,} duplicate SMILES detected.** Different reagent combinations produced identical products."),
                    kind="info"
                )
            )
            _summary_items.append(mo.md("**Duplicate Products (same SMILES, different reagents):**"))
            # Show first 30 duplicates
            _display_dup_df = sample_duplicates_df.head(30)
            _summary_items.append(_display_dup_df)

            if len(sample_duplicates_df) > 30:
                _summary_items.append(mo.md(f"*... and {len(sample_duplicates_df) - 30:,} more duplicate entries*"))

        # Add failures section if there are any
        if sample_failures_count > 0 and sample_failures_df is not None and len(sample_failures_df) > 0:
            _summary_items.append(mo.md("---"))
            _summary_items.append(
                mo.callout(
                    mo.md(f"**{sample_failures_count:,} enumeration failures detected.** See details below."),
                    kind="warn"
                )
            )
            _summary_items.append(mo.md("**Failed Reagent Combinations:**"))
            # Show first 20 failures
            _display_df = sample_failures_df.head(20).select([
                "product_name",
                "error_type",
                "reagent_names",
                "reagent_smiles",
                "message"
            ])
            _summary_items.append(_display_df)

            if len(sample_failures_df) > 20:
                _summary_items.append(mo.md(f"*... and {len(sample_failures_df) - 20:,} more failures*"))

        _output = mo.vstack(_summary_items, gap=2)

    _output
    return


@app.cell(hide_code=True)
def _(mo):
    """Full enumeration header - always visible."""
    mo.md(r"""
    ---
    ## Full Library Enumeration

    Enumerate the entire combinatorial library and view all products as a DataFrame.
    """)
    return


@app.cell
def _(mo):
    """Enumeration controls - always created."""
    n_jobs_input = mo.ui.number(
        value=1,
        start=1,
        stop=16,
        step=1,
        label="Number of Jobs (parallel workers)"
    )
    enumerate_button = mo.ui.run_button(label="Enumerate Full Library")
    return n_jobs_input, enumerate_button


@app.cell
def _(mo, n_jobs_input, enumerate_button):
    """Display enumeration controls - always visible."""
    mo.hstack([n_jobs_input, enumerate_button], gap=2)
    return


@app.cell
def _(
    mo,
    config,
    config_error,
    enumerate_button,
    n_jobs_input,
    SynthesisPipeline,
    results_to_dataframe,
    summarize_failures,
    failures_to_dataframe,
    Path,
):
    """Enumerate full library and convert to DataFrame."""
    # Check if config is built
    if config is None or config_error is not None:
        mo.stop(True, mo.callout(
            mo.md("**Build the ReactionConfig first** by clicking the 'Build ReactionConfig' button above."),
            kind="warn"
        ))
    mo.stop(not enumerate_button.value, mo.md("*Click 'Enumerate Full Library' to enumerate all products*"))

    # Resolve file paths
    try:
        _project_root = Path(__file__).parent.parent.resolve()
    except NameError:
        _project_root = Path("/Users/aakankschitnandkeolyar/Desktop/TACTICS")

    # Update reagent file paths to absolute
    _abs_files = []
    for _f in config.reagent_file_list:
        _p = Path(_f)
        if not _p.is_absolute():
            _p = _project_root / _f
        _abs_files.append(str(_p))

    # Check if files exist
    _missing_files = [_f for _f in _abs_files if not Path(_f).exists()]
    if _missing_files:
        mo.stop(True, mo.md(f"**Error:** Files not found: {_missing_files}"))

    # Create config with absolute paths
    from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig as _RC

    _config_with_abs_paths = _RC(
        reactions=config.reactions,
        reagent_file_list=_abs_files,
        step_inputs=config.step_inputs,
        step_modes=config.step_modes,
    )

    # Create pipeline and enumerate full library
    try:
        _pipeline = SynthesisPipeline(_config_with_abs_paths)
        _n_jobs = int(n_jobs_input.value)
        _results = _pipeline.enumerate_library(n_jobs=_n_jobs, show_progress=True)

        # Get success/failure summary
        _summary = summarize_failures(_results)
        enum_total = _summary["total"]
        enum_successes = _summary["successes"]
        enum_failures_count = _summary["failures"]
        enum_failure_rate = _summary["failure_rate"]

        # Convert successful results to DataFrame
        library_df = results_to_dataframe(_results)

        # Get failures DataFrame
        enum_failures_df = failures_to_dataframe(_results)

        # Detect duplicate SMILES
        if len(library_df) > 0:
            # Find SMILES that appear more than once
            _smiles_counts = library_df.group_by("product_smiles").agg(
                pl.count().alias("count"),
                pl.col("product_name").alias("reagent_combinations")
            ).filter(pl.col("count") > 1)

            if len(_smiles_counts) > 0:
                # Get all rows with duplicate SMILES
                _duplicate_smiles = _smiles_counts["product_smiles"].to_list()
                enum_duplicates_df = library_df.filter(
                    pl.col("product_smiles").is_in(_duplicate_smiles)
                ).sort("product_smiles")
                enum_unique_smiles_count = len(library_df["product_smiles"].unique())
            else:
                enum_duplicates_df = None
                enum_unique_smiles_count = len(library_df)
        else:
            enum_duplicates_df = None
            enum_unique_smiles_count = 0

        enum_error = None

    except Exception as e:
        library_df = None
        enum_total = 0
        enum_successes = 0
        enum_failures_count = 0
        enum_failure_rate = 0.0
        enum_failures_df = None
        enum_duplicates_df = None
        enum_unique_smiles_count = 0
        enum_error = str(e)

    return library_df, enum_total, enum_successes, enum_failures_count, enum_failure_rate, enum_failures_df, enum_duplicates_df, enum_unique_smiles_count, enum_error


@app.cell
def _(mo, config, config_error, enumerate_button, library_df, enum_total, enum_successes, enum_failures_count, enum_failure_rate, enum_failures_df, enum_duplicates_df, enum_unique_smiles_count, enum_error):
    """Display the enumerated library DataFrame, duplicates, and failures."""
    mo.stop(config is None or config_error is not None)
    mo.stop(not enumerate_button.value)

    if enum_error:
        _output = mo.callout(
            mo.md(f"""
            **Enumeration Error:**
            ```
            {enum_error}
            ```
            """),
            kind="danger"
        )
    else:
        _items = []

        # Calculate duplicate info
        _num_duplicates = enum_successes - enum_unique_smiles_count if enum_successes > 0 else 0

        # Summary statistics
        if enum_failures_count == 0 and _num_duplicates == 0:
            _items.append(
                mo.callout(
                    mo.md(f"**Successfully enumerated {enum_successes:,} unique products** (100% success rate)"),
                    kind="success"
                )
            )
        else:
            _items.append(
                mo.callout(
                    mo.md(f"""
                    **Enumeration Complete:**
                    - Total attempted: **{enum_total:,}**
                    - Successful: **{enum_successes:,}** ✓
                    - Failed: **{enum_failures_count:,}** ✗
                    - Success rate: **{100 - enum_failure_rate:.1f}%**
                    - Unique SMILES: **{enum_unique_smiles_count:,}** ({_num_duplicates:,} duplicates)
                    """),
                    kind="info"
                )
            )

        # Products DataFrame
        _items.append(mo.md("### Successful Products"))
        _items.append(library_df)

        # Duplicates section if there are any
        if enum_duplicates_df is not None and len(enum_duplicates_df) > 0:
            _items.append(mo.md("---"))
            _items.append(
                mo.callout(
                    mo.md(f"**{_num_duplicates:,} duplicate SMILES detected.** Different reagent combinations produced identical products."),
                    kind="info"
                )
            )
            _items.append(mo.md("### Duplicate Products (same SMILES, different reagents)"))
            _items.append(enum_duplicates_df)

        # Failures section if there are any
        if enum_failures_count > 0 and enum_failures_df is not None and len(enum_failures_df) > 0:
            _items.append(mo.md("---"))
            _items.append(
                mo.callout(
                    mo.md(f"**{enum_failures_count:,} enumeration failures detected**"),
                    kind="warn"
                )
            )
            _items.append(mo.md("### Failed Reagent Combinations"))
            _items.append(enum_failures_df)

        _output = mo.vstack(_items, gap=2)

    _output
    return


if __name__ == "__main__":
    app.run()
