# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "polars",
#     "rdkit",
# ]
# ///
"""
Library Enumeration & Thompson Sampling Demo

This notebook demonstrates the TACTICS library enumeration pipeline and its
integration with Thompson Sampling for efficient combinatorial library screening.

Topics covered:
1. SynthesisPipeline basics - single compound enumeration
2. Full library enumeration
3. Multi-step and alternative SMARTS patterns
4. Thompson Sampling integration
5. Using presets for common workflows
"""

import marimo

__generated_with = "0.18.2"
app = marimo.App(width="medium")

with app.setup:
    # Initialization code that runs before all other cells
    import marimo as mo
    import polars as pl
    from pathlib import Path
    from rdkit import Chem
    from rdkit.Chem import Draw
    # TACTICS imports
    # Synthesis Pipeline Imports
    from TACTICS.library_enumeration import SynthesisPipeline
    from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef, StepInput, InputSource, DeprotectionSpec
    # Thompson Sampler Imports
    from TACTICS.thompson_sampling import ThompsonSamplingConfig
    from TACTICS.thompson_sampling.strategies.config import GreedyConfig
    from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig
    from TACTICS.thompson_sampling.core.sampler import ThompsonSampler


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # Library Enumeration and SMARTS Toolkit Demonstration

    A tutorial showcasing the SMARTS Toolkit of TACTICS and how these toolks can be used to troubleshoot SMARTS patterns for library enumeration.
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## Overview of Linear Amide Library
    The target for this library is Thrombin, but this library can be enumerated in different ways and serves as a good example of the capabilities of this package.

    The Linear Amide library consists of 3 different reagent lists-
    * Amino Acids
    * Amino Acids
    * Carboxylic Acids

    This library consists of a dipeptide coupled to a carboxylic acid. The dipeptide is formed by the coupling of the amino acids.
    """)
    return


@app.cell
def _():
    mo.md(r"""
    # Library Enumeration & Thompson Sampling Demo

    This notebook demonstrates the **TACTICS** library enumeration pipeline and how it
    integrates with Thompson Sampling for efficient screening of combinatorial chemical libraries.

    ## Overview

    The `SynthesisPipeline` is the central class for:
    - **Single compound generation** (used by ThompsonSampler)
    - **Batch enumeration** (parallel processing)
    - **Full library enumeration** (all combinations)
    - **Reagent validation** and compatibility detection

    ## Architecture

    ```
    ReactionConfig (defines reactions + reagent files)
        |
        v
    SynthesisPipeline (executes reactions)
        |
        v
    ThompsonSamplingConfig (optimization settings)
        |
        v
    ThompsonSampler (runs the search)
    ```
    """)
    return


@app.cell
def _():
    mo.md(r"""
    ## 1. Setup: Create Sample Reagent Files

    We will load the files from the data folder.
    """)
    return


@app.function
def read_file_to_list(filepath):
    """
    Read a file and return its contents as a list.
    Each line is split by spaces if multiple values exist.

    Args:
        filepath: Path to the input file

    Returns:
        List of data (nested list if lines contain multiple values)
    """
    data = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line:  # skip empty lines
                parts = line.split()
                if len(parts) == 1:
                    data.append(parts[0])
                else:
                    data.append(parts)
    return data


@app.cell
def _():
    # Start by loading the amino acids and the carboxulic acids
    # Using absolute paths to avoid issues while executing but change these paths depending on where the data files are located

    # These amino acids are fmoc deprotected but have other protecting groups on them
    AMINO_ACIDS_PTH = "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/reagents/thrombin/amino_acids_no_fmoc.smi"
    # Carboxylic acids
    ACIDS_PTH = (
        "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/reagents/thrombin/acids.smi"
    )
    # Dipeptides with the an amine substitution to be representative of the products used in experiment
    DIPEPTIDES_PTH = "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/reagents/thrombin/coupled_aa_sub.smi"
    return ACIDS_PTH, AMINO_ACIDS_PTH, DIPEPTIDES_PTH


@app.cell
def _():
    # Reference File Path for Thompson Sampling

    # This file was generated using docking evaluator with Thrombin as the target of interest.
    SCORES_PTH = "/Users/aakankschitnandkeolyar/Desktop/TACTICS/data/scores/thrombin/product_scores.csv"
    return (SCORES_PTH,)


@app.cell
def _(ACIDS_PTH, AMINO_ACIDS_PTH, DIPEPTIDES_PTH, SCORES_PTH):
    # Laod all Data into lists
    amino_acids = read_file_to_list(AMINO_ACIDS_PTH)
    dipeptides = read_file_to_list(DIPEPTIDES_PTH)
    acids = read_file_to_list(ACIDS_PTH)
    scores = pl.read_csv(SCORES_PTH)
    return acids, amino_acids, dipeptides


@app.cell
def _():
    mo.md(r"""
    ## 2. SynthesisPipeline Basics

    The `SynthesisPipeline` wraps a `ReactionConfig` and provides methods for
    enumerating products from reagent combinations.
    """)
    return


@app.cell
def _():
    # Define the SMARTS Pattern for Amide coupling
    # This is used to couple the dipeptide and the carboxylic acid
    AMIDE_COUPLING = "[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]>>[#6:1](=[O:2])[#7:3]"
    return (AMIDE_COUPLING,)


@app.cell
def _(ACIDS_PTH, AMIDE_COUPLING, DIPEPTIDES_PTH):
    # Create the reaction configuration
    rxn_config = ReactionConfig(
        reactions=[
            ReactionDef(
                reaction_smarts=AMIDE_COUPLING,
                step_index=0,
                description="Amide coupling",
            )
        ],
        reagent_file_list=[ACIDS_PTH, DIPEPTIDES_PTH],
    )

    # Create the synthesis pipeline
    pipeline = SynthesisPipeline(rxn_config)

    print(f"Pipeline created:")
    print(f"  - Steps: {pipeline.num_steps}")
    print(f"  - Components: {pipeline.num_components}")
    print(f"  - Has alternatives: {pipeline.has_alternatives}")
    print(f"  - Reagent files: {[Path(f).name for f in pipeline.reagent_file_list]}")
    return (pipeline,)


@app.cell
def _():
    mo.md(r"""
    ### 2.1 Single Compound Enumeration
    """)
    return


@app.cell
def _(acids, dipeptides, pipeline):
    # Create molecules from SMILES
    acid_mol = Chem.MolFromSmiles(acids[0][0])  # Acetic acid
    dipeptide_mol = Chem.MolFromSmiles(dipeptides[0][0])  # Ethylamine

    # Enumerate the product
    result = pipeline.enumerate_single(
        reagent_mols=[acid_mol, dipeptide_mol],
        reagent_keys=[acids[0][1], dipeptides[0][1]],
    )

    print(f"Enumeration result:")
    print(f"  - Success: {result.success}")
    print(f"  - Product SMILES: {result.product_smiles}")
    print(f"  - Product name: {result.product_name}")
    print(f"  - Patterns used: {result.patterns_used}")
    return acid_mol, dipeptide_mol, result


@app.cell
def _(acid_mol, acids, dipeptide_mol, dipeptides, result):
    # Visualize the structure to see whether they are correct
    Draw.MolsToGridImage(
        [acid_mol, dipeptide_mol, Chem.MolFromSmiles(result.product_smiles)],
        molsPerRow=3,
        subImgSize=(500, 500),
        legends=[acids[0][1], dipeptides[0][1], result.product_name],
    )
    return


@app.cell
def _(acids, dipeptides, pipeline):
    # You can also enumerate from SMILES directly
    result2 = pipeline.enumerate_single_from_smiles(
        smiles_list=[acids[0][0], dipeptides[0][0]],  # Benzoic acid + Aniline
        reagent_keys=[acids[0][1], dipeptides[0][1]],
    )

    print(f"Benzanilide synthesis:")
    print(f"  - Success: {result2.success}")
    print(f"  - Product: {result2.product_smiles}")
    return (result2,)


@app.cell
def _(acid_mol, acids, dipeptide_mol, dipeptides, result2):
    # Visualize the structure to see whether they are correct
    Draw.MolsToGridImage(
        [acid_mol, dipeptide_mol, Chem.MolFromSmiles(result2.product_smiles)],
        molsPerRow=3,
        subImgSize=(500, 500),
        legends=[acids[0][1], dipeptides[0][1], result2.product_name],
    )
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    Both methods provide the same structures. The user can directly input a list of the SMILES from the reagent lists for enumeration. Or if the user has already generated `rdkit` molecules then that can also be input into the synthesis pipeline.
    """)
    return


@app.cell
def _():
    mo.md(r"""
    ### 2.2 Handling Failed Reactions
    """)
    return


@app.cell
def _(pipeline):
    # Try with an incompatible reagent (no NH2 group)
    acid = Chem.MolFromSmiles("CC(=O)O")
    not_amine = Chem.MolFromSmiles("CCC")  # No amine group

    failed_result = pipeline.enumerate_single(
        reagent_mols=[acid, not_amine],
        reagent_keys=["acetic_acid", "not_amine"],
    )

    print(f"Failed enumeration:")
    print(f"  - Success: {failed_result.success}")
    print(
        f"  - Error type: {failed_result.error.error_type if failed_result.error else None}"
    )
    print(
        f"  - Error message: {failed_result.error.message if failed_result.error else None}"
    )
    return


@app.cell
def _():
    mo.md(r"""
    ## 3. Full Library Enumeration

    The pipeline can enumerate all possible products from the reagent files.
    """)
    return


@app.cell
def _(pipeline):
    # Enumerate all combinations (3844 dipeptides x  130 carboxylic acids = 499720 products)
    # By changing the number of jobs > 1 multiprocessing can be used enumerate other larger libraries faster
    all_results = pipeline.enumerate_library(n_jobs=1, show_progress=False)

    # Count successes and failures
    successes = [r for r in all_results if r.success]
    failures = [r for r in all_results if not r.success]

    print(f"Library enumeration complete:")
    print(f"  - Total combinations: {len(all_results)}")
    print(f"  - Successful: {len(successes)}")
    print(f"  - Failed: {len(failures)}")

    # Show first few products
    print(f"\nFirst 5 products:")
    for r in successes[:5]:
        print(f"  - {r.product_name}: {r.product_smiles}")
    return (all_results,)


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    You can use `result.error` incase the enumeration fails, this is a custom error object that provides the user information on where the enumeration failed and why it failed.
    """)
    return


@app.cell
def _():
    mo.md(r"""
    ### 3.1 Convert Results to DataFrame
    """)
    return


@app.cell
def _(all_results):
    from TACTICS.library_enumeration import results_to_dataframe

    # Convert to Polars DataFrame
    df = results_to_dataframe(all_results)
    print(f"DataFrame shape: {df.shape}")
    df.head(10)
    return


@app.cell
def _():
    mo.md(r"""
    ## 4. Thompson Sampling Integration

    The `SynthesisPipeline` is designed to work seamlessly with Thompson Sampling.
    The pipeline is passed to `ThompsonSamplingConfig` as the single source of truth
    for reaction definitions and reagent files.
    """)
    return


@app.cell
def _(SCORES_PTH, pipeline):
    # Create the Thompson Sampling configuration
    ts_config = ThompsonSamplingConfig(
        synthesis_pipeline=pipeline,
        num_ts_iterations=50,
        num_warmup_trials=3,
        strategy_config=GreedyConfig(mode="minimize"),
        evaluator_config=LookupEvaluatorConfig(ref_filename=SCORES_PTH),
        hide_progress=True,
    )

    print("Thompson Sampling Configuration:")
    print(f"  - Reagent files: {[Path(f).name for f in ts_config.reagent_file_list]}")
    print(f"  - Num components: {ts_config.num_components}")
    print(f"  - Num steps: {ts_config.num_steps}")
    print(f"  - Iterations: {ts_config.num_ts_iterations}")
    print(f"  - Strategy: {ts_config.strategy_config.strategy_type}")
    print(f"  - Mode: {ts_config.strategy_config.mode}")
    return (ts_config,)


@app.cell
def _():
    mo.md(r"""
    ### 4.1 Running Thompson Sampling
    """)
    return


@app.cell
def _(ts_config):
    # Create and run the sampler
    sampler = ThompsonSampler.from_config(ts_config)

    print(f"Sampler created:")
    print(f"  - Reagent lists: {len(sampler.reagent_lists)}")
    print(f"  - Total possible products: {sampler.num_prods}")

    # Run warmup phase
    warmup_results = sampler.warm_up(num_warmup_trials=2)
    print(f"\nWarmup complete: {len(warmup_results)} evaluations")

    # Run search phase
    search_results = sampler.search(num_cycles=30)
    print(f"Search complete: {len(search_results)} evaluations")
    return search_results, warmup_results


@app.cell
def _(search_results, warmup_results):
    # Combine results
    all_ts_results = pl.concat([warmup_results, search_results])

    # Find top hits
    top_hits = all_ts_results.sort("score", descending=True).head(10)

    print("Top 10 hits from Thompson Sampling:")
    print(top_hits)
    return


@app.cell
def _():
    mo.md(r"""
    ## 5. Using Presets

    TACTICS provides pre-configured presets for common Thompson Sampling scenarios.
    """)
    return


@app.cell
def _(SCORES_PTH, pipeline):
    from TACTICS.thompson_sampling.presets import get_preset, ConfigPresets

    # List available presets
    print("Available presets:")
    print("  - fast_exploration: Quick screening with epsilon-greedy")
    print("  - parallel_batch: Parallel processing for slow evaluators")
    print("  - conservative_exploit: Greedy exploitation")
    print("  - balanced_sampling: UCB for exploration/exploitation balance")
    print("  - diverse_coverage: Maximum diversity with roulette wheel")
    print("  - legacy_rws_maximize: Original RWS algorithm (maximize)")
    print("  - legacy_rws_minimize: Original RWS algorithm (minimize)")

    # Get a preset configuration
    preset_config = get_preset(
        "fast_exploration",
        synthesis_pipeline=pipeline,
        evaluator_config=LookupEvaluatorConfig(ref_filename=SCORES_PTH),
        num_iterations=100,
        mode="maximize",
    )

    print(f"\nFast exploration preset:")
    print(f"  - Strategy: {preset_config.strategy_config.strategy_type}")
    print(f"  - Warmup: {preset_config.warmup_config.warmup_type}")
    print(f"  - Iterations: {preset_config.num_ts_iterations}")
    return


@app.cell
def _():
    mo.md(r"""
    ## 6. Advanced: Multi-Step Synthesis

    `SynthesisPipeline` supports multi-step synthesis where the product of one
    reaction becomes a reagent for the next step.

    ### Example: 3 Component Linear Amide Library

    Previously we had used a 2 step synthesis plan, but we can use a 3 step synthesis with the same library as well. Instead of using dipeptides, we can use amino acids directly without having to couple them together and then couple a carboxylic acid.

    **Note**: This set up uses the same amino acid reagents twice to form dipeptides. The algorithm treats them as different reagent lists but they contain the same set of amino acids.

    1. **Step 0**: amino acid + amino acid → dipeptide (with free amine)
    2. **Step 1**: dipeptide (from step 0) + carboxylic acid → linear amide product

    This requires:
    - 3 reagent files (amino acid, amino acids, carboxylic acids)
    - 2 reaction definitions with `step_index` 0 and 1
    - `step_inputs` to specify where each reactant comes from
    """)
    return


@app.cell
def _():
    mo.md(r"""
    ### 6.1 Creating the Multi-Step ReactionConfig

    The goal of this class is to provide an organized method for the user to set up a synthesis pipeline. This mimics synthesis plans used by medicinal chemists when generating synthesizing compound libraries in a lab setting.

    For multi-step synthesis, we need to specify:
    - Multiple `ReactionDef` objects with different `step_index` values. These objects are the basis for any reaction, the required inputs are a SMARTS pattern. The step index tells the `ReactionConfig` which step of the synthesis it is. The description can be used to name the reaction.
    - `step_inputs` dict mapping step indices to input sources. This maps the reactions defined to the reagent files.
        - `StepInput` objects specifying whether input comes from a reagent file or previous step. This is used to define the items in the `step_inputs` dictionary.
    - `mode` tells the reaction config whether there are alternative SMARTS patterns for a specific step in the reaction. These are cases where the reaction for that specific step are carried out by two separate SMARTS patterns to cover all of the reagents for a specific reaction component.

    By separating the reagent file inputs, SMARTS and modes
    """)
    return


@app.cell
def _():
    # Define Peptide Coupling SMARTS
    PEPTIDE_COUPLING = "[#6X3:1](=[#8X1])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):2]>>[#6X3:1](=[#8X1])[#7X3:2]"
    return (PEPTIDE_COUPLING,)


@app.cell
def _(ACIDS_PTH, AMIDE_COUPLING, AMINO_ACIDS_PTH, PEPTIDE_COUPLING):
    # Define the two-step synthesis
    multi_step_config = ReactionConfig(
        reactions=[
            # Step 0: amino acid + amino acid → dipeptide (with free amine)
            # peptide coupling to form dipeptide
            ReactionDef(
                reaction_smarts=PEPTIDE_COUPLING,
                step_index=0,
                description="peptide coupling",
            ),
            ReactionDef(
                reaction_smarts="[#6:1](=O)[OH]>>[#6:1](=O)[NH2]",
                step_index=1,
                description="Amine substitution",
            ),
            ReactionDef(
                reaction_smarts=AMIDE_COUPLING,
                step_index=2,
                description="linear amide coupling",
            ),
        ],
        reagent_file_list=[
            AMINO_ACIDS_PTH,  # Index 0: First set of Amino acids
            AMINO_ACIDS_PTH,  # Index 1: Second set of Amino acids
            ACIDS_PTH,  # Index 2: Carboxylic acids
        ],
        # Define where each step gets its inputs
        step_inputs={
            # Step 0: both inputs from reagent files
            0: [
                StepInput(source=InputSource.REAGENT_FILE, file_index=0),  # Amino Acids Set 1
                StepInput(source=InputSource.REAGENT_FILE, file_index=1),  # Amino Acids Set 2
            ],
            1: [
                StepInput(source=InputSource.PREVIOUS_STEP, step_index=0)
            ],
            # Step 1: first input from step 0's product, second from reagent file
            2: [
                StepInput(source=InputSource.REAGENT_FILE, file_index=2),
                StepInput(source=InputSource.PREVIOUS_STEP, step_index=1)
            ],
        },
    )

    print("Multi-step ReactionConfig created:")
    print(f"  - Number of reactions: {len(multi_step_config.reactions)}")
    print(f"  - Number of steps: {multi_step_config.num_steps}")
    print(f"  - Is multi-step: {multi_step_config.is_multi_step}")
    print(f"  - Reagent files: {len(multi_step_config.reagent_file_list)}")
    return (multi_step_config,)


@app.cell
def _():
    mo.md(r"""
    ### 6.2 Creating and Using the Multi-Step Pipeline
    """)
    return


@app.cell
def _(multi_step_config):
    # Create the multi-step pipeline
    multi_pipeline = SynthesisPipeline(multi_step_config)

    print(f"Multi-step pipeline created:")
    print(f"  - Steps: {multi_pipeline.num_steps}")
    print(f"  - Components (reagent files): {multi_pipeline.num_components}")
    print(f"  - Is multi-step: {multi_pipeline.is_multi_step}")
    print(f"  - Has alternatives: {multi_pipeline.has_alternatives}")
    return (multi_pipeline,)


@app.cell
def _():
    mo.md(r"""
    ### 6.3 Enumerating Multi-Step Products

    When enumerating, we provide molecules for all reagent positions.
    The pipeline automatically handles the intermediate products.
    """)
    return


@app.cell
def _(acids, amino_acids, multi_pipeline):
    # Create reagent molecules
    amino_acid_mol_1 = Chem.MolFromSmiles(amino_acids[0][0])  # Acetic acid
    amino_acid_mol_2 = Chem.MolFromSmiles(amino_acids[0][0])  # Ethylenediamine
    carboxlic_acid_mol_1 = Chem.MolFromSmiles(acids[0][0])  # Benzaldehyde

    # Enumerate with intermediate storage to see the steps
    ms_result = multi_pipeline.enumerate_single(
        reagent_mols=[amino_acid_mol_1, amino_acid_mol_2, carboxlic_acid_mol_1],
        reagent_keys=[amino_acids[0][1], amino_acids[0][1], acids[0][1]],
        store_intermediates=True,
    )

    print("Multi-step enumeration result:")
    print(f"  - Success: {ms_result.success}")
    print(f"  - Final product: {ms_result.product_smiles}")
    print(f"  - Product name: {ms_result.product_name}")
    print(f"  - Patterns used per step: {ms_result.patterns_used}")

    if ms_result.intermediates:
        print(f"\nIntermediates:")
        for step_idx, intermediate in ms_result.intermediates.items():
            intermediate_smiles = Chem.MolToSmiles(intermediate)
            print(f"  - Step {step_idx}: {intermediate_smiles}")
    return amino_acid_mol_1, amino_acid_mol_2, carboxlic_acid_mol_1, ms_result


@app.cell
def _(
    acids,
    amino_acid_mol_1,
    amino_acid_mol_2,
    amino_acids,
    carboxlic_acid_mol_1,
    dipeptides,
    ms_result,
    result2,
):
    # Visualize the structure to see whether they are correct
    Draw.MolsToGridImage(
        [amino_acid_mol_1, amino_acid_mol_2, Chem.MolFromSmiles("CC1(C)Oc2ccc(C[C@H](NC(=O)[C@@H](N)Cc3ccc4c(c3)OC(C)(C)O4)C(N)=O)cc2O1"), Chem.MolFromSmiles(dipeptides[0][0]), carboxlic_acid_mol_1, Chem.MolFromSmiles(ms_result.product_smiles), Chem.MolFromSmiles(result2.product_smiles)],
        molsPerRow=4,
        subImgSize=(500, 500),
        legends=[amino_acids[0][1], amino_acids[0][1], "intermediate Product", dipeptides[0][1], acids[0][1], ms_result.product_name, "Product from 2 component Library"],
    )
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    In the above block, we can see a comparison of the intermediate products at each reaction step. I have compared the dipeptide generated by the `SynthesisPipeline` to one that I had manually generated. They are identical. When generating the dipeptide manually I also had to write another reaction to substitute the hydroxyl on the carboxylic acid of the dipeptide to an amine, but with the `SynthesisPipeline` this has been made possible in a single step. I have also placed the final product generated from the two component library (which used manually curated list of dipeptides), the `SynthesisPipeline` is able to recreate the product in a single set of code.
    """)
    return


@app.cell
def _():
    mo.md(r"""
    ### 6.4 Multi-Step Library Enumeration

    We can enumerate all combinations across the multi-step synthesis.
    With 2 acids × 2 diamines × 2 aldehydes = 8 possible products.
    """)
    return


@app.cell
def _(multi_pipeline):
    # Enumerate the full multi-step library
    # Demonstrate Multi-processing here as well
    ms_library_results = multi_pipeline.enumerate_library(n_jobs=5, show_progress=True)

    # Analyze results
    ms_successes = [r for r in ms_library_results if r.success]
    ms_failures = [r for r in ms_library_results if not r.success]

    print(f"Multi-step library enumeration:")
    print(f"  - Total combinations: {len(ms_library_results)}")
    print(f"  - Successful: {len(ms_successes)}")
    print(f"  - Failed: {len(ms_failures)}")

    print(f"\nSuccessful products:")
    for n,res in enumerate(ms_successes):
        print(f"  - {res.product_name}: {res.product_smiles}")
        if n==50:
            break
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### 6.5 Using Alternative SMARTS Pattern in a Multi-Step Reaction Pipeline

    Previously enumerated the library with a multi-step synthesis plan using the `SynthesisPipeline`. Now we will demonstrate the enumeration with the use of alternative reaction SMARTS pattern. The `SynthesisPipeline` can also handle multiple SMARTS patterns for a single step.
    """)
    return


@app.cell
def _(ACIDS_PTH, AMIDE_COUPLING, AMINO_ACIDS_PTH):
    multi_step_alt_config = ReactionConfig(
        reactions=[
            # Step 0: amino acid + amino acid → dipeptide (with free amine)
            # peptide coupling to form dipeptide
            ReactionDef(
                reaction_smarts="[#6X3:1](=[#8X1])[OH].[#7X3;H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):2]>>[#6X3:1](=[#8X1])[#7X3:2]",
                step_index=0,
                pattern_id="primary",
                description="primary amine peptide coupling",
            ),
            ReactionDef(
                reaction_smarts="[#6X3:1](=[#8X1])[OH].[#7X3;H1;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):2]>>[#6X3:1](=[#8X1])[#7X3:2]",
                step_index=0,
                description="secondary amine peptide coupling",
            ),
            ReactionDef(
                reaction_smarts="[#6:1](=O)[OH]>>[#6:1](=O)[NH2]",
                step_index=1,
                description="Amine substitution",
            ),
            ReactionDef(
                reaction_smarts=AMIDE_COUPLING,
                step_index=2,
                description="linear amide coupling",
            ),
        ],
        reagent_file_list=[
            AMINO_ACIDS_PTH,  # Index 0: First set of Amino acids
            AMINO_ACIDS_PTH,  # Index 1: Second set of Amino acids
            ACIDS_PTH,  # Index 2: Carboxylic acids
        ],
        # Define where each step gets its inputs
        step_inputs={
            # Step 0: both inputs from reagent files
            0: [
                StepInput(source=InputSource.REAGENT_FILE, file_index=0),  # Amino Acids Set 1
                StepInput(source=InputSource.REAGENT_FILE, file_index=1),  # Amino Acids Set 2
            ],
            1: [
                StepInput(source=InputSource.PREVIOUS_STEP, step_index=0)
            ],
            # Step 1: first input from step 0's product, second from reagent file
            2: [
                StepInput(source=InputSource.REAGENT_FILE, file_index=2),
                StepInput(source=InputSource.PREVIOUS_STEP, step_index=1)
            ],
        },
        step_modes={0:"alternative"}
    )

    print("Multi-step ReactionConfig created:")
    print(f"  - Number of reactions: {len(multi_step_alt_config.reactions)}")
    print(f"  - Number of steps: {multi_step_alt_config.num_steps}")
    print(f"  - Is multi-step: {multi_step_alt_config.is_multi_step}")
    print(f"  - Reagent files: {len(multi_step_alt_config.reagent_file_list)}")
    return (multi_step_alt_config,)


@app.cell
def _(multi_step_alt_config):
    alt_multi_pipeline = SynthesisPipeline(multi_step_alt_config)

    print(f"Multi-step pipeline created:")
    print(f"  - Steps: {alt_multi_pipeline.num_steps}")
    print(f"  - Components (reagent files): {alt_multi_pipeline.num_components}")
    print(f"  - Is multi-step: {alt_multi_pipeline.is_multi_step}")
    print(f"  - Has alternatives: {alt_multi_pipeline.has_alternatives}")
    return (alt_multi_pipeline,)


@app.cell
def _(
    acids,
    alt_multi_pipeline,
    amino_acid_mol_1,
    amino_acid_mol_2,
    amino_acids,
    carboxlic_acid_mol_1,
):
    # Reacreate Results for Primary amine molecule
    alt_ms_result = alt_multi_pipeline.enumerate_single(
        reagent_mols=[amino_acid_mol_1, amino_acid_mol_2, carboxlic_acid_mol_1],
        reagent_keys=[amino_acids[0][1], amino_acids[0][1], acids[0][1]],
        store_intermediates=True,
    )

    print("Multi-step enumeration result:")
    print(f"  - Success: {alt_ms_result.success}")
    print(f"  - Final product: {alt_ms_result.product_smiles}")
    print(f"  - Product name: {alt_ms_result.product_name}")
    print(f"  - Patterns used per step: {alt_ms_result.patterns_used}")

    if alt_ms_result.intermediates:
        print(f"\nIntermediates:")
        for alt_step_idx_1, alt_intermediate_1 in alt_ms_result.intermediates.items():
            alt_intermediate_smiles_1 = Chem.MolToSmiles(alt_intermediate_1)
            print(f"  - Step {alt_step_idx_1}: {alt_intermediate_smiles_1}")
    return (alt_ms_result,)


@app.cell
def _(
    acids,
    alt_ms_result,
    amino_acid_mol_1,
    amino_acid_mol_2,
    amino_acids,
    carboxlic_acid_mol_1,
    dipeptides,
    ms_result,
):
    # Visualize the structure to see whether they are correct
    Draw.MolsToGridImage(
        [amino_acid_mol_1, amino_acid_mol_2, Chem.MolFromSmiles("CC1(C)Oc2ccc(C[C@H](NC(=O)[C@@H](N)Cc3ccc4c(c3)OC(C)(C)O4)C(N)=O)cc2O1"), Chem.MolFromSmiles(dipeptides[0][0]), carboxlic_acid_mol_1, Chem.MolFromSmiles(ms_result.product_smiles), Chem.MolFromSmiles(alt_ms_result.product_smiles)],
        molsPerRow=4,
        subImgSize=(500, 500),
        legends=[amino_acids[0][1], amino_acids[0][1], "intermediate Product", dipeptides[0][1], acids[0][1], ms_result.product_name, "Product from 2 component Library"],
    )
    return


@app.cell
def _(
    acids,
    alt_ms_result,
    alt_multi_pipeline,
    amino_acids,
    carboxlic_acid_mol_1,
):
    amino_acid_mol_3 = Chem.MolFromSmiles(amino_acids[33][0])
    amino_acid_mol_4 = Chem.MolFromSmiles(amino_acids[34][0])
    alt_ms_result_2 = alt_multi_pipeline.enumerate_single(
        reagent_mols=[amino_acid_mol_4, amino_acid_mol_3, carboxlic_acid_mol_1],
        reagent_keys=[amino_acids[34][1], amino_acids[33][1], acids[0][1]],
        store_intermediates=True,
    )

    print("Multi-step enumeration result:")
    print(f"  - Success: {alt_ms_result_2.success}")
    print(f"  - Final product: {alt_ms_result_2.product_smiles}")
    print(f"  - Product name: {alt_ms_result_2.product_name}")
    print(f"  - Patterns used per step: {alt_ms_result_2.patterns_used}")

    if alt_ms_result.intermediates:
        print(f"\nIntermediates:")
        for alt_step_idx_2, alt_intermediate_2 in alt_ms_result_2.intermediates.items():
            alt_intermediate_smiles_2 = Chem.MolToSmiles(alt_intermediate_2)
            print(f"  - Step {alt_step_idx_2}: {alt_intermediate_smiles_2}")
    return alt_ms_result_2, amino_acid_mol_3, amino_acid_mol_4


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    Notice in the patterns used here, it has used the alternative pattern for enumeration. The `SynthesisPipeline` loops through all the patterns, hence if a reagent cannot be validated using the primary pattern it will move to the next one and so on. It will pickup a SMARTS pattern that fits. The block below demonstrates the use with a amino acid that has a secondary amine.
    """)
    return


@app.cell
def _(
    acids,
    alt_ms_result_2,
    amino_acid_mol_3,
    amino_acid_mol_4,
    amino_acids,
    carboxlic_acid_mol_1,
    dipeptides,
):
    # Visualize the structure to see whether they are correct
    Draw.MolsToGridImage(
        [amino_acid_mol_4, amino_acid_mol_3, Chem.MolFromSmiles("NC(=O)C1CN(C(=O)[C@H](N)Cc2ccc(Br)cc2)C1"), Chem.MolFromSmiles(dipeptides[2141][0]), carboxlic_acid_mol_1, Chem.MolFromSmiles(alt_ms_result_2.product_smiles)],
        molsPerRow=4,
        subImgSize=(500, 500),
        legends=[amino_acids[34][1], amino_acids[33][1], "intermediate Product", dipeptides[2141][1], acids[0][1], alt_ms_result_2.product_name],
    )
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    It appears to be working correctly for primary amines as can be seen from the visualization above. The question is whether it will work correctly for secondary amines.
    """)
    return


@app.cell
def _():
    mo.md(r"""
    ### 6.5 Understanding step_inputs

    The `step_inputs` dictionary is key to multi-step synthesis:

    ```python
    step_inputs = {
        # Step 0: Both inputs from reagent files
        0: [
            StepInput(source=InputSource.REAGENT_FILE, file_index=0),  # From reagent file 0
            StepInput(source=InputSource.REAGENT_FILE, file_index=1),  # From reagent file 1
        ],
        # Step 1: First input is the product of step 0
        1: [
            StepInput(source=InputSource.PREVIOUS_STEP, step_index=0),  # Product from step 0
            StepInput(source=InputSource.REAGENT_FILE, file_index=2),   # From reagent file 2
        ],
    }
    ```

    **InputSource options:**
    - `REAGENT_FILE`: Input comes from a reagent file (specify `file_index`)
    - `PREVIOUS_STEP`: Input comes from a previous reaction's product (specify `step_index`)
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## 7. Protecting Groups

    The SMARTS Toolkit comes equipped with pre-computed deprotections. This is similar to `rdkit` but it also allows the user to add their own custom deprotection units as well.

    Default protecting groups and their respective reaction SMARTS are listed in the table below, additionally the code snippet to show the protecting groups is also provided.

    ```python

    from TACTICS.library_enumeration.smarts_toolkit.constants import (
          DEFAULT_PROTECTING_GROUPS,
          PROTECTING_GROUP_MAP,
          get_protecting_group,
          get_all_protecting_group_names,
      )

      # List all available protecting groups
      print(get_all_protecting_group_names())
      # ['Boc', 'Fmoc', 'Cbz', 'Acetamide', 'TBS', 'O-Benzyl', 'Trityl', 'tBu-ester', 'Me-ester', 'Et-ester']

      # Get a specific protecting group
      boc = get_protecting_group("Boc")
      print(f"Name: {boc.name}")
      print(f"Detection SMARTS: {boc.smarts}")
      print(f"Deprotection SMARTS: {boc.deprotection_smarts}")

    ┌───────────┬─────────────────┬────────────────────────────────────────┬─────────────────────┐
    │   Name    │    Protects     │            Detection SMARTS            │ Deprotection Result │
    ├───────────┼─────────────────┼────────────────────────────────────────┼─────────────────────┤
    │ Boc       │ Amine (N)       │ [NX3][C](=O)OC(C)(C)C                  │ Free amine          │
    ├───────────┼─────────────────┼────────────────────────────────────────┼─────────────────────┤
    │ Fmoc      │ Amine (N)       │ [NX3]C(=O)OCC1c2ccccc2-c2ccccc12       │ Free amine          │
    ├───────────┼─────────────────┼────────────────────────────────────────┼─────────────────────┤
    │ Cbz       │ Amine (N)       │ [NX3]C(=O)OCc1ccccc1                   │ Free amine          │
    ├───────────┼─────────────────┼────────────────────────────────────────┼─────────────────────┤
    │ Acetamide │ Amine (N)       │ [NX3][C](=O)[CH3]                      │ Free amine          │
    ├───────────┼─────────────────┼────────────────────────────────────────┼─────────────────────┤
    │ TBS       │ Alcohol (O)     │ [OX2][Si](C)(C)C(C)(C)C                │ Free alcohol        │
    ├───────────┼─────────────────┼────────────────────────────────────────┼─────────────────────┤
    │ O-Benzyl  │ Alcohol (O)     │ [OX2]Cc1ccccc1                         │ Free alcohol        │
    ├───────────┼─────────────────┼────────────────────────────────────────┼─────────────────────┤
    │ Trityl    │ Amine/Alcohol   │ [NX3,OX2]C(c1ccccc1)(c1ccccc1)c1ccccc1 │ N/A (complex)       │
    ├───────────┼─────────────────┼────────────────────────────────────────┼─────────────────────┤
    │ tBu-ester │ Carboxylic acid │ [CX3](=O)OC(C)(C)C                     │ Free acid (-COOH)   │
    ├───────────┼─────────────────┼────────────────────────────────────────┼─────────────────────┤
    │ Me-ester  │ Carboxylic acid │ [CX3](=O)O[CH3]                        │ Free acid (-COOH)   │
    ├───────────┼─────────────────┼────────────────────────────────────────┼─────────────────────┤
    │ Et-ester  │ Carboxylic acid │ [CX3](=O)OCC                           │ Free acid (-COOH)   │
    └───────────┴─────────────────┴────────────────────────────────────────┴─────────────────────┘
    ```
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### 7.1 Specificying Deprotection Schemes in `SynthesisPipeline`

    TACTICS utilizes `DeprotectionSpec` to specify when to deprotect and what kind of protecting groups should be deprotected from the molecule.

    ```python
      from TACTICS.library_enumeration.smarts_toolkit import DeprotectionSpec

      # Deprotect a REACTANT (before reaction runs)
      deprot_reactant = DeprotectionSpec(
          group="Boc",      # Which protecting group
          target=0          # Reactant index (0 = first reactant)
      )

      # Deprotect the PRODUCT (after reaction runs)
      deprot_product = DeprotectionSpec(
          group="Fmoc",
          target="product"  # Special keyword for product
      )

      # Check properties
      print(deprot_reactant.is_product_deprotection)  # False
      print(deprot_reactant.reactant_index)           # 0

      print(deprot_product.is_product_deprotection)   # True
      print(deprot_product.reactant_index)            # None

      Fields:
      ┌────────┬─────────────────┬───────────────────────────────────────────────┐
      │ Field  │      Type       │                  Description                  │
      ├────────┼─────────────────┼───────────────────────────────────────────────┤
      │ group  │ str             │ Name of protecting group to remove            │
      ├────────┼─────────────────┼───────────────────────────────────────────────┤
      │ target │ int | "product" │ Reactant index OR "product" for post-reaction │
      └────────┴─────────────────┴───────────────────────────────────────────────┘

    ```
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### 7.2 Deprotections with custom groups

    Use the `ProtectingGroupInfo` class to add custom protecting groups to be used in the `SynthesisPipeline` class.

    ```python
    from TACTICS.library_enumeration.smarts_toolkit.config import ProtectingGroupInfo

    # Define custom protecting group
    alloc = ProtectingGroupInfo(
      name="Alloc",
      smarts="[NX3]C(=O)OCC=C",
      deprotection_smarts="[N:1]C(=O)OCC=C>>[N:1]"
    )

    config = ReactionConfig(
      reactions=[
          ReactionDef(
              reaction_smarts="...",
              step_index=0,
              deprotections=[DeprotectionSpec(group="Alloc", target="product")],
          ),
      ],
      reagent_file_list=["..."],
      protecting_groups=[alloc],  # Add custom protecting groups here
    )


    ```
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### 7.3 All Deprotection Schemes can be used with `ReactionDef`.

    The user can specify whether the product or the reactant should be deprotected in the given reaction step.

    ```python
      from TACTICS.library_enumeration import (
          ReactionConfig,
          ReactionDef,
          DeprotectionSpec,
          get_all_protecting_group_names,  # To see available groups
      )

      # See available protecting groups
      print(get_all_protecting_group_names())
      # ['Boc', 'Fmoc', 'Cbz', 'Bn', 'TBS', 'TBDPS', 'Tr', 'PMB', 'Alloc', 'Ns']

      # Deprotection on reactant (before reaction)
      config = ReactionConfig(
          reactions=[
              ReactionDef(
                  reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                  step_index=0,
                  deprotections=[
                      DeprotectionSpec(group="Boc", target=1),  # Deprotect reactant 1
                  ],
              )
          ],
          reagent_file_list=["acids.smi", "boc_amines.smi"],
      )

      # Deprotection on product (after reaction)
      config = ReactionConfig(
          reactions=[
              ReactionDef(
                  reaction_smarts="[C:1](=O)[OH].[NH2:2]>>[C:1](=O)[NH:2]",
                  step_index=0,
                  deprotections=[
                      DeprotectionSpec(group="Fmoc", target="product"),  # Deprotect product
                  ],
              )
          ],
          reagent_file_list=["acids.smi", "amines.smi"],
      )


    ```
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## Summary

    This notebook demonstrated:

    1. **SynthesisPipeline** - The central class for compound enumeration
       - `enumerate_single()` - Single compound generation from RDKit Mol objects
       - `enumerate_single_from_smiles()` - Single compound generation from SMILES
       - `enumerate_library()` - Full library enumeration with parallelization
       - Direct construction: `SynthesisPipeline(ReactionConfig(...))`

    2. **ReactionConfig** - Configuration for reactions and reagent files
       - Single-step reactions with `ReactionDef`
       - Alternative SMARTS patterns with `step_modes`
       - Multi-step synthesis with `step_inputs`

    3. **Deprotection** - Protecting group removal during synthesis
       - Reactant deprotection: `DeprotectionSpec(group="Boc", target=0)`
       - Product deprotection: `DeprotectionSpec(group="Fmoc", target="product")`
       - Custom protecting groups via `ProtectingGroupInfo`

    4. **Thompson Sampling Integration**
       - Pipeline is passed to `ThompsonSamplingConfig`
       - Provides reagent files and compound generation
       - Works with all selection strategies

    5. **Presets** - Pre-configured settings for common workflows
       - `fast_exploration`, `parallel_batch`, `balanced_sampling`, etc.

    6. **File Writing** - Export enumerated libraries
       - `write_enumerated_library()` - Write to CSV, SMI, or SDF
       - `write_products_chunked()` - Split large libraries into multiple files
       - `results_to_dataframe()` - Convert results to Polars DataFrame

    For more details, see the [TACTICS documentation](https://tactics.readthedocs.io).
    """)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
