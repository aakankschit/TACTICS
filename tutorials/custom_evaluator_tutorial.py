"""
CustomEvaluator Demo (Merimo)

Run:
    marimo run notebooks/custom_evaluator_demo.py
Edit:
    marimo edit notebooks/custom_evaluator_demo.py
"""

import marimo

__generated_with = "0.19.4"
app = marimo.App(width="medium", app_title="CustomEvaluator Demo")


@app.cell
def _():
    """Imports and dataset paths."""
    import marimo as mo
    from rdkit import Chem
    import polars as pl
    import importlib.resources

    from TACTICS.thompson_sampling.core.evaluators import CustomEvaluator
    from TACTICS.thompson_sampling.core.evaluator_config import CustomEvaluatorConfig
    from TACTICS.thompson_sampling import ThompsonSampler
    from TACTICS.thompson_sampling.config import ThompsonSamplingConfig
    from TACTICS.thompson_sampling.strategies.config import GreedyConfig
    from TACTICS.thompson_sampling.warmup.config import StandardWarmupConfig
    from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef
    from TACTICS.library_enumeration import SynthesisPipeline

    # Thrombin dataset reagent paths
    _data_files = importlib.resources.files("TACTICS.data.thrombin")
    ACIDS_FILE = str(_data_files / "acids.smi")
    AMINES_FILE = str(_data_files / "coupled_aa_sub.smi")
    REAGENT_FILES = [ACIDS_FILE, AMINES_FILE]

    # Amide coupling SMARTS
    AMIDE_COUPLING_SMARTS = "[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]>>[#6:1](=[O:2])[#7:3]"
    return (
        AMIDE_COUPLING_SMARTS,
        Chem,
        CustomEvaluator,
        CustomEvaluatorConfig,
        GreedyConfig,
        REAGENT_FILES,
        ReactionConfig,
        ReactionDef,
        StandardWarmupConfig,
        SynthesisPipeline,
        ThompsonSampler,
        ThompsonSamplingConfig,
        mo,
        pl,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    # CustomEvaluator & CustomEvaluatorConfig Demo

    This notebook demonstrates using a **CustomEvaluator** and a **CustomEvaluatorConfig**.

    The scoring function returns the **number of heavy atoms** in an RDKit molecule.

    We will:
    1. Test on 5 small SMILES
    2. Run a minimal Thrombin TS using full reagents
    """)
    return


@app.cell
def _(Chem, CustomEvaluator):
    """Define the custom scoring function and instantiate evaluator."""

    def qed_score(mol):
        """Return the QED (drug-likeness) score of a molecule."""
        return Chem.QED.qed(mol)

    evaluator = CustomEvaluator(scoring_function=qed_score)
    return evaluator, qed_score


@app.cell
def _(Chem, evaluator, mo, pl):
    """Evaluate a small set of 5 SMILES with CustomEvaluator."""
    smiles_list = [
        "CCO",
        "c1ccccc1",
        "CC(=O)O",
        "CCN(CC)CC",
        "c1ccncc1O",
    ]

    results = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        score = evaluator.evaluate(mol)
        results.append({"SMILES": smi, "HeavyAtoms": score})

    df_eval = pl.DataFrame(results)
    mo.md("**Small SMILES evaluation complete**")
    mo.md(f"**DataFrame shape:** {df_eval.shape}")
    df_eval
    return


@app.cell
def _(
    AMIDE_COUPLING_SMARTS,
    CustomEvaluatorConfig,
    GreedyConfig,
    REAGENT_FILES,
    ReactionConfig,
    ReactionDef,
    StandardWarmupConfig,
    SynthesisPipeline,
    ThompsonSampler,
    ThompsonSamplingConfig,
    mo,
    pl,
    qed_score,
):
    """Minimal Thrombin TS run using CustomEvaluatorConfig."""
    mo.md("**Running minimal Thrombin TS using CustomEvaluatorConfig...**")

    # Setup synthesis pipeline
    reaction_config = ReactionConfig(
        reactions=[ReactionDef(reaction_smarts=AMIDE_COUPLING_SMARTS, step_index=0)],
        reagent_file_list=REAGENT_FILES,
    )
    pipeline = SynthesisPipeline(reaction_config)

    # Thompson Sampling config (minimal iterations for demo)
    ts_config = ThompsonSamplingConfig(
        synthesis_pipeline=pipeline,
        num_ts_iterations=3,
        num_warmup_trials=3,
        strategy_config=GreedyConfig(mode="maximize"),
        warmup_config=StandardWarmupConfig(),
        evaluator_config=CustomEvaluatorConfig(
            scoring_function=qed_score),
        batch_size=1,
        max_resamples=10,
        hide_progress=True,
    )

    # Run TS
    sampler = ThompsonSampler.from_config(ts_config)
    warmup_df = sampler.warm_up(num_warmup_trials=ts_config.num_warmup_trials)
    search_df = sampler.search(num_cycles=ts_config.num_ts_iterations)
    df_ts = pl.concat([warmup_df, search_df])
    sampler.close()

    mo.md("**Minimal Thrombin TS run complete.**")
    mo.md(f"**DataFrame shape:** {df_ts.shape}")
    return


if __name__ == "__main__":
    app.run()
