"""
TACTICS Library Component Comparison App

An interactive benchmark tool comparing Thompson Sampling behavior between
2-component and 3-component representations of the Thrombin Linear Amide Library.

This notebook explores whether modeling the library as a 3-component synthesis
(amino acid + amino acid + carboxylic acid) produces different optimization
behavior than the 2-component approach (dipeptide + carboxylic acid).

Features:
1. Side-by-side comparison of 2-component vs 3-component library representation
2. Selection strategy configuration
3. TS_Benchmarks visualizations for both library types
4. Hit recovery analysis comparison

Dataset: Thrombin Linear Amide Library (~500K products)

Run as app: marimo run notebooks/library_component_comparison.py
Edit mode:  marimo edit notebooks/library_component_comparison.py
"""

import marimo

__generated_with = "0.18.2"
app = marimo.App(width="full", app_title="TACTICS Library Component Comparison")


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

    # Standard imports
    import polars as pl
    import time

    # TACTICS Thompson Sampling imports
    from TACTICS.thompson_sampling import ThompsonSampler
    from TACTICS.thompson_sampling.config import ThompsonSamplingConfig
    from TACTICS.thompson_sampling.strategies.config import (
        GreedyConfig,
        RouletteWheelConfig,
        UCBConfig,
        EpsilonGreedyConfig,
        BoltzmannConfig,
    )
    from TACTICS.thompson_sampling.warmup.config import (
        BalancedWarmupConfig,
        StandardWarmupConfig,
        EnhancedWarmupConfig,
    )
    from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

    # TACTICS Library Enumeration imports
    from TACTICS.library_enumeration import SynthesisPipeline
    from TACTICS.library_enumeration.smarts_toolkit import (
        ReactionConfig,
        ReactionDef,
        StepInput,
        InputSource,
    )

    # TACTICS Library Analysis imports
    from TACTICS.library_analysis.visualization import TS_Benchmarks
    return (
        BalancedWarmupConfig,
        BoltzmannConfig,
        EnhancedWarmupConfig,
        EpsilonGreedyConfig,
        GreedyConfig,
        InputSource,
        LookupEvaluatorConfig,
        ReactionConfig,
        ReactionDef,
        RouletteWheelConfig,
        StandardWarmupConfig,
        StepInput,
        SynthesisPipeline,
        TS_Benchmarks,
        ThompsonSampler,
        ThompsonSamplingConfig,
        UCBConfig,
        mo,
        pl,
        time,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # TACTICS Library Component Comparison

    This notebook compares Thompson Sampling optimization behavior when the
    Thrombin Linear Amide Library is modeled as:

    - **2-Component Library**: Dipeptide + Carboxylic Acid (pre-coupled dipeptides)
    - **3-Component Library**: Amino Acid + Amino Acid + Carboxylic Acid (synthesized in-situ)

    Both representations describe the same ~500K product library, but the component
    structure affects how the Thompson Sampling algorithm explores the chemical space.

    **Key Question**: Does the 3-component representation provide better exploration
    or different optimization dynamics compared to the 2-component approach?
    """)
    return


@app.cell
def _():
    """Load bundled Thrombin dataset paths."""
    import importlib.resources

    _data_files = importlib.resources.files("TACTICS.data.thrombin")
    ACIDS_FILE = str(_data_files / "acids.smi")
    DIPEPTIDES_FILE = str(_data_files / "coupled_aa_sub.smi")
    AMINO_ACIDS_FILE = str(_data_files / "amino_acids_no_fmoc.smi")
    SCORES_FILE = str(_data_files / "product_scores.csv")

    # 2-Component Library Files
    REAGENT_FILES_2COMP = [ACIDS_FILE, DIPEPTIDES_FILE]

    # 3-Component Library Files
    REAGENT_FILES_3COMP = [AMINO_ACIDS_FILE, AMINO_ACIDS_FILE, ACIDS_FILE]

    # SMARTS patterns
    AMIDE_COUPLING_SMARTS = "[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]>>[#6:1](=[O:2])[#7:3]"
    PEPTIDE_COUPLING_SMARTS = "[#6X3:1](=[#8X1])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):2]>>[#6X3:1](=[#8X1])[#7X3:2]"
    AMINE_SUBSTITUTION_SMARTS = "[#6:1](=O)[OH]>>[#6:1](=O)[NH2]"
    return (
        ACIDS_FILE,
        AMIDE_COUPLING_SMARTS,
        AMINE_SUBSTITUTION_SMARTS,
        AMINO_ACIDS_FILE,
        DIPEPTIDES_FILE,
        PEPTIDE_COUPLING_SMARTS,
        REAGENT_FILES_2COMP,
        REAGENT_FILES_3COMP,
        SCORES_FILE,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Configuration

    ### Selection Strategy
    """)
    return


@app.cell
def _(mo):
    """Strategy selection with parameter customization."""
    strategy_dropdown = mo.ui.dropdown(
        options=["Epsilon-Greedy", "Greedy", "UCB", "Roulette Wheel", "Boltzmann"],
        value="Epsilon-Greedy",
        label="Selection Strategy"
    )

    # Strategy parameters
    epsilon_value = mo.ui.slider(start=0.05, stop=0.5, value=0.2, step=0.05, label="Epsilon")
    epsilon_decay = mo.ui.slider(start=0.9, stop=1.0, value=0.995, step=0.005, label="Decay")
    ucb_c = mo.ui.slider(start=0.5, stop=4.0, value=2.0, step=0.5, label="UCB c")
    rw_alpha = mo.ui.slider(start=0.05, stop=0.5, value=0.1, step=0.05, label="RW alpha")
    boltz_temp = mo.ui.slider(start=0.1, stop=2.0, value=1.0, step=0.1, label="Temperature")

    # Warmup type dropdowns for each library
    warmup_type_2comp = mo.ui.dropdown(
        options=["Balanced", "Standard", "Enhanced"],
        value="Balanced",
        label="2-Component Warmup Type"
    )
    warmup_type_3comp = mo.ui.dropdown(
        options=["Balanced", "Standard", "Enhanced"],
        value="Balanced",
        label="3-Component Warmup Type"
    )

    # K sliders for balanced warmup (shown conditionally in display cell)
    warmup_k_2comp = mo.ui.slider(
        start=1, stop=20, value=3, step=1,
        label="K (obs/reagent)"
    )
    warmup_k_3comp = mo.ui.slider(
        start=1, stop=20, value=5, step=1,
        label="K (obs/reagent)"
    )
    return (
        boltz_temp,
        epsilon_decay,
        epsilon_value,
        rw_alpha,
        strategy_dropdown,
        ucb_c,
        warmup_k_2comp,
        warmup_k_3comp,
        warmup_type_2comp,
        warmup_type_3comp,
    )


@app.cell
def _(
    boltz_temp,
    epsilon_decay,
    epsilon_value,
    mo,
    rw_alpha,
    strategy_dropdown,
    ucb_c,
    warmup_k_2comp,
    warmup_k_3comp,
    warmup_type_2comp,
    warmup_type_3comp,
):
    """Display strategy controls."""
    _strategy_params = mo.hstack([
        mo.vstack([mo.md("*Epsilon-Greedy:*"), epsilon_value, epsilon_decay]),
        mo.vstack([mo.md("*UCB:*"), ucb_c]),
        mo.vstack([mo.md("*Roulette Wheel:*"), rw_alpha]),
        mo.vstack([mo.md("*Boltzmann:*"), boltz_temp]),
    ], justify="start", gap=4)

    # Conditionally show K slider only when Balanced is selected
    _2comp_warmup_ui = (
        mo.vstack([warmup_type_2comp, warmup_k_2comp])
        if warmup_type_2comp.value == "Balanced"
        else mo.vstack([warmup_type_2comp])
    )
    _3comp_warmup_ui = (
        mo.vstack([warmup_type_3comp, warmup_k_3comp])
        if warmup_type_3comp.value == "Balanced"
        else mo.vstack([warmup_type_3comp])
    )

    _warmup_controls = mo.hstack([
        mo.vstack([mo.md("**2-Component Warmup:**"), _2comp_warmup_ui]),
        mo.vstack([mo.md("**3-Component Warmup:**"), _3comp_warmup_ui]),
    ], justify="start", gap=4)

    _controls = mo.vstack([
        strategy_dropdown,
        _warmup_controls,
        mo.md("**Strategy Parameters:**"),
        _strategy_params,
    ])
    _controls
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Search Parameters
    """)
    return


@app.cell
def _(mo):
    """Search parameters."""
    iterations_slider = mo.ui.slider(
        start=100, stop=2000, value=500, step=100,
        label="Iterations per cycle"
    )
    cycles_slider = mo.ui.slider(
        start=1, stop=10, value=3, step=1,
        label="Benchmark cycles"
    )
    top_n_slider = mo.ui.slider(
        start=50, stop=500, value=100, step=50,
        label="Top N for recovery"
    )
    return cycles_slider, iterations_slider, top_n_slider


@app.cell(hide_code=True)
def _(cycles_slider, iterations_slider, mo, top_n_slider):
    """Display search parameter controls."""
    mo.hstack([iterations_slider, cycles_slider, top_n_slider], justify="start", gap=2)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Run Comparison Benchmark
    """)
    return


@app.cell
def _(mo):
    """Run button."""
    run_button = mo.ui.run_button(label="Run Library Comparison Benchmark")
    run_button
    return (run_button,)


@app.cell
def _(
    ACIDS_FILE,
    AMIDE_COUPLING_SMARTS,
    AMINE_SUBSTITUTION_SMARTS,
    AMINO_ACIDS_FILE,
    DIPEPTIDES_FILE,
    InputSource,
    PEPTIDE_COUPLING_SMARTS,
    ReactionConfig,
    ReactionDef,
    StepInput,
    SynthesisPipeline,
):
    """Create synthesis pipelines for both library representations."""
    # 2-Component Pipeline: Dipeptide + Carboxylic Acid
    reaction_config_2comp = ReactionConfig(
        reactions=[ReactionDef(reaction_smarts=AMIDE_COUPLING_SMARTS, step_index=0)],
        reagent_file_list=[ACIDS_FILE, DIPEPTIDES_FILE],
    )
    pipeline_2comp = SynthesisPipeline(reaction_config_2comp)

    # 3-Component Pipeline: Carboxylic Acid + Amino Acid + Amino Acid
    # Note: Reagent order is [Acid, AA1, AA2] to match scores file naming (CA_AA_AA)
    reaction_config_3comp = ReactionConfig(
        reactions=[
            # Step 0: amino acid + amino acid → dipeptide
            ReactionDef(
                reaction_smarts=PEPTIDE_COUPLING_SMARTS,
                step_index=0,
                description="peptide coupling",
            ),
            # Step 1: hydroxyl → amine substitution
            ReactionDef(
                reaction_smarts=AMINE_SUBSTITUTION_SMARTS,
                step_index=1,
                description="Amine substitution",
            ),
            # Step 2: dipeptide + carboxylic acid → final product
            ReactionDef(
                reaction_smarts=AMIDE_COUPLING_SMARTS,
                step_index=2,
                description="linear amide coupling",
            ),
        ],
        reagent_file_list=[
            ACIDS_FILE,        # Index 0: Carboxylic acids (for naming: CA_...)
            AMINO_ACIDS_FILE,  # Index 1: First amino acid
            AMINO_ACIDS_FILE,  # Index 2: Second amino acid
        ],
        step_inputs={
            0: [
                StepInput(source=InputSource.REAGENT_FILE, file_index=1),  # First AA
                StepInput(source=InputSource.REAGENT_FILE, file_index=2),  # Second AA
            ],
            1: [
                StepInput(source=InputSource.PREVIOUS_STEP, step_index=0)
            ],
            2: [
                StepInput(source=InputSource.REAGENT_FILE, file_index=0),  # Acid
                StepInput(source=InputSource.PREVIOUS_STEP, step_index=1)
            ],
        },
    )
    pipeline_3comp = SynthesisPipeline(reaction_config_3comp)

    print(f"2-Component Pipeline: {pipeline_2comp.num_components} components, {pipeline_2comp.num_steps} steps")
    print(f"3-Component Pipeline: {pipeline_3comp.num_components} components, {pipeline_3comp.num_steps} steps")
    return pipeline_2comp, pipeline_3comp, reaction_config_2comp, reaction_config_3comp


@app.cell
def _(
    BalancedWarmupConfig,
    BoltzmannConfig,
    EnhancedWarmupConfig,
    EpsilonGreedyConfig,
    GreedyConfig,
    LookupEvaluatorConfig,
    RouletteWheelConfig,
    SCORES_FILE,
    StandardWarmupConfig,
    ThompsonSampler,
    ThompsonSamplingConfig,
    UCBConfig,
    boltz_temp,
    cycles_slider,
    epsilon_decay,
    epsilon_value,
    iterations_slider,
    mo,
    pipeline_2comp,
    pipeline_3comp,
    pl,
    run_button,
    rw_alpha,
    strategy_dropdown,
    time,
    ucb_c,
    warmup_k_2comp,
    warmup_k_3comp,
    warmup_type_2comp,
    warmup_type_3comp,
):
    """Execute Thompson Sampling benchmark for both library representations."""
    mo.stop(not run_button.value, mo.md("*Click 'Run Library Comparison Benchmark' to start*"))

    # Build strategy config based on selection
    strategy_name = strategy_dropdown.value
    if strategy_name == "Epsilon-Greedy":
        strategy_config = EpsilonGreedyConfig(
            mode="minimize", epsilon=epsilon_value.value, decay=epsilon_decay.value
        )
    elif strategy_name == "Greedy":
        strategy_config = GreedyConfig(mode="minimize")
    elif strategy_name == "UCB":
        strategy_config = UCBConfig(mode="minimize", c=ucb_c.value)
    elif strategy_name == "Roulette Wheel":
        strategy_config = RouletteWheelConfig(
            mode="minimize", alpha=rw_alpha.value, beta=rw_alpha.value
        )
    else:  # Boltzmann
        strategy_config = BoltzmannConfig(
            mode="minimize", temperature=boltz_temp.value
        )

    # Build warmup config for 2-component library
    _warmup_type_2comp = warmup_type_2comp.value
    _k_2comp = warmup_k_2comp.value
    if _warmup_type_2comp == "Balanced":
        warmup_config_2comp = BalancedWarmupConfig(observations_per_reagent=_k_2comp)
        _warmup_label_2comp = f"Balanced-K{_k_2comp}"
    elif _warmup_type_2comp == "Standard":
        warmup_config_2comp = StandardWarmupConfig()
        _warmup_label_2comp = "Standard"
    else:  # Enhanced
        warmup_config_2comp = EnhancedWarmupConfig()
        _warmup_label_2comp = "Enhanced"

    # Build warmup config for 3-component library
    _warmup_type_3comp = warmup_type_3comp.value
    _k_3comp = warmup_k_3comp.value
    if _warmup_type_3comp == "Balanced":
        warmup_config_3comp = BalancedWarmupConfig(observations_per_reagent=_k_3comp)
        _warmup_label_3comp = f"Balanced-K{_k_3comp}"
    elif _warmup_type_3comp == "Standard":
        warmup_config_3comp = StandardWarmupConfig()
        _warmup_label_3comp = "Standard"
    else:  # Enhanced
        warmup_config_3comp = EnhancedWarmupConfig()
        _warmup_label_3comp = "Enhanced"

    # Evaluator config
    evaluator_config = LookupEvaluatorConfig(
        ref_filename=SCORES_FILE,
        compound_col="Product_Code",
        score_col="Scores",
    )

    # Run parameters
    _num_cycles = cycles_slider.value
    _num_iterations = iterations_slider.value

    # Results storage
    results_2comp = {}
    results_3comp = {}
    start_time = time.time()

    # --- Run 2-Component Library Benchmark ---
    _cycle_results_2comp = []
    for _cycle in range(_num_cycles):
        _config = ThompsonSamplingConfig(
            synthesis_pipeline=pipeline_2comp,
            num_ts_iterations=_num_iterations,
            num_warmup_trials=3,
            strategy_config=strategy_config,
            warmup_config=warmup_config_2comp,
            evaluator_config=evaluator_config,
            batch_size=1,
            max_resamples=1000,
            hide_progress=True,
        )

        _sampler = ThompsonSampler.from_config(_config)
        _warmup_df = _sampler.warm_up(num_warmup_trials=_config.num_warmup_trials)
        _search_df = _sampler.search(num_cycles=_config.num_ts_iterations)

        _combined = pl.concat([_warmup_df, _search_df])
        _result = _combined.select(["Name", "score"])
        _cycle_results_2comp.append(_result)
        _sampler.close()

    results_2comp[f"2-Comp ({_warmup_label_2comp})"] = _cycle_results_2comp

    # --- Run 3-Component Library Benchmark ---
    _cycle_results_3comp = []
    for _cycle in range(_num_cycles):
        _config = ThompsonSamplingConfig(
            synthesis_pipeline=pipeline_3comp,
            num_ts_iterations=_num_iterations,
            num_warmup_trials=3,
            strategy_config=strategy_config,
            warmup_config=warmup_config_3comp,
            evaluator_config=evaluator_config,
            batch_size=1,
            max_resamples=1000,
            hide_progress=True,
        )

        _sampler = ThompsonSampler.from_config(_config)
        _warmup_df = _sampler.warm_up(num_warmup_trials=_config.num_warmup_trials)
        _search_df = _sampler.search(num_cycles=_config.num_ts_iterations)

        _combined = pl.concat([_warmup_df, _search_df])
        _result = _combined.select(["Name", "score"])
        _cycle_results_3comp.append(_result)
        _sampler.close()

    results_3comp[f"3-Comp ({_warmup_label_3comp})"] = _cycle_results_3comp

    total_time = time.time() - start_time

    # Compute summary statistics
    _best_2comp = [r["score"].min() for r in _cycle_results_2comp]
    _best_3comp = [r["score"].min() for r in _cycle_results_3comp]

    run_metadata = {
        "strategy": strategy_name,
        "warmup_2comp": _warmup_label_2comp,
        "warmup_3comp": _warmup_label_3comp,
        "cycles": _num_cycles,
        "iterations": _num_iterations,
        "2comp_best": min(_best_2comp),
        "2comp_avg_best": sum(_best_2comp) / len(_best_2comp),
        "3comp_best": min(_best_3comp),
        "3comp_avg_best": sum(_best_3comp) / len(_best_3comp),
    }
    return (
        evaluator_config,
        results_2comp,
        results_3comp,
        run_metadata,
        strategy_config,
        strategy_name,
        total_time,
    )


@app.cell
def _(mo, run_button, run_metadata, total_time):
    """Display results summary."""
    mo.stop(not run_button.value)

    mo.md(f"""
    ---
    ## Results Summary

    **Strategy:** {run_metadata['strategy']}

    **Total runtime:** {total_time:.1f}s ({run_metadata['cycles']} cycles × {run_metadata['iterations']} iterations)

    | Library Type | Warmup | Best Score | Avg Best Score |
    |--------------|--------|------------|----------------|
    | 2-Component  | {run_metadata['warmup_2comp']} | {run_metadata['2comp_best']:.3f} | {run_metadata['2comp_avg_best']:.3f} |
    | 3-Component  | {run_metadata['warmup_3comp']} | {run_metadata['3comp_best']:.3f} | {run_metadata['3comp_avg_best']:.3f} |

    **Difference:** {abs(run_metadata['2comp_avg_best'] - run_metadata['3comp_avg_best']):.3f}
    ({"3-Component better" if run_metadata['3comp_avg_best'] < run_metadata['2comp_avg_best'] else "2-Component better"})
    """)
    return


@app.cell
def _(SCORES_FILE, mo, pl, run_button):
    """Load reference data for visualizations."""
    mo.stop(not run_button.value)

    reference_data = (
        pl.read_csv(SCORES_FILE)
        .rename({"Product_Code": "Name", "Scores": "score"})
        .select(["Name", "score"])
    )
    return (reference_data,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Comparison Visualizations

    ### Strip Plot: Score Distributions
    """)
    return


@app.cell
def _(
    TS_Benchmarks,
    cycles_slider,
    mo,
    reference_data,
    results_2comp,
    results_3comp,
    run_button,
    top_n_slider,
):
    """Create combined TS_Benchmarks for side-by-side comparison."""
    mo.stop(not run_button.value)

    # Merge results from both libraries
    combined_results = {**results_2comp, **results_3comp}
    combined_methods = list(combined_results.keys())

    benchmarks_combined = TS_Benchmarks(
        no_of_cycles=cycles_slider.value,
        methods_list=combined_methods,
        TS_runs_data=combined_results,
        reference_data=reference_data,
        top_n=top_n_slider.value,
        sort_type="minimize",
        top_ns=[25, 50, 100, 200, 300],
    )
    return benchmarks_combined, combined_methods, combined_results


@app.cell
def _(benchmarks_combined, mo, run_button):
    """Generate combined strip plot."""
    mo.stop(not run_button.value)

    benchmarks_combined.stripplot_TS_results(
        width=900,
        height=450,
        show_plot=True,
        legend_position="bottom"
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Bar Plot: Hit Recovery by Cycle
    """)
    return


@app.cell
def _(benchmarks_combined, mo, run_button):
    """Generate combined bar plot."""
    mo.stop(not run_button.value)

    benchmarks_combined.plot_barplot_TS_results(
        width=900,
        height=450,
        show_plot=True,
        legend_position="bottom",
        dark_mode=False
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Line Plot: Performance vs Top-N
    """)
    return


@app.cell
def _(benchmarks_combined, mo, run_button):
    """Generate combined line plot."""
    mo.stop(not run_button.value)

    benchmarks_combined.plot_line_performance_with_error_bars(
        width=900,
        height=450,
        show_plot=True,
        legend_position="bottom"
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Performance Statistics
    """)
    return


@app.cell(hide_code=True)
def _(benchmarks_combined, mo, run_button):
    """Display grouped statistics."""
    mo.stop(not run_button.value)

    _output = (
        mo.vstack([
            mo.md("**Mean Performance by Top-N:**"),
            benchmarks_combined.grouped_stats,
        ])
        if benchmarks_combined.grouped_stats is not None
        else mo.md("*Statistics not available*")
    )
    _output
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Top Compounds Found
    """)
    return


@app.cell
def _(combined_results, mo, pl, run_button):
    """Display top compounds from both library representations."""
    mo.stop(not run_button.value)

    _all_dfs = []
    for _name, _results in combined_results.items():
        for _i, _df in enumerate(_results):
            _labeled = _df.with_columns([
                pl.lit(_name).alias("library_type"),
                pl.lit(_i + 1).alias("cycle"),
            ])
            _all_dfs.append(_labeled)

    combined_df = pl.concat(_all_dfs)
    top_compounds = (
        combined_df
        .sort("score")
        .unique(subset=["Name"], keep="first")
        .head(20)
        .select(["Name", "score", "library_type", "cycle"])
    )
    top_compounds
    return combined_df, top_compounds


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Additivity Analysis
    """)
    return


@app.cell
def _(SCORES_FILE, mo, pl, run_button):
    """Compute additivity metrics to explain performance differences."""
    mo.stop(not run_button.value)

    # Load and parse scores
    _scores = pl.read_csv(SCORES_FILE)
    _scores = _scores.with_columns([
        pl.col("Product_Code").str.split("_").list.get(0).alias("Acid"),
        pl.col("Product_Code").str.split("_").list.get(1).alias("AA1"),
        pl.col("Product_Code").str.split("_").list.get(2).alias("AA2"),
        (pl.col("Product_Code").str.split("_").list.get(1) + "_" +
         pl.col("Product_Code").str.split("_").list.get(2)).alias("Dipeptide"),
    ])

    # Calculate variance components
    _total_var = _scores["Scores"].var()
    _acid_var = _scores.group_by("Acid").agg(pl.col("Scores").mean())["Scores"].var()
    _aa1_var = _scores.group_by("AA1").agg(pl.col("Scores").mean())["Scores"].var()
    _aa2_var = _scores.group_by("AA2").agg(pl.col("Scores").mean())["Scores"].var()
    _dip_var = _scores.group_by("Dipeptide").agg(pl.col("Scores").mean())["Scores"].var()

    _additive_var = _aa1_var + _aa2_var
    _interaction_ratio = _dip_var / _additive_var if _additive_var > 0 else float('inf')

    # Determine recommendation
    if _interaction_ratio > 1.3:
        _recommendation = "**Strong interactions detected** → 2-Component representation recommended"
        _color = "red"
    elif _interaction_ratio > 1.0:
        _recommendation = "**Moderate interactions** → 2-Component likely better"
        _color = "orange"
    else:
        _recommendation = "**Approximately additive** → 3-Component may be viable"
        _color = "green"

    mo.md(f"""
    ### Variance Component Analysis

    | Component | Variance | % of Total |
    |-----------|----------|------------|
    | Total Score | {_total_var:.4f} | 100% |
    | Acid (mean per acid) | {_acid_var:.4f} | {100*_acid_var/_total_var:.1f}% |
    | AA1 (mean per AA1) | {_aa1_var:.4f} | {100*_aa1_var/_total_var:.1f}% |
    | AA2 (mean per AA2) | {_aa2_var:.4f} | {100*_aa2_var/_total_var:.1f}% |
    | Dipeptide (mean per dipeptide) | {_dip_var:.4f} | {100*_dip_var/_total_var:.1f}% |

    ### Additivity Test

    - **Sum of individual AA variances**: {_additive_var:.4f}
    - **Dipeptide variance**: {_dip_var:.4f}
    - **Interaction ratio** (dipeptide / sum): **{_interaction_ratio:.2f}**

    {_recommendation}

    *Ratio > 1.0 indicates interaction effects beyond what individual amino acids explain.*
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Understanding the Comparison

    ### Library Representations

    | Aspect | 2-Component | 3-Component |
    |--------|-------------|-------------|
    | **Reagent Files** | Dipeptides (3844) + Acids (130) | Amino Acids (62) × 2 + Acids (130) |
    | **Posteriors Tracked** | 130 + 3844 = **3974** | 130 + 62 + 62 = **254** |
    | **Synthesis Steps** | 1 (amide coupling) | 3 (peptide → amine sub → amide) |
    | **Implicit Assumption** | Dipeptides are atomic units | AA contributions are additive |

    ### Why 3-Component May Underperform

    Thompson Sampling with independent posteriors assumes **additivity**:
    ```
    score ≈ f(acid) + g(AA1) + h(AA2)
    ```

    But if there are strong **interaction effects** between amino acids:
    ```
    score ≈ f(acid) + interaction(AA1, AA2)
    ```

    The 3-component model fails because:

    1. **Information Loss**: When `CA5_AA32_AA61` scores well, the model credits AA32 and AA61
       independently. But AA32 may only be good *with* AA61, not generally.

    2. **Noisy Posteriors**: Each amino acid's posterior is updated based on combinations with
       many different partners, creating conflicting signals.

    3. **Fewer Parameters**: 254 posteriors must capture 499,720 product preferences, vs
       3974 posteriors in the 2-component case.

    ### Additivity Test for Your Data

    To determine which representation is appropriate, check:
    - **Dipeptide variance** vs **AA1 variance + AA2 variance**
    - If dipeptide variance >> sum of individual variances → strong interactions → use 2-component
    - If roughly equal → additive → 3-component may work

    For Thrombin data: Dipeptide variance (0.89) >> AA1 (0.25) + AA2 (0.31) = 0.56
    → **Strong interactions exist** → 2-component is more appropriate

    ### Usage Recommendations

    - **2-Component**: Use when component interactions matter (most drug discovery cases)
    - **3-Component**: Use when contributions are approximately additive, or for initial
      exploration to identify promising component classes before detailed optimization
    """)
    return


if __name__ == "__main__":
    app.run()
