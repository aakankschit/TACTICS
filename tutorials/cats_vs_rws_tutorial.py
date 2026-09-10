"""
TACTICS: CATS vs RWS Tutorial

An interactive tutorial comparing Component-Aware Thompson Sampling (CATS) with
standard Roulette Wheel Selection (RWS) on the Thrombin dataset.

CATS adds criticality modulation on top of thermal cycling — this notebook
demonstrates the difference by running the same alpha/beta configuration with
CATS enabled vs disabled (via min_observations=999999).

The tutorial covers both 2-component and 3-component library representations,
with analysis sections for each and a combined comparison at the end.

Features:
1. CATS vs RWS (thermal cycling only) comparison
2. 2-component and 3-component library benchmarks
3. TS_Benchmarks visualizations (bar, strip, line plots)
4. Custom Altair plots (cumulative recovery, exploration efficiency)
5. Combined faceted comparison across library formats

Dataset: Thrombin Linear Amide Library (~500K products)

Run as app: marimo run tutorials/cats_vs_rws_tutorial.py
Edit mode:  marimo edit tutorials/cats_vs_rws_tutorial.py
"""

import marimo

__generated_with = "0.18.2"
app = marimo.App(width="full", app_title="TACTICS: CATS vs RWS Tutorial")


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
    import altair as alt
    import time

    # TACTICS Thompson Sampling imports
    from TACTICS.thompson_sampling import ThompsonSampler
    from TACTICS.thompson_sampling.config import ThompsonSamplingConfig
    from TACTICS.thompson_sampling.strategies.config import RouletteWheelConfig
    from TACTICS.thompson_sampling.warmup.config import BalancedWarmupConfig
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
        InputSource,
        LookupEvaluatorConfig,
        ReactionConfig,
        ReactionDef,
        RouletteWheelConfig,
        StepInput,
        SynthesisPipeline,
        TS_Benchmarks,
        ThompsonSampler,
        ThompsonSamplingConfig,
        alt,
        mo,
        pl,
        time,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # CATS vs RWS: Component-Aware Thompson Sampling Tutorial

    This notebook compares two configurations of the Roulette Wheel selection strategy
    on the Thrombin Linear Amide Library:

    - **CATS** (Component-Aware Thompson Sampling): Thermal cycling **+** criticality modulation.
      CATS adjusts exploration per-component based on Shannon entropy of posterior means —
      critical components (high score variance) get exploited, flexible components get explored.

    - **RWS** (Roulette Wheel Selection): Thermal cycling **only**.
      Same alpha/beta temperature differential, but criticality modulation is disabled
      by setting `min_observations=999999` so the CATS weight never ramps up.

    By toggling only the criticality modulation, we isolate the benefit of CATS
    while keeping the base thermal cycling identical.

    ### What is Thermal Cycling?

    Thermal cycling alternates a "heated" component (higher temperature → more exploration)
    with "cooled" components (lower temperature → more exploitation). The `alpha` parameter
    sets the heated temperature, and `beta` sets the cooled temperature.

    ### What does CATS add?

    CATS measures each component's **criticality** — how concentrated the posterior means are.
    High criticality means a few reagents dominate; CATS reduces temperature (exploit).
    Low criticality means reagents are balanced; CATS increases temperature (explore).
    This per-component adaptation helps focus search effort where it matters most.
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

    Adjust the search parameters and thermal cycling temperatures below. The same
    alpha/beta values are used for both CATS and RWS to ensure a fair comparison.
    """)
    return


@app.cell
def _(mo):
    """Configuration UI elements."""
    num_cycles_slider = mo.ui.slider(
        start=1, stop=10, value=5, step=1,
        label="Benchmark cycles"
    )
    iterations_slider = mo.ui.slider(
        start=100, stop=2000, value=500, step=100,
        label="Iterations per cycle"
    )
    top_n_slider = mo.ui.slider(
        start=50, stop=500, value=100, step=50,
        label="Top N for recovery"
    )
    warmup_k_slider = mo.ui.slider(
        start=1, stop=20, value=3, step=1,
        label="Warmup K (obs/reagent)"
    )
    alpha_slider = mo.ui.slider(
        start=0.05, stop=0.5, value=0.1, step=0.05,
        label="Alpha (heated temp)"
    )
    beta_slider = mo.ui.slider(
        start=0.01, stop=0.3, value=0.05, step=0.01,
        label="Beta (cooled temp)"
    )
    return (
        alpha_slider,
        beta_slider,
        iterations_slider,
        num_cycles_slider,
        top_n_slider,
        warmup_k_slider,
    )


@app.cell(hide_code=True)
def _(
    alpha_slider,
    beta_slider,
    iterations_slider,
    mo,
    num_cycles_slider,
    top_n_slider,
    warmup_k_slider,
):
    """Display configuration controls."""
    mo.vstack([
        mo.md("### Search Parameters"),
        mo.hstack([num_cycles_slider, iterations_slider, top_n_slider], justify="start", gap=2),
        mo.md("### Warmup & Temperature"),
        mo.hstack([warmup_k_slider, alpha_slider, beta_slider], justify="start", gap=2),
    ])
    return


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
    reaction_config_3comp = ReactionConfig(
        reactions=[
            ReactionDef(
                reaction_smarts=PEPTIDE_COUPLING_SMARTS,
                step_index=0,
                description="peptide coupling",
            ),
            ReactionDef(
                reaction_smarts=AMINE_SUBSTITUTION_SMARTS,
                step_index=1,
                description="Amine substitution",
            ),
            ReactionDef(
                reaction_smarts=AMIDE_COUPLING_SMARTS,
                step_index=2,
                description="linear amide coupling",
            ),
        ],
        reagent_file_list=[
            ACIDS_FILE,
            AMINO_ACIDS_FILE,
            AMINO_ACIDS_FILE,
        ],
        step_inputs={
            0: [
                StepInput(source=InputSource.REAGENT_FILE, file_index=1),
                StepInput(source=InputSource.REAGENT_FILE, file_index=2),
            ],
            1: [
                StepInput(source=InputSource.PREVIOUS_STEP, step_index=0)
            ],
            2: [
                StepInput(source=InputSource.REAGENT_FILE, file_index=0),
                StepInput(source=InputSource.PREVIOUS_STEP, step_index=1)
            ],
        },
    )
    pipeline_3comp = SynthesisPipeline(reaction_config_3comp)

    print(f"2-Component Pipeline: {pipeline_2comp.num_components} components, {pipeline_2comp.num_steps} steps")
    print(f"3-Component Pipeline: {pipeline_3comp.num_components} components, {pipeline_3comp.num_steps} steps")
    return pipeline_2comp, pipeline_3comp


@app.cell
def _(
    BalancedWarmupConfig,
    LookupEvaluatorConfig,
    RouletteWheelConfig,
    SCORES_FILE,
    ThompsonSampler,
    ThompsonSamplingConfig,
    pl,
    time,
):
    """Define the CATS vs RWS run function."""
    def run_cats_vs_rws(pipeline, alpha, beta, warmup_k, num_cycles, num_iterations):
        """
        Run CATS and RWS (thermal cycling only) on the same pipeline.

        Returns:
            dict with keys "CATS" and "RWS", each mapping to a list of cycle DataFrames
            dict with metadata (timing, best scores)
        """
        evaluator_config = LookupEvaluatorConfig(
            ref_filename=SCORES_FILE,
            compound_col="Product_Code",
            score_col="Scores",
        )
        warmup_config = BalancedWarmupConfig(observations_per_reagent=warmup_k)

        # CATS config: default min_observations (criticality active)
        cats_strategy = RouletteWheelConfig(
            mode="minimize", alpha=alpha, beta=beta,
        )

        # RWS config: min_observations=999999 disables CATS weight ramp-up
        # weight = min_obs / (2 * 999999) ≈ 0, so criticality has no effect
        rws_strategy = RouletteWheelConfig(
            mode="minimize", alpha=alpha, beta=beta, min_observations=999999,
        )

        results = {"CATS": [], "RWS": []}
        metadata = {"CATS": [], "RWS": []}

        for method_name, strategy_config in [("CATS", cats_strategy), ("RWS", rws_strategy)]:
            _start = time.time()
            for _cycle in range(num_cycles):
                _config = ThompsonSamplingConfig(
                    synthesis_pipeline=pipeline,
                    num_ts_iterations=num_iterations,
                    num_warmup_trials=3,
                    strategy_config=strategy_config,
                    warmup_config=warmup_config,
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
                results[method_name].append(_result)
                _sampler.close()

            _elapsed = time.time() - _start
            _best_scores = [r["score"].min() for r in results[method_name]]
            metadata[method_name] = {
                "time": _elapsed,
                "best": min(_best_scores),
                "avg_best": sum(_best_scores) / len(_best_scores),
            }

        return results, metadata
    return (run_cats_vs_rws,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## 2-Component Library: CATS vs RWS

    The 2-component library models the Thrombin dataset as:
    - **130 carboxylic acids** + **3,844 dipeptides** = ~499,720 products
    - Single amide coupling reaction
    - Thompson Sampling tracks **3,974 posteriors** (130 + 3,844)

    This is the standard representation where dipeptides are treated as atomic units.
    """)
    return


@app.cell
def _(mo):
    """Run button for 2-component benchmark."""
    run_2comp_button = mo.ui.run_button(label="Run 2-Component Benchmark")
    run_2comp_button
    return (run_2comp_button,)


@app.cell
def _(
    alpha_slider,
    beta_slider,
    iterations_slider,
    mo,
    num_cycles_slider,
    pipeline_2comp,
    run_2comp_button,
    run_cats_vs_rws,
    warmup_k_slider,
):
    """Execute 2-component benchmark."""
    mo.stop(not run_2comp_button.value, mo.md("*Click 'Run 2-Component Benchmark' to start*"))

    results_2comp, metadata_2comp = run_cats_vs_rws(
        pipeline=pipeline_2comp,
        alpha=alpha_slider.value,
        beta=beta_slider.value,
        warmup_k=warmup_k_slider.value,
        num_cycles=num_cycles_slider.value,
        num_iterations=iterations_slider.value,
    )
    return metadata_2comp, results_2comp


@app.cell(hide_code=True)
def _(metadata_2comp, mo, run_2comp_button):
    """Display 2-component results summary."""
    mo.stop(not run_2comp_button.value)

    mo.md(f"""
    ### 2-Component Results Summary

    | Method | Best Score | Avg Best Score | Time (s) |
    |--------|-----------|----------------|----------|
    | CATS   | {metadata_2comp['CATS']['best']:.4f} | {metadata_2comp['CATS']['avg_best']:.4f} | {metadata_2comp['CATS']['time']:.1f} |
    | RWS    | {metadata_2comp['RWS']['best']:.4f} | {metadata_2comp['RWS']['avg_best']:.4f} | {metadata_2comp['RWS']['time']:.1f} |

    **Difference (avg best):** {abs(metadata_2comp['CATS']['avg_best'] - metadata_2comp['RWS']['avg_best']):.4f}
    ({"CATS better" if metadata_2comp['CATS']['avg_best'] < metadata_2comp['RWS']['avg_best'] else "RWS better" if metadata_2comp['RWS']['avg_best'] < metadata_2comp['CATS']['avg_best'] else "Tied"})
    """)
    return


@app.cell
def _(SCORES_FILE, pl):
    """Load reference data for visualizations."""
    reference_data = (
        pl.read_csv(SCORES_FILE)
        .rename({"Product_Code": "Name", "Scores": "score"})
        .select(["Name", "score"])
    )
    return (reference_data,)


@app.cell
def _(
    TS_Benchmarks,
    mo,
    num_cycles_slider,
    reference_data,
    results_2comp,
    run_2comp_button,
    top_n_slider,
):
    """Create TS_Benchmarks for 2-component results."""
    mo.stop(not run_2comp_button.value)

    benchmarks_2comp = TS_Benchmarks(
        no_of_cycles=num_cycles_slider.value,
        methods_list=["CATS", "RWS"],
        TS_runs_data=results_2comp,
        reference_data=reference_data,
        top_n=top_n_slider.value,
        sort_type="minimize",
        top_ns=[25, 50, 100, 200, 300],
    )
    return (benchmarks_2comp,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 2-Component: Score Distributions (Strip Plot)
    """)
    return


@app.cell
def _(benchmarks_2comp, mo, run_2comp_button):
    """2-component strip plot."""
    mo.stop(not run_2comp_button.value)

    benchmarks_2comp.stripplot_TS_results(
        width=900, height=400, show_plot=True, legend_position="bottom"
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 2-Component: Hit Recovery by Cycle (Bar Plot)
    """)
    return


@app.cell
def _(benchmarks_2comp, mo, run_2comp_button):
    """2-component bar plot."""
    mo.stop(not run_2comp_button.value)

    benchmarks_2comp.plot_barplot_TS_results(
        width=900, height=400, show_plot=True, legend_position="bottom", dark_mode=False
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 2-Component: Recovery Fraction vs Top-N (Line Plot)
    """)
    return


@app.cell
def _(benchmarks_2comp, mo, run_2comp_button):
    """2-component line plot."""
    mo.stop(not run_2comp_button.value)

    benchmarks_2comp.plot_line_performance_with_error_bars(
        width=900, height=400, show_plot=True, legend_position="bottom"
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 2-Component: Custom Analysis
    """)
    return


@app.cell
def _(alt, mo, num_cycles_slider, pl, reference_data, results_2comp, run_2comp_button, top_n_slider):
    """2-component custom Altair plots: cumulative recovery and exploration efficiency."""
    mo.stop(not run_2comp_button.value)

    _top_n = top_n_slider.value
    _num_cycles = num_cycles_slider.value

    # Get reference top-N compound names
    _ref_top = set(
        reference_data.sort("score").head(_top_n)["Name"].to_list()
    )

    # --- Cumulative Recovery Curve ---
    _cum_rows = []
    for _method in ["CATS", "RWS"]:
        _seen = set()
        for _c in range(_num_cycles):
            _cycle_names = set(results_2comp[_method][_c]["Name"].to_list())
            _hits = _cycle_names & _ref_top
            _seen = _seen | _hits
            _cum_rows.append({
                "method": _method,
                "cycle": _c + 1,
                "cumulative_fraction": len(_seen) / _top_n,
            })
    _cum_df = pl.DataFrame(_cum_rows)

    cumulative_chart_2comp = (
        alt.Chart(_cum_df.to_pandas())
        .mark_line(point=True, strokeWidth=2)
        .encode(
            x=alt.X("cycle:O", title="Cycle"),
            y=alt.Y("cumulative_fraction:Q", title=f"Cumulative Recovery (top {_top_n})", scale=alt.Scale(domain=[0, 1])),
            color=alt.Color("method:N", title="Method"),
        )
        .properties(width=420, height=300, title="Cumulative Recovery Curve")
    )

    # --- Exploration Efficiency ---
    _eff_rows = []
    for _method in ["CATS", "RWS"]:
        for _c in range(_num_cycles):
            _unique_count = results_2comp[_method][_c]["Name"].n_unique()
            _eff_rows.append({
                "method": _method,
                "cycle": _c + 1,
                "unique_compounds": _unique_count,
            })
    _eff_df = pl.DataFrame(_eff_rows)

    efficiency_chart_2comp = (
        alt.Chart(_eff_df.to_pandas())
        .mark_bar(opacity=0.8)
        .encode(
            x=alt.X("cycle:O", title="Cycle"),
            y=alt.Y("unique_compounds:Q", title="Unique Compounds Sampled"),
            color=alt.Color("method:N", title="Method"),
            xOffset="method:N",
        )
        .properties(width=420, height=300, title="Exploration Efficiency")
    )

    mo.hstack([cumulative_chart_2comp, efficiency_chart_2comp], justify="center", gap=2)
    return cumulative_chart_2comp, efficiency_chart_2comp


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## 3-Component Library: CATS vs RWS

    The 3-component library models the same Thrombin dataset as:
    - **62 amino acids** x **62 amino acids** x **130 carboxylic acids** = ~499,720 products
    - Multi-step synthesis: peptide coupling → amine substitution → amide coupling
    - Thompson Sampling tracks only **254 posteriors** (62 + 62 + 130)

    The 3-component representation assumes amino acid contributions are approximately
    additive. CATS can help here by identifying which components are most critical.
    """)
    return


@app.cell
def _(mo):
    """Run button for 3-component benchmark."""
    run_3comp_button = mo.ui.run_button(label="Run 3-Component Benchmark")
    run_3comp_button
    return (run_3comp_button,)


@app.cell
def _(
    alpha_slider,
    beta_slider,
    iterations_slider,
    mo,
    num_cycles_slider,
    pipeline_3comp,
    run_3comp_button,
    run_cats_vs_rws,
    warmup_k_slider,
):
    """Execute 3-component benchmark."""
    mo.stop(not run_3comp_button.value, mo.md("*Click 'Run 3-Component Benchmark' to start*"))

    results_3comp, metadata_3comp = run_cats_vs_rws(
        pipeline=pipeline_3comp,
        alpha=alpha_slider.value,
        beta=beta_slider.value,
        warmup_k=warmup_k_slider.value,
        num_cycles=num_cycles_slider.value,
        num_iterations=iterations_slider.value,
    )
    return metadata_3comp, results_3comp


@app.cell(hide_code=True)
def _(metadata_3comp, mo, run_3comp_button):
    """Display 3-component results summary."""
    mo.stop(not run_3comp_button.value)

    mo.md(f"""
    ### 3-Component Results Summary

    | Method | Best Score | Avg Best Score | Time (s) |
    |--------|-----------|----------------|----------|
    | CATS   | {metadata_3comp['CATS']['best']:.4f} | {metadata_3comp['CATS']['avg_best']:.4f} | {metadata_3comp['CATS']['time']:.1f} |
    | RWS    | {metadata_3comp['RWS']['best']:.4f} | {metadata_3comp['RWS']['avg_best']:.4f} | {metadata_3comp['RWS']['time']:.1f} |

    **Difference (avg best):** {abs(metadata_3comp['CATS']['avg_best'] - metadata_3comp['RWS']['avg_best']):.4f}
    ({"CATS better" if metadata_3comp['CATS']['avg_best'] < metadata_3comp['RWS']['avg_best'] else "RWS better" if metadata_3comp['RWS']['avg_best'] < metadata_3comp['CATS']['avg_best'] else "Tied"})
    """)
    return


@app.cell
def _(
    TS_Benchmarks,
    mo,
    num_cycles_slider,
    reference_data,
    results_3comp,
    run_3comp_button,
    top_n_slider,
):
    """Create TS_Benchmarks for 3-component results."""
    mo.stop(not run_3comp_button.value)

    benchmarks_3comp = TS_Benchmarks(
        no_of_cycles=num_cycles_slider.value,
        methods_list=["CATS", "RWS"],
        TS_runs_data=results_3comp,
        reference_data=reference_data,
        top_n=top_n_slider.value,
        sort_type="minimize",
        top_ns=[25, 50, 100, 200, 300],
    )
    return (benchmarks_3comp,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 3-Component: Score Distributions (Strip Plot)
    """)
    return


@app.cell
def _(benchmarks_3comp, mo, run_3comp_button):
    """3-component strip plot."""
    mo.stop(not run_3comp_button.value)

    benchmarks_3comp.stripplot_TS_results(
        width=900, height=400, show_plot=True, legend_position="bottom"
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 3-Component: Hit Recovery by Cycle (Bar Plot)
    """)
    return


@app.cell
def _(benchmarks_3comp, mo, run_3comp_button):
    """3-component bar plot."""
    mo.stop(not run_3comp_button.value)

    benchmarks_3comp.plot_barplot_TS_results(
        width=900, height=400, show_plot=True, legend_position="bottom", dark_mode=False
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 3-Component: Recovery Fraction vs Top-N (Line Plot)
    """)
    return


@app.cell
def _(benchmarks_3comp, mo, run_3comp_button):
    """3-component line plot."""
    mo.stop(not run_3comp_button.value)

    benchmarks_3comp.plot_line_performance_with_error_bars(
        width=900, height=400, show_plot=True, legend_position="bottom"
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 3-Component: Custom Analysis
    """)
    return


@app.cell
def _(alt, mo, num_cycles_slider, pl, reference_data, results_3comp, run_3comp_button, top_n_slider):
    """3-component custom Altair plots: cumulative recovery and exploration efficiency."""
    mo.stop(not run_3comp_button.value)

    _top_n = top_n_slider.value
    _num_cycles = num_cycles_slider.value

    # Get reference top-N compound names
    _ref_top = set(
        reference_data.sort("score").head(_top_n)["Name"].to_list()
    )

    # --- Cumulative Recovery Curve ---
    _cum_rows = []
    for _method in ["CATS", "RWS"]:
        _seen = set()
        for _c in range(_num_cycles):
            _cycle_names = set(results_3comp[_method][_c]["Name"].to_list())
            _hits = _cycle_names & _ref_top
            _seen = _seen | _hits
            _cum_rows.append({
                "method": _method,
                "cycle": _c + 1,
                "cumulative_fraction": len(_seen) / _top_n,
            })
    _cum_df = pl.DataFrame(_cum_rows)

    cumulative_chart_3comp = (
        alt.Chart(_cum_df.to_pandas())
        .mark_line(point=True, strokeWidth=2)
        .encode(
            x=alt.X("cycle:O", title="Cycle"),
            y=alt.Y("cumulative_fraction:Q", title=f"Cumulative Recovery (top {_top_n})", scale=alt.Scale(domain=[0, 1])),
            color=alt.Color("method:N", title="Method"),
        )
        .properties(width=420, height=300, title="Cumulative Recovery Curve")
    )

    # --- Exploration Efficiency ---
    _eff_rows = []
    for _method in ["CATS", "RWS"]:
        for _c in range(_num_cycles):
            _unique_count = results_3comp[_method][_c]["Name"].n_unique()
            _eff_rows.append({
                "method": _method,
                "cycle": _c + 1,
                "unique_compounds": _unique_count,
            })
    _eff_df = pl.DataFrame(_eff_rows)

    efficiency_chart_3comp = (
        alt.Chart(_eff_df.to_pandas())
        .mark_bar(opacity=0.8)
        .encode(
            x=alt.X("cycle:O", title="Cycle"),
            y=alt.Y("unique_compounds:Q", title="Unique Compounds Sampled"),
            color=alt.Color("method:N", title="Method"),
            xOffset="method:N",
        )
        .properties(width=420, height=300, title="Exploration Efficiency")
    )

    mo.hstack([cumulative_chart_3comp, efficiency_chart_3comp], justify="center", gap=2)
    return cumulative_chart_3comp, efficiency_chart_3comp


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Combined Comparison: 2-Component vs 3-Component

    Side-by-side comparison of CATS vs RWS behavior across both library representations.
    This section requires running both benchmarks above.
    """)
    return


@app.cell
def _(
    alt,
    cumulative_chart_2comp,
    cumulative_chart_3comp,
    efficiency_chart_2comp,
    efficiency_chart_3comp,
    metadata_2comp,
    metadata_3comp,
    mo,
    run_2comp_button,
    run_3comp_button,
):
    """Combined comparison: faceted charts and summary table."""
    mo.stop(
        not (run_2comp_button.value and run_3comp_button.value),
        mo.md("*Run both 2-Component and 3-Component benchmarks to see the combined comparison.*"),
    )

    # Faceted cumulative recovery: side by side
    _cum_combined = alt.hconcat(
        cumulative_chart_2comp.properties(title="2-Component: Cumulative Recovery"),
        cumulative_chart_3comp.properties(title="3-Component: Cumulative Recovery"),
    ).resolve_scale(color="shared")

    # Faceted exploration efficiency: side by side
    _eff_combined = alt.hconcat(
        efficiency_chart_2comp.properties(title="2-Component: Exploration Efficiency"),
        efficiency_chart_3comp.properties(title="3-Component: Exploration Efficiency"),
    ).resolve_scale(color="shared")

    # Summary comparison table
    _cats_2c = metadata_2comp["CATS"]
    _rws_2c = metadata_2comp["RWS"]
    _cats_3c = metadata_3comp["CATS"]
    _rws_3c = metadata_3comp["RWS"]

    _summary_md = mo.md(f"""
    ### Summary Comparison

    | Library | Method | Best Score | Avg Best Score | Time (s) |
    |---------|--------|-----------|----------------|----------|
    | 2-Comp  | CATS   | {_cats_2c['best']:.4f} | {_cats_2c['avg_best']:.4f} | {_cats_2c['time']:.1f} |
    | 2-Comp  | RWS    | {_rws_2c['best']:.4f} | {_rws_2c['avg_best']:.4f} | {_rws_2c['time']:.1f} |
    | 3-Comp  | CATS   | {_cats_3c['best']:.4f} | {_cats_3c['avg_best']:.4f} | {_cats_3c['time']:.1f} |
    | 3-Comp  | RWS    | {_rws_3c['best']:.4f} | {_rws_3c['avg_best']:.4f} | {_rws_3c['time']:.1f} |

    **2-Component CATS advantage:** {_rws_2c['avg_best'] - _cats_2c['avg_best']:.4f}
    {"(CATS improves avg best score)" if _cats_2c['avg_best'] < _rws_2c['avg_best'] else "(RWS performs equally or better)"}

    **3-Component CATS advantage:** {_rws_3c['avg_best'] - _cats_3c['avg_best']:.4f}
    {"(CATS improves avg best score)" if _cats_3c['avg_best'] < _rws_3c['avg_best'] else "(RWS performs equally or better)"}
    """)

    mo.vstack([_cum_combined, _eff_combined, _summary_md])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Conclusions

    ### Key Takeaways

    1. **CATS adds per-component intelligence** on top of thermal cycling. By measuring
       criticality (how concentrated posterior means are), CATS adjusts temperatures
       dynamically — exploiting critical components and exploring flexible ones.

    2. **The 2-component library** (130 acids + 3,844 dipeptides) has strong interaction
       effects between amino acids within dipeptides. CATS can detect which dipeptides
       are critical vs flexible and allocate search effort accordingly.

    3. **The 3-component library** (62 AA + 62 AA + 130 acids) assumes additivity.
       With only 254 posteriors (vs 3,974), it has less resolution but CATS can still
       identify which amino acid positions contribute most to score variation.

    4. **When to use CATS**: CATS benefits are largest when components have unequal
       criticality — i.e., when some components matter much more than others for
       the objective. If all components contribute equally, CATS reverts to neutral
       behavior (multiplier ≈ 1.0).

    5. **Thermal cycling alone** (RWS) provides exploration through temperature
       alternation but treats all components identically. CATS adds a data-driven
       layer that adapts as posteriors are updated.

    ### Recommendations

    - **Default to CATS** (`RouletteWheelConfig` with default `min_observations=5`)
      unless you have reason to believe all components are equally important.
    - **Use 2-component representation** when amino acid interactions are significant
      (most drug discovery libraries).
    - **Use 3-component representation** for initial exploration when the library
      has many building blocks and you want faster posterior convergence.
    """)
    return


if __name__ == "__main__":
    app.run()
