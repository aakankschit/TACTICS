"""
TACTICS Thompson Sampling Benchmark App

An interactive benchmark tool for comparing selection strategies and warmup methods.
Designed for use with marimo's app mode (`marimo run`).

Features:
1. Multi-strategy selection with customizable parameters
2. Warmup strategy comparison (Standard, Balanced, Enhanced)
3. Strategy x Warmup combination matrix benchmarking
4. TS_Benchmarks visualizations for hit recovery analysis
5. Interactive toggles for visualization filtering

Dataset: Thrombin (amide coupling)

Run as app: marimo run notebooks/thompson_sampling_tutorial.py
Edit mode:  marimo edit notebooks/thompson_sampling_tutorial.py
"""

import marimo

__generated_with = "0.18.2"
app = marimo.App(width="full", app_title="TACTICS Thompson Sampling Benchmark")


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
        StandardWarmupConfig,
        BalancedWarmupConfig,
        EnhancedWarmupConfig,
    )
    from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

    # TACTICS Library Enumeration imports
    from TACTICS.library_enumeration import SynthesisPipeline
    from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef

    # TACTICS Library Analysis imports
    from TACTICS.library_analysis.visualization import TS_Benchmarks
    return (
        BalancedWarmupConfig,
        BoltzmannConfig,
        EnhancedWarmupConfig,
        EpsilonGreedyConfig,
        GreedyConfig,
        LookupEvaluatorConfig,
        ReactionConfig,
        ReactionDef,
        RouletteWheelConfig,
        StandardWarmupConfig,
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
    # TACTICS Thompson Sampling Benchmark

    Compare selection strategies and warmup methods for efficient exploration of
    combinatorial chemical libraries.

    **Dataset:** Thrombin Linear Amide Library (~500K products)
    """)
    return


@app.cell
def _():
    """Load bundled Thrombin dataset paths."""
    import importlib.resources

    _data_files = importlib.resources.files("TACTICS.data.thrombin")
    ACIDS_FILE = str(_data_files / "acids.smi")
    AMINES_FILE = str(_data_files / "coupled_aa_sub.smi")
    SCORES_FILE = str(_data_files / "product_scores.csv")
    REAGENT_FILES = [ACIDS_FILE, AMINES_FILE]
    AMIDE_COUPLING_SMARTS = "[#6:1](=[O:2])[OH].[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]>>[#6:1](=[O:2])[#7:3]"
    return ACIDS_FILE, AMINES_FILE, AMIDE_COUPLING_SMARTS, REAGENT_FILES, SCORES_FILE


@app.cell
def _(mo):
    mo.md(r"""
    ---
    ## Configuration

    ### Selection Strategies
    """)
    return


@app.cell
def _(mo):
    """Strategy selection with parameter customization."""
    # Strategy checkboxes
    strategy_epsilon = mo.ui.checkbox(value=True, label="Epsilon-Greedy")
    strategy_roulette = mo.ui.checkbox(value=False, label="Roulette Wheel (CATS)")
    strategy_ucb = mo.ui.checkbox(value=False, label="UCB")
    strategy_greedy = mo.ui.checkbox(value=True, label="Greedy (baseline)")
    strategy_boltzmann = mo.ui.checkbox(value=False, label="Boltzmann")

    # Strategy parameters
    epsilon_value = mo.ui.slider(start=0.05, stop=0.5, value=0.2, step=0.05, label="Epsilon")
    epsilon_decay = mo.ui.slider(start=0.9, stop=1.0, value=0.995, step=0.005, label="Decay")
    ucb_c = mo.ui.slider(start=0.5, stop=4.0, value=2.0, step=0.5, label="UCB c")
    rw_alpha = mo.ui.slider(start=0.05, stop=0.5, value=0.1, step=0.05, label="RW alpha")
    boltz_temp = mo.ui.slider(start=0.1, stop=2.0, value=1.0, step=0.1, label="Temperature")
    return (
        boltz_temp,
        epsilon_decay,
        epsilon_value,
        rw_alpha,
        strategy_boltzmann,
        strategy_epsilon,
        strategy_greedy,
        strategy_roulette,
        strategy_ucb,
        ucb_c,
    )


@app.cell
def _(
    boltz_temp,
    epsilon_decay,
    epsilon_value,
    mo,
    rw_alpha,
    strategy_boltzmann,
    strategy_epsilon,
    strategy_greedy,
    strategy_roulette,
    strategy_ucb,
    ucb_c,
):
    """Display strategy controls."""
    _strategy_grid = mo.vstack([
        mo.hstack([
            strategy_epsilon, strategy_roulette, strategy_ucb,
            strategy_greedy, strategy_boltzmann
        ], justify="start", gap=2),
        mo.md("**Strategy Parameters:**"),
        mo.hstack([
            mo.vstack([mo.md("*Epsilon-Greedy:*"), epsilon_value, epsilon_decay]),
            mo.vstack([mo.md("*UCB:*"), ucb_c]),
            mo.vstack([mo.md("*Roulette Wheel:*"), rw_alpha]),
            mo.vstack([mo.md("*Boltzmann:*"), boltz_temp]),
        ], justify="start", gap=4),
    ])
    _strategy_grid
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Warmup Strategies
    """)
    return


@app.cell
def _(mo):
    """Warmup strategy selection."""
    warmup_standard = mo.ui.checkbox(value=False, label="Standard (random)")
    warmup_balanced_3 = mo.ui.checkbox(value=True, label="Balanced K=3")
    warmup_balanced_5 = mo.ui.checkbox(value=False, label="Balanced K=5")
    warmup_enhanced = mo.ui.checkbox(value=False, label="Enhanced (stochastic)")
    return (
        warmup_balanced_3,
        warmup_balanced_5,
        warmup_enhanced,
        warmup_standard,
    )


@app.cell
def _(
    mo,
    warmup_balanced_3,
    warmup_balanced_5,
    warmup_enhanced,
    warmup_standard,
):
    """Display warmup controls."""
    mo.hstack([
        warmup_standard, warmup_balanced_3, warmup_balanced_5, warmup_enhanced
    ], justify="start", gap=2)
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
    combination_mode = mo.ui.checkbox(
        value=False,
        label="Run all Strategy x Warmup combinations"
    )
    return combination_mode, cycles_slider, iterations_slider, top_n_slider


@app.cell(hide_code=True)
def _(combination_mode, cycles_slider, iterations_slider, mo, top_n_slider):
    """Display search parameter controls."""
    mo.vstack([
        mo.hstack([iterations_slider, cycles_slider, top_n_slider], justify="start", gap=2),
        combination_mode,
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Run Benchmark
    """)
    return


@app.cell
def _(mo):
    """Run button."""
    run_button = mo.ui.run_button(label="Run Thompson Sampling Benchmark")
    run_button
    return (run_button,)


@app.cell
def _(
    AMIDE_COUPLING_SMARTS,
    BalancedWarmupConfig,
    BoltzmannConfig,
    EnhancedWarmupConfig,
    EpsilonGreedyConfig,
    GreedyConfig,
    LookupEvaluatorConfig,
    REAGENT_FILES,
    ReactionConfig,
    ReactionDef,
    RouletteWheelConfig,
    SCORES_FILE,
    StandardWarmupConfig,
    SynthesisPipeline,
    ThompsonSampler,
    ThompsonSamplingConfig,
    UCBConfig,
    boltz_temp,
    combination_mode,
    cycles_slider,
    epsilon_decay,
    epsilon_value,
    iterations_slider,
    mo,
    pl,
    run_button,
    rw_alpha,
    strategy_boltzmann,
    strategy_epsilon,
    strategy_greedy,
    strategy_roulette,
    strategy_ucb,
    time,
    ucb_c,
    warmup_balanced_3,
    warmup_balanced_5,
    warmup_enhanced,
    warmup_standard,
):
    """Execute Thompson Sampling benchmark."""
    mo.stop(not run_button.value, mo.md("*Click 'Run Thompson Sampling Benchmark' to start*"))

    # Build selected strategies
    selected_strategies = {}
    if strategy_epsilon.value:
        selected_strategies["Epsilon-Greedy"] = lambda: EpsilonGreedyConfig(
            mode="minimize", epsilon=epsilon_value.value, decay=epsilon_decay.value
        )
    if strategy_roulette.value:
        selected_strategies["Roulette Wheel"] = lambda: RouletteWheelConfig(
            mode="minimize", alpha=rw_alpha.value, beta=rw_alpha.value
        )
    if strategy_ucb.value:
        selected_strategies["UCB"] = lambda: UCBConfig(mode="minimize", c=ucb_c.value)
    if strategy_greedy.value:
        selected_strategies["Greedy"] = lambda: GreedyConfig(mode="minimize")
    if strategy_boltzmann.value:
        selected_strategies["Boltzmann"] = lambda: BoltzmannConfig(
            mode="minimize", temperature=boltz_temp.value
        )

    if not selected_strategies:
        mo.stop(True, mo.md("**Error:** Select at least one strategy."))

    # Build selected warmups
    selected_warmups = {}
    if warmup_standard.value:
        selected_warmups["Standard"] = StandardWarmupConfig()
    if warmup_balanced_3.value:
        selected_warmups["Balanced-K3"] = BalancedWarmupConfig(observations_per_reagent=3)
    if warmup_balanced_5.value:
        selected_warmups["Balanced-K5"] = BalancedWarmupConfig(observations_per_reagent=5)
    if warmup_enhanced.value:
        selected_warmups["Enhanced"] = EnhancedWarmupConfig()

    if not selected_warmups:
        mo.stop(True, mo.md("**Error:** Select at least one warmup strategy."))

    # Create synthesis pipeline
    _reaction_config = ReactionConfig(
        reactions=[ReactionDef(reaction_smarts=AMIDE_COUPLING_SMARTS, step_index=0)],
        reagent_file_list=REAGENT_FILES,
    )
    _pipeline = SynthesisPipeline(_reaction_config)

    # Create evaluator config
    _evaluator_config = LookupEvaluatorConfig(
        ref_filename=SCORES_FILE,
        compound_col="Product_Code",
        score_col="Scores",
    )

    # Determine combinations to run
    if combination_mode.value:
        # Run all Strategy x Warmup combinations
        run_combinations = [
            (s_name, w_name, s_factory, w_config)
            for s_name, s_factory in selected_strategies.items()
            for w_name, w_config in selected_warmups.items()
        ]
    else:
        # Run each strategy with first warmup only
        first_warmup_name = list(selected_warmups.keys())[0]
        first_warmup_config = selected_warmups[first_warmup_name]
        run_combinations = [
            (s_name, first_warmup_name, s_factory, first_warmup_config)
            for s_name, s_factory in selected_strategies.items()
        ]

    # Run parameters
    _num_cycles = cycles_slider.value
    _num_iterations = iterations_slider.value

    # Store results: {combination_name: [cycle1_df, cycle2_df, ...]}
    all_results = {}
    run_metadata = []
    start_time = time.time()

    for s_name, w_name, s_factory, w_config in run_combinations:
        combo_name = f"{s_name} + {w_name}"
        _cycle_results = []

        for _cycle in range(_num_cycles):
            _config = ThompsonSamplingConfig(
                synthesis_pipeline=_pipeline,
                num_ts_iterations=_num_iterations,
                num_warmup_trials=3,
                strategy_config=s_factory(),
                warmup_config=w_config,
                evaluator_config=_evaluator_config,
                batch_size=1,
                max_resamples=1000,
                hide_progress=True,
            )

            _sampler = ThompsonSampler.from_config(_config)
            _warmup_df = _sampler.warm_up(num_warmup_trials=_config.num_warmup_trials)
            _search_df = _sampler.search(num_cycles=_config.num_ts_iterations)

            _combined = pl.concat([_warmup_df, _search_df])
            _result = _combined.select(["Name", "score"])
            _cycle_results.append(_result)
            _sampler.close()

        all_results[combo_name] = _cycle_results

        # Compute metadata
        _best_scores = [r["score"].min() for r in _cycle_results]
        run_metadata.append({
            "combination": combo_name,
            "strategy": s_name,
            "warmup": w_name,
            "cycles": _num_cycles,
            "best_score": min(_best_scores),
            "avg_best": sum(_best_scores) / len(_best_scores),
        })

    total_time = time.time() - start_time
    combination_names = list(all_results.keys())
    return all_results, combination_names, run_metadata, total_time


@app.cell
def _(mo, run_button, run_metadata, total_time):
    """Display results summary."""
    mo.stop(not run_button.value)

    _rows = []
    for m in run_metadata:
        _rows.append(
            f"| {m['combination']} | {m['strategy']} | {m['warmup']} | "
            f"{m['cycles']} | {m['best_score']:.3f} | {m['avg_best']:.3f} |"
        )
    _table = "\n".join(_rows)

    mo.md(f"""
    ---
    ## Results Summary

    **Total runtime:** {total_time:.1f}s

    | Combination | Strategy | Warmup | Cycles | Best Score | Avg Best |
    |-------------|----------|--------|--------|------------|----------|
    {_table}
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### Top Compounds Found
    """)
    return


@app.cell
def _(all_results, combination_names, mo, pl, run_button):
    """Display top compounds table."""
    mo.stop(not run_button.value)

    _all_dfs = []
    for _name in combination_names:
        _results = all_results[_name]
        for _i, _df in enumerate(_results):
            _labeled = _df.with_columns([
                pl.lit(_name).alias("combination"),
                pl.lit(_i + 1).alias("cycle"),
            ])
            _all_dfs.append(_labeled)

    combined_all = pl.concat(_all_dfs)
    top_20 = (
        combined_all
        .sort("score")
        .head(20)
        .select(["Name", "score", "combination", "cycle"])
    )
    top_20
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Benchmark Visualizations

    ### Toggle Combinations for Comparison
    """)
    return


@app.cell
def _(SCORES_FILE, mo, pl, run_button):
    """Load reference data."""
    mo.stop(not run_button.value)

    reference_data = (
        pl.read_csv(SCORES_FILE)
        .rename({"Product_Code": "Name", "Scores": "score"})
        .select(["Name", "score"])
    )
    return (reference_data,)


@app.cell
def _(combination_names, mo, run_button):
    """Create toggle checkboxes for visualization."""
    mo.stop(not run_button.value)

    viz_toggles = mo.ui.array([
        mo.ui.checkbox(value=True, label=name)
        for name in combination_names
    ])
    return (viz_toggles,)


@app.cell(hide_code=True)
def _(mo, run_button, viz_toggles):
    """Display visualization toggles."""
    mo.stop(not run_button.value)

    mo.vstack([
        mo.md("**Select combinations to visualize:**"),
        mo.hstack(viz_toggles, justify="start", gap=2, wrap=True),
    ])
    return


@app.cell
def _(all_results, combination_names, mo, run_button, viz_toggles):
    """Get selected combinations for visualization."""
    mo.stop(not run_button.value)

    viz_combinations = [
        name for i, name in enumerate(combination_names)
        if viz_toggles[i].value
    ]
    viz_results = {name: all_results[name] for name in viz_combinations}

    mo.md(f"**Visualizing:** {', '.join(viz_combinations) if viz_combinations else 'None'}")
    return viz_combinations, viz_results


@app.cell
def _(
    TS_Benchmarks,
    cycles_slider,
    mo,
    reference_data,
    run_button,
    top_n_slider,
    viz_combinations,
    viz_results,
):
    """Create TS_Benchmarks instance - reactive to toggle changes."""
    mo.stop(not run_button.value)

    # Handle empty selection gracefully without breaking reactivity
    if len(viz_combinations) == 0:
        benchmarks = None
    else:
        benchmarks = TS_Benchmarks(
            no_of_cycles=cycles_slider.value,
            methods_list=viz_combinations,
            TS_runs_data=viz_results,
            reference_data=reference_data,
            top_n=top_n_slider.value,
            sort_type="minimize",
            top_ns=[25, 50, 100, 200, 300],
        )
    return (benchmarks,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Strip Plot: Score Distributions
    """)
    return


@app.cell
def _(benchmarks, mo, run_button):
    """Generate strip plot - reactive to toggle changes."""
    mo.stop(not run_button.value)

    # Use conditional rendering to maintain reactivity
    _output = (
        mo.md("*Select at least one combination to visualize*")
        if benchmarks is None
        else benchmarks.stripplot_TS_results(
            width=900,
            height=400,
            show_plot=True,
            legend_position="bottom"
        )
    )
    _output
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Bar Plot: Hit Recovery by Cycle
    """)
    return


@app.cell
def _(benchmarks, mo, run_button):
    """Generate bar plot - reactive to toggle changes."""
    mo.stop(not run_button.value)

    # Use conditional rendering to maintain reactivity
    _output = (
        mo.md("*Select at least one combination to visualize*")
        if benchmarks is None
        else benchmarks.plot_barplot_TS_results(
            width=900,
            height=400,
            show_plot=True,
            legend_position="bottom",
            dark_mode=False
        )
    )
    _output
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Line Plot: Performance vs Top-N
    """)
    return


@app.cell
def _(benchmarks, mo, run_button):
    """Generate line plot - reactive to toggle changes."""
    mo.stop(not run_button.value)

    # Use conditional rendering to maintain reactivity
    _output = (
        mo.md("*Select at least one combination to visualize*")
        if benchmarks is None
        else benchmarks.plot_line_performance_with_error_bars(
            width=900,
            height=450,
            show_plot=True,
            legend_position="bottom"
        )
    )
    _output
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Performance Statistics
    """)
    return


@app.cell(hide_code=True)
def _(benchmarks, mo, run_button):
    """Display grouped statistics - reactive to toggle changes."""
    mo.stop(not run_button.value)

    # Use conditional rendering to maintain reactivity
    if benchmarks is None:
        _stats = mo.md("*Select at least one combination to visualize*")
    elif benchmarks.grouped_stats is not None:
        _stats = mo.vstack([
            mo.md("**Mean Performance by Top-N:**"),
            benchmarks.grouped_stats,
        ])
    else:
        _stats = mo.md("*Statistics not available*")
    _stats
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Quick Reference

    ### Selection Strategies

    | Strategy | Description | Parameters |
    |----------|-------------|------------|
    | **Greedy** | Pure exploitation - always picks best | None |
    | **Epsilon-Greedy** | Random exploration with probability epsilon | `epsilon`, `decay` |
    | **UCB** | Upper confidence bound exploration | `c` (exploration weight) |
    | **Roulette Wheel** | Thermal cycling (CATS) | `alpha`, `beta` |
    | **Boltzmann** | Softmax selection | `temperature` |

    ### Warmup Strategies

    | Warmup | Description | Observations |
    |--------|-------------|--------------|
    | **Standard** | Random partner selection | K x num_reagents per component |
    | **Balanced** | Stratified K per reagent | K x total_reagents |
    | **Enhanced** | Stochastic parallel pairing | Variable |

    ### Usage Tips

    1. **Quick comparison**: Select Epsilon-Greedy + Greedy with Balanced-K3
    2. **Comprehensive benchmark**: Enable "Run all combinations" checkbox
    3. **Fine-tuning**: Adjust strategy parameters using sliders
    4. **Analysis**: Use toggles to compare specific combinations
    """)
    return


if __name__ == "__main__":
    app.run()
