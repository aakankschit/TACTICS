"""
TACTICS Legacy vs Current Implementation Benchmark

This notebook compares the performance of:
1. Legacy StandardThompsonSampler (greedy selection)
2. Legacy EnhancedThompsonSampler (roulette wheel with CATS)
3. Current ThompsonSampler with various strategies

Dataset: Thrombin (amide coupling)

Run as app: marimo run notebooks/legacy_vs_current_benchmark.py
Edit mode:  marimo edit notebooks/legacy_vs_current_benchmark.py
"""

import marimo

__generated_with = "0.18.2"
app = marimo.App(width="full", app_title="TACTICS Legacy vs Current Benchmark")


@app.cell
def _():
    """Imports and project setup."""
    import marimo as mo
    import sys
    from pathlib import Path
    import time
    import tempfile
    import numpy as np

    # Add TACTICS project paths
    try:
        project_root = Path(__file__).parent.parent.resolve()
    except NameError:
        project_root = Path("/Users/aakankschitnandkeolyar/Desktop/TACTICS")

    sys.path.insert(0, str(project_root))
    sys.path.insert(0, str(project_root / "src"))

    import polars as pl

    # Current TACTICS imports
    from TACTICS.thompson_sampling import ThompsonSampler
    from TACTICS.thompson_sampling.config import ThompsonSamplingConfig
    from TACTICS.thompson_sampling.strategies.config import (
        GreedyConfig,
        RouletteWheelConfig,
        EpsilonGreedyConfig,
    )
    from TACTICS.thompson_sampling.warmup.config import (
        StandardWarmupConfig,
        BalancedWarmupConfig,
        EnhancedWarmupConfig,
    )
    from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

    from TACTICS.library_enumeration import SynthesisPipeline
    from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef

    from TACTICS.library_analysis.visualization import TS_Benchmarks

    # Legacy TACTICS imports
    from TACTICS.thompson_sampling.legacy.standard_thompson_sampling import StandardThompsonSampler
    from TACTICS.thompson_sampling.legacy.enhanced_thompson_sampling import EnhancedThompsonSampler
    from TACTICS.thompson_sampling.legacy.evaluators import LookupEvaluator as LegacyLookupEvaluator
    return (
        BalancedWarmupConfig,
        EnhancedThompsonSampler,
        EnhancedWarmupConfig,
        EpsilonGreedyConfig,
        GreedyConfig,
        LegacyLookupEvaluator,
        LookupEvaluatorConfig,
        ReactionConfig,
        ReactionDef,
        RouletteWheelConfig,
        StandardThompsonSampler,
        StandardWarmupConfig,
        SynthesisPipeline,
        TS_Benchmarks,
        ThompsonSampler,
        ThompsonSamplingConfig,
        mo,
        np,
        pl,
        project_root,
        tempfile,
        time,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Legacy vs Current Thompson Sampling Benchmark

    Compare the original TACTICS Thompson Sampling implementations with the current
    refactored version to validate performance parity and identify improvements.

    ## Implementations Compared

    | Implementation | Description | Selection Method |
    |----------------|-------------|------------------|
    | **Legacy Standard** | Original greedy TS | argmax/argmin of samples |
    | **Legacy Enhanced** | CATS with roulette wheel | Thermal cycling + % library |
    | **Current Greedy** | Refactored greedy | Same as legacy standard |
    | **Current Epsilon-Greedy** | Exploration-exploitation | Random with epsilon prob |
    | **Current Roulette Wheel** | Refactored CATS | Thermal cycling (iterations) |

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


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Configuration

    ### Select Implementations to Compare
    """)
    return


@app.cell
def _(mo):
    """Implementation selection checkboxes."""
    # Legacy implementations
    run_legacy_standard = mo.ui.checkbox(value=True, label="Legacy Standard (Greedy)")
    run_legacy_enhanced = mo.ui.checkbox(value=False, label="Legacy Enhanced (CATS/RWS)")

    # Current implementations
    run_current_greedy = mo.ui.checkbox(value=True, label="Current Greedy")
    run_current_epsilon = mo.ui.checkbox(value=True, label="Current Epsilon-Greedy")
    run_current_roulette = mo.ui.checkbox(value=False, label="Current Roulette Wheel")
    return (
        run_current_epsilon,
        run_current_greedy,
        run_current_roulette,
        run_legacy_enhanced,
        run_legacy_standard,
    )


@app.cell
def _(
    mo,
    run_current_epsilon,
    run_current_greedy,
    run_current_roulette,
    run_legacy_enhanced,
    run_legacy_standard,
):
    """Display implementation selection."""
    mo.vstack([
        mo.md("**Legacy Implementations:**"),
        mo.hstack([run_legacy_standard, run_legacy_enhanced], justify="start", gap=2),
        mo.md("**Current Implementations:**"),
        mo.hstack([run_current_greedy, run_current_epsilon, run_current_roulette], justify="start", gap=2),
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Warmup Strategy (for Current Implementations)
    """)
    return


@app.cell
def _(mo):
    """Warmup strategy selection for current implementations."""
    warmup_dropdown = mo.ui.dropdown(
        options={
            "Standard (random)": "standard",
            "Balanced": "balanced",
            "Enhanced (stochastic)": "enhanced",
        },
        value="Balanced",
        label="Warmup Strategy",
    )

    balanced_k_slider = mo.ui.slider(
        start=1, stop=10, value=3, step=1,
        label="Balanced K (observations per reagent)"
    )
    return balanced_k_slider, warmup_dropdown


@app.cell
def _(balanced_k_slider, mo, warmup_dropdown):
    """Display warmup controls."""
    mo.hstack([warmup_dropdown, balanced_k_slider], justify="start", gap=2)
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
    warmup_trials_slider = mo.ui.slider(
        start=1, stop=10, value=3, step=1,
        label="Warmup trials (Legacy Standard)"
    )
    top_n_slider = mo.ui.slider(
        start=50, stop=500, value=100, step=50,
        label="Top N for recovery"
    )
    return cycles_slider, iterations_slider, top_n_slider, warmup_trials_slider


@app.cell
def _(cycles_slider, iterations_slider, mo, top_n_slider, warmup_trials_slider):
    """Display search parameter controls."""
    mo.hstack([iterations_slider, cycles_slider, warmup_trials_slider, top_n_slider], justify="start", gap=2)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Legacy Enhanced (CATS) Parameters

    The legacy Enhanced sampler uses a different interface based on percent of library to screen.
    """)
    return


@app.cell
def _(mo):
    """Legacy Enhanced specific parameters."""
    legacy_percent_lib = mo.ui.slider(
        start=0.001, stop=0.02, value=0.005, step=0.001,
        label="Percent of library to screen"
    )
    legacy_scaling = mo.ui.slider(
        start=-1, stop=1, value=-1, step=1,
        label="Scaling (-1=minimize, 1=maximize)"
    )
    return legacy_percent_lib, legacy_scaling


@app.cell
def _(legacy_percent_lib, legacy_scaling, mo, run_legacy_enhanced):
    """Display legacy enhanced parameters."""
    mo.stop(not run_legacy_enhanced.value)
    mo.hstack([legacy_percent_lib, legacy_scaling], justify="start", gap=2)
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
    run_button = mo.ui.run_button(label="Run Legacy vs Current Benchmark")
    run_button
    return (run_button,)


@app.cell
def _(
    AMIDE_COUPLING_SMARTS,
    BalancedWarmupConfig,
    EnhancedThompsonSampler,
    EnhancedWarmupConfig,
    EpsilonGreedyConfig,
    GreedyConfig,
    LegacyLookupEvaluator,
    LookupEvaluatorConfig,
    REAGENT_FILES,
    ReactionConfig,
    ReactionDef,
    RouletteWheelConfig,
    SCORES_FILE,
    StandardThompsonSampler,
    StandardWarmupConfig,
    SynthesisPipeline,
    ThompsonSampler,
    ThompsonSamplingConfig,
    balanced_k_slider,
    cycles_slider,
    iterations_slider,
    legacy_percent_lib,
    legacy_scaling,
    mo,
    np,
    pl,
    run_button,
    run_current_epsilon,
    run_current_greedy,
    run_current_roulette,
    run_legacy_enhanced,
    run_legacy_standard,
    tempfile,
    time,
    warmup_dropdown,
    warmup_trials_slider,
):
    """Execute benchmark for all selected implementations."""
    mo.stop(not run_button.value, mo.md("*Click 'Run Legacy vs Current Benchmark' to start*"))

    # Check at least one is selected
    any_selected = (
        run_legacy_standard.value or
        run_legacy_enhanced.value or
        run_current_greedy.value or
        run_current_epsilon.value or
        run_current_roulette.value
    )
    if not any_selected:
        mo.stop(True, mo.md("**Error:** Select at least one implementation."))

    # Parameters
    _num_cycles = cycles_slider.value
    _num_iterations = iterations_slider.value
    _num_warmup = warmup_trials_slider.value
    _balanced_k = balanced_k_slider.value

    # Results storage
    all_results = {}
    run_metadata = []
    start_time = time.time()

    # =========================================================================
    # LEGACY STANDARD THOMPSON SAMPLER
    # =========================================================================
    if run_legacy_standard.value:
        method_name = "Legacy Standard"
        _cycle_results = []

        for _cycle in range(_num_cycles):
            # Create legacy evaluator
            _legacy_eval = LegacyLookupEvaluator({
                "ref_filename": SCORES_FILE,
                "ref_colname": "Scores",
                "compound_col": "Product_Code",
            })

            # Create and configure legacy sampler
            _legacy_sampler = StandardThompsonSampler(mode="minimize")
            _legacy_sampler.set_hide_progress(True)
            _legacy_sampler.set_evaluator(_legacy_eval)
            _legacy_sampler.read_reagents(REAGENT_FILES)
            _legacy_sampler.set_reaction(AMIDE_COUPLING_SMARTS)

            # Run warmup and search
            _warmup_results = _legacy_sampler.warm_up(num_warmup_trials=_num_warmup)
            _search_results = _legacy_sampler.search(num_cycles=_num_iterations)

            # Convert to polars DataFrame
            _all_legacy = _warmup_results + _search_results
            _legacy_df = pl.DataFrame({
                "Name": [r[2] for r in _all_legacy],
                "score": [r[0] for r in _all_legacy],
            }).filter(pl.col("score").is_not_nan())

            _cycle_results.append(_legacy_df)

        all_results[method_name] = _cycle_results
        _best_scores = [r["score"].min() for r in _cycle_results]
        run_metadata.append({
            "method": method_name,
            "type": "Legacy",
            "warmup": f"Standard (K={_num_warmup})",
            "cycles": _num_cycles,
            "best_score": min(_best_scores),
            "avg_best": sum(_best_scores) / len(_best_scores),
        })

    # =========================================================================
    # LEGACY ENHANCED THOMPSON SAMPLER (CATS/RWS)
    # =========================================================================
    if run_legacy_enhanced.value:
        method_name = "Legacy Enhanced"
        _cycle_results = []

        for _cycle in range(_num_cycles):
            # Create a temp file for results (required by legacy enhanced)
            _temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
            _temp_filename = _temp_file.name
            _temp_file.close()

            # Create legacy evaluator
            _legacy_eval = LegacyLookupEvaluator({
                "ref_filename": SCORES_FILE,
                "ref_colname": "Scores",
                "compound_col": "Product_Code",
            })

            # Create and configure legacy enhanced sampler
            _legacy_enhanced = EnhancedThompsonSampler(
                processes=1,
                scaling=legacy_scaling.value,
                percent_lib=legacy_percent_lib.value,
                search_stop=6000,
                min_cpds_per_core=50,
                results_filename=_temp_filename,
            )
            _legacy_enhanced.set_hide_progress(True)
            _legacy_enhanced.set_evaluator(_legacy_eval)
            _legacy_enhanced.read_reagents(REAGENT_FILES)
            _legacy_enhanced.set_reaction(AMIDE_COUPLING_SMARTS)

            # Run warmup and search
            _legacy_enhanced.warm_up(num_warmup_trials=_num_warmup)
            _search_results = _legacy_enhanced.search(
                percent_of_library=legacy_percent_lib.value,
                stop=6000,
                min_cpds_per_core=50,
            )

            # Read results from temp file and convert
            try:
                _legacy_df = pl.read_csv(_temp_filename, null_values=["nan"])
                _legacy_df = _legacy_df.rename({"SMILES": "smiles"}).select(["Name", "score"])
                _legacy_df = _legacy_df.filter(pl.col("score").is_not_nan())
            except Exception:
                # If file reading fails, create from search results
                _legacy_df = pl.DataFrame({
                    "Name": [r[2] for r in _search_results],
                    "score": [r[0] for r in _search_results],
                }).filter(pl.col("score").is_not_nan())

            _cycle_results.append(_legacy_df)

            # Clean up temp file
            import os
            try:
                os.unlink(_temp_filename)
            except Exception:
                pass

        all_results[method_name] = _cycle_results
        _best_scores = [r["score"].min() for r in _cycle_results if len(r) > 0]
        if _best_scores:
            run_metadata.append({
                "method": method_name,
                "type": "Legacy",
                "warmup": f"Enhanced (K={_num_warmup})",
                "cycles": _num_cycles,
                "best_score": min(_best_scores),
                "avg_best": sum(_best_scores) / len(_best_scores),
            })

    # =========================================================================
    # CURRENT IMPLEMENTATIONS
    # =========================================================================

    # Shared pipeline and evaluator config for current implementations
    _reaction_config = ReactionConfig(
        reactions=[ReactionDef(reaction_smarts=AMIDE_COUPLING_SMARTS, step_index=0)],
        reagent_file_list=REAGENT_FILES,
    )
    _pipeline = SynthesisPipeline(_reaction_config)
    _evaluator_config = LookupEvaluatorConfig(
        ref_filename=SCORES_FILE,
        compound_col="Product_Code",
        score_col="Scores",
    )

    # Build warmup config based on selection
    _warmup_name = warmup_dropdown.value
    if _warmup_name == "standard":
        _warmup_config = StandardWarmupConfig()
        _warmup_display = "Standard"
    elif _warmup_name == "balanced":
        _warmup_config = BalancedWarmupConfig(observations_per_reagent=_balanced_k)
        _warmup_display = f"Balanced (K={_balanced_k})"
    else:  # enhanced
        _warmup_config = EnhancedWarmupConfig()
        _warmup_display = "Enhanced"

    # Helper function for current implementations
    def run_current_implementation(method_name, strategy_config):
        _cycle_results = []
        for _cycle in range(_num_cycles):
            _config = ThompsonSamplingConfig(
                synthesis_pipeline=_pipeline,
                num_ts_iterations=_num_iterations,
                num_warmup_trials=_balanced_k if _warmup_name == "balanced" else _num_warmup,
                strategy_config=strategy_config,
                warmup_config=_warmup_config,
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

        all_results[method_name] = _cycle_results
        _best_scores = [r["score"].min() for r in _cycle_results]
        run_metadata.append({
            "method": method_name,
            "type": "Current",
            "warmup": _warmup_display,
            "cycles": _num_cycles,
            "best_score": min(_best_scores),
            "avg_best": sum(_best_scores) / len(_best_scores),
        })

    # Current Greedy
    if run_current_greedy.value:
        run_current_implementation("Current Greedy", GreedyConfig(mode="minimize"))

    # Current Epsilon-Greedy
    if run_current_epsilon.value:
        run_current_implementation(
            "Current Epsilon-Greedy",
            EpsilonGreedyConfig(mode="minimize", epsilon=0.2, decay=0.995)
        )

    # Current Roulette Wheel
    if run_current_roulette.value:
        run_current_implementation(
            "Current Roulette Wheel",
            RouletteWheelConfig(mode="minimize", alpha=0.1, beta=0.1)
        )

    total_time = time.time() - start_time
    method_names = list(all_results.keys())
    return all_results, method_names, run_metadata, total_time


@app.cell(hide_code=True)
def _(mo, run_button, run_metadata, total_time):
    """Display results summary."""
    mo.stop(not run_button.value)

    _rows = []
    for m in run_metadata:
        _rows.append(
            f"| {m['method']} | {m['type']} | {m['warmup']} | {m['cycles']} | "
            f"{m['best_score']:.3f} | {m['avg_best']:.3f} |"
        )
    _table = "\n".join(_rows)

    mo.md(f"""
    ---
    ## Results Summary

    **Total runtime:** {total_time:.1f}s

    | Method | Type | Warmup | Cycles | Best Score | Avg Best |
    |--------|------|--------|--------|------------|----------|
    {_table}
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Performance Comparison Analysis

    Key questions to answer:
    1. Does **Current Greedy** match **Legacy Standard** performance?
    2. Does **Current Roulette Wheel** match **Legacy Enhanced** performance?
    3. How does **Current Epsilon-Greedy** compare to both?
    4. How does **warmup strategy** affect results?
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Top Compounds Found
    """)
    return


@app.cell(hide_code=True)
def _(all_results, method_names, mo, pl, run_button):
    """Display top compounds table."""
    mo.stop(not run_button.value)

    _all_dfs = []
    for _name in method_names:
        _results = all_results[_name]
        for _i, _df in enumerate(_results):
            if len(_df) > 0:
                _labeled = _df.with_columns([
                    pl.lit(_name).alias("method"),
                    pl.lit(_i + 1).alias("cycle"),
                ])
                _all_dfs.append(_labeled)

    if _all_dfs:
        combined_all = pl.concat(_all_dfs)
        top_20 = (
            combined_all
            .sort("score")
            .head(20)
            .select(["Name", "score", "method", "cycle"])
        )
        top_20
    else:
        mo.md("*No results to display*")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Benchmark Visualizations

    ### Toggle Methods for Comparison
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
def _(method_names, mo, run_button):
    """Create toggle checkboxes for visualization."""
    mo.stop(not run_button.value)

    viz_toggles = mo.ui.array([
        mo.ui.checkbox(value=True, label=name)
        for name in method_names
    ])
    return (viz_toggles,)


@app.cell(hide_code=True)
def _(method_names, mo, run_button, viz_toggles):
    """Display visualization toggles."""
    mo.stop(not run_button.value)

    mo.vstack([
        mo.md("**Select methods to visualize:**"),
        mo.hstack(viz_toggles, justify="start", gap=2, wrap=True),
    ])
    return


@app.cell
def _(all_results, method_names, mo, run_button, viz_toggles):
    """Get selected methods for visualization."""
    mo.stop(not run_button.value)

    viz_methods = [
        name for i, name in enumerate(method_names)
        if viz_toggles[i].value
    ]
    viz_results = {name: all_results[name] for name in viz_methods}

    mo.md(f"**Visualizing:** {', '.join(viz_methods) if viz_methods else 'None'}")
    return viz_methods, viz_results


@app.cell
def _(
    TS_Benchmarks,
    cycles_slider,
    mo,
    reference_data,
    run_button,
    top_n_slider,
    viz_methods,
    viz_results,
):
    """Create TS_Benchmarks instance - reactive to toggle changes."""
    mo.stop(not run_button.value)

    # Handle empty selection gracefully without breaking reactivity
    if len(viz_methods) == 0:
        benchmarks = None
    else:
        # Filter out empty results
        _filtered_results = {
            k: v for k, v in viz_results.items()
            if all(len(df) > 0 for df in v)
        }
        _filtered_methods = list(_filtered_results.keys())

        if not _filtered_methods:
            benchmarks = None
        else:
            benchmarks = TS_Benchmarks(
                no_of_cycles=cycles_slider.value,
                methods_list=_filtered_methods,
                TS_runs_data=_filtered_results,
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
        mo.md("*Select at least one method to visualize*")
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

    Compare how legacy and current implementations recover top hits.
    """)
    return


@app.cell
def _(benchmarks, mo, run_button):
    """Generate bar plot - reactive to toggle changes."""
    mo.stop(not run_button.value)

    # Use conditional rendering to maintain reactivity
    _output = (
        mo.md("*Select at least one method to visualize*")
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

    This plot shows recovery rates at different top-N thresholds.
    Legacy and current implementations should show similar curves if they're equivalent.
    """)
    return


@app.cell
def _(benchmarks, mo, run_button):
    """Generate line plot - reactive to toggle changes."""
    mo.stop(not run_button.value)

    # Use conditional rendering to maintain reactivity
    _output = (
        mo.md("*Select at least one method to visualize*")
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
        _stats = mo.md("*Select at least one method to visualize*")
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
    ## Implementation Comparison Notes

    ### Legacy vs Current Mapping

    | Legacy | Current Equivalent | Key Differences |
    |--------|-------------------|-----------------|
    | `StandardThompsonSampler(mode="minimize")` | `ThompsonSampler` + `GreedyConfig` | Same algorithm, refactored |
    | `StandardThompsonSampler(mode="minimize_boltzmann")` | `ThompsonSampler` + `BoltzmannConfig` | Boltzmann temperature handling |
    | `EnhancedThompsonSampler` | `ThompsonSampler` + `RouletteWheelConfig` | CATS thermal cycling |

    ### Warmup Strategy Comparison

    | Warmup | Description | When to Use |
    |--------|-------------|-------------|
    | **Standard** | Random partner selection | Quick baseline |
    | **Balanced** | K observations per reagent | Best coverage (recommended) |
    | **Enhanced** | Stochastic parallel pairing | Large libraries |

    ### Expected Results

    - **Legacy Standard** and **Current Greedy** should produce statistically similar results
    - **Current Epsilon-Greedy** typically shows better exploration (more diverse hits)
    - **Current Roulette Wheel** should match **Legacy Enhanced** behavior

    ### Key Improvements in Current Version

    1. **Unified API**: Single `ThompsonSampler` class with pluggable strategies
    2. **Pydantic Configuration**: Type-safe, validated configuration objects
    3. **SynthesisPipeline Integration**: Cleaner molecule synthesis workflow
    4. **Polars DataFrames**: Better performance for large result sets
    5. **Strategy Pattern**: Easy to add new selection strategies
    6. **Warmup Flexibility**: Choose warmup strategy independently of selection strategy
    """)
    return


if __name__ == "__main__":
    app.run()
