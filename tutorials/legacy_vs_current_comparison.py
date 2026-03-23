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

__generated_with = "0.19.7"
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
        TopTwoConfig,
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
        TopTwoConfig,
        mo,
        pl,
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
    | **Current Roulette Wheel** | RWS + GMIC criticality | Thermal cycling + criticality rotation |
    | **Current TT-TS** | Top-Two TS + thermal cycling | Two-sample disagreement + criticality rotation |

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
    return AMIDE_COUPLING_SMARTS, REAGENT_FILES, SCORES_FILE


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
    run_current_roulette = mo.ui.checkbox(value=True, label="Current Roulette Wheel (GMIC)")
    run_current_ttts = mo.ui.checkbox(value=True, label="Current TT-TS")
    return (
        run_current_epsilon,
        run_current_greedy,
        run_current_roulette,
        run_current_ttts,
        run_legacy_enhanced,
        run_legacy_standard,
    )


@app.cell
def _(
    mo,
    run_current_epsilon,
    run_current_greedy,
    run_current_roulette,
    run_current_ttts,
    run_legacy_enhanced,
    run_legacy_standard,
):
    """Display implementation selection."""
    mo.vstack([
        mo.md("**Legacy Implementations:**"),
        mo.hstack([run_legacy_standard, run_legacy_enhanced], justify="start", gap=2),
        mo.md("**Current Implementations:**"),
        mo.hstack([run_current_greedy, run_current_epsilon, run_current_roulette, run_current_ttts], justify="start", gap=2),
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Shared Parameters
    """)
    return


@app.cell
def _(mo):
    """Shared parameters for all methods."""
    warmup_k_slider = mo.ui.slider(
        start=1, stop=10, value=3, step=1,
        label="Warmup K (evaluations per reagent)"
    )
    cycles_slider = mo.ui.slider(
        start=1, stop=10, value=3, step=1,
        label="Benchmark cycles"
    )
    top_n_slider = mo.ui.slider(
        start=50, stop=500, value=100, step=50,
        label="Top N for recovery"
    )
    return cycles_slider, top_n_slider, warmup_k_slider


@app.cell
def _(cycles_slider, mo, top_n_slider, warmup_k_slider):
    """Display shared parameter controls."""
    mo.hstack([warmup_k_slider, cycles_slider, top_n_slider], justify="start", gap=2)
    return


@app.cell
def _(mo):
    """Legacy Standard search parameters."""
    legacy_iterations_slider = mo.ui.slider(
        start=100, stop=20000, value=500, step=100,
        label="Search iterations (1 compound/cycle)"
    )
    return (legacy_iterations_slider,)


@app.cell
def _(legacy_iterations_slider, mo, run_legacy_standard, warmup_k_slider):
    """Display Legacy Standard parameters (conditional)."""
    mo.stop(not run_legacy_standard.value)
    _output = mo.vstack([
        mo.md(
            f"### Legacy Standard (Greedy) Parameters\n"
            f"- **Warmup:** Standard (K={warmup_k_slider.value}) — each reagent sampled K times with random partners\n"
            f"- **Search:** {legacy_iterations_slider.value:,} iterations, 1 compound per cycle"
        ),
        legacy_iterations_slider,
    ])
    _output
    return


@app.cell
def _(mo):
    """Current TACTICS search parameters — widget definitions."""
    iterations_slider = mo.ui.slider(start=1, stop=500, value=5, step=1, label="Search iterations")
    batch_size_slider = mo.ui.slider(start=10, stop=500, value=100, step=10, label="Batch size (compounds per iteration)")
    warmup_dropdown = mo.ui.dropdown(
        options={"Standard (random)": "standard", "Balanced (stratified)": "balanced", "Enhanced (stochastic)": "enhanced"},
        value="Balanced (stratified)",
        label="Warmup Strategy",
    )
    return batch_size_slider, iterations_slider, warmup_dropdown


@app.cell
def _(
    batch_size_slider,
    iterations_slider,
    mo,
    run_current_epsilon,
    run_current_greedy,
    run_current_roulette,
    run_current_ttts,
    warmup_dropdown,
    warmup_k_slider,
):
    """Display Current TACTICS parameters (conditional)."""
    _any_current = (
        run_current_greedy.value or run_current_epsilon.value
        or run_current_roulette.value or run_current_ttts.value
    )
    mo.stop(not _any_current)
    _search_evals = iterations_slider.value * batch_size_slider.value
    _output = mo.vstack([
        mo.md(
            "### Current TACTICS Parameters\n"
            "- **Warmup:** " + str(warmup_dropdown.value) + " (K=" + str(warmup_k_slider.value) + ")\n"
            "- **Search:** " + str(iterations_slider.value) + " iterations x "
            + str(batch_size_slider.value) + " batch = **" + str(f"{_search_evals:,}") + "** search evaluations"
        ),
        warmup_dropdown,
        mo.hstack([iterations_slider, batch_size_slider]),
    ])
    _output
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### Iteration Budget Breakdown

    Shows how the total evaluation budget is distributed between warmup and search
    phases for each selected method, based on the current parameter settings.
    """)
    return


@app.cell
def _(
    REAGENT_FILES,
    batch_size_slider,
    iterations_slider,
    legacy_iterations_slider,
    legacy_percent_lib,
    mo,
    run_current_epsilon,
    run_current_greedy,
    run_current_roulette,
    run_current_ttts,
    run_legacy_enhanced,
    run_legacy_standard,
    warmup_dropdown,
    warmup_k_slider,
):
    """Compute and display iteration budget breakdown per method."""
    # Read reagent file sizes
    _reagent_sizes = []
    for _f in REAGENT_FILES:
        with open(_f) as _fh:
            _reagent_sizes.append(sum(1 for _line in _fh if _line.strip()))

    _sum_sizes = sum(_reagent_sizes)
    _max_size = max(_reagent_sizes)
    _library_size = 1
    for _s in _reagent_sizes:
        _library_size *= _s

    _K = warmup_k_slider.value
    _warmup_name = warmup_dropdown.value

    _rows = []

    # Legacy Standard: Standard warmup, 1 compound per search cycle
    if run_legacy_standard.value:
        _ls_warmup = _sum_sizes * _K
        _ls_search = legacy_iterations_slider.value
        _ls_total = _ls_warmup + _ls_search
        _ls_pct = _ls_total / _library_size * 100
        _rows.append(
            f"| Legacy Standard | Standard (K={_K}) | "
            f"{_sum_sizes} x {_K} = **{_ls_warmup:,}** | "
            f"**{_ls_search:,}** (1/cycle) | **{_ls_total:,}** | {_ls_pct:.2f}% |"
        )

    # Legacy Enhanced: Enhanced warmup, batch=100 per search cycle
    # percent_of_library is the TOTAL budget (warmup + search)
    # The legacy search method computes: nsearch = pct * library_size - num_warmup
    if run_legacy_enhanced.value:
        _le_warmup = _max_size * _K
        _le_total = int(_library_size * legacy_percent_lib.value)
        _le_search = max(_le_total - _le_warmup, 0)
        _le_pct = _le_total / _library_size * 100
        _le_warn = " **[warmup > budget!]**" if _le_search == 0 else ""
        _rows.append(
            f"| Legacy Enhanced | Enhanced (K={_K}) | "
            f"{_max_size} x {_K} = **{_le_warmup:,}** | "
            f"{_le_total:,} - {_le_warmup:,} = **{_le_search:,}** (batch=100){_le_warn} | "
            f"**{_le_total:,}** | {_le_pct:.2f}% |"
        )

    # Current TACTICS: configurable warmup, batch_size per iteration
    _tactics_search = iterations_slider.value * batch_size_slider.value
    if _warmup_name == "enhanced":
        _current_warmup = _max_size * _K
        _current_warmup_label = f"Enhanced (K={_K})"
    else:  # standard or balanced — both use sum(sizes) * K
        _current_warmup = _sum_sizes * _K
        _current_warmup_label = f"{_warmup_name.title()} (K={_K})"

    _current_methods = []
    if run_current_greedy.value:
        _current_methods.append("Current Greedy")
    if run_current_epsilon.value:
        _current_methods.append("Current Epsilon-Greedy")
    if run_current_roulette.value:
        _current_methods.append("Current Roulette Wheel")
    if run_current_ttts.value:
        _current_methods.append("Current TT-TS")

    for _method in _current_methods:
        _c_total = _current_warmup + _tactics_search
        _c_pct = _c_total / _library_size * 100
        _rows.append(
            f"| {_method} | {_current_warmup_label} | "
            f"**{_current_warmup:,}** | "
            f"{iterations_slider.value} x {batch_size_slider.value} = **{_tactics_search:,}** | "
            f"**{_c_total:,}** | {_c_pct:.2f}% |"
        )

    _sizes_str = " + ".join(str(s) for s in _reagent_sizes)
    _table = "\n".join(_rows) if _rows else ""

    _output = mo.md(
        f"**Library:** {_sizes_str} = **{_library_size:,}** products "
        f"({len(_reagent_sizes)} components: {', '.join(f'{s:,}' for s in _reagent_sizes)})\n\n"
        f"| Method | Warmup Strategy | Warmup Evaluations | Search Evaluations | Total Evaluations | % Library |\n"
        f"|--------|----------------|-------------------|-------------------|-------------------|-----------|\n"
        f"{_table}\n\n"
        f"**How evaluations are calculated:**\n\n"
        f"| Phase | Method | Formula | This library |\n"
        f"|-------|--------|---------|-------------|\n"
        f"| Warmup | Standard / Balanced | `sum(component_sizes) x K` | {_sum_sizes} x K |\n"
        f"| Warmup | Enhanced | `max(component_sizes) x K` | {_max_size} x K |\n"
        f"| Search | Legacy Standard | `iterations` (1 compound/cycle) | iterations |\n"
        f"| Search | Legacy Enhanced | `pct x library_size - warmup` (batch=100/cycle) | total budget - warmup |\n"
        f"| Search | TACTICS | `iterations x batch_size` | iterations x batch_size |"
    ) if _rows else mo.md("*Select at least one implementation above to see the iteration budget.*")
    _output
    return


@app.cell
def _(mo):
    """Legacy Enhanced specific parameters."""
    legacy_percent_lib = mo.ui.slider(
        start=0.01, stop=0.10, value=0.05, step=0.01,
        label="Percent of library to screen (total budget incl. warmup)"
    )
    legacy_scaling = mo.ui.slider(
        start=-1, stop=1, value=-1, step=1,
        label="Scaling (-1=minimize, 1=maximize)"
    )
    return legacy_percent_lib, legacy_scaling


@app.cell
def _(
    legacy_percent_lib,
    legacy_scaling,
    mo,
    run_legacy_enhanced,
    warmup_k_slider,
):
    """Display Legacy Enhanced parameters (conditional)."""
    mo.stop(not run_legacy_enhanced.value)
    _output = mo.vstack([
        mo.md(
            f"### Legacy Enhanced (CATS/RWS) Parameters\n"
            f"- **Warmup:** Enhanced (K={warmup_k_slider.value}) — stochastic parallel pairing\n"
            f"- **Search:** {legacy_percent_lib.value*100:.0f}% of library (total budget), batch size = 100 per cycle"
        ),
        mo.hstack([legacy_percent_lib, legacy_scaling]),
    ])
    _output
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
    TopTwoConfig,
    batch_size_slider,
    cycles_slider,
    iterations_slider,
    legacy_iterations_slider,
    legacy_percent_lib,
    legacy_scaling,
    mo,
    pl,
    run_button,
    run_current_epsilon,
    run_current_greedy,
    run_current_roulette,
    run_current_ttts,
    run_legacy_enhanced,
    run_legacy_standard,
    tempfile,
    time,
    warmup_dropdown,
    warmup_k_slider,
):
    """Execute benchmark for all selected implementations."""
    mo.stop(not run_button.value, mo.md("*Click 'Run Legacy vs Current Benchmark' to start*"))

    # Check at least one is selected
    any_selected = (
        run_legacy_standard.value or
        run_legacy_enhanced.value or
        run_current_greedy.value or
        run_current_epsilon.value or
        run_current_roulette.value or
        run_current_ttts.value
    )
    if not any_selected:
        mo.stop(True, mo.md("**Error:** Select at least one implementation."))

    # Parameters
    _num_cycles = cycles_slider.value
    _K = warmup_k_slider.value

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
            _warmup_results = _legacy_sampler.warm_up(num_warmup_trials=_K)
            _search_results = _legacy_sampler.search(num_cycles=legacy_iterations_slider.value)

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
            "warmup": f"Standard (K={_K})",
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
            _legacy_enhanced.warm_up(num_warmup_trials=_K)
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
                "warmup": f"Enhanced (K={_K})",
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
        _warmup_display = f"Standard (K={_K})"
    elif _warmup_name == "balanced":
        _warmup_config = BalancedWarmupConfig(observations_per_reagent=_K)
        _warmup_display = f"Balanced (K={_K})"
    else:  # enhanced
        _warmup_config = EnhancedWarmupConfig()
        _warmup_display = f"Enhanced (K={_K})"

    # Helper function for current implementations
    def run_current_implementation(method_name, strategy_config):
        _cycle_results = []
        for _cycle in range(_num_cycles):
            _config = ThompsonSamplingConfig(
                synthesis_pipeline=_pipeline,
                num_ts_iterations=iterations_slider.value,
                num_warmup_trials=_K,
                strategy_config=strategy_config,
                warmup_config=_warmup_config,
                evaluator_config=_evaluator_config,
                batch_size=batch_size_slider.value,
                max_resamples=1000,
                hide_progress=True,
                use_boltzmann_weighting=True,
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

    # Current Roulette Wheel (GMIC)
    if run_current_roulette.value:
        run_current_implementation(
            "Current Roulette Wheel",
            RouletteWheelConfig(mode="minimize", alpha=0.1, beta=0.05, divergence_threshold=0.1)
        )

    # Current TT-TS
    if run_current_ttts.value:
        run_current_implementation(
            "Current TT-TS",
            TopTwoConfig(mode="minimize", beta=0.5, heated_scale=1.5, cooled_scale=0.75)
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


@app.cell
def _(mo, run_button):
    """Save benchmark results to CSV."""
    mo.stop(not run_button.value)

    import pathlib as _pathlib
    import datetime as _dt

    _timestamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    _default_dir = str(_pathlib.Path.cwd().parent / "outputs" / "benchmark_results")

    save_dir_input = mo.ui.text(
        value=_default_dir,
        label="Output directory",
        full_width=True,
    )
    save_button = mo.ui.run_button(label="Save Results to CSV")

    mo.vstack([
        mo.md("### Save Results"),
        save_dir_input,
        save_button,
    ])
    return save_button, save_dir_input


@app.cell
def _(
    all_results,
    method_names,
    mo,
    pl,
    run_button,
    run_metadata,
    save_button,
    save_dir_input,
):
    """Write CSV files when save button is clicked."""
    mo.stop(not run_button.value)
    mo.stop(not save_button.value, mo.md("*Click 'Save Results to CSV' to export.*"))

    import pathlib as _pathlib
    import datetime as _dt

    _timestamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    _out_dir = _pathlib.Path(save_dir_input.value) / _timestamp
    _out_dir.mkdir(parents=True, exist_ok=True)

    _saved = []

    # 1. Run metadata summary
    _meta_df = pl.DataFrame(run_metadata)
    _meta_path = _out_dir / "run_summary.csv"
    _meta_df.write_csv(_meta_path)
    _saved.append(f"`run_summary.csv` — {len(run_metadata)} methods")

    # 2. Per-method, per-cycle compound results
    for _name in method_names:
        _cycle_results = all_results[_name]
        _safe_name = _name.lower().replace(" ", "_").replace("-", "_")
        for _i, _df in enumerate(_cycle_results):
            if len(_df) > 0:
                _fname = f"{_safe_name}_cycle{_i + 1}.csv"
                _df.write_csv(_out_dir / _fname)
        _saved.append(f"`{_safe_name}_cycle*.csv` — {len(_cycle_results)} cycles")

    # 3. Combined results across all methods and cycles
    _all_dfs = []
    for _name in method_names:
        for _i, _df in enumerate(all_results[_name]):
            if len(_df) > 0:
                _all_dfs.append(
                    _df.with_columns(
                        pl.lit(_name).alias("method"),
                        pl.lit(_i + 1).alias("cycle"),
                    )
                )
    if _all_dfs:
        _combined = pl.concat(_all_dfs)
        _combined.write_csv(_out_dir / "all_results_combined.csv")
        _saved.append(f"`all_results_combined.csv` — {len(_combined):,} rows")

    _file_list = "\n".join(f"- {s}" for s in _saved)
    mo.md(f"""
    **Saved to** `{_out_dir}`

    {_file_list}
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
def _(mo, run_button, viz_toggles):
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
    | `EnhancedThompsonSampler` | `ThompsonSampler` + `RouletteWheelConfig` | GMIC criticality + thermal cycling |
    | — | `ThompsonSampler` + `TopTwoConfig` | TT-TS best-arm ID + thermal cycling |

    ### Warmup Strategy Comparison

    | Warmup | Description | When to Use |
    |--------|-------------|-------------|
    | **Standard** | Random partner selection | Quick baseline |
    | **Balanced** | K observations per reagent | Best coverage (recommended) |
    | **Enhanced** | Stochastic parallel pairing | Large libraries |

    ### Expected Results

    - **Legacy Standard** and **Current Greedy** should produce statistically similar results
    - **Current Epsilon-Greedy** typically shows better exploration (more diverse hits)
    - **Current Roulette Wheel** uses GMIC criticality-weighted rotation (improved over Legacy Enhanced)
    - **Current TT-TS** targets best-arm identification with thermal cycling — best recovery on manuscript benchmarks

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
