"""
TACTICS: Thrombin Recovery Benchmark

An interactive benchmark comparing TACTICS methods vs Legacy approaches on the
Thrombin docking dataset (minimize mode). Settings match
`examples/run_recovery_benchmark.py` and `examples/_recovery_worker.py`.

Methods compared:
1. Legacy TS -- greedy Thompson Sampling (random warmup, no CATS)
2. Legacy RWS -- Roulette Wheel Selection (Boltzmann-weighted, Enhanced warmup)
3. TACTICS Greedy -- modern greedy with Balanced warmup
4. TACTICS RWS (CATS) -- full CATS + adaptive temperature, batch mode

Dataset: Thrombin Linear Amide Library (~500K products, docking scores, minimize)

Run as app: marimo run tutorials/thrombin_recovery_benchmark.py
Edit mode:  marimo edit tutorials/thrombin_recovery_benchmark.py
"""

import marimo

__generated_with = "0.19.11"
app = marimo.App(
    width="full",
    app_title="TACTICS: Thrombin Recovery Benchmark",
)


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
    import numpy as np
    import time
    import matplotlib
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy import stats

    matplotlib.rcParams.update({
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "figure.dpi": 150,
    })

    # TACTICS Thompson Sampling imports
    from TACTICS.thompson_sampling import ThompsonSampler
    from TACTICS.thompson_sampling.config import ThompsonSamplingConfig
    from TACTICS.thompson_sampling.strategies.config import (
        GreedyConfig,
        RouletteWheelConfig,
    )
    from TACTICS.thompson_sampling.warmup.config import (
        StandardWarmupConfig,
        EnhancedWarmupConfig,
        BalancedWarmupConfig,
    )
    from TACTICS.thompson_sampling.core.evaluator_config import LookupEvaluatorConfig

    # TACTICS Library Enumeration imports
    from TACTICS.library_enumeration import SynthesisPipeline
    from TACTICS.library_enumeration.smarts_toolkit import (
        ReactionConfig,
        ReactionDef,
    )

    # TACTICS Library Analysis imports
    from TACTICS.library_analysis.visualization import TS_Benchmarks

    return (
        BalancedWarmupConfig,
        EnhancedWarmupConfig,
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
        alt,
        mo,
        np,
        pl,
        plt,
        sns,
        stats,
        time,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Thrombin Recovery Benchmark: TACTICS vs Legacy

    This notebook benchmarks **modern TACTICS strategies** against **legacy Thompson Sampling**
    on the Thrombin docking dataset (~500K products, minimize mode).

    **Methods compared:**
    - **Legacy TS** -- greedy selection with random warmup (no CATS)
    - **Legacy RWS** -- Roulette Wheel Selection with Boltzmann weighting, Enhanced warmup (no CATS)
    - **TACTICS Greedy** -- greedy selection with Balanced warmup
    - **TACTICS RWS (CATS)** -- Roulette Wheel with CATS + adaptive temperature, batch mode

    **Recovery metric:** fraction of the ground-truth top-N compounds found by each method.
    Settings match `examples/run_recovery_benchmark.py`.
    """)
    return


@app.cell
def _():
    """Load bundled Thrombin dataset paths."""
    import importlib.resources

    _data_files = importlib.resources.files("TACTICS.data.thrombin")
    ACIDS_FILE = str(_data_files / "acids.smi")
    DIPEPTIDES_FILE = str(_data_files / "coupled_aa_sub.smi")
    SCORES_FILE = str(_data_files / "product_scores.csv")
    return ACIDS_FILE, DIPEPTIDES_FILE, SCORES_FILE


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Configuration

    Adjust benchmark parameters below. Default settings match
    `run_recovery_benchmark.py`: 5000 iterations (~1% of library), 5 replicates, top-100 recovery.
    """)
    return


@app.cell
def _(mo):
    """Configuration UI elements."""
    replicates_slider = mo.ui.slider(
        start=3, stop=10, value=5, step=1,
        label="Replicates",
    )
    iterations_slider = mo.ui.slider(
        start=1000, stop=10000, value=5000, step=1000,
        label="Total iterations",
    )
    top_n_slider = mo.ui.slider(
        start=50, stop=500, value=100, step=50,
        label="Top N for recovery",
    )
    return iterations_slider, replicates_slider, top_n_slider


@app.cell(hide_code=True)
def _(iterations_slider, mo, replicates_slider, top_n_slider):
    """Display configuration controls."""
    mo.vstack([
        mo.md("### Benchmark Parameters"),
        mo.hstack(
            [replicates_slider, iterations_slider, top_n_slider],
            justify="start", gap=2,
        ),
    ])
    return


@app.cell
def _(
    ACIDS_FILE,
    DIPEPTIDES_FILE,
    ReactionConfig,
    ReactionDef,
    SynthesisPipeline,
):
    """Create 2-component synthesis pipeline (Thrombin amide coupling)."""
    AMIDE_COUPLING_SMARTS = (
        "[#6:1](=[O:2])[OH]."
        "[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]"
        ">>[#6:1](=[O:2])[#7:3]"
    )
    _reaction_config = ReactionConfig(
        reactions=[ReactionDef(reaction_smarts=AMIDE_COUPLING_SMARTS, step_index=0)],
        reagent_file_list=[ACIDS_FILE, DIPEPTIDES_FILE],
    )
    pipeline = SynthesisPipeline(_reaction_config)
    print(f"Pipeline: {pipeline.num_components} components, ~{pipeline.library_size:,} products")
    return (pipeline,)


@app.cell
def _(SCORES_FILE, pl):
    """Load reference data (ground truth scores)."""
    reference_data = (
        pl.read_csv(SCORES_FILE)
        .rename({"Product_Code": "Name", "Scores": "score"})
        .select(["Name", "score"])
    )
    print(f"Reference data: {len(reference_data):,} products")
    return (reference_data,)


@app.cell
def _():
    """Method definitions and color constants."""
    METHOD_ORDER = [
        "Legacy TS",
        "Legacy RWS",
        "TACTICS Greedy",
        "TACTICS RWS (CATS)",
    ]
    METHOD_COLORS = {
        "Legacy TS": "#8c8c8c",
        "Legacy RWS": "#bfbfbf",
        "TACTICS Greedy": "#2196F3",
        "TACTICS RWS (CATS)": "#E91E63",
    }
    METHOD_PALETTE = [METHOD_COLORS[m] for m in METHOD_ORDER]
    return METHOD_ORDER, METHOD_PALETTE


@app.cell
def _(
    BalancedWarmupConfig,
    EnhancedWarmupConfig,
    GreedyConfig,
    LookupEvaluatorConfig,
    RouletteWheelConfig,
    SCORES_FILE,
    StandardWarmupConfig,
    ThompsonSampler,
    ThompsonSamplingConfig,
    pl,
    time,
):
    """Define the benchmark run function."""
    def run_benchmark(pipeline, n_replicates, n_iters, top_n):
        """
        Run all 4 methods for the specified number of replicates.

        Returns:
            ts_runs_data: dict of method -> list of DataFrames (one per replicate)
            metadata: dict of method -> timing/score info
        """
        evaluator_config = LookupEvaluatorConfig(
            ref_filename=SCORES_FILE,
            compound_col="Product_Code",
            score_col="Scores",
        )

        methods = {
            "Legacy TS": {
                "strategy_config": GreedyConfig(mode="minimize"),
                "warmup_config": StandardWarmupConfig(),
                "use_boltzmann": False,
                "batch_size": 1,
            },
            "Legacy RWS": {
                "strategy_config": RouletteWheelConfig(
                    mode="minimize", alpha=0.1, beta=0.05,
                ),
                "warmup_config": EnhancedWarmupConfig(),
                "use_boltzmann": True,
                "batch_size": 1,
            },
            "TACTICS Greedy": {
                "strategy_config": GreedyConfig(mode="minimize"),
                "warmup_config": BalancedWarmupConfig(observations_per_reagent=5),
                "use_boltzmann": False,
                "batch_size": 1,
            },
            "TACTICS RWS (CATS)": {
                "strategy_config": RouletteWheelConfig(
                    mode="minimize", alpha=0.1, beta=0.05,
                    adaptive_temperature=True,
                ),
                "warmup_config": BalancedWarmupConfig(observations_per_reagent=5),
                "use_boltzmann": True,
                "batch_size": 100,
            },
        }

        ts_runs_data = {name: [] for name in methods}
        metadata = {}

        for method_name, method_cfg in methods.items():
            _start = time.time()
            _batch_size = method_cfg["batch_size"]
            _n_cycles = max(1, n_iters // _batch_size)

            for _rep in range(n_replicates):
                _config = ThompsonSamplingConfig(
                    synthesis_pipeline=pipeline,
                    num_ts_iterations=_n_cycles,
                    num_warmup_trials=5,
                    strategy_config=method_cfg["strategy_config"],
                    warmup_config=method_cfg["warmup_config"],
                    evaluator_config=evaluator_config,
                    batch_size=_batch_size,
                    use_boltzmann_weighting=method_cfg["use_boltzmann"],
                    max_resamples=1000,
                    hide_progress=True,
                )
                _sampler = ThompsonSampler.from_config(_config)
                _warmup_df = _sampler.warm_up(num_warmup_trials=_config.num_warmup_trials)
                _search_df = _sampler.search(num_cycles=_config.num_ts_iterations)
                _combined = pl.concat([_warmup_df, _search_df])
                _result = _combined.select(["Name", "score"])
                ts_runs_data[method_name].append(_result)
                _sampler.close()

            _elapsed = time.time() - _start
            _best_scores = [r["score"].min() for r in ts_runs_data[method_name]]
            metadata[method_name] = {
                "time": _elapsed,
                "best": min(_best_scores),
                "avg_best": sum(float(s) for s in _best_scores) / len(_best_scores),
                "n_evaluated": [len(r) for r in ts_runs_data[method_name]],
            }
            print(f"  {method_name}: {_elapsed:.1f}s, best={min(_best_scores):.4f}")

        return ts_runs_data, metadata

    return (run_benchmark,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Run Benchmark

    Click the button below to run all 4 methods. Each method evaluates ~5000 compounds
    per replicate for a fair comparison.
    """)
    return


@app.cell
def _(mo):
    """Run button."""
    run_button = mo.ui.run_button(label="Run Benchmark")
    run_button
    return (run_button,)


@app.cell
def _(
    iterations_slider,
    mo,
    pipeline,
    replicates_slider,
    run_benchmark,
    run_button,
    time,
    top_n_slider,
):
    """Execute the benchmark."""
    mo.stop(not run_button.value, mo.md("*Click 'Run Benchmark' to start*"))

    _start = time.time()
    print(f"Running benchmark: {replicates_slider.value} replicates, "
          f"{iterations_slider.value} iterations, top-{top_n_slider.value}")

    ts_runs_data, bench_metadata = run_benchmark(
        pipeline=pipeline,
        n_replicates=replicates_slider.value,
        n_iters=iterations_slider.value,
        top_n=top_n_slider.value,
    )

    _total_time = time.time() - _start
    print(f"\nTotal benchmark time: {_total_time:.1f}s")
    return bench_metadata, ts_runs_data


@app.cell(hide_code=True)
def _(METHOD_ORDER, bench_metadata, mo, run_button):
    """Display results summary table."""
    mo.stop(not run_button.value)

    _rows = []
    for _m in METHOD_ORDER:
        _meta = bench_metadata[_m]
        _avg_evals = sum(_meta["n_evaluated"]) / len(_meta["n_evaluated"])
        _rows.append(
            f"| {_m} | {_meta['avg_best']:.4f} | {_meta['best']:.4f} | "
            f"{_avg_evals:.0f} | {_meta['time']:.1f}s |"
        )
    _table = "\n".join(_rows)

    mo.md(f"""
    ### Results Summary

    | Method | Avg Best Score | Best Score | Avg Evaluations | Time |
    |--------|---------------|------------|-----------------|------|
    {_table}
    """)
    return


@app.cell
def _(
    TS_Benchmarks,
    mo,
    reference_data,
    replicates_slider,
    run_button,
    top_n_slider,
    ts_runs_data,
):
    """Create TS_Benchmarks instance."""
    mo.stop(not run_button.value)

    benchmarks = TS_Benchmarks(
        no_of_cycles=replicates_slider.value,
        methods_list=list(ts_runs_data.keys()),
        TS_runs_data=ts_runs_data,
        reference_data=reference_data,
        top_n=top_n_slider.value,
        sort_type="minimize",
        top_ns=[25, 50, 100, 200, 300],
    )
    return (benchmarks,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ### Per-Replicate Hit Recovery (Bar Plot)

    How many ground-truth top-N compounds each method recovered in each replicate.
    """)
    return


@app.cell
def _(benchmarks, mo, run_button):
    """Per-replicate barplot via TS_Benchmarks."""
    mo.stop(not run_button.value)

    benchmarks.plot_barplot_TS_results(
        width=900, height=400, show_plot=True,
        legend_position="bottom", dark_mode=False,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Recovery Fraction vs Top-N (Line Plot)

    Fraction of top-N ground-truth compounds found, plotted across different N values
    with error bars (mean +/- std across replicates).
    """)
    return


@app.cell
def _(benchmarks, mo, run_button):
    """Line plot via TS_Benchmarks."""
    mo.stop(not run_button.value)

    benchmarks.plot_line_performance_with_error_bars(
        width=900, height=400, show_plot=True, legend_position="bottom",
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ### Recovery Computation

    Compute per-replicate recovery counts and fractions for the summary bar and violin plots.
    """)
    return


@app.cell
def _(
    METHOD_ORDER,
    mo,
    pl,
    reference_data,
    run_button,
    top_n_slider,
    ts_runs_data,
):
    """Compute recovery data for summary plots."""
    mo.stop(not run_button.value)

    _top_n = top_n_slider.value
    _ref_top = set(
        reference_data.sort("score").head(_top_n)["Name"].to_list()
    )

    _rows = []
    for _method in METHOD_ORDER:
        for _rep, _df in enumerate(ts_runs_data[_method]):
            _found = set(_df.sort("score").head(_top_n)["Name"].to_list())
            _recovered = len(_found & _ref_top)
            _rows.append({
                "method": _method,
                "replicate": _rep,
                "recovered": _recovered,
                "recovery_frac": _recovered / _top_n,
            })

    recovery_data_df = pl.DataFrame(_rows)

    # Summary stats
    _summary = (
        recovery_data_df
        .group_by("method")
        .agg([
            pl.col("recovered").mean().round(1).alias("mean"),
            pl.col("recovered").std().round(1).alias("std"),
        ])
    )
    print("Recovery summary:")
    for _row in _summary.iter_rows(named=True):
        print(f"  {_row['method']:25s}  {_row['mean']:.1f} +/- {_row['std']:.1f}")
    return (recovery_data_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Summary Bar Plot

    Mean recovery per method with error bars (standard deviation across replicates).
    """)
    return


@app.cell
def _(METHOD_ORDER, METHOD_PALETTE, alt, mo, pl, recovery_data_df, run_button):
    """Summary barplot (Altair)."""
    mo.stop(not run_button.value)

    _summary = (
        recovery_data_df
        .group_by("method")
        .agg([
            pl.col("recovered").mean().alias("mean_recovery"),
            pl.col("recovered").std().alias("std_recovery"),
        ])
    )

    _bars = (
        alt.Chart(_summary)
        .mark_bar(cornerRadiusTopLeft=3, cornerRadiusTopRight=3)
        .encode(
            x=alt.X("method:N", title="Method",
                     sort=METHOD_ORDER,
                     axis=alt.Axis(labelAngle=-20)),
            y=alt.Y("mean_recovery:Q", title="Mean Recovery (compounds)"),
            color=alt.Color("method:N", title="Method",
                            scale=alt.Scale(domain=METHOD_ORDER, range=METHOD_PALETTE),
                            legend=None),
            tooltip=["method", "mean_recovery", "std_recovery"],
        )
    )

    _error = (
        alt.Chart(_summary)
        .mark_errorbar(extent="stdev")
        .encode(
            x=alt.X("method:N", sort=METHOD_ORDER),
            y=alt.Y("mean_recovery:Q"),
            yError=alt.YError("std_recovery:Q"),
        )
    )

    _text = (
        alt.Chart(_summary)
        .mark_text(dy=-15, fontSize=11, fontWeight="bold")
        .encode(
            x=alt.X("method:N", sort=METHOD_ORDER),
            y=alt.Y("mean_recovery:Q"),
            text=alt.Text("mean_recovery:Q", format=".1f"),
        )
    )

    summary_chart = (
        (_bars + _error + _text)
        .properties(width=500, height=350,
                    title="Mean Recovery by Method")
    )
    summary_chart
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Violin Plot with Paired t-Test Significance

    Distribution of recovery counts per method with significance brackets from paired t-tests.
    Bracket colors: **blue** = TACTICS wins, **red** = Legacy wins, **grey** = not significant.
    """)
    return


@app.function
def draw_significance_bracket(ax, x1, x2, y, h, text, color):
    """Draw a significance bracket between two x positions."""
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.2, color=color)
    ax.text(
        (x1 + x2) / 2, y + h, text,
        ha="center", va="bottom", fontsize=8, color=color, fontweight="bold",
    )


@app.cell
def _(
    METHOD_ORDER,
    METHOD_PALETTE,
    mo,
    np,
    pl,
    plt,
    recovery_data_df,
    run_button,
    sns,
    stats,
):
    """Violin plot with paired t-test significance brackets."""
    mo.stop(not run_button.value)

    # Paired t-test comparisons (matching run_recovery_benchmark.py)
    _method_pairs = [
        ("Legacy TS", "TACTICS Greedy"),
        ("Legacy RWS", "TACTICS RWS (CATS)"),
        ("Legacy TS", "TACTICS RWS (CATS)"),
        ("TACTICS Greedy", "TACTICS RWS (CATS)"),
        ("Legacy TS", "Legacy RWS"),
    ]

    # Compute paired t-tests
    _test_results = []
    for _ma, _mb in _method_pairs:
        _a_data = (
            recovery_data_df
            .filter(pl.col("method") == _ma)
            .sort("replicate")["recovered"]
            .to_numpy()
        )
        _b_data = (
            recovery_data_df
            .filter(pl.col("method") == _mb)
            .sort("replicate")["recovered"]
            .to_numpy()
        )
        _min_len = min(len(_a_data), len(_b_data))
        if _min_len < 2:
            continue
        _a_data = _a_data[:_min_len]
        _b_data = _b_data[:_min_len]
        _t, _p = stats.ttest_rel(_a_data, _b_data)
        _sig = (
            "****" if _p <= 0.0001
            else "***" if _p <= 0.001
            else "**" if _p <= 0.01
            else "*" if _p <= 0.05
            else "ns"
        )
        _test_results.append({
            "method_a": _ma,
            "method_b": _mb,
            "mean_a": float(np.mean(_a_data)),
            "mean_b": float(np.mean(_b_data)),
            "t_stat": float(_t),
            "p_value": float(_p),
            "significance": _sig,
        })

    # Plot
    _plot_data = recovery_data_df.to_pandas()

    _fig, _ax = plt.subplots(figsize=(10, 6))
    sns.violinplot(
        data=_plot_data, x="method", y="recovered", hue="method",
        order=METHOD_ORDER, hue_order=METHOD_ORDER,
        palette=METHOD_PALETTE, legend=False,
        inner="box", linewidth=0.8, saturation=0.85, cut=0,
        ax=_ax,
    )
    sns.stripplot(
        data=_plot_data, x="method", y="recovered",
        order=METHOD_ORDER,
        color="black", size=4, alpha=0.6, jitter=True,
        ax=_ax,
    )

    _ax.set_xlabel("")
    _ax.set_ylabel("Recovered Compounds")
    _ax.set_title("Recovery by Method (with paired t-test significance)", fontweight="bold")

    _short_labels = ["Legacy\nTS", "Legacy\nRWS", "TACTICS\nGreedy", "TACTICS\nRWS (CATS)"]
    _ax.set_xticks(range(len(METHOD_ORDER)))
    _ax.set_xticklabels(_short_labels, fontsize=9)

    _method_idx = {m: i for i, m in enumerate(METHOD_ORDER)}

    # Cross-category pairs where TACTICS (method_b) is compared against Legacy (method_a)
    _cross_pairs = {
        ("Legacy TS", "TACTICS Greedy"),
        ("Legacy RWS", "TACTICS RWS (CATS)"),
        ("Legacy TS", "TACTICS RWS (CATS)"),
    }

    # Draw brackets ordered by span width (narrow first)
    _bracket_order = [
        ("TACTICS Greedy", "TACTICS RWS (CATS)"),
        ("Legacy TS", "Legacy RWS"),
        ("Legacy TS", "TACTICS Greedy"),
        ("Legacy RWS", "TACTICS RWS (CATS)"),
        ("Legacy TS", "TACTICS RWS (CATS)"),
    ]

    _y_max = _plot_data["recovered"].max()
    _y_base = _y_max + 2
    _bracket_h = 0.8
    _bracket_spacing = 3.0

    for _j, (_ma, _mb) in enumerate(_bracket_order):
        _match = [t for t in _test_results if t["method_a"] == _ma and t["method_b"] == _mb]
        if not _match:
            continue
        _t_res = _match[0]
        _sig = _t_res["significance"]
        _x1 = _method_idx[_ma]
        _x2 = _method_idx[_mb]
        _y = _y_base + _j * _bracket_spacing

        if _sig == "ns":
            _color = "#999999"
            _text = "ns"
        elif (_ma, _mb) in _cross_pairs:
            _color = "#1565C0" if _t_res["mean_b"] > _t_res["mean_a"] else "#D32F2F"
            _text = _sig
        else:
            _color = "#FF9800"
            _text = _sig

        draw_significance_bracket(_ax, _x1, _x2, _y, _bracket_h, _text, _color)

    _ax.set_ylim(top=_y_base + len(_bracket_order) * _bracket_spacing + 2)
    _fig.tight_layout()

    # Print test results
    print("Paired t-test results:")
    for _t_res in _test_results:
        print(f"  {_t_res['method_a']:25s} vs {_t_res['method_b']:25s}  "
              f"p={_t_res['p_value']:.4f}  {_t_res['significance']}")

    plt.gca()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Conclusions

    ### Recovery Metric
    - **Ground truth**: top-N products by docking score from the exhaustive score matrix
    - **Recovery**: number of ground-truth products found in the top-N of each method's results

    ### Methods
    - **Legacy TS**: greedy selection, random warmup (StandardWarmup), no CATS
    - **Legacy RWS**: Roulette Wheel with thermal cycling (alpha=0.1, beta=0.05),
      Boltzmann-weighted updates, Enhanced warmup, no CATS
    - **TACTICS Greedy**: greedy selection with Balanced warmup (K=5 observations/reagent)
    - **TACTICS RWS (CATS)**: full Component-Aware Thompson Sampling with adaptive temperature,
      Balanced warmup, Boltzmann weighting, batch_size=100

    ### Statistical Testing
    - **Paired t-test** (`scipy.stats.ttest_rel`) -- replicate i of method A paired with replicate i of method B
    - Significance: \*\*\*\* (p<=0.0001), \*\*\* (p<=0.001), \*\* (p<=0.01), \* (p<=0.05), ns

    ### References
    - Zhao, H., Nittinger, E. & Tyrchan, C. Enhanced Thompson Sampling by Roulette
      Wheel Selection for Screening Ultra-Large Combinatorial Libraries.
      *J. Cheminform.* (2025)
    """)
    return


if __name__ == "__main__":
    app.run()
