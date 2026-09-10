"""
CATS Diagnostics Analysis: Online Heuristics vs Brute-Force Ground Truth

Loads pre-computed outputs from `examples/run_diagnostics_benchmark.py` and
visualizes how well CATS's real-time entropy/criticality tracking agrees with
exhaustive ground truth from the full combinatorial space.

Three data files are expected in `outputs/diagnostics_benchmark/`:
- sar_summaries.parquet       — per (lib, query, component) CATS vs GT metrics
- diagnostics_trajectories.parquet  — per-cycle enhanced diagnostics
- convergence_dynamics.parquet      — convergence timing per component

Run as app:  marimo run tutorials/cats_diagnostics_analysis.py
Edit mode:   marimo edit tutorials/cats_diagnostics_analysis.py
"""

import marimo

__generated_with = "0.19.7"
app = marimo.App(
    width="full",
    app_title="CATS Diagnostics: Heuristics vs Ground Truth",
)


@app.cell
def _():
    """Imports and project setup."""
    import marimo as mo
    import sys
    from pathlib import Path

    try:
        project_root = Path(__file__).parent.parent.resolve()
    except NameError:
        project_root = Path("/Users/aakankschitnandkeolyar/Desktop/TACTICS")

    sys.path.insert(0, str(project_root))
    sys.path.insert(0, str(project_root / "src"))

    import polars as pl
    import altair as alt
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    from scipy.stats import pearsonr, spearmanr

    matplotlib.rcParams.update({
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "figure.dpi": 150,
    })

    from TACTICS.library_analysis.diagnostic_plots import (
        plot_criticality_trajectory,
        plot_snr_trajectory,
        plot_temperature_decomposition,
    )

    OUTPUT_DIR = project_root / "outputs" / "diagnostics_benchmark"
    return (
        OUTPUT_DIR,
        alt,
        mo,
        np,
        pearsonr,
        pl,
        plot_criticality_trajectory,
        plot_snr_trajectory,
        plot_temperature_decomposition,
        plt,
        spearmanr,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # CATS Diagnostics: Online Heuristics vs Brute-Force Ground Truth

    This notebook validates whether CATS's **real-time entropy and criticality tracking**
    agrees with what we would compute from the **full combinatorial score matrix** (brute-force ground truth).

    **Key question:** Does CATS's online `concentration` metric correlate with
    the exhaustive `dominance` computed over all products?

    Data is pre-computed by `examples/run_diagnostics_benchmark.py`.
    """)
    return


@app.cell
def _(OUTPUT_DIR, mo, pl):
    """Load pre-computed benchmark outputs."""
    _sar_path = OUTPUT_DIR / "sar_summaries.parquet"
    _diag_path = OUTPUT_DIR / "diagnostics_trajectories.parquet"
    _conv_path = OUTPUT_DIR / "convergence_dynamics.parquet"

    mo.stop(
        not _sar_path.exists(),
        mo.md(
            "**Data not found.** Run `python examples/run_diagnostics_benchmark.py` first.\n\n"
            f"Expected: `{_sar_path}`"
        ),
    )

    sar_df = pl.read_parquet(_sar_path)
    diag_df = pl.read_parquet(_diag_path)
    conv_df = pl.read_parquet(_conv_path)

    print(f"SAR summaries:   {len(sar_df):,} rows")
    print(f"Diagnostics:     {len(diag_df):,} rows")
    print(f"Convergence:     {len(conv_df):,} rows")
    return conv_df, diag_df, sar_df


@app.cell
def _(mo, sar_df):
    """Library and query selectors."""
    _lib_ids = sorted(sar_df["library_id"].unique().to_list())
    library_dropdown = mo.ui.dropdown(
        options=_lib_ids,
        value=_lib_ids[0] if _lib_ids else None,
        label="Library",
    )
    return (library_dropdown,)


@app.cell
def _(library_dropdown, mo, pl, sar_df):
    """Query dropdown filtered by selected library."""
    _lib = library_dropdown.value
    _query_ids = sorted(
        sar_df.filter(pl.col("library_id") == _lib)["query_id"]
        .unique()
        .to_list()
    ) if _lib else []
    query_dropdown = mo.ui.dropdown(
        options=_query_ids,
        value=_query_ids[0] if _query_ids else None,
        label="Query",
    )
    return (query_dropdown,)


@app.cell
def _(library_dropdown, mo, pl, query_dropdown, sar_df):
    """Component slider based on selected library/query."""
    _lib = library_dropdown.value
    _query = query_dropdown.value
    _comps = sorted(
        sar_df.filter(
            (pl.col("library_id") == _lib) & (pl.col("query_id") == _query)
        )["component_idx"]
        .unique()
        .to_list()
    ) if _lib and _query else [0]
    component_slider = mo.ui.slider(
        start=min(_comps), stop=max(_comps), value=min(_comps), step=1,
        label="Component index",
    )
    return (component_slider,)


@app.cell(hide_code=True)
def _(component_slider, library_dropdown, mo, query_dropdown):
    mo.vstack([
        mo.md("---\n### Controls"),
        mo.hstack(
            [library_dropdown, query_dropdown, component_slider],
            justify="start", gap=2,
        ),
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Figure 1: CATS Concentration vs GT Dominance (Scatter)

    Does CATS's online estimate **rank** components the same way as brute-force ground truth?
    Each point is one (library, query, component).

    **Important:** These two metrics are on **different scales** by design:
    - **GT Dominance** = fraction of top-100 products sharing the most common reagent
      (range: 0–1, easily reaches 0.9+)
    - **CATS Concentration** = `1 − normalized_entropy` of the posterior after sampling ~2% of the library
      (range: 0–1, but stays low because the posterior hasn't fully converged)

    A positive correlation means CATS correctly identifies **which** components have structured SAR,
    even though the absolute values differ. The dashed line shows the linear regression fit.
    """)
    return


@app.cell
def _(alt, np, pearsonr, pl, sar_df, spearmanr):
    """Concentration vs Dominance scatter with regression line."""
    _valid = sar_df.filter(
        pl.col("cats_concentration").is_not_nan()
        & pl.col("gt_dominance").is_not_nan()
    )

    _cats = _valid["cats_concentration"].to_numpy()
    _gt = _valid["gt_dominance"].to_numpy()
    _r_pearson, _p_pearson = pearsonr(_cats, _gt) if len(_cats) > 2 else (float("nan"), float("nan"))
    _r_spearman, _p_spearman = spearmanr(_cats, _gt) if len(_cats) > 2 else (float("nan"), float("nan"))

    # Regression line fit (GT on x-axis, CATS on y-axis)
    _slope, _intercept = np.polyfit(_gt, _cats, 1) if len(_gt) > 2 else (0.0, 0.0)
    _x_range = np.array([float(_gt.min()), float(_gt.max())])
    _y_fit = _slope * _x_range + _intercept
    _reg_line = pl.DataFrame({"x": _x_range.tolist(), "y": _y_fit.tolist()})

    _scatter = (
        alt.Chart(_valid.to_pandas())
        .mark_circle(size=60, opacity=0.7)
        .encode(
            x=alt.X("gt_dominance:Q", title="GT Dominance (top-100 max reagent share)"),
            y=alt.Y("cats_concentration:Q", title="CATS Concentration (1 − norm. posterior entropy)"),
            color=alt.Color("library_id:N", title="Library"),
            tooltip=["library_id", "query_id", "component_idx",
                      "cats_concentration", "gt_dominance"],
        )
    )

    _line = (
        alt.Chart(_reg_line.to_pandas())
        .mark_line(strokeDash=[4, 4], color="gray", strokeWidth=1.5)
        .encode(x="x:Q", y="y:Q")
    )

    concentration_scatter = (
        (_scatter + _line)
        .properties(
            width=500,
            height=450,
            title=f"Pearson r={_r_pearson:.3f}, Spearman ρ={_r_spearman:.3f} (n={len(_valid)})",
        )
    )

    concentration_scatter
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Figure 2: Criticality Trajectory

    Per-component criticality evolution over search cycles for the selected library and query.
    The shaded band marks the "diffuse" zone (criticality ≤ 0.3).
    """)
    return


@app.cell
def _(
    diag_df,
    library_dropdown,
    mo,
    pl,
    plot_criticality_trajectory,
    plt,
    query_dropdown,
):
    """Criticality trajectory for selected lib+query."""
    _lib = library_dropdown.value
    _query = query_dropdown.value

    _filtered = diag_df.filter(
        (pl.col("library_id") == _lib) & (pl.col("query_id") == _query)
    )

    mo.stop(
        len(_filtered) == 0,
        mo.md(f"*No diagnostics data for {_lib} / {_query}*"),
    )

    _fig_crit = plot_criticality_trajectory(_filtered)
    _fig_crit.suptitle(f"Criticality Trajectory — {_lib} / {_query}", fontsize=12)
    plt.gca()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Figure 3: SNR Trajectory

    Signal-to-noise ratio evolution. SNR > 1 means the signal between reagent
    means exceeds the noise from posterior uncertainty.
    """)
    return


@app.cell
def _(
    diag_df,
    library_dropdown,
    mo,
    pl,
    plot_snr_trajectory,
    plt,
    query_dropdown,
):
    """SNR trajectory for selected lib+query."""
    _lib = library_dropdown.value
    _query = query_dropdown.value

    _filtered = diag_df.filter(
        (pl.col("library_id") == _lib) & (pl.col("query_id") == _query)
    )

    mo.stop(
        len(_filtered) == 0,
        mo.md(f"*No diagnostics data for {_lib} / {_query}*"),
    )

    _fig_snr = plot_snr_trajectory(_filtered)
    _fig_snr.suptitle(f"SNR Trajectory — {_lib} / {_query}", fontsize=12)
    plt.gca()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Figure 4: Convergence Timing

    How quickly does CATS identify SAR structure? Bar chart shows `cycle_stable`
    (the first cycle after which criticality stays above 0.3) per component
    across all queries.
    """)
    return


@app.cell
def _(alt, conv_df, pl):
    """Convergence timing bar chart."""
    _conv_valid = conv_df.filter(pl.col("cycle_stable").is_not_null())

    _agg = (
        _conv_valid
        .group_by(["library_id", "component_idx"])
        .agg([
            pl.col("cycle_stable").mean().round(1).alias("mean_cycle_stable"),
            pl.col("cycle_stable").median().alias("median_cycle_stable"),
            pl.len().alias("n_queries"),
        ])
        .sort(["library_id", "component_idx"])
    )

    convergence_chart = (
        alt.Chart(_agg.to_pandas())
        .mark_bar(opacity=0.8)
        .encode(
            x=alt.X("library_id:N", title="Library", sort=None),
            y=alt.Y("mean_cycle_stable:Q", title="Mean Cycle to Stable Convergence"),
            color=alt.Color("component_idx:N", title="Component"),
            xOffset="component_idx:N",
            tooltip=["library_id", "component_idx", "mean_cycle_stable",
                      "median_cycle_stable", "n_queries"],
        )
        .properties(width=600, height=350, title="Convergence Timing by Library and Component")
    )

    convergence_chart
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Figure 5: Temperature Decomposition

    Full CATS temperature pipeline decomposition for the selected library, query,
    and component. Shows how criticality, CATS multiplier, and thermal cycling
    combine to produce the final sampling temperature.
    """)
    return


@app.cell
def _(
    component_slider,
    diag_df,
    library_dropdown,
    mo,
    pl,
    plot_temperature_decomposition,
    plt,
    query_dropdown,
):
    """Temperature decomposition for selected lib+query+component."""
    _lib = library_dropdown.value
    _query = query_dropdown.value
    _comp = component_slider.value

    _filtered = diag_df.filter(
        (pl.col("library_id") == _lib) & (pl.col("query_id") == _query)
    )

    mo.stop(
        len(_filtered) == 0 or "current_cycle" not in _filtered.columns,
        mo.md(f"*No enhanced diagnostics for {_lib} / {_query}*"),
    )

    _fig_temp = plot_temperature_decomposition(_filtered, component_idx=_comp)
    _fig_temp.suptitle(
        f"Temperature Decomposition — {_lib} / {_query} / Component {_comp}",
        fontsize=12,
    )
    plt.gca()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Summary Statistics

    **Rank-order correlation** between CATS online estimates and brute-force ground truth.
    Spearman ρ is the primary metric — it tests whether CATS correctly **ranks** components
    by SAR structure, regardless of absolute scale differences.

    MAE is not reported because the two metrics are on different scales by construction
    (posterior entropy vs top-N frequency share).
    """)
    return


@app.cell
def _(mo, np, pearsonr, pl, sar_df, spearmanr):
    """Summary statistics table."""
    _valid = sar_df.filter(
        pl.col("cats_concentration").is_not_nan()
        & pl.col("gt_dominance").is_not_nan()
    )

    if len(_valid) < 3:
        mo.md("*Insufficient data for statistics.*")
    else:
        _cats = _valid["cats_concentration"].to_numpy()
        _gt = _valid["gt_dominance"].to_numpy()
        _r_p, _p_p = pearsonr(_cats, _gt)
        _r_s, _p_s = spearmanr(_cats, _gt)

        # Rank-order agreement: does the most-concentrated component in CATS
        # match the most-dominant component in GT within each (lib, query)?
        _rank_agreement = _compute_rank_agreement(_valid)

        mo.md(f"""
| Metric | Value | Interpretation |
|--------|-------|----------------|
| **Spearman ρ** | {_r_s:.4f} (p = {_p_s:.2e}) | Rank-order agreement (primary metric) |
| **Pearson r** | {_r_p:.4f} (p = {_p_p:.2e}) | Linear correlation |
| **Component rank agreement** | {_rank_agreement:.1%} | Most-structured component matches GT |
| **N components** | {len(_valid)} | |

**How to read this:** Spearman ρ > 0 means CATS assigns higher concentration to components
that genuinely have higher dominance in the full combinatorial space. A value of ~0.3–0.5
with only 2% sampling budget is encouraging — CATS identifies the SAR *direction* from
limited data.
""")
    return


@app.cell
def _(np, pl):
    """Helper: within-query rank agreement."""
    def _compute_rank_agreement(valid_df: pl.DataFrame) -> float:
        _matches = 0
        _total = 0
        for (_lib, _query), _group in valid_df.group_by(["library_id", "query_id"]):
            if len(_group) < 2:
                continue
            _cats_best = _group.sort("cats_concentration", descending=True)["component_idx"][0]
            _gt_best = _group.sort("gt_dominance", descending=True)["component_idx"][0]
            _matches += int(_cats_best == _gt_best)
            _total += 1
        return _matches / max(_total, 1)
    return (_compute_rank_agreement,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Per-Library Summary

    Mean CATS concentration and GT dominance per library, plus within-library
    Pearson correlation.
    """)
    return


@app.cell
def _(np, pearsonr, pl, sar_df):
    """Per-library summary with correlations."""
    _valid = sar_df.filter(
        pl.col("cats_concentration").is_not_nan()
        & pl.col("gt_dominance").is_not_nan()
    )

    _rows = []
    for _lib in sorted(_valid["library_id"].unique().to_list()):
        _lib_data = _valid.filter(pl.col("library_id") == _lib)
        _c = _lib_data["cats_concentration"].to_numpy()
        _g = _lib_data["gt_dominance"].to_numpy()
        _r, _p = pearsonr(_c, _g) if len(_c) > 2 else (float("nan"), float("nan"))
        _rows.append({
            "library_id": _lib,
            "mean_cats_concentration": round(float(np.mean(_c)), 3),
            "mean_gt_dominance": round(float(np.mean(_g)), 3),
            "pearson_r": round(float(_r), 3),
            "p_value": round(float(_p), 4),
            "n_components": len(_c),
        })

    per_library_summary = pl.DataFrame(_rows)
    per_library_summary
    return


if __name__ == "__main__":
    app.run()
