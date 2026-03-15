"""
TACTICS Manuscript Plots

Loads the manuscript benchmark results (87,240 trials, 4 methods, 20 libraries)
and provides DataFrames ready for plotting.

Run as app: marimo run tutorials/manuscript_plots.py
Edit mode:  marimo edit tutorials/manuscript_plots.py
"""

import marimo

__generated_with = "0.19.11"
app = marimo.App(width="full", app_title="TACTICS Manuscript Plots")


@app.cell
def _():
    """Imports."""
    import io
    import marimo as mo
    import sys
    from pathlib import Path

    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    import numpy as np
    import polars as pl
    from scipy import stats

    def fig_to_pdf_bytes(fig):
        """Convert a matplotlib figure to PDF bytes in memory."""
        buf = io.BytesIO()
        fig.savefig(buf, format="pdf", bbox_inches="tight", dpi=300)
        return buf.getvalue()

    try:
        project_root = Path(__file__).parent.parent.resolve()
    except NameError:
        project_root = Path.cwd()

    sys.path.insert(0, str(project_root))
    sys.path.insert(0, str(project_root / "src"))
    return (
        Line2D,
        Patch,
        fig_to_pdf_bytes,
        mo,
        np,
        pl,
        plt,
        project_root,
        stats,
    )


@app.cell
def _(mo, pl, project_root):
    """Load benchmark results."""
    results_dir = project_root / "outputs" / "manuscript_benchmark_array" / "manuscript_benchmark"

    benchmark_results = pl.read_parquet(results_dir / "benchmark_results.parquet")

    # Exclude thrombin (Legacy RWS crashed — scaling=-1 compatibility issue)
    df = benchmark_results.filter(pl.col("library_id") != "thrombin")

    n_trials = len(df)
    n_libraries = df["library_id"].n_unique()
    methods = sorted(df["method"].unique().to_list())

    mo.md(
        f"""
        ## Data Loaded

        - **Source**: `{results_dir}`
        - **Trials**: {n_trials:,} (thrombin excluded)
        - **Libraries**: {n_libraries}
        - **Methods**: {', '.join(methods)}
        """
    )
    return (df,)


@app.cell
def _(df, mo, pl):
    """Overall method summary."""
    method_summary = (
        df.group_by("method")
        .agg([
            pl.col("recovered").mean().round(2).alias("mean_recovery"),
            pl.col("recovered").std().round(2).alias("std_recovery"),
            pl.col("recovered").median().alias("median_recovery"),
            pl.col("n_evaluated").mean().round(0).alias("mean_evaluations"),
        ])
        .sort("mean_recovery", descending=True)
    )

    mo.md("## Overall Method Summary")
    mo.ui.table(method_summary.to_pandas())
    return


@app.cell
def _(df, mo, pl):
    """Per-library per-method mean recovery (wide format for plotting)."""
    lib_method = (
        df.group_by(["library_id", "n_components", "method"])
        .agg(pl.col("recovered").mean().round(1).alias("mean_recovery"))
        .sort(["library_id", "method"])
    )

    # Wide format: one column per method
    lib_wide = lib_method.pivot(
        on="method",
        index=["library_id", "n_components"],
        values="mean_recovery",
    ).sort("library_id")

    mo.vstack([
        mo.md("## Per-Library Recovery (Wide Format)"),
        mo.ui.table(lib_wide.to_pandas()),
    ])
    return


@app.cell
def _(df, mo, pl):
    """Per-library per-method recovery with std (for error bars)."""
    lib_method_stats = (
        df.group_by(["library_id", "n_components", "method"])
        .agg([
            pl.col("recovered").mean().round(2).alias("mean_recovery"),
            pl.col("recovered").std().round(2).alias("std_recovery"),
            pl.col("recovered").count().alias("n_trials"),
        ])
        .with_columns(
            (pl.col("std_recovery") / pl.col("n_trials").sqrt())
            .round(2)
            .alias("se_recovery")
        )
        .sort(["library_id", "method"])
    )

    mo.vstack([
        mo.md("## Per-Library Recovery with Error Bars"),
        mo.ui.table(lib_method_stats.to_pandas()),
    ])
    return


@app.cell
def _(df, mo, pl):
    """Component-count breakdown."""
    component_summary = (
        df.group_by(["n_components", "method"])
        .agg([
            pl.col("recovered").mean().round(2).alias("mean_recovery"),
            pl.col("recovered").std().round(2).alias("std_recovery"),
        ])
        .sort(["n_components", "mean_recovery"], descending=[False, True])
    )

    mo.vstack([
        mo.md("## Component-Count Breakdown"),
        mo.ui.table(component_summary.to_pandas()),
    ])
    return


@app.cell
def _(df, mo, np, pl, stats):
    """Pairwise significance tests (TT-TS vs each baseline, per library)."""
    methods_order = ["Legacy TS", "Legacy RWS", "TACTICS RWS", "TACTICS TT-TS"]
    libs = sorted(df["library_id"].unique().to_list())

    comparisons = [
        ("Legacy TS", "TACTICS TT-TS"),
        ("Legacy RWS", "TACTICS TT-TS"),
        ("Legacy RWS", "TACTICS RWS"),
        ("TACTICS RWS", "TACTICS TT-TS"),
    ]

    test_rows = []
    for baseline, challenger in comparisons:
        for _lib_id in libs:
            _lib_data = df.filter(pl.col("library_id") == _lib_id)
            _a = (
                _lib_data.filter(pl.col("method") == baseline)
                .sort(["query_id", "replicate"])["recovered"]
                .to_numpy()
            )
            _b = (
                _lib_data.filter(pl.col("method") == challenger)
                .sort(["query_id", "replicate"])["recovered"]
                .to_numpy()
            )
            _n = min(len(_a), len(_b))
            if _n < 2:
                continue
            _a, _b = _a[:_n], _b[:_n]
            _t_stat, _p_value = stats.ttest_rel(_a, _b)
            _mean_diff = float(np.mean(_b) - np.mean(_a))
            _sig = (
                "****" if _p_value <= 0.0001
                else "***" if _p_value <= 0.001
                else "**" if _p_value <= 0.01
                else "*" if _p_value <= 0.05
                else "ns"
            )
            _n_comp = _lib_data["n_components"].to_list()[0]
            test_rows.append({
                "comparison": f"{challenger} vs {baseline}",
                "library_id": _lib_id,
                "n_components": _n_comp,
                "baseline": baseline,
                "challenger": challenger,
                "mean_baseline": round(float(np.mean(_a)), 2),
                "mean_challenger": round(float(np.mean(_b)), 2),
                "mean_diff": round(_mean_diff, 2),
                "t_stat": round(float(_t_stat), 4),
                "p_value": round(float(_p_value), 6),
                "significance": _sig,
                "n_pairs": _n,
            })

    pairwise_df = pl.DataFrame(test_rows)

    mo.vstack([
        mo.md("## Pairwise Significance Tests"),
        mo.ui.table(pairwise_df.to_pandas()),
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Plot the Differences in Recovery between different Methods
    """)
    return


@app.cell
def _(df, mo):
    """Select baseline and challenger methods for difference plot."""
    _methods = sorted(df["method"].unique().to_list())

    baseline_dropdown = mo.ui.dropdown(
        options=_methods,
        value="Legacy RWS",
        label="Baseline method",
    )
    challenger_dropdown = mo.ui.dropdown(
        options=_methods,
        value="TACTICS TT-TS",
        label="Challenger method",
    )

    mo.hstack(
        [baseline_dropdown, challenger_dropdown],
        justify="start",
        gap=2,
    )
    return baseline_dropdown, challenger_dropdown


@app.cell
def _(baseline_dropdown, challenger_dropdown, df, mo, np, pl, stats):
    """Compute per-library paired recovery differences with 95% CI."""
    _baseline = baseline_dropdown.value
    _challenger = challenger_dropdown.value
    _libs = sorted(df["library_id"].unique().to_list())

    diff_rows = []
    for _lib_id in _libs:
        _lib_data = df.filter(pl.col("library_id") == _lib_id)
        _n_comp = _lib_data["n_components"].to_list()[0]

        _a = (
            _lib_data.filter(pl.col("method") == _baseline)
            .sort(["query_id", "replicate"])["recovered"]
            .to_numpy()
        )
        _b = (
            _lib_data.filter(pl.col("method") == _challenger)
            .sort(["query_id", "replicate"])["recovered"]
            .to_numpy()
        )

        _n = min(len(_a), len(_b))
        if _n < 2:
            continue
        _a, _b = _a[:_n], _b[:_n]

        _diffs = _b - _a
        _mean_diff = float(np.mean(_diffs))
        _se_diff = float(np.std(_diffs, ddof=1) / np.sqrt(_n))

        _t_crit = float(stats.t.ppf(0.975, df=_n - 1))
        _ci_low = _mean_diff - _t_crit * _se_diff
        _ci_high = _mean_diff + _t_crit * _se_diff

        _t_stat, _p_value = stats.ttest_rel(_a, _b)
        _sig = (
            "****" if _p_value <= 0.0001
            else "***" if _p_value <= 0.001
            else "**" if _p_value <= 0.01
            else "*" if _p_value <= 0.05
            else "ns"
        )

        diff_rows.append({
            "library_id": _lib_id,
            "n_components": _n_comp,
            "mean_diff": round(_mean_diff, 2),
            "se_diff": round(_se_diff, 2),
            "ci_low": round(_ci_low, 2),
            "ci_high": round(_ci_high, 2),
            "p_value": round(float(_p_value), 6),
            "significance": _sig,
            "n_pairs": _n,
        })

    diff_df = pl.DataFrame(diff_rows).sort("mean_diff", descending=True)

    mo.vstack([
        mo.md(f"## Paired Recovery Differences: {_challenger} vs {_baseline}"),
        mo.ui.table(diff_df.to_pandas()),
    ])
    return (diff_df,)


@app.cell
def _(
    Patch,
    baseline_dropdown,
    challenger_dropdown,
    diff_df,
    fig_to_pdf_bytes,
    mo,
    np,
    plt,
):
    """Per-library recovery difference bar plot."""

    _baseline = baseline_dropdown.value
    _challenger = challenger_dropdown.value
    _n_libs = len(diff_df)

    _labels = diff_df["library_id"].to_list()
    _means = diff_df["mean_diff"].to_numpy()
    _ci_low = diff_df["ci_low"].to_numpy()
    _ci_high = diff_df["ci_high"].to_numpy()
    _n_comps = diff_df["n_components"].to_list()
    _sigs = diff_df["significance"].to_list()

    _yerr = np.array([_means - _ci_low, _ci_high - _means])

    _color_map = {2: "#4C78A8", 3: "#E45756"}
    _colors = [_color_map.get(nc, "#999999") for nc in _n_comps]

    _fig_width = max(12, _n_libs * 0.7)
    _fig, _ax = plt.subplots(figsize=(_fig_width, 6), facecolor="white")
    _ax.set_facecolor("white")

    _x = np.arange(_n_libs)
    _ax.bar(
        _x, _means, yerr=_yerr,
        color=_colors, edgecolor="white", linewidth=0.5,
        capsize=4, error_kw={"linewidth": 1.2, "color": "black"},
        zorder=3,
    )

    _ax.axhline(0, color="black", linewidth=0.8, linestyle="-", zorder=2)

    # Annotation offset scaled to data range so stars stay inside the plot
    _data_range = max(abs(_ci_high.max()), abs(_ci_low.min()), 1)
    _offset = _data_range * 0.04

    for _i, (_sig, _ci_h, _ci_l, _m) in enumerate(
        zip(_sigs, _ci_high, _ci_low, _means)
    ):
        _y_pos = _ci_h + _offset if _m >= 0 else _ci_l - _offset
        _va = "bottom" if _m >= 0 else "top"
        _ax.text(
            _i, _y_pos, _sig,
            ha="center", va=_va,
            fontsize=12, fontweight="bold", color="black",
        )

    # Expand y-limits to fit annotations
    _y_min = min(_ci_low.min(), _means.min()) - _data_range * 0.15
    _y_max = max(_ci_high.max(), _means.max()) + _data_range * 0.15
    _ax.set_ylim(_y_min, _y_max)

    _ax.set_xticks(_x)
    _ax.set_xticklabels(_labels, rotation=45, ha="right", fontsize=12, color="black")
    _ax.set_ylabel(
        "\u0394 % Recovery of Top-100 Hits",
        fontsize=13, color="black",
    )
    _ax.set_title(
        f"Per-Library Recovery of Top-100 Hit Molecules: {_challenger} vs {_baseline}",
        fontsize=14, fontweight="bold", color="black",
    )
    _ax.tick_params(colors="black", labelsize=12)

    _legend_elements = [
        Patch(facecolor="#4C78A8", edgecolor="black", label="2-component"),
        Patch(facecolor="#E45756", edgecolor="black", label="3-component"),
    ]
    _ax.legend(handles=_legend_elements, loc="upper right", fontsize=12,
               facecolor="white", edgecolor="black", labelcolor="black")

    _ax.grid(axis="y", alpha=0.3, zorder=1, color="grey")
    _ax.set_xlim(-0.6, _n_libs - 0.4)
    for _spine in _ax.spines.values():
        _spine.set_color("black")
    _fig.tight_layout()

    _pdf_name = f"per_library_recovery_diff_{_challenger}_{_baseline}.pdf".replace(" ", "_")
    mo.vstack([
        plt.gca(),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename=_pdf_name,
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


@app.cell
def _(df, mo, pl):
    """Gain decomposition: criticality vs TT-TS."""
    libs_list = sorted(df["library_id"].unique().to_list())

    decomp_rows = []
    for _lib_id in libs_list:
        _lib_data = df.filter(pl.col("library_id") == _lib_id)
        _n_comp = _lib_data["n_components"].to_list()[0]

        _lrws = float(_lib_data.filter(pl.col("method") == "Legacy RWS")["recovered"].mean())
        _trws = float(_lib_data.filter(pl.col("method") == "TACTICS RWS")["recovered"].mean())
        _ttts = float(_lib_data.filter(pl.col("method") == "TACTICS TT-TS")["recovered"].mean())

        _d_crit = _trws - _lrws  # criticality-weighted rotation
        _d_ttts = _ttts - _trws  # TT-TS selection policy
        _d_total = _ttts - _lrws
        _abs_sum = abs(_d_crit) + abs(_d_ttts)

        _share_crit = round(abs(_d_crit) / _abs_sum * 100, 0) if _abs_sum > 0.01 else 0
        _share_ttts = round(abs(_d_ttts) / _abs_sum * 100, 0) if _abs_sum > 0.01 else 0

        # Signed shares: multiply share by +1/-1 based on direction
        _signed_crit = _share_crit if _d_crit >= 0 else -_share_crit
        _signed_ttts = _share_ttts if _d_ttts >= 0 else -_share_ttts
        _net_effect = _signed_crit + _signed_ttts

        decomp_rows.append({
            "library_id": _lib_id,
            "n_components": _n_comp,
            "delta_crit": round(_d_crit, 2),
            "delta_ttts": round(_d_ttts, 2),
            "delta_total": round(_d_total, 2),
            "share_crit_%": _share_crit,
            "share_ttts_%": _share_ttts,
            "signed_crit_%": _signed_crit,
            "signed_ttts_%": _signed_ttts,
            "net_effect_%": _net_effect,
        })

    gain_decomposition = pl.DataFrame(decomp_rows)

    mo.vstack([
        mo.md(
            """
            ## Gain Decomposition: Criticality vs TT-TS

            **Signed shares** express each component as a percentage of `delta_total`
            and always sum to 100%. A **positive** share means the component pushes
            in the same direction as the overall gain; a **negative** share means
            it partially offsets the other component.

            *Example*: signed\_share\_crit = +130%, signed\_share\_ttts = −30%
            means criticality drove 130% of the gain while TT-TS offset 30%,
            netting to the observed total delta.
            """
        ),
        mo.ui.table(gain_decomposition.to_pandas()),
    ])
    return (gain_decomposition,)


@app.cell
def _(
    Line2D,
    Patch,
    df,
    fig_to_pdf_bytes,
    gain_decomposition,
    mo,
    np,
    pl,
    plt,
    stats,
):
    """Diverging stacked bar: criticality vs TT-TS contributions by component count."""

    # Compute per-library significance for TT-TS vs Legacy RWS
    _sig_map = {}
    for _lib_id in gain_decomposition["library_id"].to_list():
        _lib_data = df.filter(pl.col("library_id") == _lib_id)
        _a = (
            _lib_data.filter(pl.col("method") == "Legacy RWS")
            .sort(["query_id", "replicate"])["recovered"].to_numpy()
        )
        _b = (
            _lib_data.filter(pl.col("method") == "TACTICS TT-TS")
            .sort(["query_id", "replicate"])["recovered"].to_numpy()
        )
        _n_pairs = min(len(_a), len(_b))
        if _n_pairs < 2:
            _sig_map[_lib_id] = "ns"
            continue
        _, _p = stats.ttest_rel(_a[:_n_pairs], _b[:_n_pairs])
        _sig_map[_lib_id] = (
            "****" if _p <= 0.0001
            else "***" if _p <= 0.001
            else "**" if _p <= 0.01
            else "*" if _p <= 0.05
            else "ns"
        )

    _crit_color = "#59A14F"   # green (distinct from blue/red component-count palette)
    _ttts_color = "#F28E2B"   # orange

    _groups = {
        2: gain_decomposition.filter(pl.col("n_components") == 2).sort("delta_total", descending=True),
        3: gain_decomposition.filter(pl.col("n_components") == 3).sort("delta_total", descending=True),
    }

    _fig, _axes = plt.subplots(1, 2, figsize=(16, 6), facecolor="white")

    # First pass: draw bars and collect extents for global y-limits
    _global_top = 0.0
    _global_bot = 0.0
    _panel_data = []

    for _ax, (_nc, _gdf) in zip(_axes, _groups.items()):
        _ax.set_facecolor("white")
        _n = len(_gdf)
        if _n == 0:
            _panel_data.append(None)
            continue

        _labels = _gdf["library_id"].to_list()
        _d_crit = _gdf["delta_crit"].to_numpy().astype(float)
        _d_ttts = _gdf["delta_ttts"].to_numpy().astype(float)
        _d_total = _gdf["delta_total"].to_numpy().astype(float)
        _x = np.arange(_n)

        _bar_top = np.zeros(_n)
        _bar_bot = np.zeros(_n)

        for _i in range(_n):
            _pos_bottom = 0.0
            _neg_bottom = 0.0

            _dc = _d_crit[_i]
            if _dc >= 0:
                _ax.bar(_i, _dc, bottom=_pos_bottom, color=_crit_color,
                        edgecolor="white", linewidth=0.5, width=0.7, zorder=3)
                _pos_bottom += _dc
            else:
                _ax.bar(_i, _dc, bottom=_neg_bottom, color=_crit_color,
                        edgecolor="white", linewidth=0.5, width=0.7, zorder=3)
                _neg_bottom += _dc

            _dt = _d_ttts[_i]
            if _dt >= 0:
                _ax.bar(_i, _dt, bottom=_pos_bottom, color=_ttts_color,
                        edgecolor="white", linewidth=0.5, width=0.7, zorder=3)
                _pos_bottom += _dt
            else:
                _ax.bar(_i, _dt, bottom=_neg_bottom, color=_ttts_color,
                        edgecolor="white", linewidth=0.5, width=0.7, zorder=3)
                _neg_bottom += _dt

            _bar_top[_i] = _pos_bottom
            _bar_bot[_i] = _neg_bottom

        _global_top = max(_global_top, _bar_top.max())
        _global_bot = min(_global_bot, _bar_bot.min())

        _panel_data.append((_labels, _d_total, _bar_top, _bar_bot, _x, _n))

        _ax.scatter(_x, _d_total, marker="D", color="black", s=30, zorder=5)
        _ax.axhline(0, color="black", linewidth=0.8, zorder=2)

        _ax.set_xticks(_x)
        _ax.set_xticklabels(_labels, rotation=45, ha="right", fontsize=12, color="black")
        _ax.set_title(f"{_nc}-Component Libraries", fontsize=14, fontweight="bold", color="black")
        _ax.tick_params(colors="black", labelsize=12)
        _ax.grid(axis="y", alpha=0.3, zorder=1, color="grey")
        _ax.set_xlim(-0.6, _n - 0.4)
        for _spine in _ax.spines.values():
            _spine.set_color("black")

    # Second pass: add stars and set uniform y-limits
    _data_range = max(abs(_global_top), abs(_global_bot), 1)
    _offset = _data_range * 0.05
    _y_min = _global_bot - _data_range * 0.25
    _y_max = _global_top + _data_range * 0.25

    for _ax, _pd in zip(_axes, _panel_data):
        if _pd is None:
            continue
        _labels, _d_total, _bar_top, _bar_bot, _x, _n = _pd
        for _i, _lib in enumerate(_labels):
            _sig = _sig_map.get(_lib, "")
            _y_pos = _bar_top[_i] + _offset if _d_total[_i] >= 0 else _bar_bot[_i] - _offset
            _va = "bottom" if _d_total[_i] >= 0 else "top"
            _ax.text(_i, _y_pos, _sig, ha="center", va=_va,
                     fontsize=12, fontweight="bold", color="black")
        _ax.set_ylim(_y_min, _y_max)

    _axes[0].set_ylabel("\u0394 % Recovery of Top-100 Hits\n(TACTICS TT-TS \u2212 Legacy RWS)",
                        fontsize=13, color="black")

    _legend_elements = [
        Patch(facecolor=_crit_color, edgecolor="black", label="Criticality weighting"),
        Patch(facecolor=_ttts_color, edgecolor="black", label="TT-TS selection"),
        Line2D([0], [0], marker="D", color="black", linestyle="None",
               markersize=6, label="Net total (\u0394 total)"),
    ]
    _axes[1].legend(handles=_legend_elements, loc="upper right", fontsize=12,
                    facecolor="white", edgecolor="black", labelcolor="black")

    _fig.suptitle("Gain Decomposition: Criticality Weighting vs TT-TS Selection",
                  fontsize=15, fontweight="bold", color="black", y=1.01)
    _fig.tight_layout()

    mo.vstack([
        plt.gcf(),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename="gain_decomposition.pdf",
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


@app.cell
def _(df, mo, pl):
    """Per-query recovery for detailed analysis (e.g., scatter plots)."""
    query_method_recovery = (
        df.group_by(["library_id", "n_components", "query_id", "method"])
        .agg([
            pl.col("recovered").mean().round(1).alias("mean_recovery"),
            pl.col("recovered").std().round(1).alias("std_recovery"),
        ])
    )

    mo.vstack([
        mo.md("## Per-Query Recovery"),
        mo.ui.table(query_method_recovery.to_pandas()),
    ])
    return


@app.cell
def _(mo):
    """Available DataFrames reference."""
    mo.md(
        """
        ## Available DataFrames

        | Variable | Description |
        |---|---|
        | `df` | Raw benchmark results (87,200 rows, excl. thrombin) |
        | `method_summary` | Overall method means/std/median |
        | `lib_method` | Per-library per-method mean recovery (long format) |
        | `lib_wide` | Per-library recovery, one column per method (wide format) |
        | `lib_method_stats` | Per-library per-method mean, std, SE |
        | `component_summary` | Recovery by 2-comp vs 3-comp |
        | `pairwise_df` | All pairwise significance tests per library |
        | `gain_decomposition` | Δ criticality vs Δ TT-TS per library (signed shares sum to 100%) |
        | `query_method_recovery` | Per-query per-method mean recovery |
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Sensitivity Analysis - Top N compounds
    In this section we explore the roubustness of CATS as the interval for the top N compounds increases, if CATS is interpreting the SAR landscape correctly and systematically able to uncover mid tier hits.
    """)
    return


@app.cell
def _(mo, pl, project_root):
    """Load top-N sensitivity benchmark results."""
    _topn_dir = (
        project_root / "outputs" / "topn_sensitivity_benchmark_array"
        / "topn_sensitivity_benchmark"
    )
    topn_df = pl.read_parquet(_topn_dir / "benchmark_results.parquet")

    _n_rows = len(topn_df)
    _n_libs = topn_df["library_id"].n_unique()
    _thresholds = sorted(topn_df["top_n"].unique().to_list())

    mo.vstack([
        mo.md(
            f"""
            ## Top-N Sensitivity Data Loaded

            - **Source**: `{_topn_dir}`
            - **Rows**: {_n_rows:,}
            - **Libraries**: {_n_libs}
            - **Top-N thresholds**: {_thresholds}
            """
        ),
    ])
    return (topn_df,)


@app.cell
def _(fig_to_pdf_bytes, mo, pl, plt, topn_df):
    """Overall recovery vs top-N by component count (2-comp and 3-comp panels)."""

    _method_colors = {
        "Legacy TS": "#76B7B2",
        "Legacy RWS": "#B07AA1",
        "TACTICS RWS": "#4C78A8",
        "TACTICS TT-TS": "#E45756",
    }
    _method_markers = {
        "Legacy TS": "s",
        "Legacy RWS": "^",
        "TACTICS RWS": "o",
        "TACTICS TT-TS": "D",
    }

    _fig, _axes = plt.subplots(1, 2, figsize=(16, 6), facecolor="white",
                               sharey=True)

    for _panel_idx, _nc in enumerate([2, 3]):
        _ax = _axes[_panel_idx]
        _ax.set_facecolor("white")

        _comp_data = topn_df.filter(pl.col("n_components") == _nc)

        _agg = (
            _comp_data.group_by(["method", "top_n"])
            .agg([
                pl.col("recovery_frac").mean().alias("mean_frac"),
                pl.col("recovery_frac").std().alias("std_frac"),
                pl.col("recovery_frac").count().alias("n"),
            ])
            .with_columns(
                (pl.col("std_frac") / pl.col("n").sqrt()).alias("se_frac")
            )
            .sort(["method", "top_n"])
        )

        _methods = sorted(_agg["method"].unique().to_list())

        for _method in _methods:
            _mdata = _agg.filter(pl.col("method") == _method).sort("top_n")
            _x = _mdata["top_n"].to_numpy()
            _y = _mdata["mean_frac"].to_numpy() * 100
            _se = _mdata["se_frac"].to_numpy() * 100
            _color = _method_colors.get(_method, "#999999")
            _marker = _method_markers.get(_method, "o")

            _ax.plot(_x, _y, marker=_marker, color=_color, linewidth=2,
                     markersize=7, label=_method, zorder=3)
            _ax.fill_between(_x, _y - 1.96 * _se, _y + 1.96 * _se,
                             alpha=0.15, color=_color, zorder=2)

        _ax.set_xlabel("Top-N Threshold", fontsize=13, color="black")
        _ax.set_title(f"{_nc}-Component Libraries", fontsize=14,
                      fontweight="bold", color="black")
        _ax.set_xticks(sorted(_agg["top_n"].unique().to_list()))
        _ax.tick_params(colors="black", labelsize=12)
        _ax.grid(axis="both", alpha=0.3, color="grey", zorder=1)
        for _spine in _ax.spines.values():
            _spine.set_color("black")

    _axes[0].set_ylabel("% Recovery", fontsize=13, color="black")
    _axes[1].legend(fontsize=11, facecolor="white", edgecolor="black",
                    labelcolor="black", loc="lower right")

    _fig.suptitle("Mean Recovery Rate vs Top-N Threshold (All Libraries)",
                  fontsize=15, fontweight="bold", color="black", y=1.01)
    _fig.tight_layout()

    mo.vstack([
        plt.gcf(),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename="overall_recovery_vs_topn.pdf",
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


@app.cell
def _(mo, topn_df):
    """Select library for top-N recovery curve."""
    _libs = sorted(topn_df["library_id"].unique().to_list())

    topn_library_dropdown = mo.ui.dropdown(
        options=_libs,
        value=_libs[0],
        label="Library",
    )
    topn_library_dropdown
    return (topn_library_dropdown,)


@app.cell
def _(fig_to_pdf_bytes, mo, pl, plt, topn_df, topn_library_dropdown):
    """Top-N recovery curve: recovery rate vs top-N threshold per method."""

    _lib_id = topn_library_dropdown.value
    _lib_data = topn_df.filter(pl.col("library_id") == _lib_id)

    # Aggregate: mean and SE of recovery_frac per (method, top_n)
    _agg = (
        _lib_data.group_by(["method", "top_n"])
        .agg([
            pl.col("recovery_frac").mean().alias("mean_frac"),
            pl.col("recovery_frac").std().alias("std_frac"),
            pl.col("recovery_frac").count().alias("n"),
        ])
        .with_columns(
            (pl.col("std_frac") / pl.col("n").sqrt()).alias("se_frac")
        )
        .sort(["method", "top_n"])
    )

    _methods = sorted(_agg["method"].unique().to_list())
    _method_colors = {
        "Legacy TS": "#76B7B2",
        "Legacy RWS": "#B07AA1",
        "TACTICS RWS": "#4C78A8",
        "TACTICS TT-TS": "#E45756",
    }
    _method_markers = {
        "Legacy TS": "s",
        "Legacy RWS": "^",
        "TACTICS RWS": "o",
        "TACTICS TT-TS": "D",
    }

    _fig, _ax = plt.subplots(figsize=(10, 6), facecolor="white")
    _ax.set_facecolor("white")

    for _method in _methods:
        _mdata = _agg.filter(pl.col("method") == _method).sort("top_n")
        _x = _mdata["top_n"].to_numpy()
        _y = _mdata["mean_frac"].to_numpy() * 100
        _se = _mdata["se_frac"].to_numpy() * 100
        _color = _method_colors.get(_method, "#999999")
        _marker = _method_markers.get(_method, "o")

        _ax.plot(_x, _y, marker=_marker, color=_color, linewidth=2,
                 markersize=7, label=_method, zorder=3)
        _ax.fill_between(_x, _y - 1.96 * _se, _y + 1.96 * _se,
                         alpha=0.15, color=_color, zorder=2)

    _ax.set_xlabel("Top-N Threshold", fontsize=13, color="black")
    _ax.set_ylabel("% Recovery", fontsize=13, color="black")
    _ax.set_title(
        f"Recovery Rate vs Top-N Threshold: {_lib_id}",
        fontsize=14, fontweight="bold", color="black",
    )
    _ax.set_xticks(sorted(_agg["top_n"].unique().to_list()))
    _ax.tick_params(colors="black", labelsize=12)
    _ax.legend(fontsize=12, facecolor="white", edgecolor="black", labelcolor="black")
    _ax.grid(axis="both", alpha=0.3, color="grey", zorder=1)
    for _spine in _ax.spines.values():
        _spine.set_color("black")
    _fig.tight_layout()

    _pdf_name = f"recovery_vs_topn_{_lib_id}.pdf".replace(" ", "_")
    mo.vstack([
        plt.gca(),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename=_pdf_name,
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


@app.cell
def _(mo, topn_df):
    """Select baseline and challenger methods for top-N difference plot."""
    _methods = sorted(topn_df["method"].unique().to_list())

    topn_baseline_dropdown = mo.ui.dropdown(
        options=_methods,
        value="Legacy RWS",
        label="Baseline method",
    )
    topn_challenger_dropdown = mo.ui.dropdown(
        options=_methods,
        value="TACTICS TT-TS",
        label="Challenger method",
    )

    mo.hstack(
        [topn_baseline_dropdown, topn_challenger_dropdown],
        justify="start",
        gap=2,
    )
    return topn_baseline_dropdown, topn_challenger_dropdown


@app.cell
def _(
    fig_to_pdf_bytes,
    mo,
    np,
    pl,
    plt,
    stats,
    topn_baseline_dropdown,
    topn_challenger_dropdown,
    topn_df,
    topn_library_dropdown,
):
    """Bars showing Δ % recovery at each top-N threshold with significance stars."""

    _lib_id = topn_library_dropdown.value
    _baseline = topn_baseline_dropdown.value
    _challenger = topn_challenger_dropdown.value
    _lib_data = topn_df.filter(pl.col("library_id") == _lib_id)
    _thresholds = sorted(_lib_data["top_n"].unique().to_list())

    _topn_vals = []
    _mean_diffs = []
    _ci_lows = []
    _ci_highs = []
    _sigs = []

    for _tn in _thresholds:
        _tn_data = _lib_data.filter(pl.col("top_n") == _tn)
        _a = (
            _tn_data.filter(pl.col("method") == _baseline)
            .sort(["query_id", "replicate"])["recovery_frac"]
            .to_numpy()
        )
        _b = (
            _tn_data.filter(pl.col("method") == _challenger)
            .sort(["query_id", "replicate"])["recovery_frac"]
            .to_numpy()
        )
        _n = min(len(_a), len(_b))
        if _n < 2:
            continue
        _a, _b = _a[:_n], _b[:_n]
        _diffs = (_b - _a) * 100
        _mean_d = float(np.mean(_diffs))
        _se_d = float(np.std(_diffs, ddof=1) / np.sqrt(_n))
        _t_crit = float(stats.t.ppf(0.975, df=_n - 1))
        _, _p = stats.ttest_rel(_a, _b)
        _sig = (
            "****" if _p <= 0.0001
            else "***" if _p <= 0.001
            else "**" if _p <= 0.01
            else "*" if _p <= 0.05
            else "ns"
        )

        _topn_vals.append(_tn)
        _mean_diffs.append(_mean_d)
        _ci_lows.append(_mean_d - _t_crit * _se_d)
        _ci_highs.append(_mean_d + _t_crit * _se_d)
        _sigs.append(_sig)

    _x_vals = np.array(_topn_vals)
    _y = np.array(_mean_diffs)
    _lo = np.array(_ci_lows)
    _hi = np.array(_ci_highs)

    _x = np.arange(len(_x_vals))
    _bar_width = 0.6
    _color = "#4C78A8"

    _fig, _ax = plt.subplots(figsize=(10, 6), facecolor="white")
    _ax.set_facecolor("white")

    # Bars with 95% CI error bars
    _yerr = np.array([_y - _lo, _hi - _y])
    _ax.bar(_x, _y, width=_bar_width, color=_color, alpha=0.75,
            edgecolor="white", linewidth=0.5, zorder=3,
            yerr=_yerr, capsize=4,
            error_kw={"linewidth": 1.2, "color": "black"})

    _ax.axhline(0, color="black", linewidth=0.8, zorder=2)

    # Significance stars
    _data_range = max(abs(_hi.max()), abs(_lo.min()), 1)
    _offset = _data_range * 0.04
    for _i in range(len(_x)):
        _y_pos = _hi[_i] + _offset if _y[_i] >= 0 else _lo[_i] - _offset
        _va = "bottom" if _y[_i] >= 0 else "top"
        _ax.text(_x[_i], _y_pos, _sigs[_i], ha="center", va=_va,
                 fontsize=12, fontweight="bold", color="black")

    # Expand y-limits for annotations
    _y_min = min(_lo.min(), _y.min()) - _data_range * 0.15
    _y_max = max(_hi.max(), _y.max()) + _data_range * 0.15
    _ax.set_ylim(_y_min, _y_max)

    _ax.set_xticks(_x)
    _ax.set_xticklabels([str(v) for v in _x_vals], fontsize=12, color="black")
    _ax.set_xlabel("Top-N Threshold", fontsize=13, color="black")
    _ax.set_ylabel("\u0394 % Recovery", fontsize=13, color="black")
    _ax.set_title(
        f"Recovery Difference vs Top-N: {_challenger} \u2212 {_baseline} ({_lib_id})",
        fontsize=14, fontweight="bold", color="black",
    )
    _ax.tick_params(colors="black", labelsize=12)
    _ax.grid(axis="y", alpha=0.3, color="grey", zorder=1)
    for _spine in _ax.spines.values():
        _spine.set_color("black")

    _fig.tight_layout()

    _pdf_name = f"recovery_diff_vs_topn_{_challenger}_{_baseline}_{_lib_id}.pdf".replace(" ", "_")
    mo.vstack([
        plt.gca(),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename=_pdf_name,
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


@app.cell
def _(
    fig_to_pdf_bytes,
    mo,
    np,
    pl,
    plt,
    stats,
    topn_baseline_dropdown,
    topn_challenger_dropdown,
    topn_df,
):
    """Component-aggregated Δ % recovery vs top-N (2-comp and 3-comp panels)."""

    _baseline = topn_baseline_dropdown.value
    _challenger = topn_challenger_dropdown.value

    _fig, _axes = plt.subplots(1, 2, figsize=(16, 6), facecolor="white",
                               sharey=True)

    _global_ymin = 0.0
    _global_ymax = 0.0
    _panel_annots = []  # store (sigs, y, hi, lo) per panel for second pass

    for _panel_idx, _nc in enumerate([2, 3]):
        _ax = _axes[_panel_idx]
        _ax.set_facecolor("white")

        _comp_data = topn_df.filter(pl.col("n_components") == _nc)
        _thresholds = sorted(_comp_data["top_n"].unique().to_list())

        _topn_vals = []
        _mean_diffs = []
        _ci_lows = []
        _ci_highs = []
        _sigs = []

        for _tn in _thresholds:
            _tn_data = _comp_data.filter(pl.col("top_n") == _tn)
            # Pool all libraries: pair by (library_id, query_id, replicate)
            _a = (
                _tn_data.filter(pl.col("method") == _baseline)
                .sort(["library_id", "query_id", "replicate"])["recovery_frac"]
                .to_numpy()
            )
            _b = (
                _tn_data.filter(pl.col("method") == _challenger)
                .sort(["library_id", "query_id", "replicate"])["recovery_frac"]
                .to_numpy()
            )
            _n = min(len(_a), len(_b))
            if _n < 2:
                continue
            _a, _b = _a[:_n], _b[:_n]
            _diffs = (_b - _a) * 100
            _mean_d = float(np.mean(_diffs))
            _se_d = float(np.std(_diffs, ddof=1) / np.sqrt(_n))
            _t_crit = float(stats.t.ppf(0.975, df=_n - 1))
            _, _p = stats.ttest_rel(_a, _b)
            _sig = (
                "****" if _p <= 0.0001
                else "***" if _p <= 0.001
                else "**" if _p <= 0.01
                else "*" if _p <= 0.05
                else "ns"
            )

            _topn_vals.append(_tn)
            _mean_diffs.append(_mean_d)
            _ci_lows.append(_mean_d - _t_crit * _se_d)
            _ci_highs.append(_mean_d + _t_crit * _se_d)
            _sigs.append(_sig)

        _y = np.array(_mean_diffs)
        _lo = np.array(_ci_lows)
        _hi = np.array(_ci_highs)
        _x = np.arange(len(_topn_vals))

        _color = "#4C78A8" if _nc == 2 else "#E45756"
        _yerr = np.array([_y - _lo, _hi - _y])
        _ax.bar(_x, _y, width=0.6, color=_color, alpha=0.75,
                edgecolor="white", linewidth=0.5, zorder=3,
                yerr=_yerr, capsize=4,
                error_kw={"linewidth": 1.2, "color": "black"})

        _ax.axhline(0, color="black", linewidth=0.8, zorder=2)

        _ax.set_xticks(_x)
        _ax.set_xticklabels([str(v) for v in _topn_vals], fontsize=12, color="black")
        _ax.set_xlabel("Top-N Threshold", fontsize=13, color="black")
        _ax.set_title(f"{_nc}-Component Libraries", fontsize=14,
                      fontweight="bold", color="black")
        _ax.tick_params(colors="black", labelsize=12)
        _ax.grid(axis="y", alpha=0.3, color="grey", zorder=1)
        for _spine in _ax.spines.values():
            _spine.set_color("black")

        _global_ymin = min(_global_ymin, _lo.min() if len(_lo) else 0)
        _global_ymax = max(_global_ymax, _hi.max() if len(_hi) else 0)

        _panel_annots.append((_sigs, _y, _hi, _lo))

    # Second pass: uniform y-limits and significance stars
    _data_range = max(abs(_global_ymax), abs(_global_ymin), 1)
    _offset = _data_range * 0.04
    _y_min = _global_ymin - _data_range * 0.20
    _y_max = _global_ymax + _data_range * 0.20

    for _panel_idx in range(2):
        _ax = _axes[_panel_idx]
        _ax.set_ylim(_y_min, _y_max)
        _sigs, _y_arr, _hi_arr, _lo_arr = _panel_annots[_panel_idx]
        for _i in range(len(_sigs)):
            _y_pos = _hi_arr[_i] + _offset if _y_arr[_i] >= 0 else _lo_arr[_i] - _offset
            _va = "bottom" if _y_arr[_i] >= 0 else "top"
            _ax.text(_i, _y_pos, _sigs[_i], ha="center", va=_va,
                     fontsize=12, fontweight="bold", color="black")

    _axes[0].set_ylabel("\u0394 % Recovery", fontsize=13, color="black")

    _fig.suptitle(
        f"Mean Recovery Difference vs Top-N: {_challenger} \u2212 {_baseline}",
        fontsize=15, fontweight="bold", color="black", y=1.01,
    )
    _fig.tight_layout()

    _pdf_name = f"component_agg_diff_vs_topn_{_challenger}_{_baseline}.pdf".replace(" ", "_")
    mo.vstack([
        plt.gcf(),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename=_pdf_name,
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


if __name__ == "__main__":
    app.run()
