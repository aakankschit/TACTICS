"""
TACTICS Manuscript Plots

Loads the manuscript benchmark results (109,050 trials, 5 methods, 21 libraries)
plus thrombin docking recovery data, and provides DataFrames ready for plotting.

Run as app: marimo run tutorials/manuscript_plots.py
Edit mode:  marimo edit tutorials/manuscript_plots.py
"""

import marimo

__generated_with = "0.21.1"
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
def _(mo):
    """Top-N slider for filtering pre-computed recovery data."""
    topn_slider = mo.ui.slider(
        start=25, stop=250, step=25, value=100,
        label="Top-N Recovery Threshold",
    )
    topn_slider
    return (topn_slider,)


@app.cell
def _(mo, pl, project_root, topn_slider):
    """Load pre-computed recovery and statistical analysis results."""
    analysis_dir = project_root / "outputs" / "full_benchmark" / "analysis"

    TOP_N = topn_slider.value

    # --- Recovery summary (all libraries × all top-N thresholds) ---
    recovery_all = pl.read_parquet(analysis_dir / "recovery_summary.parquet")

    # Filter to selected top-N threshold
    df = recovery_all.filter(pl.col("top_n") == TOP_N)

    # --- Statistical test results ---
    tukey_df = pl.read_parquet(analysis_dir / "tukey_hsd_results.parquet")
    tukey_at_n = tukey_df.filter(pl.col("top_n") == TOP_N)

    anova_df = pl.read_parquet(analysis_dir / "anova_results.parquet")
    anova_at_n = anova_df.filter(pl.col("top_n") == TOP_N)

    n_trials = len(df)
    n_libraries = df["library_id"].n_unique()
    libraries = sorted(df["library_id"].unique().to_list())
    methods = sorted(df["method"].unique().to_list())

    mo.md(
        f"""
        ## Data Loaded (Pre-Computed Post-Hoc Analysis)

        - **Source**: `{analysis_dir}`
        - **Top-N threshold**: {TOP_N} (use slider to change)
        - **Available thresholds**: {sorted(recovery_all['top_n'].unique().to_list())}
        - **Total trials at top-{TOP_N}**: {n_trials:,}
        - **Libraries**: {n_libraries} ({', '.join(libraries)})
        - **Methods**: {', '.join(methods)}
        - **Replicates**: {df['replicate'].n_unique()}
        - **Tukey HSD comparisons**: {len(tukey_at_n):,} ({tukey_at_n.filter(pl.col('reject')).height} significant)
        - **ANOVA tests**: {len(anova_at_n)} ({anova_at_n.filter(pl.col('significance') != 'ns').height} significant)
        """
    )
    return TOP_N, df, recovery_all, tukey_at_n


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
def _(mo, tukey_at_n):
    """Pairwise significance tests (Tukey HSD, pre-computed)."""
    mo.vstack([
        mo.md(
            """
            ## Pairwise Significance Tests (Tukey HSD)

            FWER-controlled pairwise comparisons.
            Following Ash et al. 2025 (JCIM) guidelines.
            """
        ),
        mo.ui.table(
            tukey_at_n.select([
                "library_id", "n_components", "group1", "group2",
                "mean_diff", "p_adj", "significance", "ci_low", "ci_high",
                "n_pairs",
            ]).to_pandas()
        ),
    ])
    return


@app.cell
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
        value="Legacy Enhanced-RWS",
        label="Baseline method",
    )
    challenger_dropdown = mo.ui.dropdown(
        options=_methods,
        value="TACTICS Enhanced-TT-TS (GMIC)",
        label="Challenger method",
    )

    mo.hstack(
        [baseline_dropdown, challenger_dropdown],
        justify="start",
        gap=2,
    )
    return baseline_dropdown, challenger_dropdown


@app.cell
def _(baseline_dropdown, challenger_dropdown, mo, pl, tukey_at_n):
    """Per-library paired recovery differences from Tukey HSD."""
    _baseline = baseline_dropdown.value
    _challenger = challenger_dropdown.value

    # Filter Tukey results for this comparison (either direction)
    _match = tukey_at_n.filter(
        ((pl.col("group1") == _baseline) & (pl.col("group2") == _challenger))
        | ((pl.col("group1") == _challenger) & (pl.col("group2") == _baseline))
    )

    # Normalize direction: positive mean_diff = challenger better
    diff_rows = []
    for _row in _match.iter_rows(named=True):
        if _row["group1"] == _baseline:
            diff_rows.append({
                "library_id": _row["library_id"],
                "n_components": _row["n_components"],
                "mean_diff": _row["mean_diff"],
                "ci_low": _row["ci_low"],
                "ci_high": _row["ci_high"],
                "p_value": _row["p_adj"],
                "significance": _row["significance"],
                "n_pairs": _row["n_pairs"],
            })
        else:
            diff_rows.append({
                "library_id": _row["library_id"],
                "n_components": _row["n_components"],
                "mean_diff": -_row["mean_diff"],
                "ci_low": -_row["ci_high"],
                "ci_high": -_row["ci_low"],
                "p_value": _row["p_adj"],
                "significance": _row["significance"],
                "n_pairs": _row["n_pairs"],
            })

    diff_df = pl.DataFrame(diff_rows).sort("mean_diff", descending=True)

    mo.vstack([
        mo.md(f"## Paired Recovery Differences: {_challenger} vs {_baseline} (Tukey HSD)"),
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
        "\u0394 % Recovery of Top-N Hits",
        fontsize=14, color="black",
    )
    _ax.set_title(
        f"Per-Library Recovery of Top-N Hit Molecules: {_challenger} vs {_baseline}",
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
        _lib_methods = set(_lib_data["method"].unique().to_list())

        # Skip libraries missing required methods for decomposition
        _required = {"Legacy Enhanced-RWS", "TACTICS Enhanced-RWS (GMIC)", "TACTICS Enhanced-TT-TS (GMIC)"}
        if not _required.issubset(_lib_methods):
            continue

        _lrws = float(_lib_data.filter(pl.col("method") == "Legacy Enhanced-RWS")["recovered"].mean())
        _trws = float(_lib_data.filter(pl.col("method") == "TACTICS Enhanced-RWS (GMIC)")["recovered"].mean())
        _ttts = float(_lib_data.filter(pl.col("method") == "TACTICS Enhanced-TT-TS (GMIC)")["recovered"].mean())

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
    fig_to_pdf_bytes,
    gain_decomposition,
    mo,
    np,
    pl,
    plt,
    tukey_at_n,
):
    """Diverging stacked bar: criticality vs TT-TS contributions by component count."""

    # Per-library significance for TT-TS vs Legacy Enhanced-RWS (from Tukey HSD)
    _sig_map = {}
    for _lib_id in gain_decomposition["library_id"].to_list():
        _match = tukey_at_n.filter(
            (pl.col("library_id") == _lib_id)
            & (
                ((pl.col("group1") == "Legacy Enhanced-RWS") & (pl.col("group2") == "TACTICS Enhanced-TT-TS (GMIC)"))
                | ((pl.col("group1") == "TACTICS Enhanced-TT-TS (GMIC)") & (pl.col("group2") == "Legacy Enhanced-RWS"))
            )
        )
        if len(_match) > 0:
            _sig_map[_lib_id] = _match["significance"].to_list()[0]
        else:
            _sig_map[_lib_id] = "ns"

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

    _axes[0].set_ylabel("\u0394 % Recovery of Top-N Hits\n(TACTICS TT-TS \u2212 Legacy Enhanced-RWS)",
                        fontsize=14, color="black")

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
        | `df` | Raw benchmark results (109,090 rows, 22 libraries incl. thrombin-docking) |
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


@app.cell
def _(mo):
    mo.md(r"""
    # Sensitivity Analysis - Top N compounds
    In this section we explore the roubustness of CATS as the interval for the top N compounds increases, if CATS is interpreting the SAR landscape correctly and systematically able to uncover mid tier hits.
    """)
    return


@app.cell
def _(mo, recovery_all):
    """Summary of recovery data available for sensitivity analysis."""
    _n_rows = len(recovery_all)
    _n_libs = recovery_all["library_id"].n_unique()
    _thresholds = sorted(recovery_all["top_n"].unique().to_list())
    _methods = sorted(recovery_all["method"].unique().to_list())

    mo.vstack([
        mo.md(
            f"""
            ## Top-N Sensitivity Data (from Post-Hoc Analysis)

            - **Rows**: {_n_rows:,}
            - **Libraries**: {_n_libs}
            - **Top-N thresholds**: {_thresholds}
            - **Methods**: {', '.join(_methods)}
            """
        ),
    ])
    return


@app.cell
def _(fig_to_pdf_bytes, mo, pl, plt, recovery_all):
    """Overall % recovery vs top-N by component count (ROCS libraries, 10 thresholds)."""

    _DOCK_LIBS = {"thrombin", "adenine", "amide", "quinazoline"}
    _ROCS_NCOMP = {
        "rxn101": 2, "rxn102a": 2, "rxn108b": 2, "rxn111b": 2, "rxn114b": 2,
        "rxn203": 2, "rxn205": 2, "rxn206": 2, "rxn207": 2, "rxn208": 2,
        "amide-suzuki": 3, "betti": 3, "dobener": 3,
        "groebke-blackburn-bienayme": 3, "mannich": 3, "niementowski": 3,
        "orru": 3, "passerini": 3, "petasis": 3, "poparov": 3,
    }

    _rocs_all = (
        recovery_all
        .filter(~pl.col("library_id").is_in(_DOCK_LIBS))
        .filter(pl.col("library_id").is_in(set(_ROCS_NCOMP.keys())))
        .with_columns(
            pl.col("library_id").replace(_ROCS_NCOMP).cast(pl.Int64).alias("n_components")
        )
    )

    mo.stop(len(_rocs_all) == 0, mo.md("No ROCS library data available."))

    _method_colors = {
        "TACTICS Enhanced-TT-TS (GMIC)": "#E45756",
        "TACTICS Enhanced-RWS (GMIC)": "#4C78A8",
        "TACTICS Balanced-Greedy": "#F28E2B",
        "Legacy Standard-Greedy": "#76B7B2",
        "Legacy Enhanced-RWS": "#B07AA1",
    }
    _method_markers = {
        "TACTICS Enhanced-TT-TS (GMIC)": "D",
        "TACTICS Enhanced-RWS (GMIC)": "o",
        "TACTICS Balanced-Greedy": "P",
        "Legacy Standard-Greedy": "s",
        "Legacy Enhanced-RWS": "^",
    }
    _method_order = [
        "TACTICS Enhanced-TT-TS (GMIC)",
        "TACTICS Enhanced-RWS (GMIC)",
        "TACTICS Balanced-Greedy",
        "Legacy Standard-Greedy",
        "Legacy Enhanced-RWS",
    ]

    _fig, _axes = plt.subplots(1, 2, figsize=(16, 6), facecolor="white",
                               sharey=True)

    for _panel_idx, _nc in enumerate([2, 3]):
        _ax = _axes[_panel_idx]
        _ax.set_facecolor("white")

        _comp_data = _rocs_all.filter(pl.col("n_components") == _nc)

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

        for _method in _method_order:
            _mdata = _agg.filter(pl.col("method") == _method).sort("top_n")
            if len(_mdata) == 0:
                continue
            _x = _mdata["top_n"].to_numpy()
            _y = _mdata["mean_frac"].to_numpy() * 100
            _se = _mdata["se_frac"].to_numpy() * 100
            _color = _method_colors.get(_method, "#999999")
            _marker = _method_markers.get(_method, "o")

            _ax.plot(_x, _y, marker=_marker, color=_color, linewidth=2,
                     markersize=7, label=_method, zorder=3)
            _ax.fill_between(_x, _y - 1.96 * _se, _y + 1.96 * _se,
                             alpha=0.15, color=_color, zorder=2)

        _ax.set_xlabel("Top-N Threshold", fontsize=14, color="black")
        _ax.set_title(f"{_nc}-Component Libraries", fontsize=15,
                      fontweight="bold", color="black")
        _ax.set_xticks(sorted(_agg["top_n"].unique().to_list()))
        _ax.tick_params(colors="black", labelsize=12)
        _ax.grid(axis="both", alpha=0.3, color="grey", zorder=1)
        for _spine in _ax.spines.values():
            _spine.set_color("black")

    _axes[0].set_ylabel("% Recovery", fontsize=14, color="black")

    # Legend outside the figure below the panels
    _handles, _labels = _axes[0].get_legend_handles_labels()
    _fig.legend(
        _handles, _labels,
        loc="lower center", ncol=3, fontsize=12,
        facecolor="white", edgecolor="black", labelcolor="black",
        bbox_to_anchor=(0.5, -0.08),
    )

    _fig.suptitle("Mean % Recovery vs Top-N Threshold (ROCS Libraries)",
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
def _(fig_to_pdf_bytes, mo, pl, plt, recovery_all):
    """Mean % recovery vs top-N for docking libraries (2-comp and 3-comp panels)."""

    _DOCK_LIBS = {"adenine", "amide", "quinazoline", "thrombin"}
    # n_components is 0 for thrombin/quinazoline in the data — fix with known values
    _DOCK_NCOMP = {"thrombin": 2, "amide": 2, "adenine": 3, "quinazoline": 3}

    _dock_all = (
        recovery_all
        .filter(pl.col("library_id").is_in(_DOCK_LIBS))
        .with_columns(
            pl.col("library_id").replace(_DOCK_NCOMP).cast(pl.Int64).alias("n_components")
        )
    )

    mo.stop(len(_dock_all) == 0, mo.md("No docking library data available."))

    _method_colors = {
        "TACTICS Enhanced-TT-TS (GMIC)": "#E45756",
        "TACTICS Enhanced-RWS (GMIC)": "#4C78A8",
        "TACTICS Balanced-Greedy": "#F28E2B",
        "Legacy Standard-Greedy": "#76B7B2",
        "Legacy Enhanced-RWS": "#B07AA1",
    }
    _method_markers = {
        "TACTICS Enhanced-TT-TS (GMIC)": "D",
        "TACTICS Enhanced-RWS (GMIC)": "o",
        "TACTICS Balanced-Greedy": "P",
        "Legacy Standard-Greedy": "s",
        "Legacy Enhanced-RWS": "^",
    }
    _method_order = [
        "TACTICS Enhanced-TT-TS (GMIC)",
        "TACTICS Enhanced-RWS (GMIC)",
        "TACTICS Balanced-Greedy",
        "Legacy Standard-Greedy",
        "Legacy Enhanced-RWS",
    ]

    _fig, _axes = plt.subplots(1, 2, figsize=(16, 6), facecolor="white",
                               sharey=True)

    for _panel_idx, _nc in enumerate([2, 3]):
        _ax = _axes[_panel_idx]
        _ax.set_facecolor("white")

        _comp_data = _dock_all.filter(pl.col("n_components") == _nc)

        # Aggregate: mean and SE of recovery fraction across libraries and replicates
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

        for _method in _method_order:
            _mdata = _agg.filter(pl.col("method") == _method).sort("top_n")
            if len(_mdata) == 0:
                continue
            _x = _mdata["top_n"].to_numpy()
            _y = _mdata["mean_frac"].to_numpy() * 100
            _se = _mdata["se_frac"].to_numpy() * 100
            _color = _method_colors.get(_method, "#999999")
            _marker = _method_markers.get(_method, "o")

            _ax.plot(_x, _y, marker=_marker, color=_color, linewidth=2,
                     markersize=7, label=_method, zorder=3)
            _ax.fill_between(_x, _y - 1.96 * _se, _y + 1.96 * _se,
                             alpha=0.15, color=_color, zorder=2)

        _libs_in_panel = sorted(
            _comp_data["library_id"].unique().to_list()
        )
        _ax.set_xlabel("Top-N Threshold", fontsize=14, color="black")
        _ax.set_title(
            f"{_nc}-Component Docking Libraries\n({', '.join(_libs_in_panel)})",
            fontsize=15, fontweight="bold", color="black",
        )
        _ax.set_xticks(sorted(_agg["top_n"].unique().to_list()))
        _ax.tick_params(colors="black", labelsize=12)
        _ax.grid(axis="both", alpha=0.3, color="grey", zorder=1)
        for _spine in _ax.spines.values():
            _spine.set_color("black")

    _axes[0].set_ylabel("% Recovery", fontsize=14, color="black")

    # Legend outside the figure below the panels
    _handles, _labels = _axes[0].get_legend_handles_labels()
    _fig.legend(
        _handles, _labels,
        loc="lower center", ncol=3, fontsize=12,
        facecolor="white", edgecolor="black", labelcolor="black",
        bbox_to_anchor=(0.5, -0.08),
    )

    _fig.suptitle("Mean % Recovery vs Top-N Threshold (Docking Libraries)",
                  fontsize=15, fontweight="bold", color="black", y=1.01)
    _fig.tight_layout()

    mo.vstack([
        plt.gcf(),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename="docking_recovery_vs_topn.pdf",
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


@app.cell
def _(mo, recovery_all):
    """Select library for top-N recovery curve (all libraries)."""
    _all_libs = sorted(recovery_all["library_id"].unique().to_list())

    topn_library_dropdown = mo.ui.dropdown(
        options=_all_libs,
        value=_all_libs[0],
        label="Library",
    )
    topn_library_dropdown
    return (topn_library_dropdown,)


@app.cell
def _(fig_to_pdf_bytes, mo, pl, plt, recovery_all, topn_library_dropdown):
    """Top-N recovery curve: % recovery vs top-N threshold per method."""

    _lib_id = topn_library_dropdown.value
    _lib_data = recovery_all.filter(pl.col("library_id") == _lib_id)

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
        "TACTICS Enhanced-TT-TS (GMIC)": "#E45756",
        "TACTICS Enhanced-RWS (GMIC)": "#4C78A8",
        "TACTICS Balanced-Greedy": "#F28E2B",
        "Legacy Standard-Greedy": "#76B7B2",
        "Legacy Enhanced-RWS": "#B07AA1",
    }
    _method_markers = {
        "TACTICS Enhanced-TT-TS (GMIC)": "D",
        "TACTICS Enhanced-RWS (GMIC)": "o",
        "TACTICS Balanced-Greedy": "P",
        "Legacy Standard-Greedy": "s",
        "Legacy Enhanced-RWS": "^",
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

    _ax.set_xlabel("Top-N Threshold", fontsize=14, color="black")
    _ax.set_ylabel("% Recovery", fontsize=14, color="black")
    _ax.set_title(
        f"% Recovery vs Top-N Threshold: {_lib_id}",
        fontsize=15, fontweight="bold", color="black",
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
def _(mo, recovery_all):
    """Select baseline and challenger methods for top-N difference plot."""
    _all_methods = sorted(recovery_all["method"].unique().to_list())

    _default_baseline = "Legacy Enhanced-RWS" if "Legacy Enhanced-RWS" in _all_methods else _all_methods[0]
    _default_challenger = "TACTICS Enhanced-TT-TS (GMIC)" if "TACTICS Enhanced-TT-TS (GMIC)" in _all_methods else _all_methods[-1]

    topn_baseline_dropdown = mo.ui.dropdown(
        options=_all_methods,
        value=_default_baseline,
        label="Baseline method",
    )
    topn_challenger_dropdown = mo.ui.dropdown(
        options=_all_methods,
        value=_default_challenger,
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
    recovery_all,
    stats,
    topn_baseline_dropdown,
    topn_challenger_dropdown,
    topn_library_dropdown,
):
    """Bars showing Δ % recovery at each top-N threshold with significance stars."""

    _lib_id = topn_library_dropdown.value
    _baseline = topn_baseline_dropdown.value
    _challenger = topn_challenger_dropdown.value

    _lib_data = recovery_all.filter(pl.col("library_id") == _lib_id)

    _available_methods = set(_lib_data["method"].unique().to_list())
    mo.stop(
        _baseline not in _available_methods or _challenger not in _available_methods,
        mo.md(
            f"Methods **{_baseline}** and/or **{_challenger}** not found in `{_lib_id}`. "
            f"Available: {sorted(_available_methods)}"
        ),
    )

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
    _ax.set_xlabel("Top-N Threshold", fontsize=14, color="black")
    _ax.set_ylabel("\u0394 % Recovery", fontsize=14, color="black")
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
    recovery_all,
    stats,
    topn_baseline_dropdown,
    topn_challenger_dropdown,
):
    """Component-aggregated Δ % recovery vs top-N (2-comp and 3-comp panels, ROCS libraries)."""

    _baseline = topn_baseline_dropdown.value
    _challenger = topn_challenger_dropdown.value

    _ROCS_NCOMP = {
        "rxn101": 2, "rxn102a": 2, "rxn108b": 2, "rxn111b": 2, "rxn114b": 2,
        "rxn203": 2, "rxn205": 2, "rxn206": 2, "rxn207": 2, "rxn208": 2,
        "amide-suzuki": 3, "betti": 3, "dobener": 3,
        "groebke-blackburn-bienayme": 3, "mannich": 3, "niementowski": 3,
        "orru": 3, "passerini": 3, "petasis": 3, "poparov": 3,
    }

    _rocs_data = (
        recovery_all
        .filter(pl.col("library_id").is_in(set(_ROCS_NCOMP.keys())))
        .with_columns(
            pl.col("library_id").replace(_ROCS_NCOMP).cast(pl.Int64).alias("n_components")
        )
    )

    _fig, _axes = plt.subplots(1, 2, figsize=(16, 6), facecolor="white",
                               sharey=True)

    _global_ymin = 0.0
    _global_ymax = 0.0
    _panel_annots = []  # store (sigs, y, hi, lo) per panel for second pass

    for _panel_idx, _nc in enumerate([2, 3]):
        _ax = _axes[_panel_idx]
        _ax.set_facecolor("white")

        _comp_data = _rocs_data.filter(pl.col("n_components") == _nc)
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
        _ax.set_xlabel("Top-N Threshold", fontsize=14, color="black")
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

    _axes[0].set_ylabel("\u0394 % Recovery", fontsize=14, color="black")

    _fig.suptitle(
        f"Mean Recovery Difference vs Top-N (ROCS): {_challenger} \u2212 {_baseline}",
        fontsize=15, fontweight="bold", color="black", y=1.01,
    )
    _fig.tight_layout()

    _pdf_name = f"component_agg_diff_vs_topn_rocs_{_challenger}_{_baseline}.pdf".replace(" ", "_")
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


@app.cell
def _(
    fig_to_pdf_bytes,
    mo,
    np,
    pl,
    plt,
    recovery_all,
    stats,
    topn_baseline_dropdown,
    topn_challenger_dropdown,
):
    """Component-aggregated Δ % recovery vs top-N (2-comp and 3-comp panels, docking libraries)."""

    _baseline = topn_baseline_dropdown.value
    _challenger = topn_challenger_dropdown.value

    _DOCK_NCOMP = {"thrombin": 2, "amide": 2, "adenine": 3, "quinazoline": 3}

    _dock_data = (
        recovery_all
        .filter(pl.col("library_id").is_in(set(_DOCK_NCOMP.keys())))
        .with_columns(
            pl.col("library_id").replace(_DOCK_NCOMP).cast(pl.Int64).alias("n_components")
        )
    )

    _available_methods = set(_dock_data["method"].unique().to_list())
    mo.stop(
        _baseline not in _available_methods or _challenger not in _available_methods,
        mo.md(
            f"Methods **{_baseline}** and/or **{_challenger}** not found in docking data. "
            f"Available: {sorted(_available_methods)}"
        ),
    )

    _fig, _axes = plt.subplots(1, 2, figsize=(16, 6), facecolor="white",
                               sharey=True)

    _global_ymin = 0.0
    _global_ymax = 0.0
    _panel_annots = []

    for _panel_idx, _nc in enumerate([2, 3]):
        _ax = _axes[_panel_idx]
        _ax.set_facecolor("white")

        _comp_data = _dock_data.filter(pl.col("n_components") == _nc)
        _thresholds = sorted(_comp_data["top_n"].unique().to_list())

        _libs_in_panel = sorted(_comp_data["library_id"].unique().to_list())

        _topn_vals = []
        _mean_diffs = []
        _ci_lows = []
        _ci_highs = []
        _sigs = []

        for _tn in _thresholds:
            _tn_data = _comp_data.filter(pl.col("top_n") == _tn)
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
        _ax.set_xlabel("Top-N Threshold", fontsize=14, color="black")
        _ax.set_title(
            f"{_nc}-Component Docking\n({', '.join(_libs_in_panel)})",
            fontsize=15, fontweight="bold", color="black",
        )
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

    _axes[0].set_ylabel("\u0394 % Recovery", fontsize=14, color="black")

    _fig.suptitle(
        f"Mean Recovery Difference vs Top-N (Docking): {_challenger} \u2212 {_baseline}",
        fontsize=15, fontweight="bold", color="black", y=1.01,
    )
    _fig.tight_layout()

    _pdf_name = f"component_agg_diff_vs_topn_docking_{_challenger}_{_baseline}.pdf".replace(" ", "_")
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


@app.cell
def _(mo):
    mo.md(r"""
    ## Simultaneous Confidence Interval Plots (Tukey HSD)

    Per-library: select a library to view all pairwise method comparisons.
    Overall: aggregated across ROCS libraries and docking libraries separately.
    """)
    return


@app.cell
def _(df, mo):
    """Library selector for per-library Tukey HSD CI plot."""
    _libraries = sorted(df["library_id"].unique().to_list())
    library_dropdown = mo.ui.dropdown(
        options=_libraries,
        value=_libraries[0] if _libraries else "",
        label="Library",
    )
    library_dropdown
    return (library_dropdown,)


@app.cell
def _(
    Line2D,
    TOP_N,
    fig_to_pdf_bytes,
    library_dropdown,
    mo,
    np,
    pl,
    plt,
    tukey_at_n,
):
    """Option 4 (per-library): Simultaneous CI plot for one library."""
    _lib_id = library_dropdown.value
    _lib_tukey = tukey_at_n.filter(pl.col("library_id") == _lib_id)

    mo.stop(len(_lib_tukey) == 0, mo.md(f"No Tukey HSD results for {_lib_id} at this top-N."))

    _n_comparisons = len(_lib_tukey)

    _left_labels = []
    _right_labels = []
    _means = []
    _ci_lows = []
    _ci_highs = []
    _rejects = []

    # Normalize: mean_diff = group2 - group1 from Tukey HSD.
    # Ensure the better method (higher recovery) is always the challenger (right).
    # If mean_diff > 0, group2 is better → group1=baseline(left), group2=challenger(right).
    # If mean_diff < 0, group1 is better → flip: group2=baseline(left), group1=challenger(right).
    for _row in _lib_tukey.iter_rows(named=True):
        if _row["mean_diff"] >= 0:
            _left_labels.append(_row["group1"])
            _right_labels.append(_row["group2"])
            _means.append(_row["mean_diff"])
            _ci_lows.append(_row["ci_low"])
            _ci_highs.append(_row["ci_high"])
        else:
            _left_labels.append(_row["group2"])
            _right_labels.append(_row["group1"])
            _means.append(-_row["mean_diff"])
            _ci_lows.append(-_row["ci_high"])
            _ci_highs.append(-_row["ci_low"])
        _rejects.append(_row["reject"])

    # Sort by mean_diff descending after normalization
    _order = np.argsort(_means)[::-1]
    _left_labels = [_left_labels[i] for i in _order]
    _right_labels = [_right_labels[i] for i in _order]
    _rejects = [_rejects[i] for i in _order]
    _means = np.array([_means[i] for i in _order])
    _ci_lows = np.array([_ci_lows[i] for i in _order])
    _ci_highs = np.array([_ci_highs[i] for i in _order])

    _means = np.array(_means)
    _ci_lows = np.array(_ci_lows)
    _ci_highs = np.array(_ci_highs)

    # Extract p-values for significance stars
    _p_values = []
    for _row in _lib_tukey.iter_rows(named=True):
        _p_values.append(_row["p_adj"])
    _p_values = np.array([_p_values[i] for i in _order])

    _fig, _ax = plt.subplots(
        figsize=(12, max(4, _n_comparisons * 0.55)),
        facecolor="white",
    )
    _ax.set_facecolor("white")
    _y = np.arange(_n_comparisons)

    _colors = ["#4C78A8" if r else "#999999" for r in _rejects]

    for _i in range(_n_comparisons):
        _ax.plot(
            [_ci_lows[_i], _ci_highs[_i]], [_i, _i],
            color=_colors[_i], linewidth=2.5, solid_capstyle="round", zorder=3,
        )
        _ax.plot(_means[_i], _i, "o", color=_colors[_i], markersize=8, zorder=4)
        # Significance stars
        _p = _p_values[_i]
        _sig_str = (
            "****" if _p <= 0.0001
            else "***" if _p <= 0.001
            else "**" if _p <= 0.01
            else "*" if _p <= 0.05
            else "ns"
        )
        _ax.text(
            _ci_highs[_i] + 0.15, _i, _sig_str,
            ha="left", va="center",
            fontsize=12, fontweight="bold", color="black",
        )

    _ax.axvline(0, color="black", linewidth=0.8, linestyle="-", zorder=2)

    # Expand x-limits so stars fit inside the axes
    _x_pad = max(_ci_highs.max(), abs(_ci_lows.min())) * 0.18
    _ax.set_xlim(
        min(_ci_lows.min(), 0) - _x_pad * 0.5,
        _ci_highs.max() + _x_pad,
    )

    # Define key comparisons
    _legacy_methods = {"Legacy Enhanced-RWS", "Legacy Standard-Greedy"}
    _tactics_adaptive = {"TACTICS Enhanced-RWS (GMIC)", "TACTICS Enhanced-TT-TS (GMIC)"}

    def _is_tactics_vs_legacy(_left, _right):
        return (
            (_left in _legacy_methods and _right in _tactics_adaptive)
            or (_left in _tactics_adaptive and _right in _legacy_methods)
        )

    _COLOR_KEY_TICK = "#4C78A8"

    # Left y-axis: challenger (better method)
    _ax.set_yticks(_y)
    _ax.set_yticklabels(_right_labels, fontsize=12, color="black")
    _ax.set_ylabel("Challenger", fontsize=14, color="black")

    # Pad y-limits so effect size labels above top row aren't clipped
    _ax.set_ylim(-0.6, _n_comparisons - 0.4)

    # Invert y so largest diff is at top (do this BEFORE twinx so both share orientation)
    _ax.invert_yaxis()

    # Right y-axis: baseline (worse method)
    _ax_right = _ax.twinx()
    _ax_right.set_ylim(_ax.get_ylim())
    _ax_right.set_yticks(_y)
    _ax_right.set_yticklabels(_left_labels, fontsize=12, color="black")
    _ax_right.set_ylabel("Baseline", fontsize=14, color="black")
    _ax_right.tick_params(colors="black", labelsize=12)
    for _spine in _ax_right.spines.values():
        _spine.set_color("black")

    # Per-tick styling: challenger red+bold, baseline bold for key comparisons
    for _i in range(_n_comparisons):
        if _is_tactics_vs_legacy(_left_labels[_i], _right_labels[_i]):
            _ax.get_yticklabels()[_i].set_color(_COLOR_KEY_TICK)
            _ax.get_yticklabels()[_i].set_fontweight("bold")
            _ax_right.get_yticklabels()[_i].set_fontweight("bold")

    _ax.set_xlabel(
        f"\u0394 Recovery: Challenger \u2212 Baseline (top-{TOP_N})",
        fontsize=14, color="black",
    )
    _ax.set_title(
        f"Tukey HSD Pairwise CIs: {_lib_id} (top-{TOP_N})",
        fontsize=15, fontweight="bold", color="black",
    )
    _ax.grid(axis="x", alpha=0.3, zorder=1, color="grey")
    _ax.tick_params(axis="x", colors="black", labelsize=12)
    _ax.tick_params(axis="y", color="black", labelsize=12)
    for _spine in _ax.spines.values():
        _spine.set_color("black")

    _legend_elements = [
        Line2D([0], [0], color="#4C78A8", linewidth=2.5, marker="o",
               markersize=8, label="Significant (p < 0.05)"),
        Line2D([0], [0], color="#999999", linewidth=2.5, marker="o",
               markersize=8, label="Not significant"),
    ]
    _ax.legend(handles=_legend_elements, loc="lower right", fontsize=12,
               facecolor="white", edgecolor="black", labelcolor="black")

    _fig.tight_layout()

    _pdf_name = f"tukey_ci_{_lib_id}_top{TOP_N}.pdf".replace(" ", "_")
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


@app.cell
def _(Line2D, TOP_N, df, fig_to_pdf_bytes, mo, np, pl, plt):
    """Across-library Tukey HSD pairwise CI plot (ROCS libraries only).

    Excludes docking libraries (thrombin, adenine, amide, quinazoline) which
    have different score distributions and only 1 query each. Docking results
    are shown in the separate docking Tukey HSD CI plot.

    Averages replicates within each (library, query) to get one independent
    observation per query, then runs statsmodels pairwise_tukeyhsd for
    FWER-controlled simultaneous CIs and adjusted p-values.

    Layout matches the per-library Tukey CI plot: positive values always mean
    the challenger (right y-axis) outperforms the baseline (left y-axis).
    """
    from statsmodels.stats.multicomp import pairwise_tukeyhsd

    # Exclude docking libraries — different score distributions, 1 query each
    _DOCKING_LIBS = {"thrombin", "adenine", "amide", "quinazoline"}
    _rocs_df = df.filter(~pl.col("library_id").is_in(_DOCKING_LIBS))

    # Average replicates: one observation per (library, query, method)
    _query_means = (
        _rocs_df.group_by(["library_id", "query_id", "method"])
        .agg(pl.col("recovered").mean().alias("mean_recovered"))
    )

    _values = _query_means["mean_recovered"].to_numpy()
    _groups = _query_means["method"].to_list()

    _tukey = pairwise_tukeyhsd(_values, _groups, alpha=0.05)
    _summary = _tukey.summary()

    # Parse Tukey HSD results into lists
    _left_labels = []   # baseline (worse method)
    _right_labels = []  # challenger (better method)
    _grand_means = []
    _ci_lows = []
    _ci_highs = []
    _p_values = []
    _rejects = []

    for _row in _summary.data[1:]:  # skip header row
        _g1, _g2, _meandiff, _p_adj, _lower, _upper, _reject = _row
        _meandiff = float(_meandiff)
        _p_adj = float(_p_adj)
        _lower = float(_lower)
        _upper = float(_upper)

        # Normalize direction: positive = challenger (right) is better.
        if _meandiff >= 0:
            _left_labels.append(_g1)
            _right_labels.append(_g2)
            _grand_means.append(_meandiff)
            _ci_lows.append(_lower)
            _ci_highs.append(_upper)
        else:
            _left_labels.append(_g2)
            _right_labels.append(_g1)
            _grand_means.append(-_meandiff)
            _ci_lows.append(-_upper)
            _ci_highs.append(-_lower)

        _p_values.append(_p_adj)
        _rejects.append(bool(_reject))

    _grand_means = np.array(_grand_means)
    _ci_lows = np.array(_ci_lows)
    _ci_highs = np.array(_ci_highs)
    _p_values = np.array(_p_values)

    # Sort by grand mean descending (largest positive diff at top)
    _order = np.argsort(_grand_means)[::-1]
    _left_labels = [_left_labels[i] for i in _order]
    _right_labels = [_right_labels[i] for i in _order]
    _rejects = [_rejects[i] for i in _order]
    _grand_means = _grand_means[_order]
    _ci_lows = _ci_lows[_order]
    _ci_highs = _ci_highs[_order]
    _p_values = _p_values[_order]

    _n_comparisons = len(_left_labels)
    _y = np.arange(_n_comparisons)

    _fig, _ax = plt.subplots(
        figsize=(12, max(4, _n_comparisons * 0.55)),
        facecolor="white",
    )
    _ax.set_facecolor("white")

    # Define the key TACTICS-vs-Legacy comparisons to highlight
    _legacy_methods = {"Legacy Enhanced-RWS", "Legacy Standard-Greedy"}
    _tactics_adaptive = {"TACTICS Enhanced-RWS (GMIC)", "TACTICS Enhanced-TT-TS (GMIC)"}

    def _is_tactics_vs_legacy(_left, _right):
        """True if one side is a legacy method and the other is a TACTICS adaptive method."""
        return (
            (_left in _legacy_methods and _right in _tactics_adaptive)
            or (_left in _tactics_adaptive and _right in _legacy_methods)
        )

    _COLOR_SIG = "#4C78A8"        # blue — significant pairs
    _COLOR_NS = "#999999"          # grey — not significant
    _COLOR_KEY_TICK = "#4C78A8"    # blue — TACTICS vs Legacy tick labels

    for _i in range(_n_comparisons):
        _color = _COLOR_SIG if _rejects[_i] else _COLOR_NS

        # CI line
        _ax.plot(
            [_ci_lows[_i], _ci_highs[_i]], [_i, _i],
            color=_color, linewidth=2.5, solid_capstyle="round", zorder=3,
        )
        # Point estimate
        _ax.plot(
            _grand_means[_i], _i, "o",
            color=_color, markersize=8, zorder=4,
        )
        # Significance stars (Tukey HSD adjusted p-values, FWER=0.05)
        _p = _p_values[_i]
        _sig_str = (
            "****" if _p <= 0.0001
            else "***" if _p <= 0.001
            else "**" if _p <= 0.01
            else "*" if _p <= 0.05
            else "ns"
        )
        _ax.text(
            _ci_highs[_i] + 0.15, _i, _sig_str,
            ha="left", va="center",
            fontsize=12, fontweight="bold", color="black",
        )

    _ax.axvline(0, color="black", linewidth=0.8, linestyle="-", zorder=2)

    # Expand x-limits so stars and CI lines fit inside the axes
    _x_pad = max(_ci_highs.max(), abs(_ci_lows.min())) * 0.18
    _ax.set_xlim(
        min(_ci_lows.min(), 0) - _x_pad * 0.5,
        _ci_highs.max() + _x_pad,
    )

    # Left y-axis: challenger (better method)
    _ax.set_yticks(_y)
    _ax.set_yticklabels(_right_labels, fontsize=12, color="black")
    _ax.set_ylabel("Challenger", fontsize=14, color="black")

    # Pad y-limits so effect size labels above top/bottom rows aren't clipped
    _ax.set_ylim(-0.6, _n_comparisons - 0.4)

    # Invert y so largest diff is at top (before twinx so both share orientation)
    _ax.invert_yaxis()

    # Right y-axis: baseline (worse method)
    _ax_right = _ax.twinx()
    _ax_right.set_ylim(_ax.get_ylim())
    _ax_right.set_yticks(_y)
    _ax_right.set_yticklabels(_left_labels, fontsize=12, color="black")
    _ax_right.set_ylabel("Baseline", fontsize=14, color="black")
    _ax_right.tick_params(colors="black", labelsize=12)
    for _spine in _ax_right.spines.values():
        _spine.set_color("black")

    # Per-tick styling: challenger red+bold, baseline bold for key comparisons
    # Applied after both axes are fully configured so nothing resets them
    for _i in range(_n_comparisons):
        if _is_tactics_vs_legacy(_left_labels[_i], _right_labels[_i]):
            _ax.get_yticklabels()[_i].set_color(_COLOR_KEY_TICK)
            _ax.get_yticklabels()[_i].set_fontweight("bold")
            _ax_right.get_yticklabels()[_i].set_fontweight("bold")

    _ax.set_xlabel(
        f"\u0394 Recovery: Challenger \u2212 Baseline (top-{TOP_N})",
        fontsize=14, color="black",
    )
    _n_rocs_libs = _rocs_df["library_id"].n_unique()
    _ax.set_title(
        f"Tukey HSD Pairwise CIs \u2014 {_n_rocs_libs} ROCS Libraries (top-{TOP_N}, FWER=0.05)",
        fontsize=15, fontweight="bold", color="black",
    )
    _ax.grid(axis="x", alpha=0.3, zorder=1, color="grey")
    _ax.tick_params(axis="x", colors="black", labelsize=12)
    _ax.tick_params(axis="y", color="black", labelsize=12)
    for _spine in _ax.spines.values():
        _spine.set_color("black")

    _legend_elements = [
        Line2D([0], [0], color=_COLOR_SIG, linewidth=2.5, marker="o",
               markersize=8, label="Significant (p < 0.05)"),
        Line2D([0], [0], color=_COLOR_NS, linewidth=2.5, marker="o",
               markersize=8, label="Not significant"),
    ]
    _ax.legend(handles=_legend_elements, loc="lower right", fontsize=12,
               facecolor="white", edgecolor="black", labelcolor="black")

    _fig.tight_layout()

    _pdf_name = f"tukey_hsd_ci_across_libraries_top{TOP_N}.pdf"
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


@app.cell
def _(Line2D, TOP_N, df, fig_to_pdf_bytes, mo, np, pl, plt):
    """Across-library Tukey HSD pairwise CI plot (docking libraries only).

    Docking libraries have 1 query each with 20 replicates.
    Each replicate is an independent algorithm run, so replicate is the unit of
    observation. Runs statsmodels pairwise_tukeyhsd for FWER-controlled CIs.
    """
    from statsmodels.stats.multicomp import pairwise_tukeyhsd as _pairwise_tukeyhsd

    _DOCK_LIBS = {"thrombin", "adenine", "amide", "quinazoline"}
    dock_df = df.filter(pl.col("library_id").is_in(_DOCK_LIBS))

    mo.stop(
        len(dock_df) == 0,
        mo.md("No docking library data at this top-N."),
    )

    dock_values = dock_df["recovered"].to_numpy()
    dock_groups = dock_df["method"].to_list()

    dock_tukey = _pairwise_tukeyhsd(dock_values, dock_groups, alpha=0.05)
    dock_summary = dock_tukey.summary()

    # Parse results
    dock_query_means = (
        dock_df.group_by(["library_id", "query_id", "method"])
        .agg(pl.col("recovered").mean().alias("mean_recovered"))
    )

    dock_left = []
    dock_right = []
    dock_means = []
    dock_ci_lo = []
    dock_ci_hi = []
    dock_pvals = []
    dock_reject = []

    for dock_row in dock_summary.data[1:]:
        dg1, dg2, d_meandiff, d_padj, d_lower, d_upper, d_reject = dock_row
        d_meandiff = float(d_meandiff)
        d_padj = float(d_padj)
        d_lower = float(d_lower)
        d_upper = float(d_upper)

        if d_meandiff >= 0:
            dock_left.append(dg1)
            dock_right.append(dg2)
            dock_means.append(d_meandiff)
            dock_ci_lo.append(d_lower)
            dock_ci_hi.append(d_upper)
        else:
            dock_left.append(dg2)
            dock_right.append(dg1)
            dock_means.append(-d_meandiff)
            dock_ci_lo.append(-d_upper)
            dock_ci_hi.append(-d_lower)

        dock_pvals.append(d_padj)
        dock_reject.append(bool(d_reject))

    dock_means = np.array(dock_means)
    dock_ci_lo = np.array(dock_ci_lo)
    dock_ci_hi = np.array(dock_ci_hi)
    dock_pvals = np.array(dock_pvals)

    dock_order = np.argsort(dock_means)[::-1]
    dock_left = [dock_left[i] for i in dock_order]
    dock_right = [dock_right[i] for i in dock_order]
    dock_reject = [dock_reject[i] for i in dock_order]
    dock_means = dock_means[dock_order]
    dock_ci_lo = dock_ci_lo[dock_order]
    dock_ci_hi = dock_ci_hi[dock_order]
    dock_pvals = dock_pvals[dock_order]

    dock_n_comp = len(dock_left)
    dock_y = np.arange(dock_n_comp)

    dock_fig, dock_ax = plt.subplots(
        figsize=(12, max(4, dock_n_comp * 0.55)),
        facecolor="white",
    )
    dock_ax.set_facecolor("white")

    _LEGACY = {"Legacy Enhanced-RWS", "Legacy Standard-Greedy"}
    _TACTICS = {"TACTICS Enhanced-RWS (GMIC)", "TACTICS Enhanced-TT-TS (GMIC)"}
    _CLR_SIG = "#4C78A8"
    _CLR_NS = "#999999"
    _CLR_KEY_TICK = "#4C78A8"

    def _is_tvl(left, right):
        return (
            (left in _LEGACY and right in _TACTICS)
            or (left in _TACTICS and right in _LEGACY)
        )

    for di in range(dock_n_comp):
        clr = _CLR_SIG if dock_reject[di] else _CLR_NS

        dock_ax.plot(
            [dock_ci_lo[di], dock_ci_hi[di]], [di, di],
            color=clr, linewidth=2.5, solid_capstyle="round", zorder=3,
        )
        dock_ax.plot(
            dock_means[di], di, "o",
            color=clr, markersize=8, zorder=4,
        )
        dp = dock_pvals[di]
        d_sig = (
            "****" if dp <= 0.0001
            else "***" if dp <= 0.001
            else "**" if dp <= 0.01
            else "*" if dp <= 0.05
            else "ns"
        )
        dock_ax.text(
            dock_ci_hi[di] + 0.5, di, d_sig,
            ha="left", va="center",
            fontsize=12, fontweight="bold", color="black",
        )

    dock_ax.axvline(0, color="black", linewidth=0.8, linestyle="-", zorder=2)

    dock_x_pad = max(dock_ci_hi.max(), abs(dock_ci_lo.min())) * 0.18
    dock_ax.set_xlim(
        min(dock_ci_lo.min(), 0) - dock_x_pad * 0.5,
        dock_ci_hi.max() + dock_x_pad,
    )

    dock_ax.set_yticks(dock_y)
    dock_ax.set_yticklabels(dock_right, fontsize=12, color="black")
    dock_ax.set_ylabel("Challenger", fontsize=14, color="black")

    dock_ax.set_ylim(-0.6, dock_n_comp - 0.4)
    dock_ax.invert_yaxis()

    dock_ax_right = dock_ax.twinx()
    dock_ax_right.set_ylim(dock_ax.get_ylim())
    dock_ax_right.set_yticks(dock_y)
    dock_ax_right.set_yticklabels(dock_left, fontsize=12, color="black")
    dock_ax_right.set_ylabel("Baseline", fontsize=14, color="black")
    dock_ax_right.tick_params(colors="black", labelsize=12)
    for sp in dock_ax_right.spines.values():
        sp.set_color("black")

    # Per-tick styling applied last so nothing resets them
    for di in range(dock_n_comp):
        if _is_tvl(dock_left[di], dock_right[di]):
            dock_ax.get_yticklabels()[di].set_color(_CLR_KEY_TICK)
            dock_ax.get_yticklabels()[di].set_fontweight("bold")
            dock_ax_right.get_yticklabels()[di].set_fontweight("bold")

    dock_ax.set_xlabel(
        f"\u0394 Recovery: Challenger \u2212 Baseline (top-{TOP_N})",
        fontsize=14, color="black",
    )
    dock_ax.set_title(
        f"Tukey HSD Pairwise CIs \u2014 4 Docking Datasets (top-{TOP_N}, FWER=0.05)",
        fontsize=15, fontweight="bold", color="black",
    )
    dock_ax.grid(axis="x", alpha=0.3, zorder=1, color="grey")
    dock_ax.tick_params(axis="x", colors="black", labelsize=12)
    dock_ax.tick_params(axis="y", color="black", labelsize=12)
    for sp in dock_ax.spines.values():
        sp.set_color("black")

    dock_legend = [
        Line2D([0], [0], color=_CLR_SIG, linewidth=2.5, marker="o",
               markersize=8, label="Significant (p < 0.05)"),
        Line2D([0], [0], color=_CLR_NS, linewidth=2.5, marker="o",
               markersize=8, label="Not significant"),
    ]
    dock_ax.legend(handles=dock_legend, loc="lower right", fontsize=12,
               facecolor="white", edgecolor="black", labelcolor="black")

    dock_fig.tight_layout()

    dock_pdf_name = f"tukey_hsd_ci_docking_top{TOP_N}.pdf"
    mo.vstack([
        plt.gcf(),
        mo.download(
            data=fig_to_pdf_bytes(dock_fig),
            filename=dock_pdf_name,
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## SAR Landscape Characterization — Reagent Posterior Separation

    Riley (2024) showed that for ROCS, reagents producing top compounds have high means
    and low variance (they "stand out from the crowd"), while for docking, top-hit reagents
    are statistically indistinguishable from mediocre ones. This plot reproduces that
    analysis for adenine and quinazoline docking libraries.
    """)
    return


@app.cell
def _(fig_to_pdf_bytes, mo, np, pl, plt, project_root):
    """Reagent posterior separation plot (Riley-style) for docking libraries.

    For each reagent, computes (mean score, std score) across all products containing
    that reagent. Colors reagents by whether they appear in top-100 compounds.
    Shows whether reagent-level statistics predict top-compound membership.
    """
    import re

    _SCORES_DIR = project_root / "data" / "scores"
    _TOP_N = 100

    # Library configs: (library_id, score_path, component_names, parse_function)
    def _parse_adenine(code):
        m = re.match(r"(amidine_\d+)_(isocyanide_db_\d+)_(aldehyde_\d+)", code)
        return m.groups() if m else (None, None, None)

    def _parse_quinazoline(code):
        parts = code.split("_")
        return tuple(parts) if len(parts) == 3 else (None, None, None)

    _libraries = [
        {
            "lib_id": "adenine",
            "components": ["Amidines", "Isocyanides", "Aldehydes"],
            "parse": _parse_adenine,
        },
        {
            "lib_id": "quinazoline",
            "components": ["Aminobenzoics", "Amines", "Acids"],
            "parse": _parse_quinazoline,
        },
    ]

    _fig, _axes = plt.subplots(
        2, 3, figsize=(18, 10), facecolor="white",
    )

    for _row_idx, _lib_cfg in enumerate(_libraries):
        _lib_id = _lib_cfg["lib_id"]
        _parse_fn = _lib_cfg["parse"]
        _comp_names = _lib_cfg["components"]

        _score_df = pl.read_parquet(_SCORES_DIR / f"{_lib_id}.parquet")

        # Filter NaN scores (docking failures)
        _score_df = _score_df.filter(pl.col("Scores").is_not_nan())

        # Identify top-N compounds (minimize for docking)
        _top_n_names = set(
            _score_df.sort("Scores", descending=False)
            .head(_TOP_N)["Product_Code"]
            .to_list()
        )

        # Parse product codes into component reagent IDs
        _codes = _score_df["Product_Code"].to_list()
        _scores = _score_df["Scores"].to_numpy()
        _is_top = np.array([c in _top_n_names for c in _codes])

        _parsed = [_parse_fn(c) for c in _codes]
        _n_comp = len(_comp_names)

        for _comp_idx in range(_n_comp):
            _ax = _axes[_row_idx, _comp_idx]
            _ax.set_facecolor("white")

            # Extract reagent IDs for this component
            _reagent_ids = [p[_comp_idx] for p in _parsed if p[_comp_idx] is not None]
            _reagent_scores = _scores[[i for i, p in enumerate(_parsed) if p[_comp_idx] is not None]]
            _reagent_top = _is_top[[i for i, p in enumerate(_parsed) if p[_comp_idx] is not None]]

            # Get unique reagents
            _unique_reagents = sorted(set(_reagent_ids))

            # Compute per-reagent mean and std
            _means = []
            _stds = []
            _in_top = []
            _reagent_labels = []

            # Build reagent → indices mapping
            _reagent_to_idx = {}
            for _i, _rid in enumerate(_reagent_ids):
                _reagent_to_idx.setdefault(_rid, []).append(_i)

            # Find which reagents appear in top-N products
            _top_reagent_ids = set()
            for _i, _rid in enumerate(_reagent_ids):
                if _reagent_top[_i]:
                    _top_reagent_ids.add(_rid)

            for _rid in _unique_reagents:
                _idxs = _reagent_to_idx[_rid]
                _s = _reagent_scores[_idxs]
                _means.append(float(np.mean(_s)))
                _stds.append(float(np.std(_s)))
                _in_top.append(_rid in _top_reagent_ids)
                _reagent_labels.append(_rid)

            _means = np.array(_means)
            _stds = np.array(_stds)
            _in_top = np.array(_in_top)

            # Plot non-top reagents first (grey), then top reagents (red)
            _ax.scatter(
                _means[~_in_top], _stds[~_in_top],
                c="#CCCCCC", s=40, alpha=0.6, edgecolors="none",
                label="Not in top-100", zorder=2,
            )
            _ax.scatter(
                _means[_in_top], _stds[_in_top],
                c="#E45756", s=60, alpha=0.9, edgecolors="black", linewidths=0.5,
                label="In top-100", zorder=3,
            )

            _ax.set_xlabel("Reagent Mean Score", fontsize=14, color="black")
            if _comp_idx == 0:
                _ax.set_ylabel("Reagent Std Score", fontsize=14, color="black")
            _ax.set_title(
                f"{_comp_names[_comp_idx]} ({sum(_in_top)}/{len(_in_top)} in top-100)",
                fontsize=14, fontweight="bold", color="black",
            )
            _ax.tick_params(colors="black", labelsize=12)
            _ax.grid(alpha=0.3, color="grey", zorder=1)
            for _spine in _ax.spines.values():
                _spine.set_color("black")

            if _comp_idx == 2:
                _ax.legend(fontsize=12, facecolor="white", edgecolor="black",
                           labelcolor="black", loc="upper right")

        # Row label
        _axes[_row_idx, 0].annotate(
            _lib_id.capitalize(),
            xy=(-0.35, 0.5), xycoords="axes fraction",
            fontsize=16, fontweight="bold", color="black",
            rotation=90, va="center", ha="center",
        )

    _fig.suptitle(
        f"SAR Landscape: Reagent Score Distributions vs Top-{_TOP_N} Membership (Docking)",
        fontsize=16, fontweight="bold", color="black", y=1.01,
    )
    _fig.tight_layout()

    mo.vstack([
        plt.gcf(),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename="reagent_posterior_separation_docking.pdf",
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


@app.cell
def _(fig_to_pdf_bytes, mo, np, pl, plt, project_root):
    """Reagent posterior separation plot for 2-component docking libraries (thrombin, amide)."""

    _SCORES_DIR = project_root / "data" / "scores"
    _TOP_N = 100

    def _parse_thrombin(code):
        parts = code.split("_")
        if len(parts) == 3:
            return (parts[0], parts[1] + "_" + parts[2])
        return (None, None)

    def _parse_amide(code):
        parts = code.split("_")
        return tuple(parts) if len(parts) == 2 else (None, None)

    _libraries = [
        {
            "lib_id": "thrombin",
            "components": ["Acids (130)", "Amines (3,844)"],
            "parse": _parse_thrombin,
        },
        {
            "lib_id": "amide",
            "components": ["Acids (10,000)", "Amines (1,000)"],
            "parse": _parse_amide,
        },
    ]

    _fig, _axes = plt.subplots(
        2, 2, figsize=(14, 10), facecolor="white",
    )

    for _row_idx, _lib_cfg in enumerate(_libraries):
        _lib_id = _lib_cfg["lib_id"]
        _parse_fn = _lib_cfg["parse"]
        _comp_names = _lib_cfg["components"]
        _n_comp = len(_comp_names)

        _score_df = pl.read_parquet(_SCORES_DIR / f"{_lib_id}.parquet")
        _score_df = _score_df.filter(pl.col("Scores").is_not_nan())

        # Top-N compounds (minimize for docking)
        _top_n_names = set(
            _score_df.sort("Scores", descending=False)
            .head(_TOP_N)["Product_Code"]
            .to_list()
        )

        _codes = _score_df["Product_Code"].to_list()
        _scores = _score_df["Scores"].to_numpy()
        _is_top = np.array([c in _top_n_names for c in _codes])

        _parsed = [_parse_fn(c) for c in _codes]

        for _comp_idx in range(_n_comp):
            _ax = _axes[_row_idx, _comp_idx]
            _ax.set_facecolor("white")

            _reagent_ids = [p[_comp_idx] for p in _parsed if p[_comp_idx] is not None]
            _reagent_scores = _scores[[i for i, p in enumerate(_parsed) if p[_comp_idx] is not None]]
            _reagent_top = _is_top[[i for i, p in enumerate(_parsed) if p[_comp_idx] is not None]]

            _unique_reagents = sorted(set(_reagent_ids))

            _means = []
            _stds = []
            _in_top = []

            _reagent_to_idx = {}
            for _i, _rid in enumerate(_reagent_ids):
                _reagent_to_idx.setdefault(_rid, []).append(_i)

            _top_reagent_ids = set()
            for _i, _rid in enumerate(_reagent_ids):
                if _reagent_top[_i]:
                    _top_reagent_ids.add(_rid)

            for _rid in _unique_reagents:
                _idxs = _reagent_to_idx[_rid]
                _s = _reagent_scores[_idxs]
                _means.append(float(np.mean(_s)))
                _stds.append(float(np.std(_s)))
                _in_top.append(_rid in _top_reagent_ids)

            _means = np.array(_means)
            _stds = np.array(_stds)
            _in_top = np.array(_in_top)

            _ax.scatter(
                _means[~_in_top], _stds[~_in_top],
                c="#CCCCCC", s=40, alpha=0.6, edgecolors="none",
                label="Not in top-100", zorder=2,
            )
            _ax.scatter(
                _means[_in_top], _stds[_in_top],
                c="#E45756", s=60, alpha=0.9, edgecolors="black", linewidths=0.5,
                label="In top-100", zorder=3,
            )

            _ax.set_xlabel("Reagent Mean Score", fontsize=14, color="black")
            if _comp_idx == 0:
                _ax.set_ylabel("Reagent Std Score", fontsize=14, color="black")
            _ax.set_title(
                f"{_comp_names[_comp_idx]} ({sum(_in_top)}/{len(_in_top)} in top-100)",
                fontsize=14, fontweight="bold", color="black",
            )
            _ax.tick_params(colors="black", labelsize=12)
            _ax.grid(alpha=0.3, color="grey", zorder=1)
            for _spine in _ax.spines.values():
                _spine.set_color("black")

            if _comp_idx == 1:
                _ax.legend(fontsize=12, facecolor="white", edgecolor="black",
                           labelcolor="black", loc="upper right")

        _axes[_row_idx, 0].annotate(
            _lib_cfg["lib_id"].capitalize(),
            xy=(-0.35, 0.5), xycoords="axes fraction",
            fontsize=16, fontweight="bold", color="black",
            rotation=90, va="center", ha="center",
        )

    _fig.suptitle(
        f"SAR Landscape: Reagent Score Distributions vs Top-{_TOP_N} Membership\n"
        f"2-Component Docking Libraries",
        fontsize=16, fontweight="bold", color="black", y=1.02,
    )
    _fig.tight_layout()

    mo.vstack([
        plt.gcf(),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename="reagent_posterior_separation_2comp_docking.pdf",
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


if __name__ == "__main__":
    app.run()
