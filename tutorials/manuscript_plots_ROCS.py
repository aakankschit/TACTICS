"""
TACTICS Manuscript Plots — ROCS Libraries

Loads the canonical 28-library recovery summary, filters to ROCS libraries
(20 similarity libraries × 109 queries × 20 replicates), and renders all
ROCS-related manuscript figures: aggregate Top-N recovery, Zhao Fig 5
extensions, Tukey HSD pairwise CIs, per-library breakdowns, and the ROCS
budget-sensitivity panel. The parallel docking figures live in
`tutorials/manuscript_plots.py`.

Run as app: marimo run tutorials/manuscript_plots_ROCS.py
Edit mode:  marimo edit tutorials/manuscript_plots_ROCS.py
"""

import marimo

__generated_with = "0.21.1"
app = marimo.App(
    width="full",
    app_title="TACTICS Manuscript Plots — ROCS Libraries",
)


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
    return Line2D, fig_to_pdf_bytes, mo, np, pl, plt, project_root, stats


@app.cell
def _(mo):
    """Top-N slider and method selector for filtering pre-computed recovery data."""
    topn_slider = mo.ui.slider(
        start=25, stop=250, step=25, value=100,
        label="Top-N Recovery Threshold",
    )
    _all_methods = [
        "TACTICS Enhanced-TT-TS (GMIC)",
        "TACTICS Enhanced-RWS (GMIC)",
        "Legacy Enhanced-RWS",
        "Legacy Standard-Greedy",
    ]
    method_selector = mo.ui.multiselect(
        options=_all_methods,
        value=_all_methods,
        label="Methods to include",
    )
    mo.vstack([topn_slider, method_selector])
    return method_selector, topn_slider


@app.cell
def _(method_selector, mo, pl, project_root, topn_slider):
    """Load pre-computed recovery and statistical analysis results."""
    analysis_dir = project_root / "outputs" / "full_benchmark_ema_cats" / "analysis"

    TOP_N = topn_slider.value
    selected_methods = method_selector.value

    # --- Recovery summary (all libraries × top-N thresholds ≤ 250) ---
    # Top-N range capped at 250 manuscript-wide (per 2026-05-27 scope decision);
    # wider top-N values exist in the parquet but are not used in any manuscript figure.
    recovery_all_raw = (
        pl.read_parquet(analysis_dir / "recovery_summary.parquet")
        .filter(pl.col("top_n") <= 250)
    )
    recovery_all = recovery_all_raw.filter(pl.col("method").is_in(selected_methods))

    # Filter to selected top-N threshold
    df = recovery_all.filter(pl.col("top_n") == TOP_N)

    # --- Statistical test results (filter to selected method pairs) ---
    tukey_df = pl.read_parquet(analysis_dir / "tukey_hsd_results.parquet")
    tukey_at_n = tukey_df.filter(
        (pl.col("top_n") == TOP_N)
        & pl.col("group1").is_in(selected_methods)
        & pl.col("group2").is_in(selected_methods)
    )

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
        - **Available thresholds**: {sorted(recovery_all_raw['top_n'].unique().to_list())}
        - **Selected methods** ({len(selected_methods)}/4): {', '.join(methods)}
        - **Total trials at top-{TOP_N}**: {n_trials:,}
        - **Libraries**: {n_libraries} ({', '.join(libraries)})
        - **Replicates**: {df['replicate'].n_unique()}
        - **Tukey HSD comparisons**: {len(tukey_at_n):,} ({tukey_at_n.filter(pl.col('reject')).height} significant)
        - **ANOVA tests**: {len(anova_at_n)} ({anova_at_n.filter(pl.col('significance') != 'ns').height} significant)
        """
    )
    return TOP_N, df, recovery_all, tukey_at_n


@app.cell
def _():
    """Statistical-test selection rule and star color convention.

    Rule (by scoring type):
      - Docking libraries → Welch's t-test + Bonferroni-4 on the planned
        TACTICS-vs-Legacy pairs. Docking benchmarks here use one query
        per library × 20 replicates, so within-method variance (including
        bimodal-failure modes) shows up directly. Tukey HSD's pooled-MSE
        homoscedasticity assumption is violated whenever any single method
        is bimodal, so Welch is the appropriate confirmatory test.
      - ROCS libraries → Tukey HSD per Ash et al. ROCS uses 109 queries
        × 20 replicates, and cross-query variance dominates the per-method
        variance, so the variances are approximately homogeneous across
        methods and Tukey's pooled-MSE assumption holds.

    Docking libraries: thrombin, adenine, adenine_hybrid, amide, quinazoline.
    """
    DOCKING_LIBS = {
        "thrombin", "adenine", "adenine_hybrid", "amide", "quinazoline",
        "thrombin_hybrid", "amide_hybrid", "quinazoline_hybrid",
    }
    TUKEY_STAR_COLOR = "black"
    WELCH_STAR_COLOR = "#4477AA"
    return (TUKEY_STAR_COLOR,)


@app.cell
def _(mo):
    mo.md(r"""
    # Section 1 — Overall ROCS Recovery

    Aggregate recovery across the 20 ROCS similarity libraries (109 queries × 20
    replicates per library). Methods compared: TACTICS TT-TS, TACTICS RWS,
    Legacy RWS, Legacy Standard-Greedy.
    """)
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

    _DOCK_LIBS = {
        "thrombin", "adenine", "adenine_hybrid", "amide", "quinazoline",
        "thrombin_hybrid", "amide_hybrid", "quinazoline_hybrid",
    }
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
        .filter(pl.col("top_n") <= 250)
        .with_columns(
            pl.col("library_id").replace(_ROCS_NCOMP).cast(pl.Int64).alias("n_components")
        )
    )

    mo.stop(len(_rocs_all) == 0, mo.md("No ROCS library data available."))

    _method_colors = {
        "TACTICS Enhanced-TT-TS (GMIC)": "#E45756",
        "TACTICS Enhanced-RWS (GMIC)": "#4C78A8",
        "Legacy Standard-Greedy": "#76B7B2",
        "Legacy Enhanced-RWS": "#B07AA1",
    }
    _method_markers = {
        "TACTICS Enhanced-TT-TS (GMIC)": "D",
        "TACTICS Enhanced-RWS (GMIC)": "o",
        "Legacy Standard-Greedy": "s",
        "Legacy Enhanced-RWS": "^",
    }
    _method_order = [
        "TACTICS Enhanced-TT-TS (GMIC)",
        "TACTICS Enhanced-RWS (GMIC)",
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
def _(fig_to_pdf_bytes, mo, np, pl, plt, recovery_all, stats):
    """ROCS aggregate Δ recovery bars — TACTICS TT-TS vs each Legacy method, split by component count.

    For each (library × top-N), compute mean recovery per method (across query × replicate),
    then Δ_lib = mean(TT-TS_lib) − mean(Legacy_lib). Aggregate across the 10 libraries in
    each component split: bar height = mean Δ, error bar = SE across libraries
    (std / sqrt(10)). Significance: scipy.stats.ttest_1samp on the per-library Δ vector
    (H0: Δ=0); *** p<0.001, ** p<0.01, * p<0.05, else ns. A paired across-library test
    (one-sample t on Δ) is the natural aggregate analogue of per-library Welch tests and
    avoids the pooled-variance contamination issue.

    On ROCS, TACTICS TT-TS and TACTICS RWS are practically equivalent at the per-library
    level: across 200 (library × top-N) cells the larger TACTICS method beats the other
    by at most 0.79 pp (≤ 0.2 compounds, on `poparov` at top-25). TT-TS is reported here
    as the headline TACTICS method because it's the per-library *statistical* winner on
    179/200 cells; the choice has no practical consequence on ROCS, so this figure can
    be read straightforwardly as "TACTICS vs Legacy" without method-selection caveats.
    """

    _ROCS_NCOMP_BAR = {
        "rxn101": 2, "rxn102a": 2, "rxn108b": 2, "rxn111b": 2, "rxn114b": 2,
        "rxn203": 2, "rxn205": 2, "rxn206": 2, "rxn207": 2, "rxn208": 2,
        "amide-suzuki": 3, "betti": 3, "dobener": 3,
        "groebke-blackburn-bienayme": 3, "mannich": 3, "niementowski": 3,
        "orru": 3, "passerini": 3, "petasis": 3, "poparov": 3,
    }
    _HEADLINE_TACTICS = "TACTICS Enhanced-TT-TS (GMIC)"
    _LEGACIES_BAR = [
        ("Legacy Standard-Greedy", "Legacy Greedy", "#F28E2B"),
        ("Legacy Enhanced-RWS", "Legacy RWS", "#4C78A8"),
    ]
    _ROCS_LIBS_BAR = set(_ROCS_NCOMP_BAR.keys())

    _rocs_bar_all = (
        recovery_all
        .filter(pl.col("library_id").is_in(_ROCS_LIBS_BAR))
        .filter(pl.col("top_n") <= 250)
        .filter(pl.col("method").is_in({_HEADLINE_TACTICS} | {m for m, _, _ in _LEGACIES_BAR}))
    )

    mo.stop(len(_rocs_bar_all) == 0, mo.md("No ROCS library data available."))

    # Per (library × method × top_n) mean recovery_frac (collapse query × replicate)
    _per_lib = (
        _rocs_bar_all.group_by(["library_id", "method", "top_n"])
        .agg(pl.col("recovery_frac").mean().alias("mean_frac"))
    )

    _topn_vals_bar = sorted(_rocs_bar_all["top_n"].unique().to_list())

    def _sig_stars(_p):
        if _p is None or np.isnan(_p):
            return ""
        if _p < 0.001:
            return "***"
        if _p < 0.01:
            return "**"
        if _p < 0.05:
            return "*"
        return ""

    _fig, _axes = plt.subplots(1, 2, figsize=(16, 6), facecolor="white", sharey=True)
    _bar_width = 0.35

    _panel_diffs = {2: [], 3: []}

    for _panel_idx, _nc in enumerate([2, 3]):
        _ax = _axes[_panel_idx]
        _ax.set_facecolor("white")

        _libs_nc = sorted([_lib for _lib, _v in _ROCS_NCOMP_BAR.items() if _v == _nc])
        _x_pos = np.arange(len(_topn_vals_bar))

        for _gi, (_leg_meth, _leg_lbl, _bar_color) in enumerate(_LEGACIES_BAR):
            _mean_diffs = []
            _se_diffs = []
            _sigs = []

            for _tn in _topn_vals_bar:
                _per_lib_deltas = []
                for _lib in _libs_nc:
                    _tt = _per_lib.filter(
                        (pl.col("library_id") == _lib)
                        & (pl.col("method") == _HEADLINE_TACTICS)
                        & (pl.col("top_n") == _tn)
                    )["mean_frac"]
                    _leg = _per_lib.filter(
                        (pl.col("library_id") == _lib)
                        & (pl.col("method") == _leg_meth)
                        & (pl.col("top_n") == _tn)
                    )["mean_frac"]
                    if len(_tt) == 0 or len(_leg) == 0:
                        continue
                    _delta = (float(_tt[0]) - float(_leg[0])) * 100.0
                    _per_lib_deltas.append(_delta)

                _per_lib_deltas = np.array(_per_lib_deltas, dtype=float)
                if len(_per_lib_deltas) == 0:
                    _mean_diffs.append(0.0)
                    _se_diffs.append(0.0)
                    _sigs.append("")
                    continue
                _n_libs = len(_per_lib_deltas)
                _mean_diffs.append(float(_per_lib_deltas.mean()))
                _se_diffs.append(
                    float(_per_lib_deltas.std(ddof=1) / np.sqrt(_n_libs))
                    if _n_libs > 1 else 0.0
                )
                if _n_libs >= 2 and _per_lib_deltas.std(ddof=1) > 0:
                    _tres = stats.ttest_1samp(_per_lib_deltas, popmean=0.0)
                    _sigs.append(_sig_stars(_tres.pvalue))
                else:
                    _sigs.append("")

            _mean_diffs = np.array(_mean_diffs)
            _se_diffs = np.array(_se_diffs)
            _offset = (_gi - 0.5) * _bar_width
            _ax.bar(
                _x_pos + _offset, _mean_diffs, width=_bar_width,
                color=_bar_color, edgecolor="black", linewidth=0.5,
                yerr=_se_diffs, error_kw={"ecolor": "black", "capsize": 3, "elinewidth": 1},
                zorder=3, label=f"vs {_leg_lbl}",
            )

            _dmax = max(abs(_mean_diffs).max() if len(_mean_diffs) else 1, 1e-6)
            _star_pad = _dmax * 0.05
            for _i, (_d, _se, _s) in enumerate(zip(_mean_diffs, _se_diffs, _sigs)):
                if _s:
                    _yp = _d + _se + _star_pad if _d >= 0 else _d - _se - _star_pad
                    _ax.text(
                        _x_pos[_i] + _offset, _yp, _s, ha="center",
                        va="bottom" if _d >= 0 else "top",
                        fontsize=10, fontweight="bold", color="black",
                    )

            _panel_diffs[_nc].extend(list(_mean_diffs + _se_diffs) + list(_mean_diffs - _se_diffs))

        _ax.axhline(0, color="black", linewidth=0.8, zorder=2)
        _ax.set_xticks(_x_pos)
        _ax.set_xticklabels([str(t) for t in _topn_vals_bar], fontsize=12, color="black")
        _ax.set_xlabel("Top-N Threshold", fontsize=14, color="black")
        _ax.set_title(f"{_nc}-Component Libraries", fontsize=15,
                      fontweight="bold", color="black")
        _ax.tick_params(colors="black", labelsize=12)
        _ax.grid(axis="y", alpha=0.3, color="grey", zorder=1)
        for _spine in _ax.spines.values():
            _spine.set_color("black")

    _all_extrema = np.array(_panel_diffs[2] + _panel_diffs[3], dtype=float)
    if len(_all_extrema) > 0:
        _ymin = float(_all_extrema.min())
        _ymax = float(_all_extrema.max())
        _pad = max(abs(_ymax - _ymin) * 0.20, 1.0)
        for _ax in _axes:
            _ax.set_ylim(_ymin - _pad, _ymax + _pad)

    _axes[0].set_ylabel(
        "Δ % Recovery (TACTICS TT-TS − Legacy)", fontsize=14, color="black"
    )

    _handles, _labels = _axes[0].get_legend_handles_labels()
    _fig.legend(
        _handles, _labels,
        loc="lower center", ncol=2, fontsize=12,
        facecolor="white", edgecolor="black", labelcolor="black",
        bbox_to_anchor=(0.5, -0.08),
    )

    _fig.suptitle(
        "Δ Recovery — TACTICS TT-TS vs Legacy Methods (ROCS Libraries)",
        fontsize=15, fontweight="bold", color="black", y=1.01,
    )
    _fig.tight_layout()

    mo.vstack([
        plt.gcf(),
        mo.md(
            "**Headline method.** TACTICS TT-TS is reported here as the headline "
            "TACTICS method. The two TACTICS strategies (TT-TS, RWS) are "
            "practically equivalent on ROCS at the per-library level — the "
            "largest gap on any (library × top-N) cell is 0.79 pp ≈ 0.2 "
            "compounds — so this figure is essentially \"TACTICS vs Legacy\" "
            "without a method-selection caveat."
        ),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename="rocs_tactics_vs_legacy_bars.pdf",
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


@app.cell
def _(fig_to_pdf_bytes, mo, np, pl, plt, recovery_all):
    """Chained-decomposition stacked bars — split each TACTICS method's gain over
    Legacy Greedy into its two additive components:

        Δ(TACTICS, L-Grdy) = Δ(L-RWS, L-Grdy)  +  Δ(TACTICS, L-RWS)
                              └─ "Zhao 2025 lift" ┘  └ "TACTICS algorithmic lift" ┘

    The two segments stack at the zero line. When Zhao lift is negative
    (Legacy RWS underperforms Legacy Greedy — documented on 2-comp ROCS at
    wider top-N), the Zhao segment extends below 0 and the TACTICS segment
    starts from that negative point and stacks upward; the bar's top edge
    is always Δ(TACTICS, L-Grdy). Two grouped bars per top-N (TT-TS, T-RWS);
    panels split by 2-comp vs 3-comp ROCS.
    """
    _ROCS_NCOMP = {
        "rxn101": 2, "rxn102a": 2, "rxn108b": 2, "rxn111b": 2, "rxn114b": 2,
        "rxn203": 2, "rxn205": 2, "rxn206": 2, "rxn207": 2, "rxn208": 2,
        "amide-suzuki": 3, "betti": 3, "dobener": 3,
        "groebke-blackburn-bienayme": 3, "mannich": 3, "niementowski": 3,
        "orru": 3, "passerini": 3, "petasis": 3, "poparov": 3,
    }
    _METHODS = {
        "TT-TS":   "TACTICS Enhanced-TT-TS (GMIC)",
        "T-RWS":   "TACTICS Enhanced-RWS (GMIC)",
        "L-RWS":   "Legacy Enhanced-RWS",
        "L-Grdy":  "Legacy Standard-Greedy",
    }
    _ZHAO_COLOR  = "#B07AA1"   # Zhao 2025 lift — purple, matches Legacy RWS palette
    _TACTICS_COLORS = {"TT-TS": "#E45756", "T-RWS": "#4C78A8"}

    _rocs = (
        recovery_all
        .filter(pl.col("library_id").is_in(set(_ROCS_NCOMP.keys())))
        .filter(pl.col("top_n") <= 250)
        .filter(pl.col("method").is_in(set(_METHODS.values())))
    )
    mo.stop(len(_rocs) == 0, mo.md("No ROCS data available."))

    # Per-library mean recovery (collapse query × replicate); then per-method dict per lib.
    _per_lib = (
        _rocs.group_by(["library_id", "method", "top_n"])
        .agg(pl.col("recovery_frac").mean().alias("mean_frac"))
    )

    _topn_vals = sorted(_rocs["top_n"].unique().to_list())

    _fig, _axes = plt.subplots(1, 2, figsize=(16, 6), facecolor="white", sharey=True)
    _bar_w = 0.38

    for _panel_idx, _nc in enumerate([2, 3]):
        _ax = _axes[_panel_idx]
        _ax.set_facecolor("white")
        _libs_nc = sorted([_l for _l, _v in _ROCS_NCOMP.items() if _v == _nc])
        _x_pos = np.arange(len(_topn_vals))

        for _gi, _tac_short in enumerate(["TT-TS", "T-RWS"]):
            _tac_method = _METHODS[_tac_short]
            _zhao_means = []
            _tac_means = []
            for _tn in _topn_vals:
                _zhao_per_lib = []
                _tac_per_lib = []
                for _lib in _libs_nc:
                    _tn_data = _per_lib.filter(
                        (pl.col("library_id") == _lib) & (pl.col("top_n") == _tn)
                    )
                    _by_method = {
                        _r["method"]: float(_r["mean_frac"])
                        for _r in _tn_data.iter_rows(named=True)
                    }
                    if not all(_m in _by_method for _m in _METHODS.values()):
                        continue
                    _zhao_per_lib.append(
                        (_by_method[_METHODS["L-RWS"]] - _by_method[_METHODS["L-Grdy"]]) * 100
                    )
                    _tac_per_lib.append(
                        (_by_method[_tac_method] - _by_method[_METHODS["L-RWS"]]) * 100
                    )
                _zhao_means.append(np.mean(_zhao_per_lib) if _zhao_per_lib else 0.0)
                _tac_means.append(np.mean(_tac_per_lib) if _tac_per_lib else 0.0)

            _zhao_means = np.array(_zhao_means)
            _tac_means = np.array(_tac_means)
            _offset = (_gi - 0.5) * _bar_w

            # Zhao segment anchored at 0. Hatched when negative to flag the regime
            # where Legacy RWS underperforms Legacy Greedy.
            _ax.bar(
                _x_pos + _offset, _zhao_means, width=_bar_w,
                color=_ZHAO_COLOR, edgecolor="black", linewidth=0.5,
                hatch=["////" if _z < 0 else "" for _z in _zhao_means],
                label="Zhao 2025 lift (L-RWS − L-Grdy)" if _gi == 0 else None,
                zorder=3,
            )
            # TACTICS segment stacked on top of Zhao segment.
            _ax.bar(
                _x_pos + _offset, _tac_means, width=_bar_w,
                bottom=_zhao_means,
                color=_TACTICS_COLORS[_tac_short], edgecolor="black", linewidth=0.5,
                label=f"TACTICS lift ({_tac_short} − L-RWS)",
                zorder=3,
            )

            # Net total marker (top edge of full stack = Δ(TACTICS, L-Grdy))
            _net = _zhao_means + _tac_means
            for _xi, _net_yi in zip(_x_pos + _offset, _net):
                _ax.plot([_xi - _bar_w/2, _xi + _bar_w/2], [_net_yi, _net_yi],
                         color="black", linewidth=1.2, zorder=5)

        _ax.axhline(0, color="black", linewidth=1.0, zorder=2)
        _ax.set_xticks(_x_pos)
        _ax.set_xticklabels([str(t) for t in _topn_vals], fontsize=12, color="black")
        _ax.set_xlabel("Top-N Threshold", fontsize=14, color="black")
        _ax.set_title(f"{_nc}-Component ROCS Libraries", fontsize=15,
                      fontweight="bold", color="black")
        _ax.tick_params(colors="black", labelsize=12)
        _ax.grid(axis="y", alpha=0.3, color="grey", zorder=1)
        for _sp in _ax.spines.values():
            _sp.set_color("black")

    _axes[0].set_ylabel(
        "Δ % Recovery vs Legacy Greedy (chained decomposition)",
        fontsize=14, color="black",
    )

    _handles, _labels = _axes[0].get_legend_handles_labels()
    _fig.legend(
        _handles, _labels, loc="lower center", ncol=3, fontsize=11,
        facecolor="white", edgecolor="black", labelcolor="black",
        bbox_to_anchor=(0.5, -0.10),
    )
    _fig.suptitle(
        "Chained-Decomposition Stack — Zhao 2025 Lift + TACTICS Algorithmic Lift (ROCS)",
        fontsize=14, fontweight="bold", color="black", y=1.01,
    )
    _fig.tight_layout()

    mo.vstack([
        plt.gcf(),
        mo.md(
            "**Reading guide.** Black tick at the top of each bar = net "
            "Δ(TACTICS, Legacy Greedy). Lower (purple) segment = Δ(Legacy RWS, "
            "Legacy Greedy), i.e. the lift Zhao 2025 already provided — hatched "
            "when it goes negative (Legacy RWS hurts vs Greedy, expected at "
            "wider top-N on 2-comp ROCS per STATUS line 24). Upper (colored) "
            "segment = Δ(TACTICS, Legacy RWS), the algorithmic lift this work "
            "adds on top. The total visual height of the colored region equals "
            "Δ(TACTICS, Legacy RWS), independent of how Zhao 2025 compares to "
            "Greedy — so the TACTICS contribution stays legible even where the "
            "Zhao baseline is negative."
        ),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename="rocs_chained_decomposition_stack.pdf",
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


@app.cell
def _(mo):
    mo.md(r"""
    # Section 2 — Diversity vs Recovery (Zhao Fig 5 extended)

    Per-library recovery rate plotted against top-100 chemical diversity (BMS metric,
    following Zhao et al. 2025 Fig 5). Both Legacy and TACTICS methods decline with
    diversity; TACTICS' data-driven reallocation produces gentler decay.
    """)
    return


@app.cell
def _(mo, np, pl, plt, project_root, recovery_all, stats):
    """Zhao 2025 Fig 5 — data prep + helper functions.

    Computes Zhao's diversity metric (total unique reagents in top-100, averaged
    across queries per library) for the 20 ROCS libraries (10 two-component, 10
    three-component) and joins it with our per-library top-100 recovery. The
    result is one (mean_diversity, mean_recovery) point per (library, method)
    pair — exactly Zhao's per-library aggregation.

    Defines two plotting helpers that share `zhao_data`:
    - `zhao_pairwise_plot`: 2-panel scatter (2-comp / 3-comp) for any
      (challenger, baseline) pair, with per-library vertical arrows showing the
      within-library Δ recovery (baseline → challenger).
    - `zhao_delta_plot`: 2-panel companion that plots Δ = challenger − baseline
      directly against diversity, with linear fit + slope + R².
    """
    _ROCS_NCOMP = {
        "rxn101": 2, "rxn102a": 2, "rxn108b": 2, "rxn111b": 2, "rxn114b": 2,
        "rxn203": 2, "rxn205": 2, "rxn206": 2, "rxn207": 2, "rxn208": 2,
        "amide-suzuki": 3, "betti": 3, "dobener": 3,
        "groebke-blackburn-bienayme": 3, "mannich": 3, "niementowski": 3,
        "orru": 3, "passerini": 3, "petasis": 3, "poparov": 3,
    }
    _SCORES_DIR = project_root / "data" / "scores"

    _diversity_records = []
    for _lib_id, _n_comp in _ROCS_NCOMP.items():
        _path = _SCORES_DIR / f"{_lib_id}.parquet"
        if not _path.exists():
            continue
        _df = pl.read_parquet(_path)
        _query_cols = sorted(c for c in _df.columns if c.startswith("query_"))
        for _qcol in _query_cols:
            _sub = _df.select(["Name", _qcol]).filter(pl.col(_qcol).is_not_nan())
            _top100 = _sub.sort(_qcol, descending=True).head(100)
            _names = _top100["Name"].to_list()
            _comp_sets = [set() for _ in range(_n_comp)]
            for _nm in _names:
                _parts = _nm.split("_")
                if len(_parts) != _n_comp:
                    continue
                for _i, _r in enumerate(_parts):
                    _comp_sets[_i].add(_r)
            _diversity_records.append({
                "library_id": _lib_id, "query_id": _qcol,
                "n_components": _n_comp, "diversity": sum(len(_s) for _s in _comp_sets),
            })
    _diversity_df = pl.DataFrame(_diversity_records)
    mo.stop(len(_diversity_df) == 0, mo.md("No ROCS diversity data could be computed."))

    _rec_per_query = (
        recovery_all
        .filter(pl.col("library_id").is_in(set(_ROCS_NCOMP.keys())))
        .filter(pl.col("top_n") == 100)
        .group_by(["library_id", "query_id", "method"])
        .agg(pl.col("recovery_frac").mean().alias("rec_frac"))
    )
    _joined = _rec_per_query.join(_diversity_df, on=["library_id", "query_id"], how="inner")
    zhao_data = (
        _joined.group_by(["library_id", "method", "n_components"])
        .agg([
            pl.col("rec_frac").mean().alias("mean_rec"),
            pl.col("diversity").mean().alias("mean_diversity"),
        ])
    )

    mo.md(
        f"**Zhao Fig 5 data prepared.** "
        f"{zhao_data.filter(pl.col('n_components')==2).group_by('method').len()['len'].max()} "
        f"libraries per method in 2-comp panel, "
        f"{zhao_data.filter(pl.col('n_components')==3).group_by('method').len()['len'].max()} "
        f"in 3-comp panel."
    )

    def zhao_pairwise_plot(
        challenger, baseline,
        challenger_color, baseline_color,
        challenger_label, baseline_label,
        title,
    ):
        """Pairwise scatter + linear-fit, two panels (2-comp left, 3-comp right).

        Each panel: 10 ROCS libraries as data points per method. X-axis = mean
        Zhao diversity (unique reagents in top-100, averaged across queries),
        Y-axis = mean top-100 recovery. Per-library vertical arrows connect
        baseline → challenger to make the within-library Δ visible (diversity
        is a library-level property computed from the ground-truth top-100, so
        the same x for both methods). Slopes and R² shown in legend.
        """
        _fig, _axes = plt.subplots(1, 2, figsize=(18, 8), facecolor="white", sharey=False)
        _series = [
            (challenger, challenger_color, challenger_label, "o"),
            (baseline, baseline_color, baseline_label, "o"),
        ]
        for _panel_idx, _nc in enumerate([2, 3]):
            _ax = _axes[_panel_idx]
            _ax.set_facecolor("white")
            _panel = zhao_data.filter(pl.col("n_components") == _nc)
            _n_libs = _panel.filter(pl.col("method") == challenger).height

            _ch_libs = _panel.filter(pl.col("method") == challenger).sort("library_id")
            _bs_libs = _panel.filter(pl.col("method") == baseline).sort("library_id")
            _paired = _ch_libs.join(_bs_libs, on=["library_id"], suffix="_bs")
            for _row in _paired.iter_rows(named=True):
                _x_lib = _row["mean_diversity"]
                _y_bs = _row["mean_rec_bs"] * 100
                _y_ch = _row["mean_rec"] * 100
                _ax.annotate(
                    "", xy=(_x_lib, _y_ch), xytext=(_x_lib, _y_bs),
                    arrowprops=dict(arrowstyle="-|>", color="#333333",
                                    alpha=0.65, lw=1.4,
                                    shrinkA=5, shrinkB=5,
                                    mutation_scale=20),
                    zorder=2.5,
                )

            for _method, _color, _short, _marker in _series:
                _md = _panel.filter(pl.col("method") == _method).sort("library_id")
                if len(_md) < 2:
                    continue
                _x = _md["mean_diversity"].to_numpy()
                _y = _md["mean_rec"].to_numpy() * 100
                _slope, _intercept, _r, _p, _se = stats.linregress(_x, _y)
                _x_fit = np.linspace(_x.min(), _x.max(), 50)
                _y_fit = _slope * _x_fit + _intercept
                _ax.scatter(
                    _x, _y, color=_color, marker=_marker, s=110,
                    alpha=0.85, edgecolor="black", linewidth=0.7, zorder=3,
                    label=f"{_short} (n={len(_x)}, slope={_slope:+.3f}, R²={_r**2:.2f})",
                )
                _ax.plot(_x_fit, _y_fit, color=_color, linewidth=2.5, alpha=0.7, zorder=2)
            _ax.set_xlabel(
                "Total unique reagents in top-100",
                fontsize=13, color="black",
            )
            _ax.set_ylabel("Mean Top-100 Recovery (%)", fontsize=13, color="black")
            _ax.set_title(
                f"{_nc}-Component ROCS Libraries (n={_n_libs})",
                fontsize=14, fontweight="bold", color="black",
            )
            _ax.tick_params(colors="black", labelsize=12)
            _ax.grid(alpha=0.3, color="grey", zorder=1)
            for _spine in _ax.spines.values():
                _spine.set_color("black")
            _ax.legend(
                loc="upper center", bbox_to_anchor=(0.5, -0.14),
                ncol=1, fontsize=11, facecolor="white",
                edgecolor="black", labelcolor="black", framealpha=0.95,
            )
        _fig.suptitle(title, fontsize=14, fontweight="bold", color="black", y=1.02)
        _fig.tight_layout(rect=[0, 0.08, 1, 0.96])
        return _fig

    def zhao_delta_plot(
        challenger, baseline,
        challenger_color,
        challenger_label, baseline_label,
        title,
    ):
        """Δ recovery vs diversity, two panels (2-comp left, 3-comp right).

        For each library, plots Δ = challenger recovery − baseline recovery
        (percentage points) against ground-truth top-100 diversity. Linear fit,
        slope, R², and zero-Δ reference line shown. Positive Δ = challenger
        wins on that library. Companion to `zhao_pairwise_plot`: makes the
        within-library Δ the primary visual quantity rather than a paired
        scatter.
        """
        _fig, _axes = plt.subplots(1, 2, figsize=(18, 8), facecolor="white", sharey=False)
        for _panel_idx, _nc in enumerate([2, 3]):
            _ax = _axes[_panel_idx]
            _ax.set_facecolor("white")
            _panel = zhao_data.filter(pl.col("n_components") == _nc)
            _ch = _panel.filter(pl.col("method") == challenger).sort("library_id")
            _bs = _panel.filter(pl.col("method") == baseline).sort("library_id")
            _paired = _ch.join(_bs, on=["library_id"], suffix="_bs")
            if _paired.height < 2:
                continue
            _x = _paired["mean_diversity"].to_numpy()
            _y = (_paired["mean_rec"].to_numpy() - _paired["mean_rec_bs"].to_numpy()) * 100
            _slope, _intercept, _r, _p, _se = stats.linregress(_x, _y)
            _x_fit = np.linspace(_x.min(), _x.max(), 50)
            _y_fit = _slope * _x_fit + _intercept
            _ax.axhline(0, color="black", linewidth=1, linestyle="--", alpha=0.5, zorder=1)
            _ax.scatter(
                _x, _y, color=challenger_color, marker="D", s=110,
                alpha=0.85, edgecolor="black", linewidth=0.7, zorder=3,
                label=f"Δ (n={len(_x)}, slope={_slope:+.3f}, R²={_r**2:.2f})",
            )
            _ax.plot(_x_fit, _y_fit, color=challenger_color, linewidth=2.5, alpha=0.7, zorder=2)
            _ax.set_xlabel(
                "Total unique reagents in top-100",
                fontsize=13, color="black",
            )
            _ax.set_ylabel(
                f"Δ Recovery (pp): {challenger_label} − {baseline_label}",
                fontsize=13, color="black",
            )
            _ax.set_title(
                f"{_nc}-Component ROCS Libraries (n={len(_x)})",
                fontsize=14, fontweight="bold", color="black",
            )
            _ax.tick_params(colors="black", labelsize=12)
            _ax.grid(alpha=0.3, color="grey", zorder=1)
            for _spine in _ax.spines.values():
                _spine.set_color("black")
            _ax.legend(
                loc="best", fontsize=11, facecolor="white",
                edgecolor="black", labelcolor="black", framealpha=0.95,
            )
        _fig.suptitle(title, fontsize=14, fontweight="bold", color="black", y=1.02)
        _fig.tight_layout()
        return _fig

    return zhao_delta_plot, zhao_pairwise_plot


@app.cell
def _(fig_to_pdf_bytes, mo, plt, zhao_pairwise_plot):
    """Zhao Fig 5 — TACTICS RWS vs Legacy RWS (Zhao 2025)."""
    _fig = zhao_pairwise_plot(
        challenger="TACTICS Enhanced-RWS (GMIC)",
        baseline="Legacy Enhanced-RWS",
        challenger_color="#4C78A8",
        baseline_color="#B07AA1",
        challenger_label="TACTICS RWS (GMIC)",
        baseline_label="Legacy RWS (Zhao 2025)",
        title="Top-100 Recovery vs Library Diversity — TACTICS RWS vs Legacy RWS (Zhao 2025)",
    )
    mo.vstack([
        plt.gcf(),
        mo.download(data=fig_to_pdf_bytes(_fig),
                    filename="zhao_fig5_tactics_rws_vs_legacy_rws.pdf",
                    mimetype="application/pdf", label="Download as PDF"),
    ])
    return


@app.cell
def _(fig_to_pdf_bytes, mo, plt, zhao_pairwise_plot):
    """Zhao Fig 5 — TACTICS RWS vs Legacy Greedy (Klarich 2024 TS)."""
    _fig = zhao_pairwise_plot(
        challenger="TACTICS Enhanced-RWS (GMIC)",
        baseline="Legacy Standard-Greedy",
        challenger_color="#4C78A8",
        baseline_color="#76B7B2",
        challenger_label="TACTICS RWS (GMIC)",
        baseline_label="Legacy Greedy (Klarich 2024 TS)",
        title="Top-100 Recovery vs Library Diversity — TACTICS RWS vs Legacy Greedy (Klarich TS)",
    )
    mo.vstack([
        plt.gcf(),
        mo.download(data=fig_to_pdf_bytes(_fig),
                    filename="zhao_fig5_tactics_rws_vs_legacy_greedy.pdf",
                    mimetype="application/pdf", label="Download as PDF"),
    ])
    return


@app.cell
def _(fig_to_pdf_bytes, mo, plt, zhao_pairwise_plot):
    """Zhao Fig 5 — TT-TS vs Legacy RWS (Zhao 2025)."""
    _fig = zhao_pairwise_plot(
        challenger="TACTICS Enhanced-TT-TS (GMIC)",
        baseline="Legacy Enhanced-RWS",
        challenger_color="#E45756",
        baseline_color="#B07AA1",
        challenger_label="TACTICS TT-TS (GMIC)",
        baseline_label="Legacy RWS (Zhao 2025)",
        title="Top-100 Recovery vs Library Diversity — TACTICS TT-TS vs Legacy RWS (Zhao 2025)",
    )
    mo.vstack([
        plt.gcf(),
        mo.download(data=fig_to_pdf_bytes(_fig),
                    filename="zhao_fig5_tt_ts_vs_legacy_rws.pdf",
                    mimetype="application/pdf", label="Download as PDF"),
    ])
    return


@app.cell
def _(fig_to_pdf_bytes, mo, plt, zhao_pairwise_plot):
    """Zhao Fig 5 — TT-TS vs Legacy Greedy (Klarich 2024 TS)."""
    _fig = zhao_pairwise_plot(
        challenger="TACTICS Enhanced-TT-TS (GMIC)",
        baseline="Legacy Standard-Greedy",
        challenger_color="#E45756",
        baseline_color="#76B7B2",
        challenger_label="TACTICS TT-TS (GMIC)",
        baseline_label="Legacy Greedy (Klarich 2024 TS)",
        title="Top-100 Recovery vs Library Diversity — TACTICS TT-TS vs Legacy Greedy (Klarich TS)",
    )
    mo.vstack([
        plt.gcf(),
        mo.download(data=fig_to_pdf_bytes(_fig),
                    filename="zhao_fig5_tt_ts_vs_legacy_greedy.pdf",
                    mimetype="application/pdf", label="Download as PDF"),
    ])
    return


@app.cell
def _(fig_to_pdf_bytes, mo, plt, zhao_delta_plot):
    """Zhao Fig 5 Δ — TACTICS RWS minus Legacy RWS (Zhao 2025) vs diversity."""
    _fig = zhao_delta_plot(
        challenger="TACTICS Enhanced-RWS (GMIC)",
        baseline="Legacy Enhanced-RWS",
        challenger_color="#4C78A8",
        challenger_label="TACTICS RWS (GMIC)",
        baseline_label="Legacy RWS (Zhao 2025)",
        title="Δ Top-100 Recovery vs Library Diversity — TACTICS RWS − Legacy RWS (Zhao 2025)",
    )
    mo.vstack([
        plt.gcf(),
        mo.download(data=fig_to_pdf_bytes(_fig),
                    filename="zhao_fig5_delta_tactics_rws_vs_legacy_rws.pdf",
                    mimetype="application/pdf", label="Download as PDF"),
    ])
    return


@app.cell
def _(fig_to_pdf_bytes, mo, plt, zhao_delta_plot):
    """Zhao Fig 5 Δ — TACTICS RWS minus Legacy Greedy (Klarich TS) vs diversity."""
    _fig = zhao_delta_plot(
        challenger="TACTICS Enhanced-RWS (GMIC)",
        baseline="Legacy Standard-Greedy",
        challenger_color="#4C78A8",
        challenger_label="TACTICS RWS (GMIC)",
        baseline_label="Legacy Greedy (Klarich 2024 TS)",
        title="Δ Top-100 Recovery vs Library Diversity — TACTICS RWS − Legacy Greedy (Klarich TS)",
    )
    mo.vstack([
        plt.gcf(),
        mo.download(data=fig_to_pdf_bytes(_fig),
                    filename="zhao_fig5_delta_tactics_rws_vs_legacy_greedy.pdf",
                    mimetype="application/pdf", label="Download as PDF"),
    ])
    return


@app.cell
def _(fig_to_pdf_bytes, mo, plt, zhao_delta_plot):
    """Zhao Fig 5 Δ — TT-TS minus Legacy RWS (Zhao 2025) vs diversity."""
    _fig = zhao_delta_plot(
        challenger="TACTICS Enhanced-TT-TS (GMIC)",
        baseline="Legacy Enhanced-RWS",
        challenger_color="#E45756",
        challenger_label="TACTICS TT-TS (GMIC)",
        baseline_label="Legacy RWS (Zhao 2025)",
        title="Δ Top-100 Recovery vs Library Diversity — TACTICS TT-TS − Legacy RWS (Zhao 2025)",
    )
    mo.vstack([
        plt.gcf(),
        mo.download(data=fig_to_pdf_bytes(_fig),
                    filename="zhao_fig5_delta_tt_ts_vs_legacy_rws.pdf",
                    mimetype="application/pdf", label="Download as PDF"),
    ])
    return


@app.cell
def _(fig_to_pdf_bytes, mo, plt, zhao_delta_plot):
    """Zhao Fig 5 Δ — TT-TS minus Legacy Greedy (Klarich TS) vs diversity."""
    _fig = zhao_delta_plot(
        challenger="TACTICS Enhanced-TT-TS (GMIC)",
        baseline="Legacy Standard-Greedy",
        challenger_color="#E45756",
        challenger_label="TACTICS TT-TS (GMIC)",
        baseline_label="Legacy Greedy (Klarich 2024 TS)",
        title="Δ Top-100 Recovery vs Library Diversity — TACTICS TT-TS − Legacy Greedy (Klarich TS)",
    )
    mo.vstack([
        plt.gcf(),
        mo.download(data=fig_to_pdf_bytes(_fig),
                    filename="zhao_fig5_delta_tt_ts_vs_legacy_greedy.pdf",
                    mimetype="application/pdf", label="Download as PDF"),
    ])
    return


@app.cell
def _(mo):
    mo.md(r"""
    # Section 3 — Tukey HSD (ROCS libraries) + Per-library ROCS plots

    Tukey HSD pairwise comparisons across the 20 ROCS libraries, plus per-library
    Top-N curves + best-TACTICS Δ bars with Tukey HSD significance stars (per the
    scoring-type rule — ROCS satisfies repeated-measures structure required for Tukey).
    """)
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
def _(df, mo, pl):
    """Library selector for per-library Tukey HSD CI plot (ROCS libraries only)."""
    _DOCKING_LIBS = {
        "thrombin", "adenine", "adenine_hybrid", "amide", "quinazoline",
        "thrombin_hybrid", "amide_hybrid", "quinazoline_hybrid",
    }
    _libraries = sorted(
        df.filter(~pl.col("library_id").is_in(_DOCKING_LIBS))["library_id"]
        .unique().to_list()
    )
    library_dropdown = mo.ui.dropdown(
        options=_libraries,
        value=_libraries[0] if _libraries else "",
        label="ROCS library",
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
    _DOCKING_LIBS = {
        "thrombin", "adenine", "adenine_hybrid", "amide", "quinazoline",
        "thrombin_hybrid", "amide_hybrid", "quinazoline_hybrid",
    }
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
def _(mo, recovery_all):
    """Per-library ROCS recovery — library selector (20 ROCS libraries only)."""
    _ROCS_LIBS = sorted(
        lib for lib in recovery_all["library_id"].unique().to_list()
        if lib not in {
            "thrombin", "adenine", "adenine_hybrid", "amide", "quinazoline",
            "thrombin_hybrid", "amide_hybrid", "quinazoline_hybrid",
        }
    )
    rocs_library_selector = mo.ui.dropdown(
        options=_ROCS_LIBS,
        value=_ROCS_LIBS[0] if _ROCS_LIBS else None,
        label="ROCS library",
    )
    mo.vstack([
        mo.md(r"""
        ## Per-library ROCS recovery — Top-N curves + best-TACTICS Δ bars

        Pick a ROCS library to view its 4-method Top-N recovery curves (left panel)
        plus best-TACTICS Δ-recovery bars vs Legacy Greedy / Legacy RWS (right panel).
        Significance stars use **Tukey HSD** (per the scoring-type rule: ROCS libraries
        have 109 queries × 20 replicates → repeated-measures structure → Tukey HSD applies).
        """),
        rocs_library_selector,
    ])
    return (rocs_library_selector,)


@app.cell
def _(
    TUKEY_STAR_COLOR,
    fig_to_pdf_bytes,
    mo,
    np,
    pl,
    plt,
    project_root,
    recovery_all,
    rocs_library_selector,
):
    """Per-library ROCS bar+line plot — Top-N curves + TT-TS-vs-Legacy Δ bars with Tukey HSD stars.

    Bar height = mean(TT-TS) − mean(Legacy) on this library at this top-N.
    On ROCS the two TACTICS methods are practically equivalent at the per-library
    level (max gap 0.79 pp ≈ 0.2 compounds, on `poparov` at top-25), so headlining
    TT-TS is the simpler and equally accurate story.
    """

    _LIB = rocs_library_selector.value
    _CHALLENGERS = [
        ("TACTICS Enhanced-TT-TS (GMIC)", "TT-TS", "#E45756", "D"),
        ("TACTICS Enhanced-RWS (GMIC)", "TACTICS RWS", "#4C78A8", "o"),
    ]
    _LEGACIES = [
        ("Legacy Standard-Greedy", "Legacy Greedy", "#76B7B2", "s", ""),
        ("Legacy Enhanced-RWS", "Legacy RWS", "#B07AA1", "^", "////"),
    ]

    # Tukey HSD significance results for this library (ROCS uses Tukey per the rule).
    _tukey_all = pl.read_parquet(
        project_root / "outputs" / "full_benchmark_ema_cats" / "analysis" / "tukey_hsd_results.parquet"
    ).filter(pl.col("library_id") == _LIB)

    _methods_keep = {m for m, _, _, _ in _CHALLENGERS} | {m for m, _, _, _, _ in _LEGACIES}
    _lib_all = recovery_all.filter(
        (pl.col("library_id") == _LIB) & pl.col("method").is_in(_methods_keep)
    )
    _topn_vals = sorted([_tn for _tn in _lib_all["top_n"].unique().to_list() if _tn <= 250])

    mo.stop(_lib_all.height == 0, mo.md(f"No data for {_LIB}."))

    _fig, (_ax1, _ax2) = plt.subplots(1, 2, figsize=(16, 6), facecolor="white")
    _ax1.set_facecolor("white")

    _agg = (
        _lib_all.group_by(["method", "top_n"])
        .agg([
            pl.col("recovery_frac").mean().alias("mean_frac"),
            pl.col("recovery_frac").std().alias("std_frac"),
            pl.col("recovery_frac").count().alias("n"),
        ])
        .with_columns((pl.col("std_frac") / pl.col("n").sqrt()).alias("se_frac"))
        .sort(["method", "top_n"])
    )

    for _meth, _lbl, _clr, _mkr, _ in _LEGACIES:
        _md = _agg.filter(pl.col("method") == _meth).sort("top_n")
        if len(_md) == 0:
            continue
        _x = _md["top_n"].to_numpy()
        _y = _md["mean_frac"].to_numpy() * 100
        _se = _md["se_frac"].to_numpy() * 100
        _ax1.plot(_x, _y, marker=_mkr, color=_clr, linewidth=2, markersize=7, label=_lbl, zorder=3)
        _ax1.fill_between(_x, _y - 1.96 * _se, _y + 1.96 * _se, alpha=0.15, color=_clr, zorder=2)

    for _meth, _lbl, _clr, _mkr in _CHALLENGERS:
        _md = _agg.filter(pl.col("method") == _meth).sort("top_n")
        if len(_md) == 0:
            continue
        _x = _md["top_n"].to_numpy()
        _y = _md["mean_frac"].to_numpy() * 100
        _se = _md["se_frac"].to_numpy() * 100
        _ax1.plot(_x, _y, marker=_mkr, color=_clr, linewidth=2, markersize=7, label=_lbl, zorder=4)
        _ax1.fill_between(_x, _y - 1.96 * _se, _y + 1.96 * _se, alpha=0.15, color=_clr, zorder=2)

    _ax1.set_xlabel("Top-N Threshold", fontsize=14, color="black")
    _ax1.set_ylabel("% Recovery", fontsize=14, color="black")
    _ax1.set_title(f"{_LIB} (ROCS) — Recovery vs Top-N", fontsize=14, fontweight="bold", color="black")
    _ax1.set_xticks(_topn_vals)
    _ax1.tick_params(colors="black", labelsize=12)
    _ax1.grid(alpha=0.3, color="grey", zorder=1)
    _ax1.legend(fontsize=11, facecolor="white", edgecolor="black", labelcolor="black")
    for _sp in _ax1.spines.values():
        _sp.set_color("black")

    _BAR_COLORS = ["#F28E2B", "#4E79A7"]  # orange (vs Greedy), steel blue (vs Legacy RWS)
    _ax2.set_facecolor("white")
    _bar_width = 0.35
    _x_pos = np.arange(len(_topn_vals))

    def _tukey_sig(library, top_n, method_a, method_b):
        """Look up Tukey HSD significance for (method_a, method_b) at top_n; try both directions."""
        _row = _tukey_all.filter(
            (pl.col("top_n") == top_n)
            & (
                ((pl.col("group1") == method_a) & (pl.col("group2") == method_b))
                | ((pl.col("group1") == method_b) & (pl.col("group2") == method_a))
            )
        )
        return _row["significance"][0] if _row.height > 0 else ""

    _HEADLINE_TACTICS = "TACTICS Enhanced-TT-TS (GMIC)"

    for _gi, (_leg_meth, _leg_lbl, _, _, _) in enumerate(_LEGACIES):
        _diffs = []
        _sigs = []
        for _tn in _topn_vals:
            _tt_frac = _lib_all.filter(
                (pl.col("method") == _HEADLINE_TACTICS) & (pl.col("top_n") == _tn)
            )["recovery_frac"].mean()
            _leg_frac = _lib_all.filter(
                (pl.col("method") == _leg_meth) & (pl.col("top_n") == _tn)
            )["recovery_frac"].mean()
            _diffs.append(
                float((_tt_frac - _leg_frac) * 100)
                if _tt_frac is not None and _leg_frac is not None else 0
            )
            _sigs.append(_tukey_sig(_LIB, _tn, _HEADLINE_TACTICS, _leg_meth))

        _diffs = np.array(_diffs)
        _offset = (_gi - 0.5) * _bar_width
        _ax2.bar(_x_pos + _offset, _diffs, width=_bar_width, color=_BAR_COLORS[_gi],
                 edgecolor="black", linewidth=0.5, zorder=3, label=f"vs {_leg_lbl}")

        _dmax = max(abs(_diffs)) if len(_diffs) > 0 else 1
        _star_pad = _dmax * 0.03
        for _i, (_d, _s) in enumerate(zip(_diffs, _sigs)):
            if _s and _s != "ns":
                _yp = _d + _star_pad if _d >= 0 else _d - _star_pad
                _ax2.text(_x_pos[_i] + _offset, _yp, _s, ha="center",
                          va="bottom" if _d >= 0 else "top",
                          fontsize=9, fontweight="bold", color=TUKEY_STAR_COLOR,
                          rotation=90)

    _ax2.axhline(0, color="black", linewidth=0.8, zorder=2)
    _yl, _yh = _ax2.get_ylim()
    _ax2.set_ylim(_yl - abs(_yh - _yl) * 0.10, _yh + abs(_yh - _yl) * 0.30)
    _ax2.set_xticks(_x_pos)
    _ax2.set_xticklabels([str(t) for t in _topn_vals], fontsize=12, color="black")
    _ax2.set_xlabel("Top-N Threshold", fontsize=14, color="black")
    _ax2.set_ylabel(
        "Δ % Recovery (TACTICS TT-TS − Legacy)",
        fontsize=14, color="black",
    )
    _ax2.set_title(
        f"{_LIB} (ROCS) — TACTICS TT-TS Advantage",
        fontsize=14, fontweight="bold", color="black",
    )
    _ax2.tick_params(colors="black", labelsize=12)
    _ax2.grid(axis="y", alpha=0.3, color="grey", zorder=1)
    _ax2.legend(fontsize=11, facecolor="white", edgecolor="black", labelcolor="black")
    for _sp in _ax2.spines.values():
        _sp.set_color("black")
    _fig.tight_layout()
    mo.vstack([
        _fig,
        mo.download(data=fig_to_pdf_bytes(_fig),
                    filename=f"{_LIB}_topn_tactics_vs_legacy_ROCS.pdf",
                    mimetype="application/pdf", label="Download PDF"),
    ])
    return


@app.cell
def _(mo):
    mo.md(r"""
    # Section 4 — Budget Sensitivity Analysis (1% → 5%) — ROCS

    Recovery and TACTICS-vs-Legacy gap as a function of library-fraction screening
    budget for ROCS libraries, at the canonical top-100 threshold. ROCS aggregates
    are tight (uniform per-library response within each component-count cell);
    the parallel DOCK sensitivity panel lives in `manuscript_plots.py`.

    2×2 panel grid:
    - **Top row** (line plots): per-method mean recovery vs budget fraction.
    - **Bottom row** (bar plots): Δ(best-TACTICS − best-Legacy) per budget with error bars.
    - **Left column**: 2-component libraries. **Right column**: 3-component libraries.

    Recovery is averaged per query × replicate within each library, then aggregated
    across libraries within each (n_comp) cell.
    """)
    return


@app.cell
def _(pl, project_root):
    """Load top-100 recovery from all 5 budget benchmarks into a single dataframe."""
    _budgets = [
        (1, project_root / "outputs/full_benchmark_1pct/analysis/recovery_summary.parquet"),
        (2, project_root / "outputs/full_benchmark_ema_cats/analysis/recovery_summary.parquet"),
        (3, project_root / "outputs/full_benchmark_3pct/analysis/recovery_summary.parquet"),
        (4, project_root / "outputs/full_benchmark_4pct/analysis/recovery_summary.parquet"),
        (5, project_root / "outputs/full_benchmark_5pct/analysis/recovery_summary.parquet"),
    ]
    _frames = []
    for _budget, _path in _budgets:
        if not _path.exists():
            continue
        _sub = (
            pl.read_parquet(_path)
            .filter(pl.col("top_n") == 100)
            .with_columns(pl.lit(_budget).alias("budget_pct"))
        )
        _frames.append(_sub)
    sensitivity_data = pl.concat(_frames)
    return (sensitivity_data,)


@app.cell
def _(pl, sensitivity_data):
    """Compute per-(library, method, budget) mean recovery + per-(n_comp, method, budget) aggregates.

    Pipeline:
      1. Per (library, query, replicate, method, budget) → already in recovery_summary
      2. Mean over (query, replicate) → per (library, method, budget) library-level mean
      3. Mean across libraries in each (n_comp, method, budget) cell → aggregate mean + std
      4. Per-library Δ = best-T(lib) − best-L(lib); aggregate Δ mean + std across libs

    Returns aggregate dataframes ready for plotting.
    """
    _METHOD_SHORT = {
        "TACTICS Enhanced-TT-TS (GMIC)": "TT-TS",
        "TACTICS Enhanced-RWS (GMIC)":   "T-RWS",
        "Legacy Enhanced-RWS":           "L-RWS",
        "Legacy Standard-Greedy":        "L-Grdy",
    }
    _ROCS_NCOMP = {
        "rxn101": 2, "rxn102a": 2, "rxn108b": 2, "rxn111b": 2, "rxn114b": 2,
        "rxn203": 2, "rxn205": 2, "rxn206": 2, "rxn207": 2, "rxn208": 2,
        "amide-suzuki": 3, "betti": 3, "dobener": 3, "groebke-blackburn-bienayme": 3,
        "mannich": 3, "niementowski": 3, "orru": 3, "passerini": 3, "petasis": 3, "poparov": 3,
    }
    _DOCK_NCOMP = {
        "thrombin": 2, "thrombin_hybrid": 2, "amide": 2, "amide_hybrid": 2,
        "adenine": 3, "adenine_hybrid": 3, "quinazoline": 3, "quinazoline_hybrid": 3,
    }

    # Add kind + n_comp to the data
    _ncomp_map = {**_ROCS_NCOMP, **_DOCK_NCOMP}
    _kind_map = {lib: "ROCS" for lib in _ROCS_NCOMP}
    _kind_map.update({lib: "DOCK" for lib in _DOCK_NCOMP})

    _df = (
        sensitivity_data
        .with_columns(pl.col("method").replace_strict(_METHOD_SHORT).alias("m"))
        .filter(pl.col("library_id").is_in(set(_ncomp_map.keys())))
        .with_columns([
            pl.col("library_id").replace(_ncomp_map).cast(pl.Int64).alias("n_comp"),
            pl.col("library_id").replace(_kind_map).alias("kind"),
        ])
    )

    # Step 1: per-library mean across (query, replicate)
    sens_per_lib = (
        _df.group_by(["library_id", "kind", "n_comp", "m", "budget_pct"])
        .agg(pl.col("recovery_frac").mean().alias("rec"))
    )

    # Step 2: per-method aggregate across libraries within each (kind, n_comp, budget)
    sens_per_method = (
        sens_per_lib.group_by(["kind", "n_comp", "m", "budget_pct"])
        .agg([
            pl.col("rec").mean().alias("mean_rec"),
            pl.col("rec").std().alias("std_rec"),
            pl.col("rec").count().alias("n_libs"),
        ])
        .with_columns(
            (pl.col("std_rec") / pl.col("n_libs").sqrt()).alias("se_rec")
        )
    )

    # Step 3: per-library Δ = best-T − best-L
    _wide = (
        sens_per_lib.pivot(
            values="rec", index=["library_id", "kind", "n_comp", "budget_pct"], on="m"
        )
    )
    _delta_per_lib = _wide.with_columns([
        pl.max_horizontal(["TT-TS", "T-RWS"]).alias("best_T"),
        pl.max_horizontal(["L-RWS", "L-Grdy"]).alias("best_L"),
    ]).with_columns(
        ((pl.col("best_T") - pl.col("best_L")) * 100).alias("delta_pp")
    )
    sens_delta_per_lib = _delta_per_lib

    # Step 4: aggregate Δ per (kind, n_comp, budget)
    sens_delta_agg = (
        _delta_per_lib.group_by(["kind", "n_comp", "budget_pct"])
        .agg([
            pl.col("delta_pp").mean().alias("mean_delta"),
            pl.col("delta_pp").std().alias("std_delta"),
            pl.col("delta_pp").count().alias("n_libs"),
        ])
        .with_columns(
            (pl.col("std_delta") / pl.col("n_libs").sqrt()).alias("se_delta")
        )
        .sort(["kind", "n_comp", "budget_pct"])
    )
    return sens_delta_agg, sens_per_method


@app.cell
def _(fig_to_pdf_bytes, mo, pl, plt, sens_delta_agg, sens_per_method):
    """ROCS sensitivity figure: 2×2 grid (line plots top, bar plots bottom; 2-comp left, 3-comp right)."""
    _METHOD_ORDER = ["TT-TS", "T-RWS", "L-RWS", "L-Grdy"]
    _METHOD_LABELS = {
        "TT-TS": "TACTICS TT-TS (GMIC)",
        "T-RWS": "TACTICS RWS (GMIC)",
        "L-RWS": "Legacy Enhanced-RWS",
        "L-Grdy": "Legacy Standard-Greedy",
    }
    _METHOD_COLORS = {
        "TT-TS": "#E45756", "T-RWS": "#4C78A8",
        "L-RWS": "#B07AA1", "L-Grdy": "#76B7B2",
    }
    _METHOD_MARKERS = {
        "TT-TS": "D", "T-RWS": "o", "L-RWS": "^", "L-Grdy": "s",
    }

    _fig, _axes = plt.subplots(2, 2, figsize=(14, 9), facecolor="white", sharex=True)

    for _col_idx, _nc in enumerate([2, 3]):
        _ax_line = _axes[0, _col_idx]
        _ax_bar = _axes[1, _col_idx]
        _ax_line.set_facecolor("white"); _ax_bar.set_facecolor("white")

        # --- Top row: line plot of recovery vs budget ---
        _sub_method = sens_per_method.filter(
            (pl.col("kind") == "ROCS") & (pl.col("n_comp") == _nc)
        )
        _n_libs = int(_sub_method["n_libs"][0]) if _sub_method.height > 0 else 0
        for _m in _METHOD_ORDER:
            _md = _sub_method.filter(pl.col("m") == _m).sort("budget_pct")
            if _md.height == 0: continue
            _x = _md["budget_pct"].to_numpy()
            _y = _md["mean_rec"].to_numpy() * 100
            _se = _md["se_rec"].to_numpy() * 100
            _ax_line.errorbar(
                _x, _y, yerr=1.96 * _se, color=_METHOD_COLORS[_m],
                marker=_METHOD_MARKERS[_m], markersize=7, linewidth=2,
                capsize=4, capthick=1.2, label=_METHOD_LABELS[_m], zorder=3,
            )
        _ax_line.set_ylabel("Mean Top-100 Recovery (%)", fontsize=12, color="black")
        _ax_line.set_title(f"ROCS {_nc}-Component (n = {_n_libs} libraries)",
                            fontsize=13, fontweight="bold", color="black")
        _ax_line.grid(alpha=0.3, color="grey", zorder=1)
        _ax_line.tick_params(colors="black", labelsize=11)
        for _sp in _ax_line.spines.values(): _sp.set_color("black")

        # --- Bottom row: Δ bar plot ---
        _sub_delta = sens_delta_agg.filter(
            (pl.col("kind") == "ROCS") & (pl.col("n_comp") == _nc)
        ).sort("budget_pct")
        _x = _sub_delta["budget_pct"].to_numpy()
        _y = _sub_delta["mean_delta"].to_numpy()
        _se = _sub_delta["se_delta"].to_numpy()
        _ax_bar.bar(_x, _y, yerr=1.96 * _se, color="#4C78A8", alpha=0.85,
                     edgecolor="black", linewidth=0.8, capsize=5,
                     error_kw={"linewidth": 1.2, "color": "black"}, zorder=3)
        _ax_bar.axhline(0, color="black", linewidth=0.8, zorder=2)
        _ax_bar.set_xlabel("Library Fraction (%)", fontsize=12, color="black")
        _ax_bar.set_ylabel("Δ Best-TACTICS − Best-Legacy (pp)", fontsize=12, color="black")
        _ax_bar.set_xticks([1, 2, 3, 4, 5])
        _ax_bar.grid(axis="y", alpha=0.3, color="grey", zorder=1)
        _ax_bar.tick_params(colors="black", labelsize=11)
        for _sp in _ax_bar.spines.values(): _sp.set_color("black")
        # Annotate bar values
        for _xi, _yi in zip(_x, _y):
            _ax_bar.text(_xi, _yi + 0.15, f"{_yi:+.2f}", ha="center", va="bottom",
                          fontsize=10, color="black", fontweight="bold")

    _handles, _labels = _axes[0, 0].get_legend_handles_labels()
    _fig.legend(_handles, _labels, loc="lower center", ncol=5, fontsize=11,
                 facecolor="white", edgecolor="black", labelcolor="black",
                 bbox_to_anchor=(0.5, -0.04))
    _fig.suptitle("Budget Sensitivity — ROCS Libraries (Top-100 Recovery)",
                   fontsize=14, fontweight="bold", color="black", y=1.00)
    _fig.tight_layout(rect=[0, 0.04, 1, 0.98])

    mo.vstack([
        plt.gcf(),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename="sensitivity_rocs.pdf",
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


if __name__ == "__main__":
    app.run()
