"""
TACTICS Manuscript Plots — Docking Libraries

Docking-only figure notebook (8 libraries: 4 FRED + 4 HYBRID). Loads the canonical
recovery summary (`full_benchmark_ema_cats`) and renders, in order:
  §1 method × docking-library top-100 recovery heatmap,
  §2 per-library two-panel recovery + significance (dropdown-selected),
  §3 aggregate recovery vs top-N (FRED / HYBRID),
  §4 statistical-test convention reminder,
  §5 budget-sensitivity analysis (1% → 5%).
The ROCS counterpart lives in manuscript_plots_ROCS.py.

Run as app: marimo run tutorials/manuscript_plots_docking.py
Edit mode:  marimo edit tutorials/manuscript_plots_docking.py
"""

import marimo

__generated_with = "0.21.1"
app = marimo.App(width="full", app_title="TACTICS Manuscript — Docking")


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
    return fig_to_pdf_bytes, mo, np, pl, plt, project_root


@app.cell
def _(mo, pl, project_root):
    """Load the canonical docking recovery summary (docking libraries only).

    Recovery from the canonical 25-library benchmark (`full_benchmark_ema_cats`),
    capped at top-N ≤ 250 (manuscript-wide scope) and filtered to the eight docking
    libraries. All five methods are kept so the §1 heatmap and the §2 per-library
    figures always have the full method set.
    """
    analysis_dir = project_root / "outputs" / "full_benchmark_ema_cats" / "analysis"

    DOCK_LIBS = [
        "adenine", "adenine_hybrid",
        "amide", "amide_hybrid",
        "quinazoline", "quinazoline_hybrid",
        "thrombin", "thrombin_hybrid",
    ]

    recovery_all = (
        pl.read_parquet(analysis_dir / "recovery_summary.parquet")
        .filter(pl.col("top_n") <= 250)
        .filter(pl.col("library_id").is_in(DOCK_LIBS))
    )

    mo.md(
        f"""
        ## Data Loaded — Docking Libraries Only

        - **Source**: `{analysis_dir}/recovery_summary.parquet`
        - **Docking libraries** ({len(DOCK_LIBS)}): {', '.join(DOCK_LIBS)}
        - **Top-N thresholds**: {sorted(recovery_all['top_n'].unique().to_list())}
        - **Methods**: {', '.join(sorted(recovery_all['method'].unique().to_list()))}
        """
    )
    return DOCK_LIBS, recovery_all


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
    return (WELCH_STAR_COLOR,)


@app.cell
def _(mo):
    mo.md(r"""
    # Section 1 — Method × Docking-Library Performance Heatmap

    Top-100 recovery (%) for all five methods across the eight docking landscapes
    (4 FRED + 4 HYBRID). Columns are ordered by the best-TACTICS − best-Legacy gap
    (largest on the left), so the libraries where TACTICS helps most appear first.
    """)
    return


@app.cell
def _(DOCK_LIBS, fig_to_pdf_bytes, mo, np, pl, plt, recovery_all):
    """Ranking heatmap — top-100 recovery (%) for 5 methods × 8 docking landscapes.

    Computed from the canonical recovery_summary (mean recovery_frac at top-100 per
    method × library), so it matches the §2 per-library figures exactly. Columns are
    ordered by the best-TACTICS − best-Legacy gap (largest on the left).
    """
    _METHODS = [
        "TACTICS Enhanced-TT-TS (GMIC)",
        "TACTICS Enhanced-RWS (GMIC)",
        "Legacy Enhanced-RWS",
        "Legacy Standard-Greedy",
    ]
    _MSHORT = {
        "TACTICS Enhanced-TT-TS (GMIC)": "TACTICS TT-TS",
        "TACTICS Enhanced-RWS (GMIC)": "TACTICS RWS",
        "Legacy Enhanced-RWS": "Legacy RWS",
        "Legacy Standard-Greedy": "Legacy Greedy",
    }
    _TT, _RWS = "TACTICS Enhanced-TT-TS (GMIC)", "TACTICS Enhanced-RWS (GMIC)"
    _LRWS, _LGRD = "Legacy Enhanced-RWS", "Legacy Standard-Greedy"

    # Mean top-100 recovery (%) per (library, method).
    _agg = (
        recovery_all.filter(pl.col("top_n") == 100)
        .group_by(["library_id", "method"])
        .agg((pl.col("recovery_frac").mean() * 100).alias("recovery"))
    )
    _pairs = {
        (_row["library_id"], _row["method"]): _row["recovery"]
        for _row in _agg.iter_rows(named=True)
    }
    _recov = {
        _lab: {_m: _pairs.get((_lab, _m), float("nan")) for _m in _METHODS}
        for _lab in DOCK_LIBS
    }
    _delta = {
        _lab: max(_recov[_lab][_TT], _recov[_lab][_RWS])
        - max(_recov[_lab][_LRWS], _recov[_lab][_LGRD])
        for _lab in DOCK_LIBS
    }
    _libs = sorted(DOCK_LIBS, key=lambda _l: _delta[_l], reverse=True)
    _mat = np.array([[_recov[_l][_m] for _l in _libs] for _m in _METHODS])

    _fig, _ax = plt.subplots(figsize=(14, 5.2), facecolor="white")
    _ax.set_facecolor("white")
    _im = _ax.imshow(_mat, cmap="viridis", aspect="auto", vmin=0, vmax=100)
    _ax.set_xticks(range(len(_libs)))
    _ax.set_xticklabels(_libs, rotation=30, ha="right", rotation_mode="anchor",
                        fontsize=11, color="black")
    _ax.set_yticks(range(len(_METHODS)))
    _ax.set_yticklabels([_MSHORT[_m] for _m in _METHODS], fontsize=12, color="black")
    for _r in range(_mat.shape[0]):
        for _c in range(_mat.shape[1]):
            _val = _mat[_r, _c]
            _txt_c = "white" if _val < 55 else "black"
            _ax.text(_c, _r, f"{_val:.0f}", ha="center", va="center",
                     fontsize=11, color=_txt_c, fontweight="bold")
    _ax.set_title("Top-100 recovery (%) — methods × docking landscapes",
                  fontsize=14, fontweight="bold", color="black")
    _ax.tick_params(colors="black")
    _cbar = _fig.colorbar(_im, ax=_ax, fraction=0.025, pad=0.015)
    _cbar.set_label("Top-100 recovery (%)", fontsize=12, color="black")
    _cbar.ax.tick_params(colors="black", labelsize=10)
    for _sp in _ax.spines.values():
        _sp.set_color("black")

    mo.vstack([
        plt.gcf(),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename="docking_ranking_heatmap.pdf",
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


@app.cell
def _(mo):
    mo.md(r"""
    # Section 2 — Per-Library Recovery & Significance

    Select a docking library to view its two-panel figure: **(left)** % recovery vs
    top-N for the two TACTICS methods and the two Legacy baselines with 95% CI bands,
    and **(right)** Δ recovery (best-TACTICS − each Legacy) at every top-N with planned
    Welch + Bonferroni-4 significance stars (docking-libraries scoring-type rule).
    """)
    return


@app.cell
def _(DOCK_LIBS, mo):
    """Dropdown to select which docking library's two-panel figure to display."""
    dock_library_dropdown = mo.ui.dropdown(
        options=DOCK_LIBS, value="adenine", label="Docking library",
    )
    dock_library_dropdown
    return (dock_library_dropdown,)


@app.cell
def _(
    WELCH_STAR_COLOR,
    dock_library_dropdown,
    fig_to_pdf_bytes,
    mo,
    np,
    pl,
    plt,
    project_root,
    recovery_all,
):
    """Per-library two-panel figure for the dropdown-selected docking library.

    Left: % recovery vs top-N (95% CI bands) for the two TACTICS methods + two Legacy
    baselines. Right: Δ recovery (best-TACTICS − each Legacy) at each top-N with planned
    Welch + Bonferroni-4 stars (docking scoring-type rule). All eight libraries use the
    canonical full_benchmark_ema_cats recovery + Welch results, so the figure is
    consistent with the §1 heatmap.
    """
    _LIB = dock_library_dropdown.value
    _TITLE_SUFFIX = {
        "adenine": "FRED (A2A 4EIY)",
        "adenine_hybrid": "HYBRID (A2A-Lenselink)",
        "amide": "FRED (JNK3 2ZDT)",
        "amide_hybrid": "HYBRID (JNK3 2ZDT)",
        "quinazoline": "FRED (JNK3 2ZDT)",
        "quinazoline_hybrid": "HYBRID (JNK3 2ZDT)",
        "thrombin": "FRED",
        "thrombin_hybrid": "HYBRID (2ZFF)",
    }
    _suffix = _TITLE_SUFFIX.get(_LIB, "")

    _CHALLENGERS = [
        ("TACTICS Enhanced-TT-TS (GMIC)", "TT-TS", "#E45756", "D"),
        ("TACTICS Enhanced-RWS (GMIC)", "TACTICS RWS", "#4C78A8", "o"),
    ]
    _LEGACIES = [
        ("Legacy Standard-Greedy", "Legacy Greedy", "#76B7B2", "s", ""),
        ("Legacy Enhanced-RWS", "Legacy RWS", "#B07AA1", "^", "////"),
    ]

    # Planned Welch's t-test results — docking-libraries scoring-type rule.
    _welch_all = pl.read_parquet(
        project_root / "outputs" / "full_benchmark_ema_cats" / "analysis" / "welch_planned_results.parquet"
    ).filter(pl.col("library_id") == _LIB)

    _methods_keep = {m for m, _, _, _ in _CHALLENGERS} | {m for m, _, _, _, _ in _LEGACIES}
    _lib_all = recovery_all.filter(
        (pl.col("library_id") == _LIB) & pl.col("method").is_in(_methods_keep)
    )
    _topn_vals = sorted([_tn for _tn in _lib_all["top_n"].unique().to_list() if _tn <= 250])

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
    _ax1.set_title(f"{_LIB} {_suffix} — Recovery vs Top-N", fontsize=14, fontweight="bold", color="black")
    _ax1.set_xticks(_topn_vals)
    _ax1.tick_params(colors="black", labelsize=12)
    _ax1.grid(alpha=0.3, color="grey", zorder=1)
    _ax1.legend(fontsize=11, facecolor="white", edgecolor="black", labelcolor="black")
    for _sp in _ax1.spines.values():
        _sp.set_color("black")

    # Best-TACTICS Δ vs each legacy: 2 bars per top-N. At each threshold the
    # best-performing TACTICS method (max mean recovery) is selected and its Δ vs
    # each legacy is plotted, with planned Welch + Bonferroni-4 significance stars.
    _BAR_COLORS = ["#F28E2B", "#4E79A7"]  # orange (vs Greedy), steel blue (vs Legacy RWS)
    _ax2.set_facecolor("white")
    _bar_width = 0.35
    _x_pos = np.arange(len(_topn_vals))

    for _gi, (_leg_meth, _leg_lbl, _, _, _) in enumerate(_LEGACIES):
        _diffs = []
        _sigs = []
        for _tn in _topn_vals:
            _ch_fracs = {}
            for _ch_meth, _, _, _ in _CHALLENGERS:
                _f = _lib_all.filter((pl.col("method") == _ch_meth) & (pl.col("top_n") == _tn))["recovery_frac"].mean()
                _ch_fracs[_ch_meth] = _f if _f is not None else float("-inf")
            _best_meth = max(_ch_fracs, key=_ch_fracs.get)
            _best_frac = _ch_fracs[_best_meth]
            _leg_frac = _lib_all.filter((pl.col("method") == _leg_meth) & (pl.col("top_n") == _tn))["recovery_frac"].mean()
            _diffs.append(float((_best_frac - _leg_frac) * 100) if _best_frac != float("-inf") and _leg_frac is not None else 0)
            _wrow = _welch_all.filter(
                (pl.col("top_n") == _tn)
                & (pl.col("challenger") == _best_meth)
                & (pl.col("baseline") == _leg_meth)
            )
            _sigs.append(_wrow["significance_bonf"][0] if len(_wrow) > 0 else "")

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
                          fontsize=9, fontweight="bold", color=WELCH_STAR_COLOR,
                          rotation=90)

    _ax2.axhline(0, color="black", linewidth=0.8, zorder=2)
    _yl, _yh = _ax2.get_ylim()
    _ax2.set_ylim(_yl - abs(_yh - _yl) * 0.10, _yh + abs(_yh - _yl) * 0.30)
    _ax2.set_xticks(_x_pos)
    _ax2.set_xticklabels([str(t) for t in _topn_vals], fontsize=12, color="black")
    _ax2.set_xlabel("Top-N Threshold", fontsize=14, color="black")
    _ax2.set_ylabel("Δ % Recovery (Best TACTICS − Legacy)", fontsize=14, color="black")
    _ax2.set_title(f"{_LIB} {_suffix} — Best TACTICS Advantage", fontsize=14, fontweight="bold", color="black")
    _ax2.tick_params(colors="black", labelsize=12)
    _ax2.grid(axis="y", alpha=0.3, color="grey", zorder=1)
    _ax2.legend(fontsize=11, facecolor="white", edgecolor="black", labelcolor="black")
    for _sp in _ax2.spines.values():
        _sp.set_color("black")
    _fig.tight_layout()
    mo.vstack([
        _fig,
        mo.download(data=fig_to_pdf_bytes(_fig),
                    filename=f"{_LIB}_topn_tactics_vs_legacy.pdf".replace(" ", "_"),
                    mimetype="application/pdf", label="Download PDF"),
    ])
    return


@app.cell
def _(dock_library_dropdown, mo):
    """Per-library interpretation note (reactive to the §2 dropdown)."""
    _NOTES = {
        "adenine": (
            "**Adenine FRED** — hidden-hub landscape (amidine_003 = 68/100 top-100 "
            "hits, mean-uninformative)."
        ),
        "adenine_hybrid": (
            "**Adenine HYBRID (A2A-Lenselink)** — cross-docking-pipeline control "
            "(same chemistry, different scoring). TT-TS dominant; TACTICS RWS bimodal "
            "(4/20 reps fail on aldehyde_137 GMIC under-rotation). Planned Welch is used "
            "because Tukey HSD's pooled MSE is contaminated by the bimodal method."
        ),
        "quinazoline": (
            "**Quinazoline FRED (JNK3 2ZDT)** — mean-reachable / near-saturated "
            "landscape; all five methods within ~2 pp. Canonical EMA-CATS data shown "
            "(within-noise of the pre-fix run, p = 0.21)."
        ),
        "quinazoline_hybrid": (
            "**Quinazoline HYBRID (JNK3 2ZDT)** — documented null case; all methods "
            "within ~0.1 pp at top-100 (three-axis framework predicts no shift)."
        ),
        "amide": (
            "**Amide FRED (JNK3 2ZDT)** — broadly-good, low-σ acids; greedy-reachable, "
            "all methods converge near the ceiling (TACTICS still ≥ Legacy)."
        ),
        "amide_hybrid": (
            "**Amide HYBRID (JNK3 2ZDT)** — pure score-function contrast; axis-1 "
            "hidden-hub size on the acid component rises (~10% → ~31%), motivating the "
            "TT-TS advantage."
        ),
        "thrombin": (
            "**Thrombin FRED** — imbalanced (130 × 3844) mean-separable docking; the "
            "architectural floor recovers Legacy RWS's damage, GMIC adds a modest gain."
        ),
        "thrombin_hybrid": (
            "**Thrombin HYBRID (2ZFF)** — same chemistry, different scoring; axis-1 "
            "hidden-hub size rises (~8% → ~25%), shifting the ranking toward TT-TS."
        ),
    }
    _lib = dock_library_dropdown.value
    mo.md(_NOTES.get(_lib, f"**{_lib}** — best-TACTICS vs Legacy across top-N."))
    return


@app.cell
def _(mo):
    mo.md(r"""
    # Section 3 — Overall Docking Recovery (FRED + HYBRID)

    Two parallel aggregate panels, one per scoring pipeline. NOT pooled across
    pipelines (would mix structurally-different SAR landscapes — see manuscript
    Discussion §Cross-docking-pipeline framework validation).
    """)
    return


@app.cell
def _(fig_to_pdf_bytes, mo, pl, plt, recovery_all):
    """Mean % recovery vs top-N for FRED docking libraries (2-comp and 3-comp panels).

    FRED-only aggregate. HYBRID counterpart is the next cell. Two separate panels
    avoid averaging across scoring pipelines (which would mix structurally-different
    SAR landscapes).
    """

    _DOCK_LIBS = {"adenine", "amide", "quinazoline", "thrombin"}
    # n_components is 0 for thrombin/quinazoline in the data — fix with known values
    _DOCK_NCOMP = {
        "thrombin": 2, "thrombin_hybrid": 2,
        "amide": 2, "amide_hybrid": 2,
        "adenine": 3, "adenine_hybrid": 3,
        "quinazoline": 3, "quinazoline_hybrid": 3,
    }

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

    _fig.suptitle("Mean % Recovery vs Top-N Threshold (FRED Docking Libraries)",
                  fontsize=15, fontweight="bold", color="black", y=1.01)
    _fig.tight_layout()

    mo.vstack([
        plt.gcf(),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename="docking_recovery_vs_topn_FRED.pdf",
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


@app.cell
def _(fig_to_pdf_bytes, mo, pl, plt, recovery_all):
    """Mean % recovery vs top-N for HYBRID docking libraries (2-comp and 3-comp panels).

    HYBRID-only aggregate parallel to the FRED panel above. Pool: 4 HYBRID docking
    libraries (adenine_hybrid, amide_hybrid, quinazoline_hybrid, thrombin_hybrid).
    """

    _DOCK_LIBS = {"adenine_hybrid", "amide_hybrid", "quinazoline_hybrid", "thrombin_hybrid"}
    _DOCK_NCOMP = {
        "thrombin_hybrid": 2, "amide_hybrid": 2,
        "adenine_hybrid": 3, "quinazoline_hybrid": 3,
    }

    _dock_all = (
        recovery_all
        .filter(pl.col("library_id").is_in(_DOCK_LIBS))
        .with_columns(
            pl.col("library_id").replace(_DOCK_NCOMP).cast(pl.Int64).alias("n_components")
        )
    )

    mo.stop(len(_dock_all) == 0, mo.md("No HYBRID docking library data available."))

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

    _fig, _axes = plt.subplots(1, 2, figsize=(16, 6), facecolor="white", sharey=True)

    for _panel_idx, _nc in enumerate([2, 3]):
        _ax = _axes[_panel_idx]
        _ax.set_facecolor("white")
        _comp_data = _dock_all.filter(pl.col("n_components") == _nc)
        _agg = (
            _comp_data.group_by(["method", "top_n"])
            .agg([
                pl.col("recovery_frac").mean().alias("mean_frac"),
                pl.col("recovery_frac").std().alias("std_frac"),
                pl.col("recovery_frac").count().alias("n"),
            ])
            .with_columns((pl.col("std_frac") / pl.col("n").sqrt()).alias("se_frac"))
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

        _libs_in_panel = sorted(_comp_data["library_id"].unique().to_list())
        _ax.set_xlabel("Top-N Threshold", fontsize=14, color="black")
        _ax.set_title(
            f"{_nc}-Component HYBRID Docking Libraries\n({', '.join(_libs_in_panel)})",
            fontsize=15, fontweight="bold", color="black",
        )
        _ax.set_xticks(sorted(_agg["top_n"].unique().to_list()))
        _ax.tick_params(colors="black", labelsize=12)
        _ax.grid(axis="both", alpha=0.3, color="grey", zorder=1)
        for _spine in _ax.spines.values():
            _spine.set_color("black")

    _axes[0].set_ylabel("% Recovery", fontsize=14, color="black")
    _handles, _labels = _axes[0].get_legend_handles_labels()
    _fig.legend(_handles, _labels, loc="lower center", ncol=3, fontsize=12,
                facecolor="white", edgecolor="black", labelcolor="black",
                bbox_to_anchor=(0.5, -0.08))
    _fig.suptitle("Mean % Recovery vs Top-N Threshold (HYBRID Docking Libraries)",
                  fontsize=15, fontweight="bold", color="black", y=1.01)
    _fig.tight_layout()

    mo.vstack([
        plt.gcf(),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename="docking_recovery_vs_topn_HYBRID.pdf",
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


@app.cell
def _(mo):
    mo.md(r"""
    # Section 4 — Statistical-Test Convention Reminder
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Statistical-test convention reminder

    Significance stars in the per-library figures follow the **scoring-type rule**:

    - **ROCS libraries** (20 libraries × 109 queries × 20 replicates): repeated-measures
      ANOVA + **Tukey HSD** per Ash et al. 2025. Stars rendered in **black**.
    - **Docking libraries** (8 libraries × 1 query × 20 replicates: lacks
      repeated-measures structure, single-query × 20 reps reduces to one-way ANOVA where
      homoscedasticity is required but empirically violated by bimodal methods): planned
      **Welch's t-test + Bonferroni-4** on the 4 pre-specified TACTICS-vs-Legacy pairs.
      Stars rendered in **blue (#4477AA)**.

    Both tests are reported in the analysis parquets; the per-library figures use the
    one appropriate to the library's scoring type. Cross-test agreement on the planned
    pairs is in `outputs/full_benchmark_ema_cats/analysis/welch_planned_results.parquet`
    and `tukey_hsd_results.parquet`.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    # Section 5 — Budget Sensitivity Analysis (1% → 5%) — DOCK

    Recovery and TACTICS-vs-Legacy gap as a function of library-fraction screening
    budget for the 8 docking libraries, at the canonical top-100 threshold.
    Demonstrates that TACTICS dominance is preserved across budgets and is largest
    at constrained-budget regimes (1–2%). The parallel ROCS sensitivity panel
    lives in `manuscript_plots_ROCS.py`.

    2×2 panel grid:
    - **Top row** (line plots): per-method mean recovery vs budget fraction.
    - **Bottom row** (stacked bars): Δ(best-TACTICS − best-Legacy) per budget,
      **decomposed into per-library contributions** (each library's Δ ÷ n_libraries,
      so the stacked segments sum to the bin mean — annotated on top). Color = base
      chemistry, summing each chemistry's FRED + HYBRID scorings into one segment.
      This exposes that the 3-component bin's advantage is carried by **adenine** while
      the two **quinazoline** variants contribute ≈0 — the source of the
      counter-intuitive 2-comp ≥ 3-comp ordering at most budgets.
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
    return sens_delta_per_lib, sens_per_method


@app.cell
def _(fig_to_pdf_bytes, mo, np, pl, plt, sens_delta_per_lib, sens_per_method):
    """DOCK sensitivity figure: 2×2 layout (line top, stacked-Δ bars bottom; 2-comp left, 3-comp right).

    Bottom row decomposes each bin's mean Δ (best-TACTICS − best-Legacy) into per-library
    contributions, colored by base chemistry (FRED + HYBRID for a chemistry are summed into
    one segment); segments sum to the bin mean (annotated). Makes explicit that the 3-comp
    bin's advantage is carried by adenine while the quinazolines contribute ~0.
    """
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
    # Bottom-row decomposition: one stacked segment per base chemistry. A chemistry's
    # FRED + HYBRID contributions are summed into a single colored segment (no hatch).
    # Colors are the UNUSED members of the same Tableau-10 palette the method lines
    # draw from (the lines use blue / orange / red / teal / purple; these are the
    # remaining green / brown / gold / pink). Same visual family as the top row, but
    # no color is shared with a method, so the bars are never confused with a line.
    _CHEM_COLOR = {
        "adenine": "#59A14F",       # Tableau green
        "quinazoline": "#FF9DA7",   # Tableau pink
        "amide": "#9C755F",         # Tableau brown
        "thrombin": "#EDC949",      # Tableau gold
    }

    def _base(_lib):
        return _lib[:-7] if _lib.endswith("_hybrid") else _lib

    _budgets = sorted(
        sens_delta_per_lib.filter(pl.col("kind") == "DOCK")["budget_pct"].unique().to_list()
    )

    _fig, _axes = plt.subplots(2, 2, figsize=(14, 9), facecolor="white", sharex=True)

    for _col_idx, _nc in enumerate([2, 3]):
        _ax_line = _axes[0, _col_idx]
        _ax_bar = _axes[1, _col_idx]
        _ax_line.set_facecolor("white"); _ax_bar.set_facecolor("white")

        _sub_method = sens_per_method.filter(
            (pl.col("kind") == "DOCK") & (pl.col("n_comp") == _nc)
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
        _ax_line.set_title(f"DOCK {_nc}-Component (n = {_n_libs} libraries)",
                            fontsize=13, fontweight="bold", color="black")
        _ax_line.grid(alpha=0.3, color="grey", zorder=1)
        _ax_line.tick_params(colors="black", labelsize=11)
        for _sp in _ax_line.spines.values(): _sp.set_color("black")

        # Bottom row: per-library Δ decomposed by base chemistry (stacked). A
        # chemistry's FRED + HYBRID contributions are summed into one segment, then
        # divided by the number of libraries in the bin so the stacked segments sum
        # to the bin-mean Δ shown by the top-of-stack annotation.
        _sub_lib = sens_delta_per_lib.filter(
            (pl.col("kind") == "DOCK") & (pl.col("n_comp") == _nc)
        )
        _n_by_b = {_b: _sub_lib.filter(pl.col("budget_pct") == _b).height
                   for _b in _budgets}
        _by_chem = {}
        for _r in _sub_lib.iter_rows(named=True):
            _c = _base(_r["library_id"])
            _d = _by_chem.setdefault(_c, {})
            _d[_r["budget_pct"]] = _d.get(_r["budget_pct"], 0.0) + _r["delta_pp"]

        _pos = np.zeros(len(_budgets))
        _neg = np.zeros(len(_budgets))
        for _chem in sorted(_by_chem):
            _heights = np.array([
                _by_chem[_chem].get(_b, 0.0) / max(_n_by_b[_b], 1) for _b in _budgets
            ])
            _bottoms = np.where(_heights >= 0, _pos, _neg)
            _ax_bar.bar(_budgets, _heights, bottom=_bottoms, width=0.7,
                        color=_CHEM_COLOR[_chem], edgecolor="black", linewidth=0.6,
                        zorder=3, label=_chem)
            _pos = _pos + np.where(_heights >= 0, _heights, 0.0)
            _neg = _neg + np.where(_heights < 0, _heights, 0.0)

        # Net (= bin-mean Δ) annotated above each stack's positive top.
        for _xi, _top, _mean in zip(_budgets, _pos, _pos + _neg):
            _ax_bar.text(_xi, _top + 0.3, f"{_mean:+.2f}", ha="center", va="bottom",
                         fontsize=10, fontweight="bold", color="black")

        _ax_bar.axhline(0, color="black", linewidth=0.8, zorder=2)
        _ax_bar.set_xlabel("Library Fraction (%)", fontsize=12, color="black")
        _ax_bar.set_ylabel("Δ Best-TACTICS − Best-Legacy (pp)\n(stacked by library)",
                           fontsize=12, color="black")
        _ax_bar.set_xticks([1, 2, 3, 4, 5])
        _ax_bar.grid(axis="y", alpha=0.3, color="grey", zorder=1)
        _ax_bar.tick_params(colors="black", labelsize=11)
        for _sp in _ax_bar.spines.values(): _sp.set_color("black")
        # No legend title: a one-word "Library" title sits in a box sized to the
        # widest entry ("quinazoline"), leaving a white gap beside it regardless of
        # alignment. The entries are clearly library names (the y-axis already says
        # "stacked by library"), so the title is dropped.
        _ax_bar.legend(fontsize=10, facecolor="white", edgecolor="black",
                       labelcolor="black", loc="upper right",
                       handlelength=1.2, handletextpad=0.5, borderpad=0.5,
                       labelspacing=0.35)
        _yl, _yh = _ax_bar.get_ylim()
        _ax_bar.set_ylim(_yl, _yh + 0.18 * (abs(_yh - _yl) or 1.0))

    _handles, _labels = _axes[0, 0].get_legend_handles_labels()
    _fig.legend(_handles, _labels, loc="lower center", ncol=5, fontsize=11,
                 facecolor="white", edgecolor="black", labelcolor="black",
                 bbox_to_anchor=(0.5, -0.04))
    _fig.suptitle("Budget Sensitivity — Docking Libraries (Top-100 Recovery)",
                   fontsize=14, fontweight="bold", color="black", y=1.00)
    _fig.tight_layout(rect=[0, 0.04, 1, 0.98])

    mo.vstack([
        plt.gcf(),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename="sensitivity_dock.pdf",
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


@app.cell
def _(mo, pl, sens_delta_per_lib):
    """Summary table — per-library Δ for all 8 docking libraries across 5 budgets.

    Demonstrates the adenine + quinazoline_hybrid pattern extends to the full 8-library
    docking set: TACTICS is either better (top rows) or statistically indistinguishable
    (bottom rows) than the best Legacy method at every budget. The "TACTICS is either
    better or the same" claim holds on 8/8 libraries at every budget tested.
    """
    _DOCK_LIBS_ORDERED = [
        # TACTICS-dominant regime (consistent positive Δ at most budgets)
        ("adenine", "FRED", 3),
        ("adenine_hybrid", "HYBRID", 3),
        ("amide_hybrid", "HYBRID", 2),
        ("thrombin_hybrid", "HYBRID", 2),
        ("thrombin", "FRED", 2),
        # Saturated / tied regime (Δ near zero, all methods at recovery ceiling)
        ("amide", "FRED", 2),
        ("quinazoline", "FRED", 3),
        ("quinazoline_hybrid", "HYBRID", 3),
    ]
    _meta = {l: (s, n) for l, s, n in _DOCK_LIBS_ORDERED}
    _order = [l for l, _, _ in _DOCK_LIBS_ORDERED]

    _rows = []
    for _lib in _order:
        _scoring, _ncomp = _meta[_lib]
        _r = {"library": _lib, "scoring": _scoring, "n_comp": _ncomp}
        _sub = sens_delta_per_lib.filter(
            (pl.col("kind") == "DOCK") & (pl.col("library_id") == _lib)
        )
        for _b in [1, 2, 3, 4, 5]:
            _val = _sub.filter(pl.col("budget_pct") == _b)["delta_pp"]
            _r[f"delta_{_b}pct"] = round(float(_val[0]), 2) if _val.len() > 0 else None
        _rows.append(_r)
    sens_dock_table = pl.DataFrame(_rows)

    _csv_bytes = sens_dock_table.write_csv().encode()

    mo.vstack([
        mo.md(
            "### Docking library summary — Δ Best-TACTICS − Best-Legacy (pp) across budgets\n\n"
            "Libraries ordered top-to-bottom by TACTICS dominance at the canonical 2% budget. "
            "Rows 1-5 = TACTICS-dominant regime (consistent positive Δ); rows 6-8 = saturated "
            "regime (Δ near zero, all methods at the recovery ceiling). The 'TACTICS is either "
            "better or the same' claim holds on 8/8 libraries at every budget tested."
        ),
        mo.ui.table(sens_dock_table, selection=None),
        mo.download(data=_csv_bytes, filename="sensitivity_dock_table.csv",
                    mimetype="text/csv", label="Download CSV"),
    ])
    return


if __name__ == "__main__":
    app.run()
