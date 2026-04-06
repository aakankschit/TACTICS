"""
Matplotlib plotting functions for CATS diagnostics.

All functions accept Polars DataFrames and return ``matplotlib.figure.Figure``
objects so callers can save, show, or compose them freely.
"""

from typing import Optional
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure


def plot_criticality_trajectory(
    diagnostics_df: pl.DataFrame,
    figsize: tuple = (10, 4),
) -> Figure:
    """Plot per-component criticality over cycles.

    A horizontal band at criticality <= 0.3 is shaded as the "diffuse" zone.

    Args:
        diagnostics_df: Enhanced or legacy diagnostics DataFrame.
        figsize: Figure size.

    Returns:
        Matplotlib Figure.
    """
    cycle_col = "current_cycle" if "current_cycle" in diagnostics_df.columns else "cycle"

    fig, ax = plt.subplots(figsize=figsize)
    ax.axhspan(0, 0.3, alpha=0.1, color="blue", label="Diffuse zone")

    for comp_idx in sorted(diagnostics_df["component_idx"].unique().to_list()):
        comp = diagnostics_df.filter(
            pl.col("component_idx") == comp_idx
        ).sort(cycle_col)
        ax.plot(
            comp[cycle_col].to_numpy(),
            comp["criticality"].to_numpy(),
            label=f"Component {comp_idx}",
            linewidth=1.5,
        )

    ax.set_xlabel("Cycle")
    ax.set_ylabel("Criticality")
    ax.set_title("CATS Criticality Trajectory")
    ax.set_ylim(-0.05, 1.05)
    ax.legend()
    fig.tight_layout()
    return fig


def plot_snr_trajectory(
    diagnostics_df: pl.DataFrame,
    figsize: tuple = (10, 4),
) -> Figure:
    """Plot SNR evolution with a threshold line at SNR = 1.

    Requires enhanced diagnostics (``snr`` column).

    Args:
        diagnostics_df: Enhanced diagnostics DataFrame.
        figsize: Figure size.

    Returns:
        Matplotlib Figure.
    """
    cycle_col = "current_cycle" if "current_cycle" in diagnostics_df.columns else "cycle"

    fig, ax = plt.subplots(figsize=figsize)
    ax.axhline(1.0, color="red", linestyle="--", alpha=0.6, label="SNR = 1 (noise threshold)")

    for comp_idx in sorted(diagnostics_df["component_idx"].unique().to_list()):
        comp = diagnostics_df.filter(
            pl.col("component_idx") == comp_idx
        ).sort(cycle_col)
        snr = comp["snr"].to_numpy()
        cycles = comp[cycle_col].to_numpy()
        valid = np.isfinite(snr)
        ax.plot(
            cycles[valid],
            snr[valid],
            label=f"Component {comp_idx}",
            linewidth=1.5,
        )

    ax.set_xlabel("Cycle")
    ax.set_ylabel("Signal-to-Noise Ratio")
    ax.set_title("SNR Trajectory")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_temperature_decomposition(
    diagnostics_df: pl.DataFrame,
    component_idx: int,
    figsize: tuple = (10, 6),
) -> Figure:
    """Plot full temperature pipeline decomposition for one component.

    Shows base_temp, cats_multiplier, criticality_weight, effective_multiplier,
    and final_temperature on separate subplots.

    Requires enhanced diagnostics.

    Args:
        diagnostics_df: Enhanced diagnostics DataFrame.
        component_idx: Which component to plot.
        figsize: Figure size.

    Returns:
        Matplotlib Figure.
    """
    cycle_col = "current_cycle"
    comp = diagnostics_df.filter(
        pl.col("component_idx") == component_idx
    ).sort(cycle_col)
    cycles = comp[cycle_col].to_numpy()

    fields = [
        ("criticality", "Criticality"),
        ("criticality_weight", "Criticality Weight"),
        ("cats_multiplier", "CATS Multiplier"),
        ("effective_multiplier", "Effective Multiplier"),
        ("final_temperature", "Final Temperature"),
    ]

    fig, axes = plt.subplots(len(fields), 1, figsize=figsize, sharex=True)

    for ax, (col, title) in zip(axes, fields):
        vals = comp[col].to_numpy()
        ax.plot(cycles, vals, linewidth=1.5)
        ax.set_ylabel(title, fontsize=9)
        ax.grid(alpha=0.3)

    axes[-1].set_xlabel("Cycle")
    fig.suptitle(f"Temperature Decomposition — Component {component_idx}", fontsize=12)
    fig.tight_layout()
    return fig


def plot_disagreement_trajectory(
    diagnostics_df: pl.DataFrame,
    figsize: tuple = (10, 4),
) -> Figure:
    """Plot per-component TT-TS disagreement EMA over cycles.

    Shows how disagreement rate evolves as posteriors sharpen.  High
    disagreement means the two posterior samples frequently pick different
    reagents (exploration); low disagreement means they agree (exploitation).

    Horizontal bands mark the adaptive thresholds: above 0.8 (saturated,
    scale decays toward 1.0) and below 0.3 (resolved, scale grows).

    Requires TT-TS diagnostics (``disagreement_ema`` column).

    Args:
        diagnostics_df: TT-TS diagnostics DataFrame.
        figsize: Figure size.

    Returns:
        Matplotlib Figure.
    """
    cycle_col = "current_cycle" if "current_cycle" in diagnostics_df.columns else "cycle"

    fig, ax = plt.subplots(figsize=figsize)
    ax.axhspan(0.8, 1.0, alpha=0.1, color="red", label="Saturated zone (>0.8)")
    ax.axhspan(0, 0.3, alpha=0.1, color="green", label="Resolved zone (<0.3)")

    for comp_idx in sorted(diagnostics_df["component_idx"].unique().to_list()):
        comp = diagnostics_df.filter(
            pl.col("component_idx") == comp_idx
        ).sort(cycle_col)
        ema = comp["disagreement_ema"].to_numpy()
        cycles = comp[cycle_col].to_numpy()
        valid = np.isfinite(ema)
        ax.plot(
            cycles[valid],
            ema[valid],
            label=f"Component {comp_idx}",
            linewidth=1.5,
        )

    ax.set_xlabel("Cycle")
    ax.set_ylabel("Disagreement Rate (EMA)")
    ax.set_title("TT-TS Disagreement Trajectory")
    ax.set_ylim(-0.05, 1.05)
    ax.legend()
    fig.tight_layout()
    return fig


def plot_scale_adaptation(
    diagnostics_df: pl.DataFrame,
    figsize: tuple = (10, 4),
) -> Figure:
    """Plot per-component adaptive heated_scale evolution over cycles.

    Shows how the bidirectional EMA adapter adjusts posterior widening
    per component.  Components with saturated disagreement decay toward
    scale=1.0; under-explored components grow.

    Requires TT-TS diagnostics (``heated_scale`` column).

    Args:
        diagnostics_df: TT-TS diagnostics DataFrame.
        figsize: Figure size.

    Returns:
        Matplotlib Figure.
    """
    cycle_col = "current_cycle" if "current_cycle" in diagnostics_df.columns else "cycle"

    fig, ax = plt.subplots(figsize=figsize)
    ax.axhline(1.0, color="gray", linestyle="--", alpha=0.5, label="No widening (scale=1.0)")

    for comp_idx in sorted(diagnostics_df["component_idx"].unique().to_list()):
        comp = diagnostics_df.filter(
            pl.col("component_idx") == comp_idx
        ).sort(cycle_col)
        scales = comp["heated_scale"].to_numpy()
        cycles = comp[cycle_col].to_numpy()
        valid = np.isfinite(scales)
        ax.plot(
            cycles[valid],
            scales[valid],
            label=f"Component {comp_idx}",
            linewidth=1.5,
        )

    ax.set_xlabel("Cycle")
    ax.set_ylabel("Heated Scale")
    ax.set_title("TT-TS Adaptive Scale Evolution")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_ttts_mechanism_decomposition(
    diagnostics_df: pl.DataFrame,
    component_idx: int,
    figsize: tuple = (10, 8),
) -> Figure:
    """Plot full TT-TS mechanism decomposition for one component.

    Four-panel view: GMIC criticality, disagreement EMA, adaptive
    heated_scale, and effective scale (heated vs cooled).

    Requires TT-TS diagnostics.

    Args:
        diagnostics_df: TT-TS diagnostics DataFrame.
        component_idx: Which component to plot.
        figsize: Figure size.

    Returns:
        Matplotlib Figure.
    """
    cycle_col = "current_cycle"
    comp = diagnostics_df.filter(
        pl.col("component_idx") == component_idx
    ).sort(cycle_col)
    cycles = comp[cycle_col].to_numpy()

    fields = [
        ("criticality", "GMIC Criticality"),
        ("disagreement_ema", "Disagreement Rate (EMA)"),
        ("heated_scale", "Heated Scale"),
        ("effective_scale", "Effective Scale (heated or cooled)"),
    ]

    fig, axes = plt.subplots(len(fields), 1, figsize=figsize, sharex=True)

    for ax, (col, title) in zip(axes, fields):
        vals = comp[col].to_numpy()
        valid = np.isfinite(vals)
        ax.plot(cycles[valid], vals[valid], linewidth=1.5)
        ax.set_ylabel(title, fontsize=9)
        ax.grid(alpha=0.3)

    # Add reference lines
    axes[0].axhline(0.3, color="red", linestyle="--", alpha=0.4)
    axes[1].axhline(0.8, color="red", linestyle="--", alpha=0.4, label="Saturated")
    axes[1].axhline(0.3, color="green", linestyle="--", alpha=0.4, label="Resolved")
    axes[1].legend(fontsize=7)
    axes[2].axhline(1.0, color="gray", linestyle="--", alpha=0.4)

    axes[-1].set_xlabel("Cycle")
    fig.suptitle(
        f"TT-TS Mechanism Decomposition — Component {component_idx}",
        fontsize=12,
    )
    fig.tight_layout()
    return fig


def plot_posterior_heatmap(
    landscape_df: pl.DataFrame,
    component_idx: int,
    top_n: int = 20,
    mode: str = "maximize",
    figsize: tuple = (10, 6),
) -> Figure:
    """Plot reagent posterior mean heatmap for one component.

    Shows the top-N reagents by posterior mean as a horizontal bar chart.

    Args:
        landscape_df: DataFrame from ``sampler.get_posterior_landscape()``.
        component_idx: Which component.
        top_n: Number of top reagents to show.
        mode: "maximize" or "minimize".
        figsize: Figure size.

    Returns:
        Matplotlib Figure.
    """
    comp = landscape_df.filter(
        (pl.col("component_idx") == component_idx) & (pl.col("n_samples") > 0)
    )

    ascending = mode == "minimize"
    comp = comp.sort("mean", descending=not ascending).head(top_n)

    names = comp["reagent_name"].to_list()[::-1]
    means = comp["mean"].to_numpy()[::-1]
    stds = comp["std"].to_numpy()[::-1]

    fig, ax = plt.subplots(figsize=figsize)
    y_pos = np.arange(len(names))
    ax.barh(y_pos, means, xerr=stds, align="center", alpha=0.8, capsize=3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=7)
    ax.set_xlabel("Posterior Mean")
    ax.set_title(f"Top-{top_n} Reagents — Component {component_idx}")
    fig.tight_layout()
    return fig


def plot_convergence_comparison(
    diagnostics_df: pl.DataFrame,
    landscape_df: pl.DataFrame,
    mode: str = "maximize",
    figsize: tuple = (12, 5),
) -> Figure:
    """Side-by-side: CATS trajectory vs post-hoc snapshot.

    Left panel: criticality trajectory over cycles (dynamic view).
    Right panel: concentration from posterior snapshot (static view).

    This is the key manuscript figure showing what CATS adds beyond
    a simple post-hoc reagent frequency analysis.

    Args:
        diagnostics_df: Enhanced diagnostics DataFrame.
        landscape_df: Posterior landscape DataFrame.
        mode: "maximize" or "minimize".
        figsize: Figure size.

    Returns:
        Matplotlib Figure.
    """
    from ..thompson_sampling.diagnostics import compute_posterior_entropy

    snapshot = compute_posterior_entropy(landscape_df, mode=mode)
    cycle_col = "current_cycle" if "current_cycle" in diagnostics_df.columns else "cycle"

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Left: trajectory
    ax1.axhspan(0, 0.3, alpha=0.1, color="blue")
    for comp_idx in sorted(diagnostics_df["component_idx"].unique().to_list()):
        comp = diagnostics_df.filter(
            pl.col("component_idx") == comp_idx
        ).sort(cycle_col)
        ax1.plot(
            comp[cycle_col].to_numpy(),
            comp["criticality"].to_numpy(),
            label=f"Comp {comp_idx}",
            linewidth=1.5,
        )
    ax1.set_xlabel("Cycle")
    ax1.set_ylabel("Criticality")
    ax1.set_title("CATS Trajectory (dynamic)")
    ax1.set_ylim(-0.05, 1.05)
    ax1.legend()

    # Right: snapshot
    comp_indices = snapshot["component_idx"].to_list()
    concentrations = snapshot["concentration"].to_numpy()
    bars = ax2.bar(
        [str(c) for c in comp_indices],
        concentrations,
        color=["#2196F3" if c > 0.3 else "#90CAF9" for c in concentrations],
    )
    ax2.axhline(0.3, color="red", linestyle="--", alpha=0.6, label="Structured threshold")
    ax2.set_xlabel("Component")
    ax2.set_ylabel("Concentration (1 − norm. entropy)")
    ax2.set_title("Post-hoc Snapshot (static)")
    ax2.set_ylim(-0.05, 1.05)
    ax2.legend()

    fig.suptitle("Convergence Comparison: Trajectory vs Snapshot", fontsize=13)
    fig.tight_layout()
    return fig


# ── Manuscript-quality diagnostic plots ──────────────────────────────


def _component_palette(n: int) -> list[str]:
    """Return *n* distinguishable colours from the tab10 colourmap."""
    cmap = plt.cm.get_cmap("tab10")
    return [cmap(i) for i in range(n)]


def _rolling_mean(arr: np.ndarray, window: int) -> np.ndarray:
    """Centred rolling mean, NaN-padded at edges."""
    if window < 2 or len(arr) < window:
        return arr.copy()
    kernel = np.ones(window) / window
    padded = np.pad(arr, (window // 2, window - 1 - window // 2), mode="edge")
    return np.convolve(padded, kernel, mode="valid")[: len(arr)]


def _ci95(std: np.ndarray, n: int) -> np.ndarray:
    """Convert std to 95 % confidence interval half-width."""
    if n <= 1:
        return std
    from scipy.stats import t as t_dist
    return t_dist.ppf(0.975, df=n - 1) * std / np.sqrt(n)


def plot_rws_diagnostic(
    diagnostics_df: pl.DataFrame,
    title: str = "",
    figsize: tuple = (12, 5),
) -> Figure:
    """RWS diagnostic: GMIC + CATS multiplier + divergence gate.

    Dual y-axis plot averaged across replicates with 95 % confidence
    interval bands.  Vertical background shading indicates batches
    where the divergence gate blocks GMIC modulation (diversity mode).

    Left y-axis: GMIC per component (solid lines with CI bands).
    Right y-axis: CATS multiplier per component (dashed lines).
        Shows the GMIC-driven temperature modulation directly, without
        the noisy heated/cooled oscillation that final_temperature has.
        Values >1.0 = amplified exploration, <1.0 = amplified
        exploitation, =1.0 = base temperature (no modulation).
    Background: light pink shading with red boundary lines when any
        component is in diversity mode in >30% of replicates.

    Args:
        diagnostics_df: RWS diagnostics for a single query, all
            replicates.  Works with any number of components.
        title: Plot title (e.g. "rxn206 — query_008").
        figsize: Figure size.

    Returns:
        Matplotlib Figure.
    """
    cycle_col = "current_cycle"
    comps = sorted(diagnostics_df["component_idx"].unique().to_list())
    n_reps = diagnostics_df["replicate"].n_unique()
    colors = _component_palette(len(comps))
    cycles = sorted(diagnostics_df[cycle_col].unique().to_list())

    fig, ax_gmic = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor("white")
    ax_gmic.set_facecolor("white")
    ax_mult = ax_gmic.twinx()

    # ── Background shading for diversity mode ──────────────────────
    div_fracs = []
    for cycle in cycles:
        cd = diagnostics_df.filter(pl.col(cycle_col) == cycle)
        n_div = (cd["cats_mode"] == "diversity").sum()
        div_fracs.append(n_div / len(cd) if len(cd) > 0 else 0)

    _gate_spans = []
    in_span = False
    span_start = cycles[0]
    for cycle, frac in zip(cycles, div_fracs):
        if frac > 0.3 and not in_span:
            span_start = cycle
            in_span = True
        elif frac <= 0.3 and in_span:
            _gate_spans.append((span_start - 0.5, cycle - 0.5))
            in_span = False
    if in_span:
        _gate_spans.append((span_start - 0.5, cycles[-1] + 0.5))

    for x0, x1 in _gate_spans:
        ax_gmic.axvspan(x0, x1, alpha=0.15, color="#ffe0e0", zorder=0)
        ax_gmic.axvline(x0, color="red", linewidth=0.8, alpha=0.6, zorder=0)
        ax_gmic.axvline(x1, color="red", linewidth=0.8, alpha=0.6, zorder=0)

    # ── Per-component traces ──────────────────────────────────────
    for i, comp in enumerate(comps):
        color = colors[i]
        comp_data = diagnostics_df.filter(pl.col("component_idx") == comp)

        agg = (
            comp_data.group_by(cycle_col)
            .agg([
                pl.col("gmic").mean().alias("gmic_mean"),
                pl.col("gmic").std().alias("gmic_std"),
                pl.col("cats_multiplier").mean().alias("mult_mean"),
                pl.col("cats_multiplier").std().alias("mult_std"),
            ])
            .sort(cycle_col)
        )
        c = agg[cycle_col].to_numpy()
        gmic_mean = agg["gmic_mean"].to_numpy()
        gmic_ci = _ci95(agg["gmic_std"].fill_null(0).to_numpy(), n_reps)
        mult_mean = agg["mult_mean"].to_numpy()
        mult_ci = _ci95(agg["mult_std"].fill_null(0).to_numpy(), n_reps)

        # GMIC on left axis (solid line + 95% CI)
        ax_gmic.plot(c, gmic_mean, color=color, linewidth=1.8, zorder=3)
        ax_gmic.fill_between(c, gmic_mean - gmic_ci, gmic_mean + gmic_ci,
                             color=color, alpha=0.20, zorder=2)

        # CATS multiplier on right axis (dashed line + 95% CI)
        ax_mult.plot(c, mult_mean, color=color, linewidth=1.4,
                     linestyle="--", alpha=0.85, zorder=3)
        ax_mult.fill_between(c, mult_mean - mult_ci, mult_mean + mult_ci,
                             color=color, alpha=0.10, zorder=1)

    ax_gmic.set_xlabel("Batch")
    ax_gmic.set_ylabel("GMIC (solid)")
    ax_mult.set_ylabel("CATS Multiplier (dashed)")
    ax_gmic.set_title(title or "RWS Diagnostic")

    # Reference line: multiplier = 1.0 (no modulation)
    ax_mult.axhline(1.0, color="gray", linestyle=":", alpha=0.5, linewidth=0.8)

    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    _handles = [
        Line2D([0], [0], color=colors[i], linewidth=1.8, label=f"Component {c}")
        for i, c in enumerate(comps)
    ]
    _handles.append(Patch(
        facecolor="#ffe0e0", edgecolor="red", linewidth=0.8, alpha=0.7,
        label="Divergence gate active",
    ))
    ax_gmic.legend(handles=_handles, loc="upper left", fontsize=7)
    fig.tight_layout()
    return fig


def plot_ttts_diagnostic(
    diagnostics_df: pl.DataFrame,
    title: str = "",
    figsize: tuple = (12, 5),
) -> Figure:
    """TT-TS diagnostic: disagreement EMA + heated_scale.

    Dual y-axis plot averaged across replicates with 95 % confidence
    interval bands.

    Left y-axis: disagreement EMA per component (solid lines).
    Right y-axis: heated_scale per component (dashed lines).
    Horizontal bands mark the adaptive thresholds (0.3 and 0.8).

    Works with any number of components.

    Args:
        diagnostics_df: TT-TS diagnostics for a single query, all
            replicates.
        title: Plot title.
        figsize: Figure size.

    Returns:
        Matplotlib Figure.
    """
    cycle_col = "current_cycle"
    comps = sorted(diagnostics_df["component_idx"].unique().to_list())
    n_reps = diagnostics_df["replicate"].n_unique()
    colors = _component_palette(len(comps))

    fig, ax_ema = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor("white")
    ax_ema.set_facecolor("white")
    ax_scale = ax_ema.twinx()

    # Threshold bands
    ax_ema.axhspan(0.8, 1.0, alpha=0.06, color="red")
    ax_ema.axhspan(0.0, 0.3, alpha=0.06, color="green")

    for i, comp in enumerate(comps):
        color = colors[i]
        comp_data = diagnostics_df.filter(pl.col("component_idx") == comp)

        agg = (
            comp_data.group_by(cycle_col)
            .agg([
                pl.col("disagreement_ema").mean().alias("ema_mean"),
                pl.col("disagreement_ema").std().alias("ema_std"),
                pl.col("heated_scale").mean().alias("hs_mean"),
                pl.col("heated_scale").std().alias("hs_std"),
            ])
            .sort(cycle_col)
        )
        c = agg[cycle_col].to_numpy()
        ema_mean = agg["ema_mean"].to_numpy()
        ema_ci = _ci95(agg["ema_std"].fill_null(0).to_numpy(), n_reps)
        hs_mean = agg["hs_mean"].to_numpy()
        hs_ci = _ci95(agg["hs_std"].fill_null(0).to_numpy(), n_reps)

        valid_ema = np.isfinite(ema_mean)
        valid_hs = np.isfinite(hs_mean)

        # Disagreement EMA (solid + 95% CI)
        ax_ema.plot(c[valid_ema], ema_mean[valid_ema], color=color,
                    linewidth=1.5, label=f"Component {comp}")
        ax_ema.fill_between(
            c[valid_ema],
            (ema_mean - ema_ci)[valid_ema],
            (ema_mean + ema_ci)[valid_ema],
            color=color, alpha=0.15,
        )

        # Heated scale (dashed + 95% CI)
        ax_scale.plot(c[valid_hs], hs_mean[valid_hs], color=color,
                      linewidth=1.2, linestyle="--", alpha=0.85)
        ax_scale.fill_between(
            c[valid_hs],
            (hs_mean - hs_ci)[valid_hs],
            (hs_mean + hs_ci)[valid_hs],
            color=color, alpha=0.07,
        )

    ax_ema.set_xlabel("Batch")
    ax_ema.set_ylabel("Disagreement EMA (solid)")
    ax_scale.set_ylabel("Heated Scale (dashed)")
    ax_ema.set_ylim(-0.05, 1.05)
    ax_ema.set_title(title or "TT-TS Diagnostic")

    ax_scale.axhline(1.0, color="gray", linestyle=":", alpha=0.4, linewidth=0.8)
    ax_ema.legend(loc="upper right", fontsize=8)

    ax_ema.annotate("saturated (>0.8)", xy=(0.02, 0.88),
                    xycoords="axes fraction", fontsize=7, color="red", alpha=0.5)
    ax_ema.annotate("resolved (<0.3)", xy=(0.02, 0.22),
                    xycoords="axes fraction", fontsize=7, color="green", alpha=0.5)

    fig.tight_layout()
    return fig


def plot_min_gmic_vs_advantage(
    query_data: list[dict],
    figsize: tuple = (7, 5),
) -> Figure:
    """Scatter: min-component GMIC vs TACTICS advantage per query.

    Each point is one query.  Marker shape encodes scoring type (ROCS
    vs docking), colour encodes component count.  Background zones
    shade the hard / moderate / easy GMIC regimes.

    Args:
        query_data: List of dicts with keys *library*, *query*,
            *n_components*, *scoring_type*, *min_gmic*, *advantage*.
        figsize: Figure size.

    Returns:
        Matplotlib Figure.
    """
    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    comp_colors = {2: "#1f77b4", 3: "#ff7f0e", 4: "#2ca02c"}
    for d in query_data:
        marker = "o" if d["scoring_type"] == "ROCS" else "D"
        color = comp_colors.get(d["n_components"], "#7f7f7f")
        edge = "black" if d["n_components"] >= 3 else "none"
        ax.scatter(
            d["min_gmic"], d["advantage"],
            marker=marker, c=color, edgecolors=edge,
            s=80, linewidths=0.8, zorder=3,
        )
        ax.annotate(
            f"{d['library']}\n{d['query']}",
            (d["min_gmic"], d["advantage"]),
            fontsize=6, ha="center", va="bottom",
            xytext=(0, 6), textcoords="offset points",
        )

    # Threshold zones
    x_max = max(d["min_gmic"] for d in query_data) + 0.2
    ax.axvspan(0, 1.0, alpha=0.05, color="red", label="min_GMIC < 1.0 (hard)")
    ax.axvspan(1.0, 1.5, alpha=0.05, color="orange", label="1.0\u20131.5 (moderate)")
    ax.axvspan(1.5, x_max, alpha=0.05, color="green", label="min_GMIC > 1.5 (easy)")

    ax.set_xlabel("Minimum Component GMIC (late search)")
    ax.set_ylabel("TACTICS Advantage (recovery pts)")
    ax.set_title("SAR Difficulty Predicts TACTICS Advantage")
    ax.legend(fontsize=7, loc="upper right")

    # Marker legend
    unique_n = sorted(set(d["n_components"] for d in query_data))
    unique_scoring = sorted(set(d["scoring_type"] for d in query_data))
    elems = []
    for n in unique_n:
        for sc in unique_scoring:
            marker = "o" if sc == "ROCS" else "D"
            c = comp_colors.get(n, "#7f7f7f")
            edge = "black" if n >= 3 else "none"
            elems.append(Line2D(
                [0], [0], marker=marker, color="w", markerfacecolor=c,
                markeredgecolor=edge, markersize=8,
                label=f"{n}-comp ({sc})",
            ))
    ax2 = ax.twinx()
    ax2.set_yticks([])
    ax2.legend(handles=elems, loc="center right", fontsize=7,
               title="Library type", title_fontsize=8)

    fig.tight_layout()
    return fig
