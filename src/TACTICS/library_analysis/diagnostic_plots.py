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
