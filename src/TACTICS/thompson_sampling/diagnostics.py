"""
Pure analysis functions for CATS diagnostics.

All functions operate on Polars DataFrames and have no dependency on the
ThompsonSampler or strategy objects.  This makes them usable in notebooks,
scripts, and post-hoc analysis pipelines.
"""

from typing import Dict, Any, Optional
import numpy as np
import polars as pl


def compute_posterior_entropy(
    landscape_df: pl.DataFrame,
    mode: str = "maximize",
    criticality_metric: str = "ipr",
    n_adaptive_sharpening: bool = True,
) -> pl.DataFrame:
    """Compute per-component entropy metrics from a posterior landscape.

    Uses z-score softmax (same as CATS ``_calculate_criticality()``) to
    convert posterior means into a probability distribution, then measures
    concentration via IPR (default) or Shannon entropy.

    Args:
        landscape_df: DataFrame from ``sampler.get_posterior_landscape()``
            with columns ``component_idx``, ``reagent_name``, ``mean``,
            ``std``, ``n_samples``.
        mode: "maximize" or "minimize" — determines sign convention.
        criticality_metric: "ipr" or "shannon" — metric for concentration.
        n_adaptive_sharpening: Apply sqrt(log(N)) sharpening when using IPR.

    Returns:
        DataFrame with one row per component:
        ``component_idx``, ``n_active``, ``entropy``,
        ``normalized_entropy``, ``concentration``, ``snr``.
    """
    records = []

    for comp_idx in landscape_df["component_idx"].unique().sort().to_list():
        comp = landscape_df.filter(
            (pl.col("component_idx") == comp_idx) & (pl.col("n_samples") > 0)
        )
        n_active = len(comp)

        if n_active < 2:
            records.append({
                "component_idx": comp_idx,
                "n_active": n_active,
                "entropy": float("nan"),
                "normalized_entropy": float("nan"),
                "concentration": float("nan"),
                "snr": float("nan"),
            })
            continue

        means = comp["mean"].to_numpy().copy()
        stds = comp["std"].to_numpy()
        n_samples = comp["n_samples"].to_numpy()

        if mode == "minimize":
            means = -means

        mean_std = np.std(means)
        if mean_std < 1e-10:
            records.append({
                "component_idx": comp_idx,
                "n_active": n_active,
                "entropy": np.log(n_active),
                "normalized_entropy": 1.0,
                "concentration": 0.0,
                "snr": 0.0,
            })
            continue

        # SNR
        se_sq = stds ** 2 / np.maximum(n_samples, 1)
        noise_std = np.sqrt(np.mean(se_sq))
        snr = float(mean_std / max(noise_std, 1e-10))

        imbalance = 1.0 - np.exp(-max(snr - 1.0, 0.0))
        z = (means - np.mean(means)) / mean_std * imbalance

        # N-adaptive sharpening
        if criticality_metric == "ipr" and n_adaptive_sharpening and n_active > 2:
            sharpening = max(1.0, np.sqrt(np.log(n_active)))
            z *= sharpening

        exp_z = np.exp(z - z.max())  # numerical stability
        probs = exp_z / exp_z.sum()

        entropy = -np.sum(probs * np.log(probs + 1e-10))
        max_entropy = np.log(n_active)
        norm_ent = entropy / max_entropy if max_entropy > 1e-10 else 1.0

        if criticality_metric == "ipr":
            ipr = float(np.sum(probs ** 2))
            effective_n = 1.0 / ipr
            concentration = 1.0 - (effective_n / n_active)
        else:
            concentration = 1.0 - norm_ent

        records.append({
            "component_idx": comp_idx,
            "n_active": n_active,
            "entropy": float(entropy),
            "normalized_entropy": float(norm_ent),
            "concentration": float(concentration),
            "snr": snr,
        })

    return pl.DataFrame(records)


def compute_convergence_point(
    diagnostics_df: pl.DataFrame,
    threshold: float = 0.3,
    stability_window: int = 5,
) -> pl.DataFrame:
    """Find the cycle at which each component first reaches stable convergence.

    A component is considered converged at cycle *c* if criticality > threshold
    for all cycles from *c* through *c + stability_window - 1*.

    Args:
        diagnostics_df: Enhanced diagnostics DataFrame (must have
            ``current_cycle``, ``component_idx``, ``criticality`` columns).
        threshold: Criticality threshold for "structured" (default 0.3).
        stability_window: Number of consecutive cycles above threshold
            required for stability (default 5).

    Returns:
        DataFrame with ``component_idx``, ``cycle_first_above``,
        ``cycle_stable``.
    """
    # Support both enhanced (current_cycle) and legacy (cycle) column names
    cycle_col = "current_cycle" if "current_cycle" in diagnostics_df.columns else "cycle"

    records = []
    for comp_idx in diagnostics_df["component_idx"].unique().sort().to_list():
        comp = diagnostics_df.filter(
            pl.col("component_idx") == comp_idx
        ).sort(cycle_col)

        crits = comp["criticality"].to_numpy()
        cycles = comp[cycle_col].to_numpy()

        above = crits > threshold
        first_above = int(cycles[above][0]) if above.any() else None

        stable = None
        if above.any():
            for i in range(len(crits)):
                end = i + stability_window
                if end <= len(crits) and np.all(crits[i:end] > threshold):
                    stable = int(cycles[i])
                    break

        records.append({
            "component_idx": comp_idx,
            "cycle_first_above": first_above,
            "cycle_stable": stable,
        })

    return pl.DataFrame(records)


def compare_trajectory_vs_snapshot(
    diagnostics_df: pl.DataFrame,
    landscape_df: pl.DataFrame,
    mode: str = "maximize",
) -> pl.DataFrame:
    """Compare CATS trajectory information with a post-hoc snapshot.

    Shows what the trajectory reveals that the static snapshot doesn't:
    convergence timing, stability, and whether the final state was reached
    early or late.

    Args:
        diagnostics_df: Enhanced diagnostics DataFrame from
            ``sampler.get_diagnostics()``.
        landscape_df: Posterior landscape from
            ``sampler.get_posterior_landscape()``.
        mode: "maximize" or "minimize".

    Returns:
        DataFrame with per-component comparison: ``component_idx``,
        ``snapshot_concentration``, ``trajectory_final_criticality``,
        ``cycle_first_above``, ``cycle_stable``,
        ``trajectory_adds_info`` (bool).
    """
    snapshot = compute_posterior_entropy(landscape_df, mode=mode)
    convergence = compute_convergence_point(diagnostics_df)

    # Get final criticality from trajectory
    cycle_col = "current_cycle" if "current_cycle" in diagnostics_df.columns else "cycle"
    final_crits = {}
    for comp_idx in diagnostics_df["component_idx"].unique().sort().to_list():
        comp = diagnostics_df.filter(
            pl.col("component_idx") == comp_idx
        ).sort(cycle_col)
        final_crits[comp_idx] = float(comp["criticality"].to_list()[-1])

    records = []
    for comp_idx in snapshot["component_idx"].to_list():
        snap_conc = snapshot.filter(
            pl.col("component_idx") == comp_idx
        )["concentration"].to_list()[0]

        traj_crit = final_crits.get(comp_idx, float("nan"))

        conv_row = convergence.filter(pl.col("component_idx") == comp_idx)
        first_above = conv_row["cycle_first_above"].to_list()[0]
        stable = conv_row["cycle_stable"].to_list()[0]

        # Trajectory adds information if convergence timing is available
        adds_info = first_above is not None

        records.append({
            "component_idx": comp_idx,
            "snapshot_concentration": snap_conc,
            "trajectory_final_criticality": traj_crit,
            "cycle_first_above": first_above,
            "cycle_stable": stable,
            "trajectory_adds_info": adds_info,
        })

    return pl.DataFrame(records)


def compute_disagreement_convergence(
    diagnostics_df: pl.DataFrame,
    stability_window: int = 10,
    low_threshold: float = 0.3,
    high_threshold: float = 0.8,
) -> pl.DataFrame:
    """Analyze per-component disagreement EMA convergence from TT-TS diagnostics.

    Identifies when each component's disagreement rate stabilizes, indicating
    that TT-TS has resolved uncertainty about reagent rankings.

    Args:
        diagnostics_df: TT-TS diagnostics DataFrame (must have
            ``current_cycle``, ``component_idx``, ``disagreement_ema``).
        stability_window: Consecutive cycles within threshold band for
            stability (default 10).
        low_threshold: Disagreement below this means component is
            well-resolved (default 0.3).
        high_threshold: Disagreement above this means component is
            saturated / exploration is random (default 0.8).

    Returns:
        DataFrame with ``component_idx``, ``final_disagreement``,
        ``cycle_below_low``, ``cycle_stable_low``, ``mean_disagreement``,
        ``regime`` ("resolved" / "saturated" / "exploring").
    """
    cycle_col = "current_cycle" if "current_cycle" in diagnostics_df.columns else "cycle"

    records = []
    for comp_idx in diagnostics_df["component_idx"].unique().sort().to_list():
        comp = diagnostics_df.filter(
            pl.col("component_idx") == comp_idx
        ).sort(cycle_col)

        ema_vals = comp["disagreement_ema"].to_numpy()
        cycles = comp[cycle_col].to_numpy()
        valid = np.isfinite(ema_vals)

        if not valid.any():
            records.append({
                "component_idx": comp_idx,
                "final_disagreement": float("nan"),
                "cycle_below_low": None,
                "cycle_stable_low": None,
                "mean_disagreement": float("nan"),
                "regime": "insufficient_data",
            })
            continue

        final_val = float(ema_vals[valid][-1])
        mean_val = float(np.mean(ema_vals[valid]))

        # First cycle below low threshold
        below_low = ema_vals < low_threshold
        cycle_below = int(cycles[below_low][0]) if below_low.any() else None

        # Stable below low: stays below for stability_window
        cycle_stable = None
        if below_low.any():
            for i in range(len(ema_vals)):
                end = i + stability_window
                if end <= len(ema_vals) and np.all(ema_vals[i:end] < low_threshold):
                    cycle_stable = int(cycles[i])
                    break

        # Classify regime
        if final_val < low_threshold:
            regime = "resolved"
        elif final_val > high_threshold:
            regime = "saturated"
        else:
            regime = "exploring"

        records.append({
            "component_idx": comp_idx,
            "final_disagreement": final_val,
            "cycle_below_low": cycle_below,
            "cycle_stable_low": cycle_stable,
            "mean_disagreement": mean_val,
            "regime": regime,
        })

    return pl.DataFrame(records)


def compute_scale_adaptation(
    diagnostics_df: pl.DataFrame,
) -> pl.DataFrame:
    """Summarize per-component adaptive heated_scale evolution from TT-TS diagnostics.

    Args:
        diagnostics_df: TT-TS diagnostics DataFrame (must have
            ``current_cycle``, ``component_idx``, ``heated_scale``).

    Returns:
        DataFrame with ``component_idx``, ``initial_scale``, ``final_scale``,
        ``min_scale``, ``max_scale``, ``scale_range``, ``adaptation_direction``
        ("decayed" / "grew" / "stable").
    """
    cycle_col = "current_cycle" if "current_cycle" in diagnostics_df.columns else "cycle"

    records = []
    for comp_idx in diagnostics_df["component_idx"].unique().sort().to_list():
        comp = diagnostics_df.filter(
            pl.col("component_idx") == comp_idx
        ).sort(cycle_col)

        scales = comp["heated_scale"].to_numpy()
        valid = np.isfinite(scales)

        if not valid.any():
            records.append({
                "component_idx": comp_idx,
                "initial_scale": float("nan"),
                "final_scale": float("nan"),
                "min_scale": float("nan"),
                "max_scale": float("nan"),
                "scale_range": float("nan"),
                "adaptation_direction": "insufficient_data",
            })
            continue

        valid_scales = scales[valid]
        initial = float(valid_scales[0])
        final = float(valid_scales[-1])
        delta = final - initial

        if abs(delta) < 0.01:
            direction = "stable"
        elif delta < 0:
            direction = "decayed"
        else:
            direction = "grew"

        records.append({
            "component_idx": comp_idx,
            "initial_scale": initial,
            "final_scale": final,
            "min_scale": float(valid_scales.min()),
            "max_scale": float(valid_scales.max()),
            "scale_range": float(valid_scales.max() - valid_scales.min()),
            "adaptation_direction": direction,
        })

    return pl.DataFrame(records)


def format_sar_report(summary: Dict[str, Any]) -> str:
    """Format a ``get_sar_summary()`` dict into a publication-ready text block.

    Args:
        summary: Dict returned by ``ThompsonSampler.get_sar_summary()``.

    Returns:
        Multi-line string suitable for inclusion in a manuscript or notebook.
    """
    lines = []
    lines.append("SAR Landscape Assessment")
    lines.append("=" * 40)
    lines.append("")
    lines.append(f"Overall: {summary['landscape_type']}")
    lines.append(summary["summary_text"])
    lines.append("")

    lines.append("Per-Component Details")
    lines.append("-" * 40)

    for comp in summary["components"]:
        idx = comp["component_idx"]
        lines.append(f"\nComponent {idx}:")
        lines.append(f"  Assessment:    {comp['convergence_assessment']}")

        if np.isfinite(comp.get("concentration", float("nan"))):
            lines.append(f"  Concentration: {comp['concentration']:.3f}")
            lines.append(f"  SNR:           {comp['snr']:.2f}")
            lines.append(f"  Gini:          {comp['gini_coefficient']:.3f}")
            lines.append(f"  Top-1 share:   {comp['top1_share']:.1%}")
            lines.append(f"  Top-5 share:   {comp['top5_share']:.1%}")
            lines.append(f"  Dominant:      {comp['dominant_reagent']}")

    if "convergence_dynamics" in summary:
        lines.append("")
        lines.append("Convergence Dynamics (CATS trajectory)")
        lines.append("-" * 40)
        for dyn in summary["convergence_dynamics"]:
            idx = dyn["component_idx"]
            first = dyn["cycle_first_structured"]
            stable = dyn["cycle_stable"]
            slope = dyn["signal_trajectory_slope"]
            lines.append(f"\nComponent {idx}:")
            lines.append(f"  First structured:  cycle {first}")
            lines.append(f"  Stable from:       cycle {stable}")
            if np.isfinite(slope):
                lines.append(f"  Signal trend:      {slope:+.4f}/cycle")

    return "\n".join(lines)
