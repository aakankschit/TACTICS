"""
Matplotlib plotting functions for CATS diagnostics.

All functions accept Polars DataFrames and return ``matplotlib.figure.Figure``
objects so callers can save, show, or compose them freely.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Optional
import numpy as np
import polars as pl

if TYPE_CHECKING:
    from matplotlib.figure import Figure

_VIZ_HINT = "matplotlib is required for plotting: pip install 'chem-tactics[viz]'"


def _pyplot():
    """Import matplotlib.pyplot on first use (optional ``viz`` extra)."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError(_VIZ_HINT) from exc
    return plt


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

    fig, ax = _pyplot().subplots(figsize=figsize)
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

    fig, ax = _pyplot().subplots(figsize=figsize)
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

    fig, axes = _pyplot().subplots(len(fields), 1, figsize=figsize, sharex=True)

    for ax, (col, title) in zip(axes, fields):
        vals = comp[col].to_numpy()
        ax.plot(cycles, vals, linewidth=1.5)
        ax.set_ylabel(title, fontsize=13)
        ax.grid(alpha=0.3)

    axes[-1].set_xlabel("Cycle")
    fig.suptitle(f"Temperature Decomposition — Component {component_idx}", fontsize=15)
    fig.tight_layout()
    return fig


# ── Manuscript-quality diagnostic plots ──────────────────────────────


def _component_palette(n: int) -> list[str]:
    """Return *n* distinguishable colours from the tab10 colourmap."""
    _pyplot()  # ensures matplotlib is importable, with the [viz] hint on failure
    from matplotlib import colormaps

    cmap = colormaps["tab10"]
    return [cmap(i) for i in range(n)]


def _smooth_trim(arr: np.ndarray, window: int) -> np.ndarray:
    """Centred rolling mean with an edge-shrinking window (no padding).

    Averages only the points that actually exist within the window, shrinking
    the window at the array ends rather than edge-padding.  This avoids the
    boundary bias that a repeated heated/cooled extreme would otherwise inject
    into the last few cycles of an oscillating per-cycle signal.  NaNs are
    ignored.
    """
    n = len(arr)
    if window < 2 or n == 0:
        return arr.copy()
    half = window // 2
    out = np.full(n, np.nan)
    for i in range(n):
        seg = arr[max(0, i - half): min(n, i + half + 1)]
        seg = seg[~np.isnan(seg)]
        if seg.size:
            out[i] = seg.mean()
    return out


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

    fig, ax_gmic = _pyplot().subplots(figsize=figsize)
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
    ax_gmic.legend(handles=_handles, loc="upper left", fontsize=12)
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

    fig, ax_ema = _pyplot().subplots(figsize=figsize)
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
    ax_ema.legend(loc="upper right", fontsize=13)

    ax_ema.annotate("saturated (>0.8)", xy=(0.02, 0.88),
                    xycoords="axes fraction", fontsize=12, color="red", alpha=0.5)
    ax_ema.annotate("resolved (<0.3)", xy=(0.02, 0.22),
                    xycoords="axes fraction", fontsize=12, color="green", alpha=0.5)

    fig.tight_layout()
    return fig


# ── Reagent-usage action panel ───────────────────────────────────────


def _parse_name(
    name: str,
    reagent_id_sets: list[set[str]],
) -> Optional[tuple[str, ...]]:
    """Parse a product ``Name`` into per-component reagent IDs.

    Reagent IDs may themselves contain underscores (e.g. adenine's
    ``amidine_034`` or ``isocyanide_db_511``), so a simple split by
    ``_`` does not recover the components.  Instead we greedily consume
    the leftmost tokens until they form a known ID for the next
    component.

    Returns ``None`` if the name cannot be parsed.
    """
    tokens = name.split("_")
    n_comps = len(reagent_id_sets)
    out: list[str] = []
    i = 0
    for c in range(n_comps):
        remaining_comps = n_comps - c - 1
        found = False
        # Try shortest prefix first (greedy).  Leave at least one token
        # per remaining component.
        max_take = len(tokens) - i - remaining_comps
        for take in range(1, max_take + 1):
            candidate = "_".join(tokens[i : i + take])
            if candidate in reagent_id_sets[c]:
                out.append(candidate)
                i += take
                found = True
                break
        if not found:
            return None
    if i != len(tokens):
        return None
    return tuple(out)


def _load_reagent_ids(reagent_files: list[Path]) -> list[set[str]]:
    """Load reagent IDs from .smi files — second column of each row."""
    out: list[set[str]] = []
    for rf in reagent_files:
        ids: set[str] = set()
        with open(rf, encoding="latin-1") as fh:
            for line in fh:
                parts = line.strip().split()
                if len(parts) >= 2:
                    ids.add(parts[1])
        out.append(ids)
    return out


def _infer_component_name(reagent_ids: set[str], idx: int) -> str:
    """Infer a human-readable name from the common prefix of reagent IDs.

    Examples:
        {"amidine_001", "amidine_034"}      → "Amidine"
        {"isocyanide_db_511", "isocyanide_db_042"} → "Isocyanide"
        {"aldehyde_305", "aldehyde_112"}    → "Aldehyde"
        {"3043690", "1526634"}              → "Component {idx}"

    Returns the first common non-numeric token, capitalised.  Falls
    back to ``"Component {idx}"`` when no chemistry-suggesting prefix
    can be found (e.g. purely numeric IDs).
    """
    if not reagent_ids:
        return f"Component {idx}"
    tokens_per_id = [rid.split("_") for rid in reagent_ids]
    first_tokens = {t[0] for t in tokens_per_id if t}
    if len(first_tokens) == 1:
        token = next(iter(first_tokens))
        if not token.isdigit() and token.isalpha():
            return token.capitalize()
    return f"Component {idx}"


def _per_batch_reagent_usage(
    names: list[str],
    reagent_id_sets: list[set[str]],
    batch_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (unique_count, new_count) arrays of shape (n_batches, n_components).

    ``unique_count[b, c]`` = number of distinct reagent IDs for component
    *c* seen in batch *b*.  ``new_count[b, c]`` = subset of those that had
    never appeared in batches 0..b-1.
    """
    n_components = len(reagent_id_sets)
    n_batches = len(names) // batch_size
    unique_count = np.zeros((n_batches, n_components), dtype=int)
    new_count = np.zeros((n_batches, n_components), dtype=int)
    seen: list[set[str]] = [set() for _ in range(n_components)]

    for b in range(n_batches):
        batch_names = names[b * batch_size : (b + 1) * batch_size]
        per_comp_ids: list[set[str]] = [set() for _ in range(n_components)]
        for name in batch_names:
            parsed = _parse_name(name, reagent_id_sets)
            if parsed is None:
                continue
            for c in range(n_components):
                per_comp_ids[c].add(parsed[c])
        for c in range(n_components):
            unique_count[b, c] = len(per_comp_ids[c])
            new_count[b, c] = len(per_comp_ids[c] - seen[c])
            seen[c] |= per_comp_ids[c]
    return unique_count, new_count


def _find_activity_cutoff(
    new_count: np.ndarray,
    heated_by_cycle: dict[int, int],
    window: int = 20,
    peak_fraction: float = 0.15,
    buffer: int = 10,
) -> int:
    """Find the last batch where the mechanism is still 'active'.

    Active = rolling-window new-reagent discoveries are at least
    ``peak_fraction`` of the peak rolling total over the full search.
    Once discoveries drop permanently below the threshold the search
    has entered an exploitation regime where later batches look
    essentially identical to each other.

    The peak is taken over batches from index ``window`` onward so the
    initial warmup-adjacent burst of novel reagents does not inflate
    the threshold.  ``heated_by_cycle`` is reserved for future tuning
    and is currently unused.
    """
    del heated_by_cycle  # reserved for future use
    n_batches = new_count.shape[0]
    if n_batches == 0:
        return 0
    total_new_per_batch = new_count.sum(axis=1)
    rolling_new = np.zeros(n_batches)
    for b in range(n_batches):
        lo = max(0, b - window + 1)
        rolling_new[b] = total_new_per_batch[lo : b + 1].sum()
    # Peak over post-startup portion so initial all-new burst does not
    # dominate the threshold.
    peak = float(rolling_new[window:].max() if n_batches > window else rolling_new.max())
    threshold_new = peak * peak_fraction if peak > 0 else 0.0
    if threshold_new == 0:
        return n_batches

    last_active = -1
    for b in range(n_batches):
        if rolling_new[b] >= threshold_new:
            last_active = b
    if last_active < 0:
        return n_batches
    return min(n_batches, last_active + 1 + buffer)


def plot_reagent_usage_action_panel(
    search_df: pl.DataFrame,
    diag_df: pl.DataFrame,
    reagent_files: list[Path],
    batch_size: int = 100,
    replicate: Optional[int] = None,
    aggregate_window: int = 1,
    normalize: str = "library",
    auto_trim: bool = True,
    max_batches: Optional[int] = None,
    component_names: Optional[list[str]] = None,
    title: str = "",
    figsize: tuple = (14, 5.5),
    heated_label: str = "Heated",
) -> Figure:
    """Two-panel reagent usage diagnostic: heating timeline + new-reagent lines.

    **Top panel (action strip).**  A thin categorical timeline showing
    which component was marked ``is_heated`` at each batch — one
    coloured block per batch, coloured by the component's palette
    entry.  Reads left-to-right like a sequence diagram.  Batches with
    no heated record appear white.  The strip's y-label is set by
    ``heated_label`` so the figure can use method-appropriate
    terminology (e.g. ``"Heated"`` for RWS, ``"Explore"`` for TT-TS).

    **Bottom panel (lines).**  Per-component lines showing the
    fraction of that component's reagents first-seen in this batch
    (``new_count / n_reagents_per_component``) — always
    library-normalized so the scree-plot decay reads the same across
    components.  Uses the same component palette as the heating strip,
    so the reader can match each line to the heated-block colours in
    the strip above.  Interpretation: the decay rate summarizes when
    the algorithm has exhausted the informative reagents per
    component.

    Note: the ``normalize`` parameter is retained for API
    compatibility but no longer affects the plot (the bar panel it
    controlled has been replaced by the heating strip).  The value is
    still validated so stale callers passing unknown strings still
    raise.

    Args:
        search_df: Search-phase rows for a single trial (library, query,
            method, replicate), in evaluation order.  Must have ``Name``
            (underscore-joined reagent IDs) and ``phase`` columns.  If
            ``replicate`` is None the DataFrame must already be filtered
            to a single replicate.
        diag_df: Per-cycle diagnostics for the same trial.  Must have
            ``component_idx``, ``current_cycle``, ``is_heated``.
        reagent_files: Ordered list of ``.smi`` file paths, one per
            component.  Used to count total reagents (for ``library``
            normalization) and to correctly parse product names where
            reagent IDs themselves may contain underscores
            (e.g. adenine's ``isocyanide_db_511``).
        batch_size: Evaluations per batch (default 100).
        replicate: If provided, filter both DataFrames to this replicate.
        aggregate_window: If >1, bin every ``aggregate_window`` batches
            into a single bar.  Useful for long searches (e.g., 1000+
            batches on docking libraries).  Diversity is averaged
            across batches in the bin; *new* is recomputed against all
            prior bins; heated component is the majority within the bin.
        auto_trim: If True (default), automatically crop to the active
            region of the search — batches where the mechanism is still
            switching heated components or discovering new reagents.
            Trailing exploitation-only batches (same component heated
            for many batches with no new reagents) are dropped and a
            note is annotated on the figure.
        max_batches: Hard upper limit on the number of batches shown.
            Takes precedence over ``auto_trim`` if smaller.
        component_names: Optional list of display names (one per
            component) to use in the legend instead of "Component N".
            If ``None``, names are inferred from the common prefix of
            each component's reagent IDs (e.g., adenine automatically
            becomes ["Amidine", "Isocyanide", "Aldehyde"]).  Components
            with purely numeric reagent IDs fall back to "Component N".
        normalize: ``"library"`` (default) normalizes by total reagents
            in each component — interpretation: fraction of the library
            sampled this batch.  ``"batch"`` normalizes by
            ``min(batch_size, n_reagents)`` — interpretation: fraction
            of the maximum achievable diversity.  ``"share"`` normalizes
            each batch's per-component count by the batch's *total*
            unique-reagent count across components, producing a
            stacked-to-1 view.  Use ``"share"`` on imbalanced libraries
            where components differ in size by many-fold, so small
            components don't disappear visually when heated.
        title: Plot title.
        figsize: Figure size.
        heated_label: Y-axis label for the action strip.  Use
            ``"Heated"`` for RWS-family methods (default) and
            ``"Explore"`` for TT-TS, where the strip marks the
            exploring component each batch (cooled / exploit
            components are the unfilled cells).

    Returns:
        Matplotlib Figure.
    """
    search = search_df
    diag = diag_df
    if replicate is not None:
        search = search.filter(pl.col("replicate") == replicate)
        diag = diag.filter(pl.col("replicate") == replicate)

    search = search.filter(pl.col("phase") == "search")
    reagent_id_sets = _load_reagent_ids(reagent_files)
    n_comps = len(reagent_id_sets)
    n_reagents_per_component = [len(s) for s in reagent_id_sets]
    if component_names is None:
        component_names = [
            _infer_component_name(reagent_id_sets[c], c) for c in range(n_comps)
        ]
    elif len(component_names) != n_comps:
        raise ValueError(
            f"component_names has {len(component_names)} entries but "
            f"there are {n_comps} components"
        )

    names = search["Name"].to_list()
    unique_count, new_count = _per_batch_reagent_usage(
        names, reagent_id_sets, batch_size
    )
    n_batches_raw = unique_count.shape[0]
    old_count = unique_count - new_count

    heated_by_cycle: dict[int, int] = {}
    for row in diag.iter_rows(named=True):
        if row["is_heated"]:
            heated_by_cycle[row["current_cycle"]] = row["component_idx"]

    if aggregate_window > 1:
        w = aggregate_window
        n_bins = n_batches_raw // w
        unique_b = np.zeros((n_bins, n_comps))
        new_b = np.zeros((n_bins, n_comps))
        old_b = np.zeros((n_bins, n_comps))
        heated_b: dict[int, int] = {}
        for i in range(n_bins):
            slc = slice(i * w, (i + 1) * w)
            unique_b[i] = unique_count[slc].mean(axis=0)
            new_b[i] = new_count[slc].sum(axis=0)
            old_b[i] = unique_b[i] - new_b[i] / w
            batch_votes: list[int] = []
            for j in range(i * w, (i + 1) * w):
                if j in heated_by_cycle:
                    batch_votes.append(heated_by_cycle[j])
            if batch_votes:
                heated_b[i] = max(set(batch_votes), key=batch_votes.count)
        unique_count = unique_b
        new_count = new_b / w
        old_count = old_b
        heated_by_cycle = heated_b
        n_batches = n_bins
    else:
        n_batches = n_batches_raw

    # Determine plotting cutoff — auto-trim trailing exploitation-only
    # batches, respecting any user-supplied max_batches hard cap.
    n_total = n_batches
    cutoff = n_total
    if auto_trim:
        cutoff = _find_activity_cutoff(new_count, heated_by_cycle)
    if max_batches is not None:
        cutoff = min(cutoff, max_batches)
    cutoff = max(1, min(cutoff, n_total))
    n_trimmed = n_total - cutoff

    if cutoff < n_total:
        unique_count = unique_count[:cutoff]
        new_count = new_count[:cutoff]
        old_count = old_count[:cutoff]
        heated_by_cycle = {
            b: c for b, c in heated_by_cycle.items() if b < cutoff
        }
        n_batches = cutoff

    # normalize is retained for API compatibility but no longer affects
    # the plotted data — the bar panel was replaced with a heating
    # timeline strip.  Validate the value so stale callers still error.
    if normalize not in {"library", "batch", "share"}:
        raise ValueError(
            f"normalize must be 'library', 'batch', or 'share', got {normalize!r}"
        )
    del unique_count, old_count  # no longer used after bar panel removal

    # Line panel: fraction of a component's reagents first-seen in this
    # batch — a scree-plot decay, comparable across components.
    n_reagents_arr = np.asarray(n_reagents_per_component, dtype=float)
    new_frac_lib = new_count / n_reagents_arr[None, :]

    fig, (ax_heat, ax_line) = _pyplot().subplots(
        2, 1, figsize=figsize, sharex=True,
        gridspec_kw={"height_ratios": [0.35, 3], "hspace": 0.08},
    )
    fig.patch.set_facecolor("white")
    ax_heat.set_facecolor("white")
    ax_line.set_facecolor("white")
    colors = _component_palette(n_comps)

    # Top panel: heating timeline — one colored block per batch, colored
    # by the heated component.  Batches with no heated record (e.g. all
    # cooled or missing diagnostics) are left white.
    for b in range(n_batches):
        hc = heated_by_cycle.get(b)
        if hc is None:
            continue
        ax_heat.axvspan(b - 0.5, b + 0.5, color=colors[hc], linewidth=0)
    ax_heat.set_ylim(0, 1)
    ax_heat.set_yticks([])
    ax_heat.set_ylabel(
        heated_label, fontsize=11, rotation=0, ha="right", va="center", labelpad=10,
    )
    for spine in ("top", "right", "left"):
        ax_heat.spines[spine].set_visible(False)
    ax_heat.tick_params(axis="x", length=0)
    ax_heat.set_title(title or "Reagent Usage per Batch", fontsize=15, loc="left")

    if n_trimmed > 0:
        _unit = "batches" if aggregate_window == 1 else "batch bins"
        ax_heat.annotate(
            f"algorithm converged after {cutoff} {_unit}\n"
            f"(remaining {n_trimmed} {_unit} were exploitation-only "
            f"— no heated-component switching or new reagents)",
            xy=(0.995, 2.2),
            xycoords="axes fraction",
            ha="right", va="top",
            fontsize=11, color="#333333",
            bbox=dict(boxstyle="round,pad=0.35",
                      facecolor="white", edgecolor="#bbbbbb", linewidth=0.6),
        )

    # Bottom panel: per-component lines, fraction of reagents first-seen
    # this batch.  Shares x-axis with the heating strip.
    x = np.arange(n_batches)
    for c in range(n_comps):
        ax_line.plot(
            x, new_frac_lib[:, c],
            color=colors[c], linewidth=1.8,
            marker="o", markersize=3.0,
            label=component_names[c], zorder=3,
        )
    ax_line.set_xlabel(
        f"Batch bin (window={aggregate_window})" if aggregate_window > 1 else "Batch",
        fontsize=13,
    )
    ax_line.set_ylabel(
        "Fraction of new reagents\nper reaction component", fontsize=13,
    )
    ax_line.tick_params(axis="both", labelsize=12)
    _line_ymax = (
        max(1e-3, float(new_frac_lib.max()) * 1.12) if new_frac_lib.size else 0.05
    )
    ax_line.set_ylim(0, _line_ymax)
    ax_line.set_xlim(-0.5, n_batches - 0.5)
    ax_line.grid(axis="y", color="#dddddd", linewidth=0.5, zorder=0)

    # Legend: component colors (shared between strip and lines).
    from matplotlib.patches import Patch

    handles = [
        Patch(facecolor=colors[c], edgecolor="none", label=component_names[c])
        for c in range(n_comps)
    ]
    ax_line.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.22),
        ncol=len(handles),
        fontsize=12,
        frameon=False,
        handlelength=1.6,
        columnspacing=1.6,
    )
    fig.tight_layout()
    return fig


def plot_gmic_directed_exploration(
    search_df: pl.DataFrame,
    diag_df: pl.DataFrame,
    reagent_files: list[Path],
    methods: list[str],
    *,
    method_labels=None,
    replicate=None,
    batch_size: int = 100,
    component_names=None,
    component_colors=None,
    gmic_ref: float = 1.0,
    title=None,
    figsize=None,
):
    """Layer 1: how GMIC directs WHICH component to explore (signal -> decision -> result).

    One column per method; three stacked rows sharing the search-cycle x-axis:

      SIGNAL    -- per-component GMIC (criticality).  High = the component's reagent
                   means are well separated relative to noise (resolved -> exploit);
                   low = unresolved (explore).
      DECISION  -- the per-component "explored" strip (``is_heated``): which component
                   the GMIC-weighted rotation sends exploration to each cycle.
      RESULT    -- cumulative per-component reagent coverage (fraction of *that*
                   component's reagents discovered), each normalised to its own size.

    On an imbalanced library (thrombin: 130 acids x 3844 dipeptides) the small
    component resolves immediately (high GMIC -> exploited -> coverage saturates to
    100%) while the large one stays unresolved (low GMIC -> explored -> coverage
    climbs), making the GMIC-directed explore/exploit reallocation visible end to end.

    Args:
        search_df: Per-evaluation search log (``Name``, ``method``, ``replicate``,
            ``phase``); only ``phase == "search"`` rows are used.
        diag_df: Per-cycle diagnostics (``method``, ``replicate``, ``current_cycle``,
            ``component_idx``, ``gmic``, ``is_heated``).
        reagent_files: One ``.smi`` per component (defines each component's reagent IDs).
        methods: Method names to show as columns (e.g. the RWS and TT-TS labels).
        method_labels: Optional ``{method: display label}`` for the column titles.
        replicate: Replicate to plot; defaults to the smallest in ``diag_df``.
        batch_size: Evaluations per cycle/batch (default 100).
        component_names: Display names per component; inferred from reagent-ID
            prefixes when ``None``.
        component_colors: Line colours per component.
        gmic_ref: Reference line on the GMIC row (default 1.0 = critical/flexible boundary).
        title: Figure suptitle.
        figsize: Defaults to ``(7.5 * n_methods, 9.2)``.

    Returns:
        Matplotlib Figure.
    """
    rid = _load_reagent_ids(reagent_files)
    ncomp = len(rid)
    nreag = [len(s) for s in rid]
    if component_names is None:
        component_names = [_infer_component_name(rid[c], c) for c in range(ncomp)]
    if component_colors is None:
        _palette = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#AA3377"]
        component_colors = [_palette[c % len(_palette)] for c in range(ncomp)]
    if method_labels is None:
        method_labels = {}
    if replicate is None:
        replicate = int(diag_df["replicate"].min())

    ncols = len(methods)
    if figsize is None:
        figsize = (7.5 * ncols, 9.2)
    fig, axes = _pyplot().subplots(
        3, ncols, figsize=figsize, squeeze=False,
        gridspec_kw={"height_ratios": [1.1, 0.35, 1.1]},
    )

    for col, m in enumerate(methods):
        d = diag_df.filter(
            (pl.col("method") == m) & (pl.col("replicate") == replicate)
        )
        ncyc = int(d["current_cycle"].max()) + 1
        gmic = np.full((ncyc, ncomp), np.nan)
        heat = np.zeros((ncomp, ncyc))
        for r in d.iter_rows(named=True):
            gmic[r["current_cycle"], r["component_idx"]] = r["gmic"]
            if r["is_heated"]:
                heat[r["component_idx"], r["current_cycle"]] = 1.0

        s = search_df.filter(
            (pl.col("method") == m)
            & (pl.col("replicate") == replicate)
            & (pl.col("phase") == "search")
        )
        _, new = _per_batch_reagent_usage(s["Name"].to_list(), rid, batch_size)
        cov = np.cumsum(new, axis=0) / np.array(nreag) * 100.0
        nb = cov.shape[0]

        a0, a1, a2 = axes[0, col], axes[1, col], axes[2, col]
        a0.set_facecolor("white")
        for c in range(ncomp):
            a0.plot(np.arange(ncyc), gmic[:, c], color=component_colors[c],
                    lw=2, label=f"{component_names[c]} ({nreag[c]})")
        a0.axhline(gmic_ref, color="grey", ls=":", lw=0.8)
        a0.set_title(method_labels.get(m, m), fontsize=13, fontweight="bold",
                     color="black")
        a0.legend(fontsize=9, loc="center right", framealpha=0.9)
        a0.grid(alpha=0.3, color="grey")
        a0.tick_params(colors="black")

        _im = a1.imshow(heat, aspect="auto", cmap="Oranges", vmin=0, vmax=1,
                        extent=[0, ncyc, ncomp - 0.5, -0.5], interpolation="nearest")
        a1.set_yticks(range(ncomp))
        a1.set_yticklabels([cn[:9] for cn in component_names], fontsize=9,
                           color="black")
        a1.tick_params(colors="black")

        a2.set_facecolor("white")
        for c in range(ncomp):
            a2.plot(np.arange(nb), cov[:, c], color=component_colors[c], lw=2)
        a2.set_xlabel("search cycle (batch)", fontsize=11, color="black")
        a2.set_ylim(-3, 103)
        a2.grid(alpha=0.3, color="grey")
        a2.tick_params(colors="black")

        if col == 0:
            a0.set_ylabel("SIGNAL\nGMIC (criticality)", fontsize=11, color="black")
            a1.set_ylabel("DECISION\nexplored", fontsize=11, color="black")
            a2.set_ylabel("RESULT\ncomponent coverage (%)", fontsize=11,
                          color="black")

    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold", color="black", y=1.0)
    fig.tight_layout(rect=(0, 0, 0.93, 1))

    # Colorbar key for the DECISION (explored) strip -- placed in the reserved
    # right margin, vertically aligned with the (middle) decision row, so it
    # does not disturb the shared-x alignment of the three rows.
    _pos = axes[1, -1].get_position()
    _cax = fig.add_axes([0.945, _pos.y0 - 0.02, 0.012, _pos.height + 0.04])
    _cb = fig.colorbar(_im, cax=_cax)
    _cb.set_ticks([0.05, 0.95])
    _cb.ax.set_yticklabels(
        ["exploit\n(cooled)", "explore\n(heated)"], fontsize=8, color="black"
    )
    _cb.ax.tick_params(length=0)
    _cb.set_label("exploration\ndecision", fontsize=8.5, color="black")
    return fig


# Per-method within-component intensity knob.  RWS modulates the Boltzmann
# temperature; TT-TS modulates the posterior-σ inflation scale.  Both are the
# "how hard do I explore *within* the chosen component" lever, one level below
# the GMIC "which component" decision shown by plot_gmic_directed_exploration.
#
# Both knobs are reported as dimensionless multipliers centred on 1.0, so the
# two methods are directly comparable on a shared interpretation (1.0 = no
# adjustment; >1 = amplify exploration):
#   RWS  -> ``cats_multiplier``: the GMIC-driven factor that scales the base
#           thermal-cycling temperature (final_temperature = base_temp x mult).
#   TT-TS -> ``heated_scale``: the disagreement-driven factor that scales the
#           posterior σ.  This is the smooth per-component state the controller
#           sets, NOT the per-cycle ``effective_scale`` (= heated_scale on heated
#           cycles, cooled_scale otherwise), which just re-encodes the heating
#           rotation already shown in the DECISION row.
# ``effective_scale`` / ``final_temperature`` are kept only as last-resort
# fallbacks for older diagnostics that lack the multiplier columns.
_INTENSITY_SPECS = (
    # (column, display label, reference value)
    ("cats_multiplier", "CATS multiplier (×T)", 1.0),
    ("heated_scale", "σ-inflation multiplier (×σ)", 1.0),
    ("effective_scale", "σ-inflation scale (effective)", 1.0),
    ("final_temperature", "Boltzmann temperature", None),
)


def _detect_intensity_col(d: pl.DataFrame) -> tuple[str, str, Optional[float]]:
    """Pick the within-component intensity column present for these rows.

    RWS diagnostics carry ``final_temperature`` (TT-TS columns null) and TT-TS
    diagnostics carry ``effective_scale`` (RWS columns null); the canonical
    benchmark stores both methods in one union-schema parquet.  Returns the
    ``(column, label, reference)`` of the first spec whose column exists and has
    at least one non-null value.
    """
    for col, label, ref in _INTENSITY_SPECS:
        if col in d.columns and d[col].drop_nulls().len() > 0:
            return col, label, ref
    raise ValueError(
        "No within-component intensity column "
        f"({', '.join(s[0] for s in _INTENSITY_SPECS)}) found in diagnostics."
    )


def plot_adaptive_intensity(
    search_df: pl.DataFrame,
    diag_df: pl.DataFrame,
    reagent_files: list[Path],
    methods: list[str],
    *,
    method_labels=None,
    replicate=None,
    batch_size: int = 100,
    component_names=None,
    component_colors=None,
    smooth_window: int = 9,
    show_decision: bool = True,
    result: str = "new",
    title=None,
    figsize=None,
):
    """Layer 2: how each method tunes exploration *within* a component.

    Where :func:`plot_gmic_directed_exploration` (layer 1) shows WHICH component
    GMIC sends exploration to, this shows HOW HARD each method explores the
    chosen component, and what that yields in reagent discovery.  One column per
    method; rows share the search-cycle x-axis:

      INTENSITY -- the per-component within-component knob, smoothed.  Both knobs
                   are reported as **dimensionless multipliers centred on 1.0**
                   (1.0 = no adjustment; >1 = amplify exploration), so the two
                   methods read on a common interpretation: RWS modulates the
                   **CATS temperature multiplier** (``cats_multiplier``, ×base
                   temperature -- the GMIC-driven factor); TT-TS modulates the
                   **σ-inflation multiplier** (``heated_scale``, ×posterior σ --
                   the disagreement-driven factor).  Note the asymmetry the knobs
                   reveal: RWS's multiplier is GMIC-driven, so it *amplifies above
                   1.0* on the unresolved (low-GMIC) component; TT-TS's heated_scale
                   is disagreement-driven, and because two-sample disagreement
                   saturates (>0.8) on both components of a large library it
                   *decays toward 1.0* -- most on the highest-disagreement
                   (unresolved) component -- i.e. the TT-TS adaptation is largely
                   one-directional and cannot manufacture exploration beyond what
                   the raw posteriors already provide.
      DECISION  -- per-component **heating probability** (rolling fraction of
                   cycles in which the component was the explored/heated one).
                   Omitted when ``show_decision=False`` (leaving INTENSITY +
                   RESULT, i.e. "just the temperature").
      RESULT    -- per-component within-component discovery.  ``result="new"``
                   (default) plots **new reagents discovered per batch** (the
                   per-batch *rate* -- the discrete increment of coverage);
                   ``result="coverage"`` plots **cumulative coverage** (% of the
                   component's reagents seen at least once -- the running total,
                   identical to the layer-1 RESULT row).  New-per-batch is the
                   derivative of coverage: it falls toward zero as the reagent
                   pool is exhausted even while coverage is still rising.

    On thrombin (130 acids × 3844 dipeptides) both methods keep discovering new
    dipeptides long after the small acid component is exhausted -- the same
    reallocation, reached through two different within-component mechanisms.

    Args:
        search_df: Per-evaluation search log (``Name``, ``method``,
            ``replicate``, ``phase``); only ``phase == "search"`` rows are used.
        diag_df: Per-cycle diagnostics (``method``, ``replicate``,
            ``current_cycle``, ``component_idx``, ``is_heated`` plus the
            method-specific intensity column ``cats_multiplier`` for RWS or
            ``heated_scale`` for TT-TS).
        reagent_files: One ``.smi`` per component (defines each component's
            reagent IDs).
        methods: Method names to show as columns.
        method_labels: Optional ``{method: display label}`` for column titles.
        replicate: Replicate to plot; defaults to the smallest in ``diag_df``.
        batch_size: Evaluations per cycle/batch (default 100).
        component_names: Display names per component; inferred from reagent-ID
            prefixes when ``None``.
        component_colors: Line colours per component.
        smooth_window: Centred (edge-shrinking) rolling-mean window (cycles)
            applied to every row so the per-cycle heated/cooled square-wave reads
            as an effective average (default 9).
        show_decision: Include the heating-probability row (default True).
        result: ``"new"`` (new reagents per batch, default) or ``"coverage"``
            (cumulative % coverage, as in layer 1).
        title: Figure suptitle.
        figsize: Defaults to ``(7.5 * n_methods, 9.2 or 7.0)``.

    Returns:
        Matplotlib Figure.
    """
    if result not in ("new", "coverage"):
        raise ValueError(f"result must be 'new' or 'coverage', got {result!r}")
    rid = _load_reagent_ids(reagent_files)
    ncomp = len(rid)
    nreag = [len(s) for s in rid]
    if component_names is None:
        component_names = [_infer_component_name(rid[c], c) for c in range(ncomp)]
    if component_colors is None:
        _palette = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#AA3377"]
        component_colors = [_palette[c % len(_palette)] for c in range(ncomp)]
    if method_labels is None:
        method_labels = {}
    if replicate is None:
        replicate = int(diag_df["replicate"].min())

    ncols = len(methods)
    nrows = 3 if show_decision else 2
    if figsize is None:
        figsize = (7.5 * ncols, 9.2 if show_decision else 7.0)
    height_ratios = [1.1, 0.35, 1.1] if show_decision else [1.1, 1.1]
    fig, axes = _pyplot().subplots(
        nrows, ncols, figsize=figsize, squeeze=False,
        gridspec_kw={"height_ratios": height_ratios},
    )
    decision_im = None

    for col, m in enumerate(methods):
        d = diag_df.filter(
            (pl.col("method") == m) & (pl.col("replicate") == replicate)
        )
        icol, ilabel, iref = _detect_intensity_col(d)
        ncyc = int(d["current_cycle"].max()) + 1
        intensity = np.full((ncyc, ncomp), np.nan)
        heat = np.zeros((ncyc, ncomp))
        for r in d.iter_rows(named=True):
            intensity[r["current_cycle"], r["component_idx"]] = r[icol]
            if r["is_heated"]:
                heat[r["current_cycle"], r["component_idx"]] = 1.0

        s = search_df.filter(
            (pl.col("method") == m)
            & (pl.col("replicate") == replicate)
            & (pl.col("phase") == "search")
        )
        _, new = _per_batch_reagent_usage(s["Name"].to_list(), rid, batch_size)
        nb = new.shape[0]
        if result == "coverage":
            res = np.cumsum(new, axis=0) / np.array(nreag) * 100.0
            res_label = "RESULT\ncomponent coverage (%)"
        else:
            res = new.astype(float)
            res_label = "RESULT\nnew reagents / batch"

        a_int = axes[0, col]
        a_res = axes[nrows - 1, col]

        # --- INTENSITY row -------------------------------------------------
        a_int.set_facecolor("white")
        for c in range(ncomp):
            a_int.plot(
                np.arange(ncyc), _smooth_trim(intensity[:, c], smooth_window),
                color=component_colors[c], lw=2, label=component_names[c],
            )
        if iref is not None:
            a_int.axhline(iref, color="grey", ls=":", lw=0.8)
        a_int.set_title(method_labels.get(m, m), fontsize=13, fontweight="bold",
                        color="black")
        a_int.legend(fontsize=9, loc="best", framealpha=0.9)
        a_int.grid(alpha=0.3, color="grey")
        a_int.tick_params(colors="black")
        a_int.set_ylabel(f"INTENSITY\n{ilabel}", fontsize=11, color="black")

        # --- DECISION row (heating probability) ----------------------------
        if show_decision:
            a_dec = axes[1, col]
            heat_prob = np.column_stack(
                [_smooth_trim(heat[:, c], smooth_window) for c in range(ncomp)]
            ).T  # shape (ncomp, ncyc)
            decision_im = a_dec.imshow(
                heat_prob, aspect="auto", cmap="Oranges", vmin=0, vmax=1,
                extent=[0, ncyc, ncomp - 0.5, -0.5], interpolation="nearest",
            )
            a_dec.set_yticks(range(ncomp))
            a_dec.set_yticklabels([cn[:9] for cn in component_names], fontsize=9,
                                  color="black")
            a_dec.tick_params(colors="black")
            if col == 0:
                a_dec.set_ylabel("DECISION\nP(explore)", fontsize=11,
                                 color="black")

        # --- RESULT row (within-component discovery) -----------------------
        a_res.set_facecolor("white")
        for c in range(ncomp):
            a_res.plot(
                np.arange(nb), _smooth_trim(res[:, c], smooth_window),
                color=component_colors[c], lw=2,
            )
        a_res.set_xlabel("search cycle (batch)", fontsize=11, color="black")
        if result == "coverage":
            a_res.set_ylim(-3, 103)
        else:
            a_res.set_ylim(bottom=-0.5)
        a_res.grid(alpha=0.3, color="grey")
        a_res.tick_params(colors="black")
        if col == 0:
            a_res.set_ylabel(res_label, fontsize=11, color="black")

    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold", color="black", y=1.0)
    fig.tight_layout(rect=(0, 0, 0.93, 1) if show_decision else None)

    if show_decision and decision_im is not None:
        _pos = axes[1, -1].get_position()
        _cax = fig.add_axes([0.945, _pos.y0 - 0.02, 0.012, _pos.height + 0.04])
        _cb = fig.colorbar(decision_im, cax=_cax)
        _cb.set_ticks([0.05, 0.95])
        _cb.ax.set_yticklabels(["never\nexplore", "always\nexplore"], fontsize=8,
                               color="black")
        _cb.ax.tick_params(length=0)
        _cb.set_label("heating\nprobability", fontsize=8.5, color="black")
    return fig

