"""
TACTICS Manuscript — Reagent Score Landscape Plots (Docking only)

Heavy plots that read multi-million-row brute-force ground-truth parquets and
compute per-reagent posterior statistics. Separated from manuscript_plots.py
to avoid blocking iteration on the recovery/diversity plots in that notebook.

Coverage: 8 docking libraries (4 FRED + 4 HYBRID) for both the single-library
explorer and the multi-top-N grid.

Run as app: marimo run tutorials/manuscript_sar_plots.py
Edit mode:  marimo edit tutorials/manuscript_sar_plots.py
"""

import marimo

__generated_with = "0.21.1"
app = marimo.App(
    width="full",
    app_title="TACTICS Manuscript — Reagent Score Landscape",
)


@app.cell
def _():
    """Imports."""
    import io
    import marimo as mo
    import sys
    from pathlib import Path

    import matplotlib.pyplot as plt
    import numpy as np
    import polars as pl

    def fig_to_pdf_bytes(fig):
        """Convert a matplotlib figure to PDF bytes in memory.

        Includes any per-axes legends as bbox_extra_artists so legends placed
        OUTSIDE the axes (e.g. the M-S map's right-gutter legends) are not
        clipped by the tight bbox.
        """
        from matplotlib.legend import Legend as _Legend

        buf = io.BytesIO()
        # Collect ALL Legend artists (an axes can hold several via add_artist),
        # plus any figure-level legend.
        _extra = []
        for _ax in fig.axes:
            _extra.extend(_c for _c in _ax.get_children() if isinstance(_c, _Legend))
        _extra.extend(_lg for _lg in fig.legends if isinstance(_lg, _Legend))
        fig.savefig(buf, format="pdf", bbox_inches="tight", dpi=300,
                    bbox_extra_artists=_extra or None)
        return buf.getvalue()

    try:
        project_root = Path(__file__).parent.parent.resolve()
    except NameError:
        project_root = Path.cwd()

    sys.path.insert(0, str(project_root))
    sys.path.insert(0, str(project_root / "src"))
    return fig_to_pdf_bytes, mo, np, pl, plt, project_root


@app.cell
def _(mo):
    """SAR section header + per-library dropdown (docking libraries only)."""
    _docking_libs = [
        "adenine", "adenine_hybrid",
        "amide", "amide_hybrid",
        "quinazoline", "quinazoline_hybrid",
        "thrombin", "thrombin_hybrid",
    ]
    sar_library_selector = mo.ui.dropdown(
        options=_docking_libs,
        value="adenine",
        label="Library",
    )
    sar_topn = mo.ui.slider(start=25, stop=250, step=25, value=100, label="Top-N")
    mo.vstack([
        mo.md(r"""
        # SAR Landscape Characterization — Reagent Posterior Separation

        Riley (2024) showed that for ROCS, reagents producing top compounds have high
        means and low variance ("stand out from the crowd"), while for docking, top-hit
        reagents are statistically indistinguishable from mediocre ones (mean-uninformative).
        Below: per-reagent (z(mean), z(std)) scatter for the selected docking library,
        colored by hit count in the chosen top-N.
        """),
        mo.hstack([sar_library_selector, sar_topn], gap=2),
    ])
    return sar_library_selector, sar_topn


@app.cell
def _(
    fig_to_pdf_bytes,
    mo,
    np,
    pl,
    plt,
    project_root,
    sar_library_selector,
    sar_topn,
):
    """Reagent posterior separation plot (Riley-style) for the selected docking library.

    For each reagent, computes (mean score, std score) across all products containing
    that reagent. Colors reagents by how many times they appear in top-N compounds.
    """
    from matplotlib.colors import Normalize as _Normalize
    from matplotlib import cm as _cm
    import re as _re

    _SCORES_DIR = project_root / "data" / "scores"
    _lib_id = sar_library_selector.value
    _TOP_N = sar_topn.value

    def _parse_adenine(code):
        _m = _re.match(r"(amidine_\d+)_(isocyanide_db_\d+)_(aldehyde_\d+)", code)
        return _m.groups() if _m else None

    def _parse_thrombin(code):
        _pts = code.split("_")
        return (_pts[0], _pts[1] + "_" + _pts[2]) if len(_pts) == 3 else None

    def _parse_split(code, n):
        _pts = code.split("_")
        return tuple(_pts) if len(_pts) == n else None

    _CONFIGS = {
        "adenine":           {"components": ["Amidines", "Isocyanides", "Aldehydes"],
                              "parse": _parse_adenine},
        "adenine_hybrid":    {"components": ["Amidines", "Isocyanides", "Aldehydes"],
                              "parse": _parse_adenine},
        "quinazoline":       {"components": ["Aminobenzoics", "Amines", "Acids"],
                              "parse": lambda c: _parse_split(c, 3)},
        "quinazoline_hybrid": {"components": ["Aminobenzoics", "Amines", "Acids"],
                              "parse": lambda c: _parse_split(c, 3)},
        "thrombin":          {"components": ["Acids", "Dipeptides"],
                              "parse": _parse_thrombin},
        "thrombin_hybrid":   {"components": ["Acids", "Dipeptides"],
                              "parse": _parse_thrombin},
        "amide":             {"components": ["Amines", "Acids"],
                              "parse": lambda c: _parse_split(c, 2)},
        "amide_hybrid":      {"components": ["Amines", "Acids"],
                              "parse": lambda c: _parse_split(c, 2)},
    }

    _cfg = _CONFIGS[_lib_id]
    _comp_names = _cfg["components"]
    _n_comp = len(_comp_names)
    _parse_fn = _cfg["parse"]

    _score_df = pl.read_parquet(_SCORES_DIR / f"{_lib_id}.parquet")
    _score_df = _score_df.select(["Product_Code", "Scores"]).filter(pl.col("Scores").is_not_nan())

    _codes = _score_df["Product_Code"].to_list()
    _scores_arr = _score_df["Scores"].to_numpy()
    _parsed = [_parse_fn(c) for c in _codes]

    _sorted_idx = np.argsort(_scores_arr)  # ascending = best for docking minimize
    _top_codes = [_codes[i] for i in _sorted_idx[:_TOP_N]]
    _top_parsed = [_parse_fn(c) for c in _top_codes]

    _cmap = _cm.YlOrRd
    _global_max_count = 0
    _plot_data = []

    for _ci in range(_n_comp):
        _reagent_to_idx = {}
        for _i, _p in enumerate(_parsed):
            if _p is not None:
                _reagent_to_idx.setdefault(_p[_ci], []).append(_i)

        _top_counts = {}
        for _p in _top_parsed:
            if _p is not None:
                _top_counts[_p[_ci]] = _top_counts.get(_p[_ci], 0) + 1

        _unique = sorted(_reagent_to_idx.keys())
        _means = np.array([float(np.mean(_scores_arr[_reagent_to_idx[r]])) for r in _unique])
        _stds = np.array([float(np.std(_scores_arr[_reagent_to_idx[r]])) for r in _unique])
        _counts = np.array([_top_counts.get(r, 0) for r in _unique])

        _hc = _counts[_counts > 0]
        if len(_hc) > 0:
            _global_max_count = max(_global_max_count, int(_hc.max()))

        _sig = float(np.var(_means))
        _noise = float(np.mean(_stds**2)) if len(_stds) > 0 else 1.0
        _gmic = 0.5 * np.log(1 + _sig / _noise) if _noise > 0 else 0.0
        _plot_data.append((_ci, _means, _stds, _counts, _unique, _gmic))

    _gs = plt.GridSpec(1, _n_comp + 1, width_ratios=[1] * _n_comp + [0.05], wspace=0.35)
    _fig = plt.figure(figsize=(7 * _n_comp, 6), facecolor="white")
    _axes = [_fig.add_subplot(_gs[0, c]) for c in range(_n_comp)]
    _norm = _Normalize(vmin=1, vmax=max(_global_max_count, 2))

    for _ci, _means, _stds, _counts, _unique, _gmic in _plot_data:
        _ax = _axes[_ci]
        _ax.set_facecolor("white")

        _m_mu, _m_sig = np.mean(_means), np.std(_means)
        _s_mu, _s_sig = np.mean(_stds), np.std(_stds)
        _zm = (_means - _m_mu) / _m_sig if _m_sig > 0 else _means * 0
        _zs = (_stds - _s_mu) / _s_sig if _s_sig > 0 else _stds * 0

        _zero = _counts == 0
        _ax.scatter(_zm[_zero], _zs[_zero], c="#CCCCCC", s=40, alpha=0.6,
                    edgecolors="none", zorder=2)

        _hit = _counts > 0
        if _hit.any():
            _ax.scatter(_zm[_hit], _zs[_hit], c=_counts[_hit], cmap=_cmap,
                        norm=_norm, s=60, alpha=0.9, edgecolors="black",
                        linewidths=0.5, zorder=3)

        _ax.set_xlabel("Reagent Mean Score (z)", fontsize=14, color="black")
        if _ci == 0:
            _ax.set_ylabel("Reagent Stdev (z)", fontsize=14, color="black")
        _n_hit = int(_hit.sum())
        _ax.set_title(
            f"{_comp_names[_ci]} ({_n_hit}/{len(_unique)}) — GMIC {_gmic:.3f}",
            fontsize=13, fontweight="bold", color="black",
        )
        _ax.tick_params(colors="black", labelsize=12)
        _ax.grid(alpha=0.3, color="grey", zorder=1)
        for _sp in _ax.spines.values():
            _sp.set_color("black")

    _cbar_ax = _fig.add_subplot(_gs[0, _n_comp])
    _cbar = _fig.colorbar(_cm.ScalarMappable(norm=_norm, cmap=_cmap), cax=_cbar_ax)
    _cbar.set_label(f"Count in top-{_TOP_N}", fontsize=14, color="black")
    _cbar.ax.tick_params(colors="black", labelsize=12)

    _fig.suptitle(
        f"{_lib_id.capitalize()} — Reagent Score Distributions vs Top-{_TOP_N} (Docking)",
        fontsize=16, fontweight="bold", color="black", y=1.02,
    )
    _fig.tight_layout()

    mo.vstack([
        plt.gcf(),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename=f"reagent_posterior_separation_{_lib_id}.pdf",
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


@app.cell
def _(mo):
    """ROCS SAR section header + per-library dropdown + top-N slider."""
    _rocs_libs = [
        "rxn101", "rxn102a", "rxn108b", "rxn111b", "rxn114b",
        "rxn203", "rxn205", "rxn206", "rxn207", "rxn208",
        "amide-suzuki", "betti", "dobener", "groebke-blackburn-bienayme",
        "mannich", "niementowski", "orru", "passerini", "petasis", "poparov",
    ]
    sar_rocs_library_selector = mo.ui.dropdown(
        options=_rocs_libs,
        value="rxn108b",
        label="ROCS Library",
    )
    sar_rocs_topn = mo.ui.slider(
        start=25, stop=2500, step=25, value=2500, label="Top-N"
    )
    mo.vstack([
        mo.md(r"""
        # SAR Landscape — ROCS Libraries (Mean Across All Queries)

        For ROCS-similarity libraries, each library has 109 query molecules. The reagent
        posterior used here is the per-reagent (mean similarity, std similarity)
        computed by **pooling all (product × query) similarity scores** across every
        product containing the reagent and every query.

        This mean-SAR view is paired with the per-library mean diversity used in the
        Zhao Fig 5 analysis — both averaged across queries — to ensure structural
        consistency between the diversity-vs-recovery story and the per-library SAR
        characterization.

        Reagents are coloured by **total hit count in top-N across all queries** (sum
        over the 109 queries of how many times the reagent appears in that query's
        top-N ground truth). The hit-count coloring serves as a complementary
        check: if the reagents flagged as hidden hubs or mean-separable by the
        mean-SAR view are also consistently hit across queries, the library-level
        aggregation is a representative summary of the actual enriched building
        blocks.
        """),
        mo.hstack([sar_rocs_library_selector, sar_rocs_topn], gap=2),
    ])
    return sar_rocs_library_selector, sar_rocs_topn


@app.cell
def _(
    fig_to_pdf_bytes,
    mo,
    np,
    pl,
    plt,
    project_root,
    sar_rocs_library_selector,
    sar_rocs_topn,
):
    """Mean SAR landscape for the selected ROCS library.

    Per reagent at each component: pool all (product × query) similarity scores
    across all products containing that reagent and across all 109 queries.
    Mean and std define the Riley-style (z(mean), z(std)) point. Colour = sum
    across queries of per-query top-N hit count.
    """
    from matplotlib.colors import Normalize as _Normalize
    from matplotlib import cm as _cm

    _ROCS_NCOMP = {
        "rxn101": 2, "rxn102a": 2, "rxn108b": 2, "rxn111b": 2, "rxn114b": 2,
        "rxn203": 2, "rxn205": 2, "rxn206": 2, "rxn207": 2, "rxn208": 2,
        "amide-suzuki": 3, "betti": 3, "dobener": 3,
        "groebke-blackburn-bienayme": 3, "mannich": 3, "niementowski": 3,
        "orru": 3, "passerini": 3, "petasis": 3, "poparov": 3,
    }

    _SCORES_DIR = project_root / "data" / "scores"
    _lib_id = sar_rocs_library_selector.value
    _TOP_N = sar_rocs_topn.value
    _n_comp = _ROCS_NCOMP[_lib_id]
    _comp_names = [f"Component {i + 1}" for i in range(_n_comp)]

    _df = pl.read_parquet(_SCORES_DIR / f"{_lib_id}.parquet")
    _name_col = "Product_Code" if "Product_Code" in _df.columns else "Name"
    _query_cols = sorted(c for c in _df.columns if c.startswith("query_"))

    _codes = _df[_name_col].to_list()
    _parsed = []
    for _c in _codes:
        _pts = _c.split("_")
        _parsed.append(tuple(_pts) if len(_pts) == _n_comp else None)

    # Map reagent → list of product row indices, per component
    _reagent_to_idx = [{} for _ in range(_n_comp)]
    for _i, _p in enumerate(_parsed):
        if _p is None:
            continue
        for _ci in range(_n_comp):
            _reagent_to_idx[_ci].setdefault(_p[_ci], []).append(_i)

    # Materialize a (n_products × n_queries) score matrix
    _qarrs = [_df[_q].to_numpy() for _q in _query_cols]
    _qstack = np.column_stack(_qarrs)  # shape (n_products, n_queries)

    # Per-query top-N hit counts per reagent per component
    _hit_per_comp = [{} for _ in range(_n_comp)]
    for _qi in range(len(_query_cols)):
        _scores = _qstack[:, _qi]
        _valid = ~np.isnan(_scores)
        _valid_idx = np.where(_valid)[0]
        # Descending sort for similarity (higher is better)
        _sorted_idx = _valid_idx[np.argsort(-_scores[_valid_idx])][:_TOP_N]
        for _idx in _sorted_idx:
            _p = _parsed[_idx]
            if _p is None:
                continue
            for _ci in range(_n_comp):
                _hit_per_comp[_ci][_p[_ci]] = _hit_per_comp[_ci].get(_p[_ci], 0) + 1

    _cmap = _cm.YlOrRd
    _global_max_count = 0
    _plot_data = []

    for _ci in range(_n_comp):
        _unique = sorted(_reagent_to_idx[_ci].keys())
        _means = np.empty(len(_unique))
        _stds = np.empty(len(_unique))
        for _ri, _r in enumerate(_unique):
            _idx = _reagent_to_idx[_ci][_r]
            _block = _qstack[_idx, :]
            _flat = _block[~np.isnan(_block)]
            if len(_flat) == 0:
                _means[_ri] = np.nan
                _stds[_ri] = np.nan
            else:
                _means[_ri] = float(np.mean(_flat))
                _stds[_ri] = float(np.std(_flat))
        # Mean hits per query (sum across queries / n_queries). Aligns the
        # coloring convention with the diversity-vs-recovery analysis, which
        # also uses per-query means.
        _counts = np.array([_hit_per_comp[_ci].get(r, 0) for r in _unique], dtype=float) / max(len(_query_cols), 1)

        _hc = _counts[_counts > 0]
        if len(_hc) > 0:
            _global_max_count = max(_global_max_count, float(_hc.max()))

        _vmask = ~np.isnan(_means) & ~np.isnan(_stds)
        if _vmask.sum() > 0:
            _sig = float(np.var(_means[_vmask]))
            _noise = float(np.mean(_stds[_vmask] ** 2))
            _gmic = 0.5 * np.log(1 + _sig / _noise) if _noise > 0 else 0.0
        else:
            _gmic = 0.0
        _plot_data.append((_ci, _means, _stds, _counts, _unique, _gmic))

    _gs = plt.GridSpec(1, _n_comp + 1, width_ratios=[1] * _n_comp + [0.05], wspace=0.35)
    _fig = plt.figure(figsize=(7 * _n_comp, 6), facecolor="white")
    _axes = [_fig.add_subplot(_gs[0, c]) for c in range(_n_comp)]
    _norm = _Normalize(vmin=0.0, vmax=max(_global_max_count, 0.1))

    for _ci, _means, _stds, _counts, _unique, _gmic in _plot_data:
        _ax = _axes[_ci]
        _ax.set_facecolor("white")

        _vmask = ~np.isnan(_means) & ~np.isnan(_stds)
        _vmeans = _means[_vmask]
        _vstds = _stds[_vmask]
        _vcounts = _counts[_vmask]

        _m_mu, _m_sig = float(np.mean(_vmeans)), float(np.std(_vmeans))
        _s_mu, _s_sig = float(np.mean(_vstds)), float(np.std(_vstds))
        _zm = (_vmeans - _m_mu) / _m_sig if _m_sig > 0 else _vmeans * 0
        _zs = (_vstds - _s_mu) / _s_sig if _s_sig > 0 else _vstds * 0

        _zero = _vcounts == 0
        _ax.scatter(_zm[_zero], _zs[_zero], c="#CCCCCC", s=40, alpha=0.6,
                    edgecolors="none", zorder=2)
        _hit = _vcounts > 0
        if _hit.any():
            _ax.scatter(_zm[_hit], _zs[_hit], c=_vcounts[_hit], cmap=_cmap,
                        norm=_norm, s=60, alpha=0.9, edgecolors="black",
                        linewidths=0.5, zorder=3)

        _ax.set_xlabel("Reagent Mean Similarity (z)", fontsize=14, color="black")
        if _ci == 0:
            _ax.set_ylabel("Reagent Stdev (z)", fontsize=14, color="black")
        _n_hit = int(_hit.sum())
        _ax.set_title(
            f"{_comp_names[_ci]} ({_n_hit}/{int(_vmask.sum())}) — GMIC {_gmic:.3f}",
            fontsize=13, fontweight="bold", color="black",
        )
        _ax.tick_params(colors="black", labelsize=12)
        _ax.grid(alpha=0.3, color="grey", zorder=1)
        for _sp in _ax.spines.values():
            _sp.set_color("black")

    _cbar_ax = _fig.add_subplot(_gs[0, _n_comp])
    _cbar = _fig.colorbar(_cm.ScalarMappable(norm=_norm, cmap=_cmap), cax=_cbar_ax)
    _cbar.set_label(
        f"Mean hits per query in top-{_TOP_N} (averaged across 109 queries)",
        fontsize=13, color="black",
    )
    _cbar.ax.tick_params(colors="black", labelsize=12)

    _fig.suptitle(
        f"{_lib_id} — Mean SAR Across Queries vs Top-{_TOP_N} Hit Frequency (ROCS)",
        fontsize=16, fontweight="bold", color="black", y=1.02,
    )
    _fig.tight_layout()

    mo.vstack([
        plt.gcf(),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename=f"reagent_posterior_separation_rocs_{_lib_id}_top{_TOP_N}.pdf",
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


@app.cell
def _(fig_to_pdf_bytes, mo, np, plt, project_root):
    """z* sensitivity of the COVERAGE axis — why the reachability reference is 1.5, not 0.9/1.0.

    The coverage coordinate is x = max(0, z* + z(mean)); a library leaves the origin
    (is flagged "coverage-failing") once z* exceeds its activation threshold t = -z(mean).
    Two things behave differently as z* moves:
      - the RANK ORDER of flagged libraries by severity is invariant (severity is monotone
        in z(mean)) -- this is what stayed put when z* was swept up from 1.0;
      - the SET of flagged libraries (those off the origin) GROWS with z* -- this is what
        actually moves points on the landscape.
    Lowering z* to 0.9 or 1.0 drops adenine_hybrid back onto the "easy" origin even though
    greedy genuinely fails it (Legacy Greedy ~ 30 vs TACTICS ~ 95 top-100). The window that
    flags BOTH coverage libraries (adenine, adenine_hybrid) while keeping the mean-separable
    libraries (quinazoline ...) at the origin is z* in (1.04, 1.67); 1.5 sits inside it.

    NB: x is the COVERAGE coordinate ONLY (one of three failure modes). A high coverage-
    reachability does NOT imply greedy-easy -- both amide_hybrid (allocation, high-sigma hub)
    and thrombin (imbalance, 30:1) sit far-right (mean-reachable hubs) yet are greedy-hard via
    OTHER modes that live on the main landscape (y-axis / edge). Recovery is explained by the
    COMBINATION of modes, never this one axis. Colour encodes the OPERATIVE mode (red coverage
    / blue allocation / orange imbalance / grey none), computed from the three severities -- so
    a non-grey point far-right is coherent (coverage-reachable, hard via another mode).
    Source: dominant-hub z(mean) + greedy/TACTICS top-100 recovery, landscape probe JSONs
    (experiments/_debug/docking_landscape_probe.py).
    """
    import json as _json
    from matplotlib.lines import Line2D as _Line2D
    from matplotlib.patches import Patch as _Patch

    _PROBE_DIR = project_root / "outputs" / "debug" / "docking_landscape_probe"
    _MAP_LIBS = [
        "adenine", "adenine_hybrid", "quinazoline", "quinazoline_hybrid",
        "thrombin", "thrombin_hybrid", "amide", "amide_hybrid",
    ]
    _TT = "TACTICS Enhanced-TT-TS (GMIC)"
    _RWS = "TACTICS Enhanced-RWS (GMIC)"
    _LGRD = "Legacy Standard-Greedy"
    _Z_CANDIDATES = [0.9, 1.0, 1.5]

    _rows = []
    for _lab in _MAP_LIBS:
        with open(_PROBE_DIR / f"{_lab}.json") as _fh:
            _j = _json.load(_fh)
        _h = _j["dominant_hub"]
        _zm = float(_h["z_mean"])
        _m2 = _j["M2_concentration_by_topn"]["100"][f"comp_{_h['component_idx']}"]
        _counts = list(_j["n_reagents"].values())
        _pm = _j["M6_search"]["per_method"]
        _grd = _pm[_LGRD]["topn_recovered_mean"]
        _tac = max(_pm[_TT]["topn_recovered_mean"], _pm[_RWS]["topn_recovered_mean"])
        _rows.append({
            "lab": _lab, "zmean": _zm, "thr": -_zm,
            "grd": _grd, "tac": _tac, "delta": _tac - _grd,
            # the three landscape severities (same computation as the landscape cell)
            "cov": max(0.0, 1.5 + _zm),
            "alloc": _h["z_std"] * (1.0 - _m2["top1_share"]),
            "imb": max(_counts) / min(_counts),
        })
    _rows.sort(key=lambda _r: _r["thr"])  # ascending activation threshold == descending z(mean)

    # Rank-order invariance: the off-origin libraries at any z* are a prefix of the same
    # z(mean)-sorted master list -- their relative order never inverts.
    _master = [_r["lab"] for _r in _rows]

    def _flagged_prefix(_zs):
        return [_r["lab"] for _r in _rows if _zs + _r["zmean"] > 0]

    _rank_invariant = all(
        _flagged_prefix(_zs) == _master[:len(_flagged_prefix(_zs))]
        for _zs in [1.0, 1.25, 1.5, 1.75, 2.0]
    )

    # Colour = OPERATIVE failure mode, computed from the three landscape severities (NOT
    # from the recovery magnitude -- that mis-binned thrombin, a 30:1 imbalance library, as
    # "competitive" merely because its greedy gap is modest). Priority coverage > allocation
    # > imbalance picks the dominant mode when several co-occur (e.g. adenine is coverage AND
    # 24:1 imbalanced -> coverage). Thresholds are display choices; the continuous severities
    # live on the main landscape. Palette matches the landscape: coverage=red, allocation=
    # blue, imbalance=orange.
    _COV_T, _ALLOC_T, _IMB_T = 0.30, 4.0, 15.0
    _C_COV, _C_ALLOC, _C_IMB, _C_EASY = "#B23030", "#2A5DA0", "#E8820E", "#888888"

    def _mech(_r):
        if _r["cov"] > _COV_T:
            return "coverage"
        if _r["alloc"] > _ALLOC_T:
            return "allocation"
        if _r["imb"] > _IMB_T:
            return "imbalance"
        return "easy"

    _MECH_COLOR = {"coverage": _C_COV, "allocation": _C_ALLOC,
                   "imbalance": _C_IMB, "easy": _C_EASY}

    def _color(_r):
        return _MECH_COLOR[_mech(_r)]

    # valid z* window: flag the coverage libraries but not the no-failure (easy) ones
    _win_lo = max(_r["thr"] for _r in _rows if _mech(_r) == "coverage")   # 1.04 (adenine_hybrid)
    _win_hi = min(_r["thr"] for _r in _rows if _mech(_r) == "easy")       # 1.67 (quinazoline)

    _fig, _ax = plt.subplots(figsize=(11.5, 6.2), facecolor="white")
    _ax.set_facecolor("white")
    for _i, _r in enumerate(_rows):
        _c = _color(_r)
        _ax.plot([0, _r["thr"]], [_i, _i], color=_c, lw=2.0, alpha=0.55, zorder=2)
        _ax.scatter(_r["thr"], _i, s=130, color=_c, edgecolors="black",
                    linewidths=0.8, zorder=3)
        _ax.annotate(f'{_r["lab"]}   (greedy {_r["grd"]:.0f} -> TACTICS {_r["tac"]:.0f})',
                     (_r["thr"], _i), xytext=(8, 0), textcoords="offset points",
                     va="center", ha="left", fontsize=9.5, color="black", zorder=4)

    _ax.axvspan(_win_lo, _win_hi, color="#9CC16A", alpha=0.20, zorder=0)
    _ax.annotate(f"valid window\nz* in ({_win_lo:.2f}, {_win_hi:.2f})",
                 ((_win_lo + _win_hi) / 2.0, len(_rows) - 0.4), ha="center", va="top",
                 fontsize=9.5, color="#2F6B2F", style="italic", zorder=1)

    # stagger the 0.9/1.0 labels vertically (they are only 0.1 apart on the x-axis)
    _zlab_y = {0.9: -0.65, 1.0: -1.15, 1.5: -0.65}
    for _zs in _Z_CANDIDATES:
        _in_win = _win_lo < _zs < _win_hi
        _ax.axvline(_zs, ls="--", lw=1.6,
                    color=("#2F6B2F" if _in_win else "#444444"), alpha=0.9, zorder=1)
        _ax.annotate(f"z*={_zs}", (_zs, _zlab_y.get(_zs, -0.65)), ha="center", va="top",
                     fontsize=10, color=("#2F6B2F" if _in_win else "#444444"),
                     fontweight="bold")

    _ax.set_yticks(np.arange(len(_rows)))
    _ax.set_yticklabels([])
    _ax.set_ylim(-1.6, len(_rows) - 0.2)
    _ax.set_xlim(0, max(_r["thr"] for _r in _rows) * 1.15)
    _ax.set_xlabel("activation threshold  t = -z(mean) of dominant hub   "
                   "(library flagged coverage-failing once z* > t)",
                   fontsize=11.5, color="black")
    _ax.set_title("Coverage-axis z* sensitivity - a library leaves the origin when z* > -z(mean)\n"
                  "0.9 / 1.0 exclude adenine_hybrid (a real greedy-failure); 1.5 is in the valid window",
                  fontsize=12.5, fontweight="bold", color="black")
    _ax.tick_params(colors="black", labelsize=11)
    _ax.grid(axis="x", alpha=0.3, color="grey")
    for _sp in _ax.spines.values():
        _sp.set_color("black")

    # Legend: colour = OPERATIVE failure mode (red coverage / blue allocation / orange
    # imbalance / grey none). A far-right point that is NOT grey is coherent: its hub is
    # coverage-reachable, but the library is hard via allocation (amide_hybrid) or imbalance
    # (thrombin). Only x (coverage) is on this axis; allocation/imbalance live on the landscape.
    _leg_handles = [
        _Line2D([0], [0], marker="o", ls="", ms=11, mfc=_C_COV, mec="black",
                label="coverage failure — mean-blind hub (this axis)"),
        _Line2D([0], [0], marker="o", ls="", ms=11, mfc=_C_ALLOC, mec="black",
                label="allocation — partner-specific / high-σ hub (landscape y)"),
        _Line2D([0], [0], marker="o", ls="", ms=11, mfc=_C_IMB, mec="black",
                label="imbalance — under-observed large component (landscape edge)"),
        _Line2D([0], [0], marker="o", ls="", ms=11, mfc=_C_EASY, mec="black",
                label="no failure mode — greedy-competitive"),
        _Patch(facecolor="#9CC16A", alpha=0.45, edgecolor="none",
               label=f"valid z* window ({_win_lo:.2f}–{_win_hi:.2f})"),
        _Line2D([0], [0], ls="--", lw=1.6, color="#444444", label="candidate z* cut"),
    ]
    _ax.legend(handles=_leg_handles, loc="lower right", fontsize=8.5, framealpha=0.96,
               edgecolor="black", borderpad=0.7, labelspacing=0.6)
    _fig.tight_layout()

    _flag_09 = ", ".join(_r["lab"] for _r in _rows if 0.9 + _r["zmean"] > 0) or "(none)"
    _flag_10 = ", ".join(_r["lab"] for _r in _rows if 1.0 + _r["zmean"] > 0) or "(none)"
    _flag_15 = ", ".join(_r["lab"] for _r in _rows if 1.5 + _r["zmean"] > 0) or "(none)"
    _md = mo.md(
        f"""
        **Does swapping z* change the figure?** The *rank order* of flagged libraries is
        invariant (rank-invariant across z* in [1.0, 2.0]: **{_rank_invariant}**) — that is
        what stayed fixed when z* was swept. But the *set* of off-origin libraries grows with
        z*, and that moves points:

        | z* | libraries flagged coverage-failing (off origin) |
        |----|--------------------------------------------------|
        | 0.9 | {_flag_09} |
        | 1.0 | {_flag_10} |
        | 1.5 | {_flag_15} |

        At **0.9 / 1.0** only adenine is flagged — **adenine_hybrid collapses onto the 'easy'
        origin** despite greedy failing it (greedy ≈ 30 vs TACTICS ≈ 95). The window that
        catches both coverage libraries yet excludes the mean-separable ones is
        **z* ∈ ({_win_lo:.2f}, {_win_hi:.2f})**; 1.5 sits inside it, 0.9 / 1.0 do not.

        **Why are some far-right points still hard?** This axis is *coverage only*, and colour
        marks the **operative** mode. A mean-reachable hub (high −z(mean), far-right) means
        greedy is not starved on coverage — but the library can still be hard via the other two
        modes: **amide_hybrid** (blue) is allocation-hard (partner-specific / high-σ hub);
        **thrombin** (orange) is imbalance-hard (30:1 → noisy estimates, Legacy RWS breaks).
        Both sit far-right yet greedy-hard, and TACTICS wins on hub uncertainty / GMIC routing.
        **Post-hoc statement:** greedy underperforms when the dominant hub triggers **any one —
        or a combination — of the three modes** (coverage, allocation, imbalance); none
        triggered ⇒ all methods converge (quinazoline). Outcome = the *combination*, never one
        coordinate.
        """
    )

    mo.vstack([
        plt.gcf(),
        _md,
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename="coverage_zstar_sensitivity.pdf",
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


@app.cell
def _(fig_to_pdf_bytes, mo, np, plt, project_root):
    """Docking reagent score landscape: the greedy-failure modes behind a TACTICS advantage.

    POST-HOC ground-truth characterization, one point per docking library, placed by how
    much it suffers each of the THREE greedy-failure modes that make TACTICS beat greedy
    (a point at the origin with a thin edge has none -> easy to screen). The COMBINATION of
    modes determines WHICH method wins (fill colour).
      x = COVERAGE-failure severity = max(0, 1.5 + z(mean) of the dominant hub).
          0 = the hit-holding hub is greedy-reachable (mean-separated past the reachability
          reference z* = -1.5, ~ top 7% by mean, Phi(-1.5)=0.067); large = a mean-blind hub
          holds the hits and greedy starves it (adenine). The max() is the severity floor and
          z* sets where it bites -- the clamp and the constant are ONE mechanism, not two
          (max(0, z(mean)) with no offset is degenerate: every hit-hub has z(mean) < 0).
      y = ALLOCATION-failure severity = z(std) × (1 − top-1 share).
          0 = greedy can cover the hit mass; large = the hits are partner-specific AND
          dispersed, so greedy's budget can't cover the spread (amide_hybrid).
      EDGE (orange, width) = reagent-count IMBALANCE = max/min component reagent count --
          a search-dynamics mode: a heavily imbalanced library (thrombin 30:1) breaks Legacy
          RWS's symmetric thermal cycling, so TACTICS wins even with NO hub failure.
      FILL = which method wins (the combination of failure modes determines it):
        GREY easy/method-agnostic (Δ_TL<1); PURPLE TACTICS-either; GREEN TT-TS; BLUE RWS.
    Reading: x away from origin = coverage failure (adenine); up = allocation failure
    (amide_hybrid); thick edge at origin = imbalance (thrombin); origin + thin edge = easy
    (quinazoline). FRED→HYBRID arrows = reagent-score shift, same chemistry. Descriptive
    structural summary (n=8, constructed thresholds), NOT a prospective classifier.
    """
    import json as _json
    from matplotlib.lines import Line2D as _Line2D

    _PROBE_DIR = project_root / "outputs" / "debug" / "docking_landscape_probe"
    _MAP_LIBS = [
        "thrombin", "thrombin_hybrid",
        "adenine", "adenine_hybrid",
        "amide", "amide_hybrid",
        "quinazoline", "quinazoline_hybrid",
    ]
    _TT = "TACTICS Enhanced-TT-TS (GMIC)"
    _RWS = "TACTICS Enhanced-RWS (GMIC)"
    _LRWS = "Legacy Enhanced-RWS"
    _LGRD = "Legacy Standard-Greedy"
    _TIE_DELTA = 1.0  # delta_TL below this => easy / method-agnostic (grey)
    _GAP_DELTA = 1.0  # |TT-TS - RWS| below this => TACTICS-either (purple)
    _C_TIE, _C_TT, _C_RWS = "#343838", "#349E41", "#4477AA"
    _C_EITHER = "#471A6E"  # purple: TACTICS-either
    _CAT_COLOR = {"easy": _C_TIE, "either": _C_EITHER, "ttts": _C_TT, "rws": _C_RWS}
    _CAT_LABEL = {
        "easy": "easy / method-agnostic  (Δ$_{T-L}$ < 1)",
        "either": "TACTICS (either)  (|TT-TS − RWS| < 1)",
        "ttts": "TT-TS wins",
        "rws": "RWS wins",
    }

    def _categorize(_delta_tl, _gap_tr):
        if _delta_tl < _TIE_DELTA:
            return "easy"
        if abs(_gap_tr) < _GAP_DELTA:
            return "either"
        return "ttts" if _gap_tr > 0 else "rws"

    _EDGE = "#E8820E"  # imbalance edge colour

    def _imb_lw(_ratio):
        # edge width grows with log imbalance ratio (2:1 thin -> 30:1 thick)
        _t = (np.log10(_ratio) - np.log10(2.0)) / (np.log10(30.0) - np.log10(2.0))
        return 1.0 + 3.2 * float(np.clip(_t, 0.0, 1.0))

    _rows = {}
    for _lab in _MAP_LIBS:
        with open(_PROBE_DIR / f"{_lab}.json") as _fh:
            _j = _json.load(_fh)
        _h = _j["dominant_hub"]
        _cidx = _h["component_idx"]
        _m2 = _j["M2_concentration_by_topn"]["100"][f"comp_{_cidx}"]
        _pm = _j["M6_search"]["per_method"]
        _tt = _pm[_TT]["topn_recovered_mean"]
        _rws = _pm[_RWS]["topn_recovered_mean"]
        _delta = max(_tt, _rws) - max(_pm[_LRWS]["topn_recovered_mean"],
                                      _pm[_LGRD]["topn_recovered_mean"])
        _gap = _tt - _rws
        _counts = list(_j["n_reagents"].values())
        _rows[_lab] = {
            "cov": max(0.0, 1.5 + _h["z_mean"]),               # coverage-failure severity
            "alloc": _h["z_std"] * (1.0 - _m2["top1_share"]),  # allocation-failure severity
            "imb": max(_counts) / min(_counts),                # reagent-count imbalance
            "category": _categorize(_delta, _gap),
        }

    # Reserve a right-hand gutter for the legends so they never collide with the
    # plot region.
    _fig, _ax = plt.subplots(figsize=(13.6, 8.5), facecolor="white")
    _ax.set_facecolor("white")

    _xs = np.array([_r["cov"] for _r in _rows.values()])
    _ys = np.array([_r["alloc"] for _r in _rows.values()])
    _x_hi = float(_xs.max())
    _y_hi = float(_ys.max())
    _ax.set_xlim(-0.10 * _x_hi, 1.22 * _x_hi)
    _ax.set_ylim(-0.06 * _y_hi, 1.18 * _y_hi)

    # Faint failure-mode region guides (axes are severities: origin = no failure).
    _ax.axhline(0, color="#CCCCCC", lw=0.8, zorder=0)
    _ax.axvline(0, color="#CCCCCC", lw=0.8, zorder=0)
    _ax.annotate("easy\n(reachable, coverable, balanced)", (0.02 * _x_hi, 0.45),
                 fontsize=9.5, color="#888888", style="italic", ha="left",
                 va="bottom", zorder=1)
    _ax.annotate("COVERAGE failure →\n(mean-blind hub)", (0.50 * _x_hi, 0.28),
                 fontsize=10, color="#B23030", style="italic", ha="left",
                 va="bottom", zorder=1)
    _ax.annotate("↑ ALLOCATION failure\n(partner-specific + dispersed)",
                 (0.0, 0.97 * _y_hi),
                 fontsize=10, color="#2A5DA0", style="italic", ha="left",
                 va="top", zorder=1)

    # Per-label offsets (several libraries cluster at cov=0; spread their labels).
    _label_off = {
        "amide_hybrid": (13, 0, "left", "center"),
        "thrombin": (13, 11, "left", "bottom"),
        "quinazoline": (13, 1, "left", "center"),
        "quinazoline_hybrid": (13, -12, "left", "top"),
        "amide": (13, -1, "left", "center"),
        "thrombin_hybrid": (10, -11, "left", "top"),
        "adenine": (13, 2, "left", "center"),
        "adenine_hybrid": (13, 9, "left", "bottom"),
    }
    for _lab, _r in _rows.items():
        _ax.scatter(_r["cov"], _r["alloc"], s=320,
                    c=_CAT_COLOR[_r["category"]], edgecolors=_EDGE,
                    linewidths=_imb_lw(_r["imb"]), alpha=0.95, zorder=3)
        _dx, _dy, _hha, _vva = _label_off.get(_lab, (10, 7, "left", "bottom"))
        _ax.annotate(_lab, (_r["cov"], _r["alloc"]), fontsize=10.5, color="black",
                     xytext=(_dx, _dy), textcoords="offset points", ha=_hha,
                     va=_vva, zorder=4)

    # FRED -> HYBRID reagent-score-shift arrows (neutral SHIFT in the ground-truth
    # reagent scores; same chemistry, NOT a quality judgment about a scoring function).
    _pairs = [("amide", "amide_hybrid"), ("thrombin", "thrombin_hybrid"),
              ("quinazoline", "quinazoline_hybrid"), ("adenine", "adenine_hybrid")]
    for _a, _b in _pairs:
        if _a in _rows and _b in _rows:
            _ax.annotate("", xy=(_rows[_b]["cov"], _rows[_b]["alloc"]),
                         xytext=(_rows[_a]["cov"], _rows[_a]["alloc"]),
                         arrowprops=dict(arrowstyle="->", color="#999999", lw=1.1,
                                         alpha=0.5, connectionstyle="arc3,rad=0.18"),
                         zorder=2)

    _ax.set_xlabel("COVERAGE-failure severity   max(0, 1.5 + z(mean))\n"
                   "(0 = greedy-reachable hub past z*=−1.5;  large = mean-blind hidden hub)",
                   fontsize=12.5, color="black")
    _ax.set_ylabel("ALLOCATION-failure severity   z(std) × (1 − top-1 share)\n"
                   "(0 = greedy can cover;  large = partner-specific + dispersed)",
                   fontsize=12, color="black")
    _ax.set_title("Docking reagent score landscape — three greedy-failure modes behind a TACTICS advantage:\n"
                  "coverage (x)  ×  allocation (y)  ×  imbalance (edge width)",
                  fontsize=13.5, fontweight="bold", color="black")
    _ax.tick_params(colors="black", labelsize=12)
    _ax.grid(alpha=0.3, color="grey", zorder=1)
    for _sp in _ax.spines.values():
        _sp.set_color("black")

    # Outcome (fill) + FRED->HYBRID legend, upper-right gutter.
    _cat_handles = [
        _Line2D([0], [0], marker="o", ls="", ms=11, mfc=_CAT_COLOR[_c],
                mec="#888888", label=_CAT_LABEL[_c])
        for _c in ("easy", "either", "ttts", "rws")
    ]
    _cat_handles.append(
        _Line2D([0], [0], color="#999999", lw=1.3,
                label="FRED → HYBRID (reagent-score\nshift, same chemistry)")
    )
    _leg1 = _ax.legend(handles=_cat_handles, loc="upper left",
                       bbox_to_anchor=(1.02, 1.0), fontsize=10, framealpha=0.95,
                       edgecolor="black", title="Benchmark outcome (fill)",
                       title_fontsize=11, borderaxespad=0.0)
    _ax.add_artist(_leg1)

    _imb_refs = [2, 10, 30]
    _imb_handles = [
        _Line2D([0], [0], marker="o", ls="", mfc="#DDDDDD", mec=_EDGE,
                mew=_imb_lw(_r), markersize=14, label=f"{_r}:1")
        for _r in _imb_refs
    ]
    _leg2 = _ax.legend(handles=_imb_handles, loc="lower left",
                       bbox_to_anchor=(1.02, 0.0), fontsize=10, framealpha=0.95,
                       edgecolor="black", labelspacing=1.5, borderpad=1.0,
                       borderaxespad=0.0,
                       title="Reagent imbalance (edge)\nmax / min reagent count",
                       title_fontsize=10)

    # Reserve the right gutter via subplots_adjust (NOT tight_layout(rect=...),
    # which clips out-of-axes legends).
    _fig.subplots_adjust(left=0.07, right=0.74, top=0.90, bottom=0.10)

    # Display a TIGHT-bbox PNG (legends included as extra artists) so the two
    # right-gutter legends are never clipped by marimo's canvas crop at this
    # figure size — plt.gcf() crops to the figure box and cut them off.
    import io as _io
    _png = _io.BytesIO()
    _fig.savefig(_png, format="png", dpi=150, bbox_inches="tight",
                 bbox_extra_artists=[_leg1, _leg2])
    mo.vstack([
        mo.image(_png.getvalue()),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename="reagent_score_landscape_map.pdf",
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


@app.cell
def _(fig_to_pdf_bytes, mo, np, pl, plt, project_root):
    """Greedy reach mechanism — coverage failure IS greedy starving the mean-blind hub.

    Companion to the coverage axis, measured from the search shards (Legacy Greedy, search
    phase only) vs ground-truth reagent z(mean) -- NO recovery/outcome is read, so it is a
    non-circular mechanism measurement. enrichment = greedy's selection rate of a reagent
    divided by uniform 1/N_component; >1 = greedy over-samples (reaches) it, <1 = starves it.
    Greedy's enrichment is monotone in z(mean) and its over/under crossover is z(mean) ~ -0.9.
    The adenine hub (amidine_003, z(mean) -0.58, holds 49/100 hits) sits at ~0.06x uniform --
    greedy picks it ~18x LESS than blind chance. The coverage reference z*=-1.5 sits deep in
    greedy's strong-exploitation region (~2.5x uniform), i.e. where greedy concentrates enough
    budget to RECOVER a hub -- which is why -1.5 is the recovery-relevant frontier rather than
    the bare over/under crossover at -0.9.
    Source: outputs/debug/greedy_reach_calibration/ (experiments/_debug/greedy_reach_calibration.py).
    """
    import json as _json

    _CAL = project_root / "outputs" / "debug" / "greedy_reach_calibration"
    mo.stop(
        not (_CAL / "reagents.csv").exists(),
        mo.md("Run `python experiments/_debug/greedy_reach_calibration.py` first "
              "(writes `outputs/debug/greedy_reach_calibration/`)."),
    )

    _reag = pl.read_csv(_CAL / "reagents.csv")
    with open(_CAL / "summary.json") as _fh:
        _summary = _json.load(_fh)
    _hubs = _summary["hubs"]
    _binned = _summary["binned_enrichment"]
    _cross = _summary["pooled_crossover_mean"]

    _floor = 0.02
    _z = _reag["z_mean"].to_numpy()
    _Ec = np.clip(_reag["enrichment"].to_numpy(), _floor, None)

    _fig, _ax = plt.subplots(figsize=(11.5, 6.6), facecolor="white")
    _ax.set_facecolor("white")
    _ax.set_yscale("log")
    _ax.scatter(_z, _Ec, s=6, color="#BBBBBB", alpha=0.25, zorder=1,
                label="reagents (all components, 8 libraries)")
    _bz = [_b["z_center"] for _b in _binned]
    _bE = np.clip([_b["mean_E"] for _b in _binned], _floor, None)
    _ax.plot(_bz, _bE, color="#222222", lw=2.0, zorder=3, label="binned mean enrichment")

    _ax.axhline(1.0, color="#888888", ls=":", lw=1.4, zorder=2)
    _ax.annotate("enrichment = 1  (reach / starve boundary)", (2.7, 1.15), fontsize=9,
                 color="#666666", va="bottom", ha="right")
    _ax.axvline(_cross, color="#E8820E", ls="--", lw=1.6, zorder=2)
    _ax.annotate(f"greedy crossover\nz(mean) ≈ {_cross:.1f}", (_cross, _floor * 1.5),
                 fontsize=9, color="#B5650E", va="bottom", ha="center")
    _ax.axvline(-1.5, color="#B23030", ls="-", lw=1.6, alpha=0.85, zorder=2)
    _ax.annotate("coverage\nreference z*=−1.5", (-1.5, 28.0), fontsize=9,
                 color="#B23030", va="top", ha="center")

    _hub_off = {  # fan out the labels of the hubs that cluster near z(mean) ~ -1.4..-1.9
        "quinazoline_hybrid": (-6, 11, "right"),
        "quinazoline": (4, -13, "left"),
        "thrombin_hybrid": (8, 7, "left"),
    }
    for _h in _hubs:
        _hy = max(_h["hub_enrichment"], _floor)
        _ax.scatter(_h["hub_z_mean"], _hy, s=150, color="#349E41", edgecolors="black",
                    linewidths=0.9, zorder=5)
        _dx, _dy, _hha = _hub_off.get(_h["library"], (6, 4, "left"))
        _ax.annotate(_h["library"], (_h["hub_z_mean"], _hy), xytext=(_dx, _dy),
                     textcoords="offset points", fontsize=8.5, color="black",
                     ha=_hha, zorder=6)

    _ax.set_xlabel("ground-truth reagent z(mean)   "
                   "(more negative = better mean; docking minimizes)",
                   fontsize=11.5, color="black")
    _ax.set_ylabel("greedy selection enrichment\n(selection rate ÷ uniform 1/N)",
                   fontsize=11.5, color="black")
    _ax.set_title("Greedy reach: enrichment is monotone in hub mean — greedy starves the mean-blind hub\n"
                  "(dominant hubs in green; adenine sits at ~0.06× uniform — actively avoided)",
                  fontsize=12.5, fontweight="bold", color="black")
    _ax.tick_params(colors="black", labelsize=11)
    _ax.grid(alpha=0.3, color="grey", zorder=0)
    for _sp in _ax.spines.values():
        _sp.set_color("black")
    _ax.legend(loc="upper right", fontsize=9, framealpha=0.95, edgecolor="black")
    _fig.tight_layout()

    mo.vstack([
        plt.gcf(),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename="greedy_reach_mechanism.pdf",
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


@app.cell
def _(fig_to_pdf_bytes, mo, np, plt, project_root):
    """Ranking heatmap: top-100 recovery for 5 methods × 8 docking landscapes.

    Cell = M6 topn_recovered_mean (annotated). Columns ordered by the
    best-TACTICS minus best-Legacy gap (largest left), so the libraries where
    TACTICS helps most appear first (ordering not stated in the title -- the
    manuscript text covers it). Single heatmap, no rank strip.
    """
    import json as _json

    _PROBE_DIR = project_root / "outputs" / "debug" / "docking_landscape_probe"
    _MAP_LIBS = [
        "thrombin", "thrombin_hybrid",
        "adenine", "adenine_hybrid",
        "amide", "amide_hybrid",
        "quinazoline", "quinazoline_hybrid",
    ]
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

    _recov, _delta = {}, {}
    for _lab in _MAP_LIBS:
        with open(_PROBE_DIR / f"{_lab}.json") as _fh:
            _j = _json.load(_fh)
        _pm = _j["M6_search"]["per_method"]
        _recov[_lab] = {_m: _pm[_m]["topn_recovered_mean"] for _m in _METHODS}
        _bt = max(_pm[_TT]["topn_recovered_mean"], _pm[_RWS]["topn_recovered_mean"])
        _bl = max(_pm[_LRWS]["topn_recovered_mean"], _pm[_LGRD]["topn_recovered_mean"])
        _delta[_lab] = _bt - _bl

    _libs = sorted(_MAP_LIBS, key=lambda l: _delta[l], reverse=True)
    _mat = np.array([[_recov[l][_m] for l in _libs] for _m in _METHODS])

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

    # Footnote: quinazoline's within-noise TT-TS caveat (data unchanged).
    _fig.text(0.5, -0.06,
              "quinazoline is method-agnostic (mean-reachable landscape); its "
              "TT-TS value reflects a single 2% replicate set and is within-noise "
              "of the others.",
              ha="center", va="top", fontsize=9, color="#444444", style="italic",
              wrap=True)

    mo.vstack([
        plt.gcf(),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename="sar_ranking_heatmap.pdf",
            mimetype="application/pdf",
            label="Download as PDF",
        ),
    ])
    return


if __name__ == "__main__":
    app.run()
