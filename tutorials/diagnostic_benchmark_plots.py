"""
TACTICS Diagnostic Benchmark Plots

Visualises the per-query GMIC mechanism from the diagnostic benchmark
(110 trials, 5 libraries, 11 queries, 5 replicates).

Figure types:
  - **RWS Diagnostic**: GMIC trajectory + effective temperature + divergence gate
  - **TT-TS Diagnostic**: disagreement EMA + heated_scale adaptation
  - **Layers 1 & 2**: GMIC component-routing + within-component intensity,
    with per-library dropdowns for docking and ROCS landscapes

Run as app:  marimo run tutorials/diagnostic_benchmark_plots.py
Edit mode:   marimo edit tutorials/diagnostic_benchmark_plots.py
"""

import marimo

__generated_with = "0.21.1"
app = marimo.App(width="full", app_title="Diagnostic Benchmark Plots")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Diagnostic Plots and Mechanism Plots
    This notebook is designed to showcase the inner workings of TACTICS methods.
    """)
    return


@app.cell(hide_code=True)
def _():
    import marimo as mo
    import sys
    from pathlib import Path

    import matplotlib
    matplotlib.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": "black",
        "axes.labelcolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
        "text.color": "black",
        "savefig.facecolor": "white",
    })
    import matplotlib.pyplot as plt
    import polars as pl

    try:
        project_root = Path(__file__).parent.parent.resolve()
    except NameError:
        project_root = Path.cwd()

    sys.path.insert(0, str(project_root / "src"))

    from TACTICS.library_analysis.diagnostic_plots import (
        plot_reagent_usage_action_panel,
        plot_rws_diagnostic,
        plot_ttts_diagnostic,
    )

    diag_dir = project_root / "outputs" / "diagnostic_benchmark" / "diagnostics"
    search_dir = project_root / "outputs" / "diagnostic_benchmark" / "search"
    scores_dir = project_root / "data" / "scores"
    reagents_dir = project_root / "data" / "reagents"
    plot_dir = project_root / "plots"
    plot_dir.mkdir(exist_ok=True)
    return (
        diag_dir,
        mo,
        pl,
        plot_dir,
        plot_reagent_usage_action_panel,
        plot_rws_diagnostic,
        plot_ttts_diagnostic,
        plt,
        reagents_dir,
        search_dir,
    )


@app.cell
def _(diag_dir, mo, pl):
    """Load diagnostic data and build metadata index."""
    bench_data = {"libraries": {}, "meta": {}}

    for _lib_dir in sorted(diag_dir.iterdir()):
        if not _lib_dir.is_dir():
            continue
        for _f in _lib_dir.iterdir():
            if _f.suffix == ".parquet":
                bench_data["libraries"][_lib_dir.name] = pl.read_parquet(_f)
                break

    bench_data["meta"] = {
        "rxn206":      {"gain": 6.1,  "scoring": "ROCS",    "n_comps": 2},
        "mannich":     {"gain": 9.1,  "scoring": "ROCS",    "n_comps": 3},
        "orru":        {"gain": 1.7,  "scoring": "ROCS",    "n_comps": 3},
        "adenine":     {"gain": 64.6, "scoring": "docking", "n_comps": 3},
        "quinazoline": {"gain": -0.9, "scoring": "docking", "n_comps": 3},
    }

    mo.md(
        f"**Loaded {len(bench_data['libraries'])} libraries**: "
        f"{', '.join(bench_data['libraries'].keys())}\n\n"
        + "\n".join(
            f"- **{_name}**: {_df.shape[0]} rows, "
            f"{_df['component_idx'].n_unique()} components, "
            f"{_df['query_id'].n_unique()} queries, "
            f"{_df['replicate'].n_unique()} replicates"
            for _name, _df in bench_data["libraries"].items()
        )
    )
    return (bench_data,)


@app.cell
def _(mo):
    mo.md("""
    # RWS Diagnostic: GMIC + Temperature + Divergence Gate
    """)
    return


@app.cell
def _(bench_data, plot_dir, plot_rws_diagnostic, plt):
    for _lib, _df in bench_data["libraries"].items():
        _rws = _df.filter(_df["method"].str.contains("RWS"))
        if len(_rws) == 0:
            continue
        for _q in _rws["query_id"].unique().sort().to_list():
            _q_data = _rws.filter(_rws["query_id"] == _q)
            _fig = plot_rws_diagnostic(_q_data, title=f"{_lib} — {_q} (RWS)")
            _fig.savefig(plot_dir / f"rws_{_lib}_{_q}.pdf", bbox_inches="tight")
            plt.show()
            plt.close(_fig)
    return


@app.cell
def _(mo):
    mo.md("""
    # TT-TS Diagnostic: Disagreement EMA + Heated Scale
    """)
    return


@app.cell
def _(bench_data, plot_dir, plot_ttts_diagnostic, plt):
    for _lib, _df in bench_data["libraries"].items():
        _ttts = _df.filter(_df["method"].str.contains("TT-TS"))
        if len(_ttts) == 0:
            continue
        for _q in _ttts["query_id"].unique().sort().to_list():
            _q_data = _ttts.filter(_ttts["query_id"] == _q)
            _fig = plot_ttts_diagnostic(_q_data, title=f"{_lib} — {_q} (TT-TS)")
            _fig.savefig(plot_dir / f"ttts_{_lib}_{_q}.pdf", bbox_inches="tight")
            plt.show()
            plt.close(_fig)
    return


@app.cell
def _(mo):
    mo.md("""
    # Action Panel: Per-Batch Reagent Diversity + Heated Component

    Shows the GMIC/TT-TS rotation mechanism in action batch by batch.

    - **Grouped bars**: one bar per component per batch.  Bar height
      = fraction of that component's reagents sampled in the batch.
      Each bar is split into a *revisit* segment (bottom, dimmer) and
      a *new* segment (top, brighter) — "new" means first-time-seen
      this search.
    - **Heated encoding**: the component marked ``is_heated`` that
      batch is rendered at full opacity; cooled components' bars are
      muted.  Mechanism working → heated bars stand out batch by
      batch.
    """)
    return


@app.cell
def _(
    bench_data,
    pl,
    plot_dir,
    plot_reagent_usage_action_panel,
    plt,
    reagents_dir,
    search_dir,
):
    _FOLDER_MAP = {
        "amide-suzuki": "amide_suzuki",
        "groebke-blackburn-bienayme": "groebke",
    }
    _WINDOW = {
        "mannich": 1, "orru": 1, "rxn206": 1,
        "adenine": 1, "quinazoline": 1,
    }
    _FIGSIZE = {
        "mannich": (16, 5.5), "orru": (16, 5.5), "rxn206": (16, 5.5),
        "adenine": (16, 5.5), "quinazoline": (16, 5.5),
    }
    # Imbalanced libraries — normalize per batch (unique / min(batch, n_reagents))
    # so bar heights are comparable across differently-sized components.
    _NORMALIZE = {
        "mannich": "library", "orru": "library", "rxn206": "library",
        "adenine": "batch", "quinazoline": "batch",
    }
    # Component labels derived from reaction SMARTS in data/reactions.csv.
    # Order matches reagent_0, reagent_1, ... Adenine uses auto-inferred prefix
    # names (amidine/isocyanide/aldehyde) and is left to the default path.
    _COMPONENT_NAMES = {
        "mannich":     ["Amine", "Aldehyde", "Ketone"],
        "orru":        ["Aldehyde", "Cyanoester", "Amine"],
        "rxn206":      ["Alcohol", "Sulfonamide"],
        "quinazoline": ["Aminobenzoic acid", "Amine", "Carboxylic acid"],
    }

    # Quinazoline has too many batches (20M library, ~1000+ batches) for the
    # action-panel diagnostic to be legible — skip it here.
    _SKIP_LIBS = {"quinazoline"}

    for _lib, _df in bench_data["libraries"].items():
        if _lib in _SKIP_LIBS:
            continue
        _folder = _FOLDER_MAP.get(_lib, _lib)
        _rfiles = sorted((reagents_dir / _folder).glob(f"{_folder}_reagent_*.smi"))
        if not _rfiles:
            _rfiles = sorted((reagents_dir / _folder).glob(f"{_lib}_reagent_*.smi"))
        if not _rfiles:
            continue

        _search_path = search_dir / _lib / f"{_lib}.parquet"
        if not _search_path.exists():
            continue
        _search_df = pl.read_parquet(_search_path)
        _win = _WINDOW.get(_lib, 1)

        for _method_pat, _method_tag in [("RWS", "rws"), ("TT-TS", "ttts")]:
            _m_df = _df.filter(_df["method"].str.contains(_method_pat))
            if len(_m_df) == 0:
                continue
            for _q in _m_df["query_id"].unique().sort().to_list():
                _q_diag = _m_df.filter(
                    (_m_df["query_id"] == _q) & (pl.col("replicate") == 0)
                )
                _q_search = _search_df.filter(
                    (pl.col("query_id") == _q)
                    & (pl.col("method").str.contains(_method_pat))
                    & (pl.col("replicate") == 0)
                )
                if len(_q_search) == 0 or len(_q_diag) == 0:
                    continue

                _heated_label = "Explore" if _method_pat == "TT-TS" else "Heated"
                _fig = plot_reagent_usage_action_panel(
                    _q_search, _q_diag,
                    reagent_files=_rfiles,
                    batch_size=100,
                    aggregate_window=_win,
                    normalize=_NORMALIZE.get(_lib, "library"),
                    figsize=_FIGSIZE.get(_lib, (16, 5)),
                    component_names=_COMPONENT_NAMES.get(_lib),
                    title=f"{_lib} — {_q} — {_method_pat} (rep 0)",
                    heated_label=_heated_label,
                )
                _fig.savefig(
                    plot_dir / f"action_{_lib}_{_q}_{_method_tag}.pdf",
                    bbox_inches="tight",
                )
                plt.show()
                plt.close(_fig)
    return


@app.cell
def _():
    """Shared deps for the MECH cells: numpy, a fig→PDF-bytes helper (this
    notebook's original import cell defines neither), the canonical-benchmark
    paths, and the heterogeneous-schema diagnostics loader.

    These are kept separate from the original import cell so none of the
    existing cells are modified.
    """
    import io as _io
    from pathlib import Path as _Path

    import numpy as np

    def fig_to_pdf_bytes(fig):
        """Convert a matplotlib figure to PDF bytes in memory (tight bbox)."""
        _buf = _io.BytesIO()
        fig.savefig(_buf, format="pdf", bbox_inches="tight", dpi=300)
        return _buf.getvalue()

    try:
        _root = _Path(__file__).parent.parent.resolve()
    except NameError:
        _root = _Path.cwd()

    # Diagnostics: the gate-fixed 2% TT-TS GMIC/disagreement trajectories live in
    # full_benchmark_nogate (merged per-library, all 28 libs incl. hybrids).
    # ema_cats' adenine TT-TS diagnostics are STALE (pre-merge, gated) and were
    # never regenerated. Recovery stays on the canonical ema_cats: its TT-TS is
    # the spliced nogate TT-TS (identical), only non-TT-TS RNG differs.
    mech_bench = _root / "outputs" / "full_benchmark_nogate"
    mech_diag_dir = mech_bench / "diagnostics"
    mech_recovery_path = (
        _root / "outputs" / "full_benchmark_ema_cats" / "analysis" / "recovery_summary.parquet"
    )
    return fig_to_pdf_bytes, mech_diag_dir, np


@app.cell
def _(mo):
    mo.md(r"""
    # Layers 1 & 2 — Exploration Mechanism for Any Docking Library

    Pick a docking library from the dropdown; **both** mechanism layers below
    update to that library.

    - **Layer 1 — *which* component?**  Signal → Decision → Result.  GMIC reads
      each component as resolved (high GMIC → exploit) or unresolved (low GMIC →
      explore); the *explored* strip follows GMIC, and **cumulative coverage**
      (running total of reagents sampled) shows the consequence.
    - **Layer 2 — *how hard* within that component?**  Intensity → Decision →
      Result, one level down.  Both knobs are **dimensionless ×-multipliers
      centred on 1.0** (1.0 = no adjustment; >1 = amplify exploration): **RWS**
      tunes the **CATS temperature multiplier** (`cats_multiplier`, ×base
      temperature, GMIC-driven); **TT-TS** tunes the **σ-inflation multiplier**
      (`heated_scale`, ×posterior σ, disagreement-driven).  The Result row is
      **new reagents per batch** — the *per-batch increment* of coverage (it falls
      to zero as a component's reagent pool is exhausted, even while coverage is
      still climbing).

    **Data source (per library, labelled on each figure):** the **2% canonical**
    benchmark (`full_benchmark_ema_cats`) — the manuscript operating point — is
    used for the 4 base docking libraries (thrombin, amide, adenine, quinazoline).
    The 4 docking *hybrids* have no 2% diagnostics saved (search only), so they
    currently fall back to the **5%** sweep (`full_benchmark_5pct`) — flagged in
    the title — until 2% diagnostics are regenerated.  Replicate 0.  Reagent files
    resolve via the benchmark folder convention (hybrids reuse the base library's
    reagents).
    """)
    return


@app.cell
def _(mo):
    """Dropdown selector over the 8 docking libraries."""
    _dock_libs = [
        "thrombin", "thrombin_hybrid", "amide", "amide_hybrid",
        "adenine", "adenine_hybrid", "quinazoline", "quinazoline_hybrid",
    ]
    dock_lib_dd = mo.ui.dropdown(
        options=_dock_libs, value="thrombin", label="Docking library"
    )
    dock_lib_dd
    return (dock_lib_dd,)


@app.cell
def _(dock_lib_dd, fig_to_pdf_bytes, mo, pl, plt, reagents_dir):
    """Render both layer plots for the selected docking library."""
    import glob as _glob3
    from pathlib import Path as _Path

    from TACTICS.library_analysis.diagnostic_plots import (
        plot_adaptive_intensity as _ai,
        plot_gmic_directed_exploration as _gde,
    )

    _lib = dock_lib_dd.value
    # Reagent-folder convention (Local Norm #11): hybrids reuse the base folder.
    _FOLDER = {
        "adenine_hybrid": "adenine", "thrombin_hybrid": "thrombin",
        "amide_hybrid": "amide", "quinazoline_hybrid": "quinazoline",
    }
    # Friendly component names where reagent IDs are numeric (else auto-inferred).
    _NAMES = {
        "thrombin": ["Acids", "Dipeptides"],
        "thrombin_hybrid": ["Acids", "Dipeptides"],
        "amide": ["Acids", "Amines"],
        "amide_hybrid": ["Acids", "Amines"],
        "quinazoline": ["Aminobenzoic acid", "Amine", "Carboxylic acid"],
        "quinazoline_hybrid": ["Aminobenzoic acid", "Amine", "Carboxylic acid"],
        "adenine": ["Amidine", "Isocyanide", "Aldehyde"],
        "adenine_hybrid": ["Amidine", "Isocyanide", "Aldehyde"],
    }
    _folder = _FOLDER.get(_lib, _lib)
    _key = _folder if _lib.endswith("_hybrid") else _lib
    _rf = sorted(
        _glob3.glob(str(reagents_dir / _folder / f"{_key}_reagent_*.smi"))
    )

    # Budget resolution: prefer any 2% source — the CANONICAL benchmark
    # (full_benchmark_ema_cats, base docking libs) or the regenerated hybrid
    # diagnostics (hybrid_diag_2pct, docking hybrids) — both the manuscript 2%
    # operating point.  Fall back to the 5% sweep only if no 2% diagnostics
    # exist yet for this library.
    _outputs = reagents_dir.parent.parent / "outputs"

    def _parqs(_d):
        # Glob .parquet but skip hidden/AppleDouble (._*) sidecar files.
        return [
            f for f in sorted(_glob3.glob(str(_d / "*.parquet")))
            if not _Path(f).name.startswith(".")
        ]

    def _shards(_bdir):
        return _parqs(_bdir / "diagnostics" / _lib)

    _candidates = [
        (_outputs / "full_benchmark_nogate", "2% gate-fixed"),
        (_outputs / "full_benchmark_ema_cats", "2% canonical (gated TT-TS — fallback)"),
        (_outputs / "hybrid_diag_2pct", "2% canonical"),
        (_outputs / "full_benchmark_5pct", "5% (2% diagnostics not yet generated)"),
    ]
    _bench, _budget = next(
        ((b, lbl) for b, lbl in _candidates if _shards(b)), (None, None)
    )
    _dfiles = _shards(_bench) if _bench else []
    _sfiles = _parqs(_bench / "search" / _lib) if _bench else []

    if not _dfiles or not _sfiles or not _rf:
        _out3 = mo.md(f"⚠️ diagnostics/search/reagents not found for **{_lib}**.")
    else:
        # ema_cats shards are per-method heterogeneous-schema; union-concat them
        # (null-fill missing columns) so RWS + TT-TS rows coexist.
        _cols = set()
        for _f in _dfiles:
            _cols |= set(pl.read_parquet_schema(_f).keys())
        _cols = sorted(_cols)
        _parts = []
        for _f in _dfiles:
            _d = pl.read_parquet(_f)
            for _c in _cols:
                if _c not in _d.columns:
                    _d = _d.with_columns(pl.lit(None).alias(_c))
            _parts.append(_d.select(_cols))
        _diag = pl.concat(_parts, how="vertical_relaxed")
        _search = pl.concat(
            [pl.read_parquet(_sf, columns=["Name", "method", "replicate", "phase"])
             for _sf in _sfiles],
            how="vertical_relaxed",
        )
        _cn = _NAMES.get(_lib)  # None -> auto-infer from reagent-ID prefixes
        _methods = ["TACTICS Enhanced-RWS (GMIC)", "TACTICS Enhanced-TT-TS (GMIC)"]
        _mlabels = {
            "TACTICS Enhanced-RWS (GMIC)": "TACTICS RWS",
            "TACTICS Enhanced-TT-TS (GMIC)": "TACTICS TT-TS",
        }

        _f1 = _gde(
            _search, _diag, _rf, methods=_methods, method_labels=_mlabels,
            component_names=_cn,
            title=f"Layer 1 — GMIC directs WHICH component to explore  ({_lib}, {_budget})",
        )
        _f2 = _ai(
            _search, _diag, _rf, methods=_methods, method_labels=_mlabels,
            component_names=_cn,
            title=(
                f"Layer 2 — within-component exploration intensity  ({_lib}, {_budget})\n"
                "both knobs are ×-multipliers centred on 1.0:  "
                "RWS = CATS temperature multiplier;  TT-TS = σ-inflation multiplier"
            ),
        )
        # Rasterise each figure explicitly.  Embedding live matplotlib Figure
        # objects in a single marimo cell is unreliable (only the last tends to
        # show), so we convert both to PNG, then close them.
        import io as _io3

        def _png(_f):
            _b = _io3.BytesIO()
            _f.savefig(_b, format="png", dpi=110, bbox_inches="tight")
            return _b.getvalue()

        _png1, _pdf1 = _png(_f1), fig_to_pdf_bytes(_f1)
        _png2, _pdf2 = _png(_f2), fig_to_pdf_bytes(_f2)
        plt.close(_f1)
        plt.close(_f2)
        _out3 = mo.vstack([
            mo.md(f"### Layer 1 — GMIC directs *which* component  ({_lib} · {_budget})"),
            mo.image(_png1),
            mo.download(
                data=_pdf1, filename=f"diag_gmic_directed_{_lib}.pdf",
                mimetype="application/pdf", label="Download Layer 1 PDF",
            ),
            mo.md(f"### Layer 2 — within-component intensity  ({_lib} · {_budget})"),
            mo.image(_png2),
            mo.download(
                data=_pdf2, filename=f"diag_adaptive_intensity_{_lib}.pdf",
                mimetype="application/pdf", label="Download Layer 2 PDF",
            ),
        ])
    _out3
    return


@app.cell
def _(mo):
    mo.md(r"""
    # Layers 1 & 2 — Exploration Mechanism for Any ROCS Library

    The same two mechanism layers as the docking section, for the ROCS
    (similarity) libraries.  On ROCS the GMIC *direction* is query-specific
    (different query molecules flag different components as flexible), so each
    library shows a **single representative query** — picked for the clearest
    sustained between-component GMIC separation with cross-replicate consensus
    = 1.0 (`experiments/_debug/pick_rocs_demo_queries.py`; `rxn206`/`groebke`
    reproduce the earlier manual demo-query picks).  Pick a library from the dropdown;
    replicate 0.

    - **Layer 1 — *which* component?**  GMIC reads each component as resolved
      (high GMIC → exploit) or unresolved (low GMIC → explore); the *explored*
      strip follows GMIC and cumulative coverage shows the consequence.
    - **Layer 2 — *how hard* within that component?**  RWS tunes the CATS
      temperature multiplier; TT-TS tunes the σ-inflation multiplier (both
      dimensionless ×-multipliers centred on 1.0).

    **Data source:** the 2% canonical benchmark (`full_benchmark_nogate`,
    gate-fixed; falls back to `full_benchmark_ema_cats`).
    """)
    return


@app.cell
def _(mo):
    """Dropdown selector over the 20 ROCS libraries."""
    _rocs_libs = [
        "amide-suzuki", "betti", "dobener", "groebke-blackburn-bienayme",
        "mannich", "niementowski", "orru", "passerini", "petasis", "poparov",
        "rxn101", "rxn102a", "rxn108b", "rxn111b", "rxn114b", "rxn203",
        "rxn205", "rxn206", "rxn207", "rxn208",
    ]
    rocs_lib_dd = mo.ui.dropdown(
        options=_rocs_libs, value="rxn206", label="ROCS library"
    )
    rocs_lib_dd
    return (rocs_lib_dd,)


@app.cell
def _(fig_to_pdf_bytes, mo, pl, plt, reagents_dir, rocs_lib_dd):
    """Render both layer plots for the selected ROCS library at its chosen
    representative query (clearest GMIC component separation; replicate 0)."""
    import glob as _glob4
    import io as _io4
    from pathlib import Path as _Path4

    from TACTICS.library_analysis.diagnostic_plots import (
        plot_adaptive_intensity as _ai_r,
        plot_gmic_directed_exploration as _gde_r,
    )

    # Representative query per library (experiments/_debug/pick_rocs_demo_queries.py):
    # clearest sustained between-component GMIC separation, cross-rep consensus 1.0.
    _RQUERY = {
        "amide-suzuki": "query_066", "betti": "query_098", "dobener": "query_018",
        "groebke-blackburn-bienayme": "query_035", "mannich": "query_035",
        "niementowski": "query_004", "orru": "query_103", "passerini": "query_037",
        "petasis": "query_020", "poparov": "query_012", "rxn101": "query_018",
        "rxn102a": "query_081", "rxn108b": "query_104", "rxn111b": "query_060",
        "rxn114b": "query_018", "rxn203": "query_036", "rxn205": "query_024",
        "rxn206": "query_047", "rxn207": "query_038", "rxn208": "query_024",
    }
    # Chemical reactant name per component (reagent_0, reagent_1, [reagent_2]),
    # derived from each reaction's SMARTS in data/reactions.csv. Order follows the
    # reagent-file order (reagent_0 = first SMARTS reactant). Edit freely — these
    # replace the default "Component {idx}" labels on the layer-plot components.
    _ROCS_COMPONENT_NAMES = {
        # 2-component (rxn*) libraries
        "rxn102a": ["Aryl halide", "Amine"],
        "rxn108b": ["Aryl/vinyl halide", "Alkyne"],
        "rxn111b": ["Alkene", "Aryl/vinyl halide"],
        "rxn114b": ["Amine", "Amine"],
        "rxn203": ["Amine", "Amine"],
        "rxn205": ["Alcohol", "Phenol"],
        "rxn206": ["Alcohol", "Sulfonamide"],
        "rxn207": ["Carbonyl", "Amine"],
        "rxn208": ["Aryl halide", "Amine"],
        # rxn101: SMARTS absent from data/reactions.csv -> falls back to "Component {idx}"
        # 3-component (named multicomponent) reactions
        "amide-suzuki": ["Boronic acid", "Aryl halide", "Carboxylic acid"],
        "betti": ["Naphthol", "Aldehyde", "Amine"],
        "dobener": ["Aniline", "Aldehyde", "Ketoester"],
        "groebke-blackburn-bienayme": ["Aminoazine", "Aldehyde", "Isocyanide"],
        "mannich": ["Secondary amine", "Aldehyde", "Ketone"],
        "niementowski": ["Anthranilic acid", "Amine", "Carboxylic acid"],
        "orru": ["Aldehyde", "Cyanoacetate", "Amine"],
        "passerini": ["Carboxylic acid", "Carbonyl", "Isocyanide"],
        "petasis": ["Amine", "Carbonyl", "Boronic acid"],
        "poparov": ["Aniline", "Alkene", "Aldehyde"],
    }
    # Reagent-folder convention (Local Norm #11): a couple of ROCS libs live in a
    # differently-named folder; the reagent-file prefix stays the lib_id.
    _FOLDER = {"amide-suzuki": "amide_suzuki",
               "groebke-blackburn-bienayme": "groebke"}

    _lib = rocs_lib_dd.value
    _rcn = _ROCS_COMPONENT_NAMES.get(_lib)  # chemical names; None -> "Component {idx}"
    _q = _RQUERY[_lib]
    _folder = _FOLDER.get(_lib, _lib)
    _rf = sorted(
        f for f in _glob4.glob(str(reagents_dir / _folder / f"{_lib}_reagent_*.smi"))
        if not _Path4(f).name.startswith(".")
    )
    _outputs = reagents_dir.parent.parent / "outputs"

    def _parqs(_d):
        # Glob .parquet but skip hidden/AppleDouble (._*) sidecar files.
        return [
            f for f in sorted(_glob4.glob(str(_d / "*.parquet")))
            if not _Path4(f).name.startswith(".")
        ]

    _candidates = [
        (_outputs / "full_benchmark_nogate", "2% gate-fixed"),
        (_outputs / "full_benchmark_ema_cats", "2% canonical"),
    ]
    _bench, _budget = next(
        ((b, lbl) for b, lbl in _candidates if _parqs(b / "diagnostics" / _lib)),
        (None, None),
    )
    _dfiles = _parqs(_bench / "diagnostics" / _lib) if _bench else []
    _sfiles = _parqs(_bench / "search" / _lib) if _bench else []

    if not _dfiles or not _sfiles or not _rf:
        _outR = mo.md(f"⚠️ diagnostics/search/reagents not found for **{_lib}**.")
    else:
        # Per-shard column-union (RWS + TT-TS schemas differ), then filter to the
        # single representative query so the layer plots see one trial per method.
        _cols = set()
        for _f in _dfiles:
            _cols |= set(pl.read_parquet_schema(_f).keys())
        _cols = sorted(_cols)
        _parts = []
        for _f in _dfiles:
            _d = pl.read_parquet(_f)
            for _c in _cols:
                if _c not in _d.columns:
                    _d = _d.with_columns(pl.lit(None).alias(_c))
            _parts.append(_d.select(_cols))
        _diag = (
            pl.concat(_parts, how="vertical_relaxed")
            .filter(pl.col("query_id") == _q)
        )
        _search = (
            pl.concat(
                [pl.read_parquet(
                    _sf, columns=["Name", "method", "replicate", "phase", "query_id"])
                 for _sf in _sfiles],
                how="vertical_relaxed",
            )
            .filter(pl.col("query_id") == _q)
            .select(["Name", "method", "replicate", "phase"])
        )
        _methods = ["TACTICS Enhanced-RWS (GMIC)", "TACTICS Enhanced-TT-TS (GMIC)"]
        _mlabels = {
            "TACTICS Enhanced-RWS (GMIC)": "TACTICS RWS",
            "TACTICS Enhanced-TT-TS (GMIC)": "TACTICS TT-TS",
        }

        _f1 = _gde_r(
            _search, _diag, _rf, methods=_methods, method_labels=_mlabels,
            component_names=_rcn,
            title=f"Layer 1 — GMIC directs WHICH component to explore  ({_lib} · {_q}, {_budget})",
        )
        _f2 = _ai_r(
            _search, _diag, _rf, methods=_methods, method_labels=_mlabels,
            component_names=_rcn,
            title=(
                f"Layer 2 — within-component exploration intensity  ({_lib} · {_q}, {_budget})\n"
                "both knobs are ×-multipliers centred on 1.0:  "
                "RWS = CATS temperature multiplier;  TT-TS = σ-inflation multiplier"
            ),
        )

        def _png(_f):
            _b = _io4.BytesIO()
            _f.savefig(_b, format="png", dpi=110, bbox_inches="tight")
            return _b.getvalue()

        _png1, _pdf1 = _png(_f1), fig_to_pdf_bytes(_f1)
        _png2, _pdf2 = _png(_f2), fig_to_pdf_bytes(_f2)
        plt.close(_f1)
        plt.close(_f2)
        _outR = mo.vstack([
            mo.md(f"### Layer 1 — GMIC directs *which* component  ({_lib} · {_q} · {_budget})"),
            mo.image(_png1),
            mo.download(
                data=_pdf1, filename=f"diag_gmic_directed_{_lib}_{_q}.pdf",
                mimetype="application/pdf", label="Download Layer 1 PDF",
            ),
            mo.md(f"### Layer 2 — within-component intensity  ({_lib} · {_q} · {_budget})"),
            mo.image(_png2),
            mo.download(
                data=_pdf2, filename=f"diag_adaptive_intensity_{_lib}_{_q}.pdf",
                mimetype="application/pdf", label="Download Layer 2 PDF",
            ),
        ])
    _outR
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## SI · GMIC routing is query-dependent on ROCS (scales with separation)

    Backing figure for the one-line main-text claim that GMIC selects a
    *different* component to explore depending on the query.  One point per
    (ROCS library x query): **x** = between-component GMIC separation
    (max - min late-search GMIC); **y** = fraction of search cycles the
    *lowest-GMIC* (explored) component was heated.  Dashed line = round-robin
    (1/n_components).  GMIC tilts exploration toward the unresolved component
    **in proportion to how separated the components are** -- real adaptive
    routing, but modest, and near round-robin when the components are similar
    (the typical ROCS query).  The Layer-1 demo queries sit in the
    high-separation tail.  Loads all 20 ROCS libraries (~1-2 min).
    """)
    return


@app.cell
def _(fig_to_pdf_bytes, mech_diag_dir, mo, np, pl, plt):
    """SI: per-query GMIC routing tilt vs between-component separation, pooled
    over the 20 ROCS libraries and split by component count.  Reads the merged
    per-library RWS diagnostics from full_benchmark_nogate."""
    _RWS = "TACTICS Enhanced-RWS (GMIC)"
    _LAST_N = 20
    _ROCS = [
        "rxn101", "rxn102a", "rxn108b", "rxn111b", "rxn114b", "rxn203", "rxn205",
        "rxn206", "rxn207", "rxn208", "amide-suzuki", "betti", "dobener",
        "groebke-blackburn-bienayme", "mannich", "niementowski", "orru",
        "passerini", "petasis", "poparov",
    ]
    _cols = ["method", "query_id", "component_idx", "replicate",
             "current_cycle", "total_cycles", "gmic", "is_heated"]

    _rows = []
    for _lib in _ROCS:
        _p = mech_diag_dir / _lib / f"{_lib}.parquet"
        if not _p.exists():
            continue
        _df = pl.read_parquet(_p, columns=_cols).filter(pl.col("method") == _RWS)
        _nc = _df["component_idx"].n_unique()
        _late = _df.filter(pl.col("current_cycle") > (pl.col("total_cycles") - _LAST_N))
        _g = (
            _late.group_by(["query_id", "component_idx", "replicate"])
            .agg(pl.col("gmic").mean().alias("g"))
            .group_by(["query_id", "component_idx"])
            .agg(pl.col("g").mean().alias("lg"))
        )
        _hf = _df.group_by(["query_id", "component_idx"]).agg(
            pl.col("is_heated").cast(pl.Float64).mean().alias("hf")
        )
        _gq = _g.group_by("query_id").agg(
            (pl.col("lg").max() - pl.col("lg").min()).alias("gap"),
            pl.col("component_idx").sort_by("lg").first().alias("chosen"),
        )
        _j = _gq.join(_hf, left_on=["query_id", "chosen"],
                      right_on=["query_id", "component_idx"])
        for _r in _j.iter_rows(named=True):
            _rows.append(dict(lib=_lib, n_comp=_nc, gap=_r["gap"],
                              chosen_share=_r["hf"]))

    _data = pl.DataFrame(_rows)

    _fig, _axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for _ax, _ncomp, _base, _title in [
        (_axes[0], 2, 0.50, "2-component ROCS"),
        (_axes[1], 3, 1.0 / 3.0, "3-component ROCS"),
    ]:
        _d = _data.filter(pl.col("n_comp") == _ncomp)
        _x, _y = _d["gap"].to_numpy(), _d["chosen_share"].to_numpy()
        _ax.scatter(_x, _y, s=12, alpha=0.30, color="#3b6fb6", edgecolors="none")
        _ax.axhline(_base, color="crimson", ls="--", lw=1.4,
                    label=f"round-robin = {_base:.2f}")
        _m, _b = np.polyfit(_x, _y, 1)
        _xs = np.linspace(_x.min(), _x.max(), 50)
        _ax.plot(_xs, _m * _xs + _b, color="black", lw=1.6)
        _ax.text(0.96, 0.07,
                 f"r = {np.corrcoef(_x, _y)[0, 1]:.2f}\n"
                 f"{len(_x)} queries, {_d['lib'].n_unique()} libs",
                 transform=_ax.transAxes, ha="right", va="bottom", fontsize=9)
        _ax.set_title(_title, fontsize=11)
        _ax.set_xlabel("between-component GMIC separation (max - min)")
        _ax.set_ylabel("heating share of lowest-GMIC component")
        _ax.legend(loc="upper left", fontsize=8, frameon=False)
    _fig.suptitle("GMIC exploration routing on ROCS is query-dependent and "
                  "scales with component separation", fontsize=12)
    _fig.tight_layout(rect=[0, 0, 1, 0.95])

    mo.vstack([
        plt.gcf(),
        mo.download(
            data=fig_to_pdf_bytes(_fig),
            filename="SI_gmic_routing_query_dependence.pdf",
            mimetype="application/pdf", label="Download as PDF",
        ),
    ])
    return


if __name__ == "__main__":
    app.run()
