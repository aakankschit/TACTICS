"""
Interactive SAR Explorer — Hover-to-Structure with molplotly

Hover over reagent dots to see molecular structures rendered by RDKit.
Per-component oracle GMIC shows whether a component needs exploration
(low GMIC) or can be exploited (high GMIC).

Run as app:  marimo run tutorials/interactive_sar_explorer.py
Edit mode:   marimo edit tutorials/interactive_sar_explorer.py
"""

import marimo

__generated_with = "0.21.1"
app = marimo.App(width="full", app_title="Interactive SAR Explorer")


@app.cell
def _():
    import marimo as mo
    import socket
    import threading
    from pathlib import Path
    from io import BytesIO

    import numpy as np
    import pandas as pd
    import polars as pl
    import plotly.express as px
    import molplotly
    import re

    try:
        project_root = Path(__file__).parent.parent.resolve()
    except NameError:
        project_root = Path.cwd()

    SCORES_DIR = project_root / "data" / "scores"
    REAGENTS_DIR = project_root / "data" / "reagents"

    return (
        SCORES_DIR, REAGENTS_DIR, BytesIO, mo, molplotly,
        np, pd, pl, project_root, px, re, socket, threading,
    )


@app.cell
def _(mo):
    library_selector = mo.ui.dropdown(
        options={
            "Adenine (3-comp docking, +62.1)": "adenine",
            "Quinazoline (3-comp docking, +0.1 ns)": "quinazoline",
            "Thrombin (2-comp docking, +3.2)": "thrombin",
            "Amide (2-comp docking, +5.1)": "amide",
            "rxn206 (2-comp ROCS, +7.6)": "rxn206",
        },
        value="Adenine (3-comp docking, +62.1)",
        label="Library",
    )
    top_n_selector = mo.ui.slider(
        start=25, stop=250, step=25, value=100, label="Top-N threshold",
    )
    mo.vstack([
        mo.md(r"""
# Interactive SAR Explorer

**Hover** over dots to see reagent structures.
Each panel is one reaction component. Oracle GMIC quantifies
how informative per-reagent mean scores are: **low GMIC** = high
flexibility (reagent means don't distinguish hits) /
**high GMIC** = high criticality (reagent means separate hits).
        """),
        mo.hstack([library_selector, top_n_selector], gap=2),
    ])
    return library_selector, top_n_selector


@app.cell
def _(library_selector, mo, pl, re, SCORES_DIR, REAGENTS_DIR):
    """Load scores and reagent SMILES; create query selector for ROCS."""
    _lib_id = library_selector.value

    def _parse_adenine(code):
        _m = re.match(r"(amidine_\d+)_(isocyanide_db_\d+)_(aldehyde_\d+)", code)
        return list(_m.groups()) if _m else None

    def _parse_thrombin(code):
        _parts = code.split("_")
        return [_parts[0], _parts[1] + "_" + _parts[2]] if len(_parts) == 3 else None

    def _parse_split(code, n):
        _parts = code.split("_")
        return _parts if len(_parts) == n else None

    _CONFIGS = {
        "adenine": {
            "components": ["Amidines", "Isocyanides", "Aldehydes"],
            "reagent_files": [
                "adenine_reagent_0.smi", "adenine_reagent_1.smi",
                "adenine_reagent_2.smi",
            ],
            "code_col": "Product_Code", "mode": "minimize",
            "parse": _parse_adenine,
        },
        "quinazoline": {
            "components": ["Aminobenzoics", "Amines", "Acids"],
            "reagent_files": [
                "quinazoline_reagent_0.smi", "quinazoline_reagent_1.smi",
                "quinazoline_reagent_2.smi",
            ],
            "code_col": "Product_Code", "mode": "minimize",
            "parse": lambda c: _parse_split(c, 3),
        },
        "thrombin": {
            "components": ["Acids", "Dipeptides"],
            "reagent_files": [
                "thrombin_reagent_0.smi", "thrombin_reagent_1.smi",
            ],
            "code_col": "Product_Code", "mode": "minimize",
            "parse": _parse_thrombin,
        },
        "amide": {
            "components": ["Amines", "Acids"],
            "reagent_files": [
                "amide_reagent_0.smi", "amide_reagent_1.smi",
            ],
            "code_col": "Product_Code", "mode": "minimize",
            "parse": lambda c: _parse_split(c, 2),
        },
        "rxn206": {
            "components": ["Component 0", "Component 1"],
            "reagent_files": [
                "rxn206_reagent_0.smi", "rxn206_reagent_1.smi",
            ],
            "code_col": "Name", "mode": "maximize",
            "parse": lambda c: _parse_split(c, 2),
        },
    }

    lib_cfg = _CONFIGS[_lib_id]
    score_df = pl.read_parquet(SCORES_DIR / f"{_lib_id}.parquet")

    if lib_cfg["mode"] == "maximize":
        _qcols = sorted(c for c in score_df.columns if c.startswith("query_"))
        query_selector = mo.ui.dropdown(
            options=_qcols, value=_qcols[0], label="Query molecule",
        )
    else:
        query_selector = mo.ui.dropdown(
            options=["Scores"], value="Scores", label="Score column",
        )

    smiles_maps = []
    for _fname in lib_cfg["reagent_files"]:
        _smap = {}
        with open(REAGENTS_DIR / _lib_id / _fname, errors="replace") as _f:
            for _line in _f:
                _line = _line.strip()
                if not _line:
                    continue
                _parts = _line.split("\t", 1) if "\t" in _line else _line.split(None, 1)
                if len(_parts) >= 2:
                    _smap[_parts[1].strip()] = _parts[0].strip()
        smiles_maps.append(_smap)

    query_selector if lib_cfg["mode"] == "maximize" else mo.md("")
    return lib_cfg, query_selector, score_df, smiles_maps


@app.cell
def _(lib_cfg, query_selector, top_n_selector, score_df, smiles_maps, np, pd, pl):
    """Compute per-reagent stats and build chart DataFrame."""
    _score_col = query_selector.value
    _top_n = top_n_selector.value
    _mode = lib_cfg["mode"]
    _parse = lib_cfg["parse"]
    _code_col = lib_cfg["code_col"]
    _comp_names = lib_cfg["components"]

    _df = score_df.filter(pl.col(_score_col).is_not_nan())
    _ascending = _mode == "minimize"
    _top_codes = set(
        _df.sort(_score_col, descending=not _ascending)
        .head(_top_n)[_code_col]
        .to_list()
    )

    _codes = _df[_code_col].to_list()
    _scores = _df[_score_col].to_numpy()
    _parsed = [_parse(c) for c in _codes]
    _parsed_top = {c: _parse(c) for c in _top_codes}

    _rows = []
    for _ci, _cname in enumerate(_comp_names):
        _r2s = {}
        for _p, _s in zip(_parsed, _scores):
            if _p is not None:
                _r2s.setdefault(_p[_ci], []).append(_s)

        _top_counts = {}
        for _c, _p in _parsed_top.items():
            if _p is not None:
                _top_counts[_p[_ci]] = _top_counts.get(_p[_ci], 0) + 1

        _means = {r: float(np.mean(s)) for r, s in _r2s.items()}
        _stds = {r: float(np.std(s)) for r, s in _r2s.items()}

        _all_means = np.array(list(_means.values()))
        _all_stds = np.array(list(_stds.values()))
        _sig_var = float(np.var(_all_means))
        _noise_var = float(np.mean(_all_stds**2))
        _oracle_gmic = 0.5 * np.log(1 + _sig_var / _noise_var) if _noise_var > 0 else 0.0

        # Min-max normalization to [0, 1] for cross-component comparability
        _m_min, _m_max = float(np.min(_all_means)), float(np.max(_all_means))
        _s_min, _s_max = float(np.min(_all_stds)), float(np.max(_all_stds))
        _m_range = _m_max - _m_min if _m_max > _m_min else 1.0
        _s_range = _s_max - _s_min if _s_max > _s_min else 1.0

        _comp_title = f"{_cname} \u2014 GMIC {_oracle_gmic:.3f}"

        _smap = smiles_maps[_ci]
        for _rid in sorted(_r2s.keys()):
            _nm = (_means[_rid] - _m_min) / _m_range
            _ns = (_stds[_rid] - _s_min) / _s_range
            _rows.append({
                "component": _comp_title,
                "reagent_id": str(_rid),
                "smiles": _smap.get(str(_rid), ""),
                "mean": round(_means[_rid], 3),
                "std": round(_stds[_rid], 3),
                "norm_mean": round(_nm, 3),
                "norm_std": round(_ns, 3),
                "count": _top_counts.get(_rid, 0),
                "n_obs": len(_r2s[_rid]),
                "oracle_gmic": round(_oracle_gmic, 3),
            })

    chart_df = pd.DataFrame(_rows)
    return (chart_df,)


@app.cell
def _(chart_df, px, top_n_selector):
    """Build Plotly figure with faceted scatter per component."""
    _top_n = top_n_selector.value

    plotly_fig = px.scatter(
        chart_df,
        x="norm_mean",
        y="norm_std",
        color="count",
        facet_col="component",
        facet_col_wrap=3,
        hover_data=["reagent_id", "mean", "std", "n_obs", "oracle_gmic"],
        color_continuous_scale="YlOrRd",
        labels={
            "norm_mean": "Reagent Mean (normalized)",
            "norm_std": "Reagent Stdev (normalized)",
            "count": f"Top-{_top_n}",
        },
    )
    plotly_fig.update_traces(
        marker=dict(
            size=8,
            line=dict(width=0.5, color="#666"),
        ),
    )
    plotly_fig.update_layout(
        template="plotly_white",
        width=400 * chart_df["component"].nunique(),
        height=450,
        font=dict(size=13),
    )
    # Clean up facet titles (remove "component = " prefix)
    plotly_fig.for_each_annotation(
        lambda a: a.update(text=a.text.split("=", 1)[-1].strip())
    )
    return (plotly_fig,)


@app.cell
def _(plotly_fig, chart_df, molplotly, mo, socket, threading, top_n_selector, BytesIO):
    """Launch molplotly Dash app for hover-to-structure and embed via iframe."""
    _top_n = top_n_selector.value

    # Find a free port
    _sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    _sock.bind(("127.0.0.1", 0))
    _dash_port = _sock.getsockname()[1]
    _sock.close()

    _dash_app = molplotly.add_molecules(
        fig=plotly_fig,
        df=chart_df,
        smiles_col="smiles",
        title_col="reagent_id",
        show_coords=False,
        caption_cols=["count", "mean", "std", "oracle_gmic"],
        caption_transform={
            "count": lambda x: f"Top-{_top_n}: {x}",
            "oracle_gmic": lambda x: f"GMIC: {x}",
        },
        svg_size=250,
    )

    _thread = threading.Thread(
        target=lambda: _dash_app.run(
            host="127.0.0.1", port=_dash_port,
            debug=False, use_reloader=False,
        ),
        daemon=True,
    )
    _thread.start()

    # PDF export via kaleido
    _pdf_buf = BytesIO()
    plotly_fig.write_image(_pdf_buf, format="pdf", scale=2)
    _pdf_bytes = _pdf_buf.getvalue()

    mo.vstack([
        mo.md(
            f"**Interactive viewer** (hover for structures) — "
            f"[open in new tab](http://localhost:{_dash_port})"
        ),
        mo.Html(
            f'<iframe src="http://localhost:{_dash_port}" '
            f'width="100%" height="550" '
            f'style="border:1px solid #ddd; border-radius:8px;"></iframe>'
        ),
        mo.download(
            data=_pdf_bytes,
            filename="sar_scatter.pdf",
            mimetype="application/pdf",
            label="Download PDF",
        ),
    ])
    return


if __name__ == "__main__":
    app.run()
