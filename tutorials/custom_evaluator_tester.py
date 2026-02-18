"""
Custom Evaluator Tester (Marimo)

Run:
    marimo run tutorials/custom_evaluator_tester.py
"""

import marimo

__generated_with = "0.19.4"
app = marimo.App(width="medium", app_title="Custom Evaluator Tester")


@app.cell
def _():
    """Imports, fallback path setup, and constants."""
    import ast
    import importlib.resources
    import inspect
    import logging
    import sys
    import traceback
    from io import StringIO
    from pathlib import Path

    import marimo as mo
    import polars as pl
    from rdkit import Chem

    try:
        from TACTICS.library_enumeration import SynthesisPipeline
        from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef
        from TACTICS.thompson_sampling import ThompsonSampler
        from TACTICS.thompson_sampling.config import ThompsonSamplingConfig
        from TACTICS.thompson_sampling.core.evaluator_config import CustomEvaluatorConfig
        from TACTICS.thompson_sampling.core.evaluators import CustomEvaluator
        from TACTICS.thompson_sampling.strategies.config import GreedyConfig
        from TACTICS.thompson_sampling.warmup.config import StandardWarmupConfig
    except ModuleNotFoundError:
        project_root = Path(__file__).resolve().parents[1]
        src_root = project_root / "src"
        if str(src_root) not in sys.path:
            sys.path.insert(0, str(src_root))
        from TACTICS.library_enumeration import SynthesisPipeline
        from TACTICS.library_enumeration.smarts_toolkit import ReactionConfig, ReactionDef
        from TACTICS.thompson_sampling import ThompsonSampler
        from TACTICS.thompson_sampling.config import ThompsonSamplingConfig
        from TACTICS.thompson_sampling.core.evaluator_config import CustomEvaluatorConfig
        from TACTICS.thompson_sampling.core.evaluators import CustomEvaluator
        from TACTICS.thompson_sampling.strategies.config import GreedyConfig
        from TACTICS.thompson_sampling.warmup.config import StandardWarmupConfig

    _data_files = importlib.resources.files("TACTICS.data.thrombin")
    acids_file = str(_data_files / "acids.smi")
    amines_file = str(_data_files / "coupled_aa_sub.smi")
    reagent_files = [acids_file, amines_file]

    amide_coupling_smarts = (
        "[#6:1](=[O:2])[OH]."
        "[#7X3;H1,H2;!$(N[!#6]);!$(N[#6]=[O]);!$(N[#6]~[!#6;!#16]):3]"
        ">>[#6:1](=[O:2])[#7:3]"
    )

    smiles_shortlist = [
        "CCO",
        "c1ccccc1",
        "CC(=O)O",
        "CCN(CC)CC",
        "c1ccncc1O",
    ]

    return (
        ast,
        Chem,
        CustomEvaluator,
        CustomEvaluatorConfig,
        GreedyConfig,
        Path,
        ReactionConfig,
        ReactionDef,
        StandardWarmupConfig,
        SynthesisPipeline,
        ThompsonSampler,
        ThompsonSamplingConfig,
        amide_coupling_smarts,
        inspect,
        logging,
        mo,
        pl,
        reagent_files,
        smiles_shortlist,
        StringIO,
        traceback,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
        # Custom Evaluator Tester

        Validate a user-provided scoring function before using it in Thompson Sampling.

        ## Workflow
        1. Select code source mode with the checkbox.
        2. Provide code snippet or file path (based on selected mode).
        3. Set function name (snippet mode) or select it from dropdown (file mode).
        4. Click **Run SMILES shortlist**.
        5. Click **Run minimal TS (Thrombin)**.

        ## Required function template

        ```python
        def your_function_name(mol):
            ...
            return <float>
        ```

        - Function must accept one positional argument (`mol`).
        - Function must return a numeric value (`int` or `float`). Non-numeric returns or exceptions will be treated as errors.
        - Include all required imports in your snippet/file (for example, `from rdkit import Chem` if needed).
        
        """
    )
    return


@app.cell
def _(mo):
    """Code source and function-name widgets."""
    use_snippet_checkbox = mo.ui.checkbox(value=True, label="Toggle Code Snippet input / File path input")

    _default_snippet = (
        "def user_scoring_function(mol):\n"
        "    # TODO: implement your score\n"
        "    # Example:\n"
        "    #return mol.GetNumHeavyAtoms()\n"
        "    raise NotImplementedError('Implement your scoring function')\n"
    )

    if hasattr(mo.ui, "text_area"):
        code_snippet_input = mo.ui.text_area(
            value=_default_snippet,
            label="Code snippet (full function definition)",
            rows=10,
            full_width=True,
        )
    else:
        code_snippet_input = mo.ui.text(
            value=_default_snippet,
            label="Code snippet (full function definition)",
            full_width=True,
        )

    code_file_input = mo.ui.text(
        value="",
        label="Code file path (.py)",
        placeholder="/path/to/custom_score.py",
        full_width=True,
    )

    function_name_input = mo.ui.text(
        value="user_scoring_function",
        label="Scoring function name",
        full_width=True,
    )

    return (
        code_file_input,
        code_snippet_input,
        function_name_input,
        use_snippet_checkbox,
    )


@app.cell
def _(Path, ast, code_file_input, mo, use_snippet_checkbox):
    """Read/parse file mode source and expose function dropdown options."""
    file_functions = []
    file_function_error = None

    if not use_snippet_checkbox.value:
        _file_path = (code_file_input.value or "").strip()
        if not _file_path:
            file_function_error = "File-path mode selected, but file path is empty."
        else:
            try:
                _source_text = Path(_file_path).read_text(encoding="utf-8")
                _tree = ast.parse(_source_text)
                file_functions = [
                    _node.name
                    for _node in _tree.body
                    if isinstance(_node, (ast.FunctionDef, ast.AsyncFunctionDef))
                ]
                if not file_functions:
                    file_function_error = "No top-level functions found in the selected file."
            except Exception as _exc:
                file_function_error = f"Failed to parse file: {_exc}"

    _dropdown_options = file_functions if file_functions else ["<no functions available>"]
    _dropdown_value = _dropdown_options[0]
    file_function_dropdown = mo.ui.dropdown(
        options=_dropdown_options,
        value=_dropdown_value,
        label="Scoring function (from file)",
    )

    return file_function_dropdown, file_function_error, file_functions


@app.cell
def _(
    code_file_input,
    code_snippet_input,
    file_function_dropdown,
    file_function_error,
    file_functions,
    function_name_input,
    mo,
    use_snippet_checkbox,
):
    _source_mode_label = "snippet" if use_snippet_checkbox.value else "file path"

    if use_snippet_checkbox.value:
        _source_panel = mo.vstack(
            [
                mo.md(f"## Function Input\nSelected input mode: **{_source_mode_label}**"),
                use_snippet_checkbox,
                function_name_input,
                code_snippet_input,
            ]
        )
    else:
        _children = [
            mo.md(f"## Function Input\nSelected input mode: **{_source_mode_label}**"),
            use_snippet_checkbox,
            code_file_input,
        ]
        if file_function_error is None and len(file_functions) > 0:
            _children.append(file_function_dropdown)
        else:
            _children.append(
                mo.md(
                    f"**Error:** {file_function_error or 'No functions available from file.'}"
                )
            )
        _source_panel = mo.vstack(_children)

    _source_panel
    return


@app.cell
def _(
    Path,
    code_file_input,
    code_snippet_input,
    file_function_dropdown,
    file_function_error,
    file_functions,
    function_name_input,
    inspect,
    traceback,
    use_snippet_checkbox,
):
    """Load and validate scoring function from selected source mode."""
    _snippet = (code_snippet_input.value or "").strip()
    _file_path_raw = (code_file_input.value or "").strip()
    _snippet_function_name = (function_name_input.value or "").strip()
    _file_function_name = (file_function_dropdown.value or "").strip()

    resolved_source = ""
    resolved_source_kind = None
    resolved_function_name = ""
    load_error = None
    loaded_function = None

    if use_snippet_checkbox.value:
        resolved_source_kind = "snippet"
        resolved_function_name = _snippet_function_name
        if not resolved_function_name:
            load_error = "Function name cannot be empty in snippet mode."
        elif not _snippet:
            load_error = "Snippet mode selected, but snippet is empty."
        else:
            resolved_source = _snippet
    else:
        resolved_source_kind = "file"
        if file_function_error is not None:
            load_error = file_function_error
        elif len(file_functions) == 0:
            load_error = "No functions available to select from file."
        elif _file_function_name == "<no functions available>" or not _file_function_name:
            load_error = "Select a scoring function from the file dropdown."
        elif not _file_path_raw:
            load_error = "File-path mode selected, but file path is empty."
        else:
            resolved_function_name = _file_function_name
            try:
                _path_obj = Path(_file_path_raw)
                resolved_source = _path_obj.read_text(encoding="utf-8")
            except Exception:
                load_error = (
                    f"Failed to read file '{_file_path_raw}'.\n"
                    f"{traceback.format_exc()}"
                )

    if load_error is None:
        _namespace = {}
        try:
            exec(resolved_source, _namespace, _namespace)
            _candidate = _namespace.get(resolved_function_name)

            if _candidate is None or not callable(_candidate):
                load_error = (
                    f"Code must define callable '{resolved_function_name}(mol)'."
                )
            else:
                _signature = inspect.signature(_candidate)
                _params = list(_signature.parameters.values())
                _has_varargs = any(
                    _p.kind == inspect.Parameter.VAR_POSITIONAL for _p in _params
                )
                _positional_like = [
                    _p
                    for _p in _params
                    if _p.kind
                    in (
                        inspect.Parameter.POSITIONAL_ONLY,
                        inspect.Parameter.POSITIONAL_OR_KEYWORD,
                    )
                ]
                if not _has_varargs and len(_positional_like) < 1:
                    load_error = (
                        f"'{resolved_function_name}' must accept at least one "
                        "positional argument (mol)."
                    )
                else:
                    loaded_function = _candidate
        except Exception:
            load_error = (
                "Error while executing provided code.\n"
                f"{traceback.format_exc()}"
            )

    return (
        load_error,
        loaded_function,
        resolved_function_name,
        resolved_source,
        resolved_source_kind,
    )


@app.cell
def _(load_error, mo, resolved_function_name, resolved_source, resolved_source_kind):
    if load_error is not None:
        _loader_display = mo.vstack(
            [
                mo.md("### Loader status: ❌ Error"),
                mo.md(f"```text\n{load_error}\n```")
            ]
        )
    else:
        _preview = resolved_source.strip().splitlines()
        _preview_text = "\n".join(_preview[:12])
        if len(_preview) > 12:
            _preview_text += "\n..."
        _loader_display = mo.vstack(
            [
                mo.md(
                    f"### Loader status: ✅ Loaded from **{resolved_source_kind}** as "
                    f"`{resolved_function_name}`"
                ),
                mo.md(f"```python\n{_preview_text}\n```")
            ]
        )

    _loader_display
    return


@app.cell
def _(mo):
    run_shortlist_button = mo.ui.run_button(label="Run SMILES shortlist")
    run_shortlist_button
    return (run_shortlist_button,)


@app.cell
def _(mo, run_shortlist_button):
    if run_shortlist_button.value:
        _shortlist_status_display = mo.md("**Running shortlist ...**")
    else:
        _shortlist_status_display = mo.md("*Click 'Run SMILES shortlist' to start.*")
    _shortlist_status_display
    return


@app.cell
def _(Chem, load_error, loaded_function, mo, pl, run_shortlist_button, smiles_shortlist):
    """Run shortlist validation and always render output in marimo run."""
    _empty_df = pl.DataFrame(
        {
            "SMILES": [],
            "is_valid_smiles": [],
            "status": [],
            "score_raw": [],
            "score_float": [],
            "error_message": [],
        }
    )

    shortlist_df = _empty_df
    shortlist_summary = {
        "available": False,
        "executed": False,
        "total_tested": 0,
        "valid_mol_count": 0,
        "success_count": 0,
        "failure_count": 0,
    }

    if not run_shortlist_button.value:
        _shortlist_display = mo.md("## 1) Shortlist Validation\nWaiting for run button.")
    elif load_error is not None or loaded_function is None:
        _shortlist_display = mo.vstack(
            [
                mo.md("## 1) Shortlist Validation"),
                mo.md("Cannot run shortlist: function loader failed."),
                mo.md(f"```text\n{load_error}\n```")
            ]
        )
    else:
        _records = []
        for _smi in smiles_shortlist:
            _mol = Chem.MolFromSmiles(_smi)
            _row = {
                "SMILES": _smi,
                "is_valid_smiles": _mol is not None,
                "status": "error",
                "score_raw": None,
                "score_float": None,
                "error_message": None,
            }
            if _mol is None:
                _row["error_message"] = "Invalid SMILES"
                _records.append(_row)
                continue

            try:
                _raw_score = loaded_function(_mol)
                _row["score_raw"] = repr(_raw_score)
                _row["score_float"] = float(_raw_score)
                _row["status"] = "ok"
            except Exception as _exc:
                _row["error_message"] = f"{type(_exc).__name__}: {_exc}"
            _records.append(_row)

        shortlist_df = pl.DataFrame(_records)
        _total_tested = len(_records)
        _valid_mol_count = sum(1 for _r in _records if _r["is_valid_smiles"])
        _success_count = sum(1 for _r in _records if _r["status"] == "ok")
        _failure_count = _total_tested - _success_count

        shortlist_summary = {
            "available": True,
            "executed": True,
            "total_tested": _total_tested,
            "valid_mol_count": _valid_mol_count,
            "success_count": _success_count,
            "failure_count": _failure_count,
        }

        _shortlist_display = mo.vstack(
            [
                mo.md("## 1) Shortlist Validation"),
                mo.md(
                    f"- total tested: **{_total_tested}**\n"
                    f"- valid mol count: **{_valid_mol_count}**\n"
                    f"- success count: **{_success_count}**\n"
                    f"- failure count: **{_failure_count}**"
                ),
                shortlist_df,
            ]
        )

    _shortlist_display
    return shortlist_df, shortlist_summary


@app.cell
def _(Chem, CustomEvaluator, load_error, loaded_function, mo, pl, run_shortlist_button, smiles_shortlist):
    """CustomEvaluator wrapper check; always returns display in run mode."""
    _empty_df = pl.DataFrame(
        {
            "SMILES": [],
            "evaluator_score": [],
            "status": [],
            "error_message": [],
        }
    )

    integration_df = _empty_df

    if not run_shortlist_button.value:
        _integration_display = mo.md("## 2) CustomEvaluator Integration Check\nWaiting for shortlist run.")
    elif load_error is not None or loaded_function is None:
        _integration_display = mo.vstack(
            [
                mo.md("## 2) CustomEvaluator Integration Check"),
                mo.md("Skipped: function loader failed."),
            ]
        )
    else:
        _evaluator = CustomEvaluator(scoring_function=loaded_function)
        _records_eval = []

        for _smi in smiles_shortlist:
            _mol = Chem.MolFromSmiles(_smi)
            if _mol is None:
                _records_eval.append(
                    {
                        "SMILES": _smi,
                        "evaluator_score": None,
                        "status": "error",
                        "error_message": "Invalid SMILES",
                    }
                )
                continue

            try:
                _score = _evaluator.evaluate(_mol)
                _is_nan = _score != _score
                _status = "ok" if not _is_nan else "error"
                _error_message = (
                    None
                    if not _is_nan
                    else "Evaluator returned NaN (function raised or returned non-numeric)"
                )
                _records_eval.append(
                    {
                        "SMILES": _smi,
                        "evaluator_score": None if _is_nan else float(_score),
                        "status": _status,
                        "error_message": _error_message,
                    }
                )
            except Exception as _exc:
                _records_eval.append(
                    {
                        "SMILES": _smi,
                        "evaluator_score": None,
                        "status": "error",
                        "error_message": f"{type(_exc).__name__}: {_exc}",
                    }
                )

        integration_df = pl.DataFrame(_records_eval)
        _integration_display = mo.vstack(
            [
                mo.md("## 2) CustomEvaluator Integration Check"),
                integration_df,
            ]
        )

    _integration_display
    return (integration_df,)


@app.cell
def _(mo):
    run_ts_button = mo.ui.run_button(label="Run minimal TS (Thrombin)")
    run_ts_button
    return (run_ts_button,)


@app.cell
def _(mo, run_ts_button):
    if run_ts_button.value:
        _ts_status_display = mo.md("**Running minimal TS ...**")
    else:
        _ts_status_display = mo.md("*Click 'Run minimal TS (Thrombin)' to start.*")
    _ts_status_display
    return


@app.cell
def _(
    Chem,
    CustomEvaluatorConfig,
    GreedyConfig,
    ReactionConfig,
    ReactionDef,
    StringIO,
    StandardWarmupConfig,
    SynthesisPipeline,
    ThompsonSampler,
    ThompsonSamplingConfig,
    amide_coupling_smarts,
    load_error,
    loaded_function,
    logging,
    mo,
    pl,
    reagent_files,
    run_ts_button,
    shortlist_summary,
    smiles_shortlist,
    traceback,
):
    """Run minimal TS and always render output in marimo run."""
    combined = None
    search_df = None
    warmup_df = None
    _captured_logs = ""

    if not run_ts_button.value:
        _ts_display = mo.md("## 3) Minimal Thrombin TS Run\nWaiting for run button.")
    elif load_error is not None or loaded_function is None:
        _ts_display = mo.vstack(
            [
                mo.md("## 3) Minimal Thrombin TS Run"),
                mo.md("Cannot run TS: function loader failed."),
                mo.md(f"```text\n{load_error}\n```")
            ]
        )
    elif not shortlist_summary.get("executed", False):
        _ts_display = mo.vstack(
            [
                mo.md("## 3) Minimal Thrombin TS Run"),
                mo.md("Run-order guard: click **Run SMILES shortlist** first."),
            ]
        )
    elif shortlist_summary.get("success_count", 0) == 0:
        _ts_display = mo.vstack(
            [
                mo.md("## 3) Minimal Thrombin TS Run"),
                mo.md(
                    "Shortlist shows **0 successful numeric scores**. Fix function and rerun shortlist."
                ),
            ]
        )
    else:
        _precheck_success = 0
        for _smi in smiles_shortlist:
            _mol = Chem.MolFromSmiles(_smi)
            if _mol is None:
                continue
            try:
                _value = float(loaded_function(_mol))
                if _value == _value:
                    _precheck_success += 1
            except Exception:
                continue

        if _precheck_success == 0:
            _ts_display = mo.vstack(
                [
                    mo.md("## 3) Minimal Thrombin TS Run"),
                    mo.md("Precheck indicates 0 valid numeric outputs. Fix function first."),
                ]
            )
        else:
            _sampler = None
            _log_buffer = StringIO()
            _capture_handler = logging.StreamHandler(_log_buffer)
            _capture_handler.setLevel(logging.INFO)
            _capture_handler.setFormatter(
                logging.Formatter("%(asctime)s %(levelname)s %(message)s")
            )
            _sampler_logger = logging.getLogger("TACTICS.thompson_sampling.core.sampler")
            _sampler_logger.addHandler(_capture_handler)
            try:
                _reaction_config = ReactionConfig(
                    reactions=[
                        ReactionDef(
                            reaction_smarts=amide_coupling_smarts,
                            step_index=0,
                            description="Amide coupling",
                        )
                    ],
                    reagent_file_list=reagent_files,
                )
                _pipeline = SynthesisPipeline(_reaction_config)

                _ts_config = ThompsonSamplingConfig(
                    synthesis_pipeline=_pipeline,
                    num_ts_iterations=3,
                    num_warmup_trials=3,
                    strategy_config=GreedyConfig(mode="maximize"),
                    warmup_config=StandardWarmupConfig(),
                    evaluator_config=CustomEvaluatorConfig(scoring_function=loaded_function),
                    batch_size=1,
                    max_resamples=10,
                    hide_progress=True,
                )

                _sampler = ThompsonSampler.from_config(_ts_config)
                warmup_df = _sampler.warm_up(num_warmup_trials=_ts_config.num_warmup_trials)
                search_df = _sampler.search(num_cycles=_ts_config.num_ts_iterations)
                combined = pl.concat([warmup_df, search_df])
                _captured_logs = _log_buffer.getvalue().strip()
                if not _captured_logs:
                    _captured_logs = "No logs captured."

                _ts_display = mo.vstack(
                    [
                        mo.md("## 3) Minimal Thrombin TS Run"),
                        mo.md(
                            f"Run complete. Rows: warmup={len(warmup_df)}, "
                            f"search={len(search_df)}, combined={len(combined)}"
                        ),
                        mo.md("### Logged output"),
                        mo.md(f"```text\n{_captured_logs}\n```"),
                        combined.head(10),
                    ]
                )
            except Exception:
                _err = traceback.format_exc()
                _captured_logs = _log_buffer.getvalue().strip()
                if not _captured_logs:
                    _captured_logs = "No logs captured."
                _ts_display = mo.vstack(
                    [
                        mo.md("## 3) Minimal Thrombin TS Run"),
                        mo.md(
                            """
                            TS run failed. Common causes:
                            - scoring function raises exceptions on some molecules
                            - scoring function returns non-numeric outputs
                            - environment/dependency issues
                            """
                        ),
                        mo.md(f"```text\n{_err}\n```"),
                        mo.md("### Logged output"),
                        mo.md(f"```text\n{_captured_logs}\n```"),
                    ]
                )
            finally:
                _sampler_logger.removeHandler(_capture_handler)
                _capture_handler.close()
                if _sampler is not None:
                    _sampler.close()

    _ts_display
    return combined, search_df, warmup_df


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
        When custom scoring function works as intended, use it for your own TS pipelines in `thompson_sampling_tutorial.py`.
        """
    )
    return


if __name__ == "__main__":
    app.run()
