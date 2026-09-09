"""Guard the lazy-import contract.

Each check runs in a fresh subprocess so nothing is pre-cached in
``sys.modules``. These tests pin the 2.0 import-weight guarantees: a config
import must not load the chemistry/plotting stack, the package hubs resolve
names lazily, and the plotting modules do not import matplotlib/altair until
a plot is actually requested.
"""

import subprocess
import sys

import pytest


def _run(code: str) -> str:
    proc = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert proc.returncode == 0, proc.stderr
    return proc.stdout


HEAVY = ("rdkit", "scipy", "matplotlib", "altair", "sqlitedict", "pandas", "seaborn")


def test_config_import_does_not_load_heavy_stack():
    out = _run(
        "import sys\n"
        "from TACTICS.thompson_sampling.config import ThompsonSamplingConfig\n"
        f"print(sorted(m for m in {HEAVY!r} if m in sys.modules))\n"
    )
    assert out.strip() == "[]", f"config import pulled in: {out.strip()}"


def test_config_import_is_fast():
    out = _run(
        "import time\n"
        "t = time.perf_counter()\n"
        "from TACTICS.thompson_sampling.config import ThompsonSamplingConfig\n"
        "print(time.perf_counter() - t)\n"
    )
    assert float(out) < 1.0, f"config import took {float(out):.2f}s"


def test_top_level_hub_is_lazy():
    out = _run(
        "import sys, TACTICS\n"
        "before = 'rdkit' in sys.modules\n"
        "from TACTICS import ThompsonSampler, get_preset\n"
        "after = 'rdkit' in sys.modules\n"
        "print(before, after, ThompsonSampler.__name__, get_preset.__name__)\n"
    )
    assert out.split() == ["False", "True", "ThompsonSampler", "get_preset"]


def test_hub_unknown_attribute_raises():
    out = _run(
        "import TACTICS\n"
        "try:\n"
        "    TACTICS.does_not_exist\n"
        "except AttributeError as e:\n"
        "    print('AttributeError', 'does_not_exist' in str(e))\n"
    )
    assert out.split() == ["AttributeError", "True"]


def test_diagnostic_plots_import_does_not_load_plotting():
    out = _run(
        "import sys\n"
        "import TACTICS.library_analysis.diagnostic_plots\n"
        "print(sorted(m for m in ('matplotlib', 'altair', 'seaborn', 'pandas') if m in sys.modules))\n"
    )
    assert out.strip() == "[]", f"diagnostic_plots import pulled in: {out.strip()}"


def test_submodule_paths_still_work():
    _run(
        "from TACTICS.thompson_sampling.core.evaluators import LookupEvaluator\n"
        "from TACTICS.thompson_sampling.core import Reagent\n"
        "from TACTICS.thompson_sampling import TopTwoSelection, EnhancedWarmup\n"
        "from TACTICS.library_enumeration import SynthesisPipeline\n"
    )


@pytest.mark.skipif(
    subprocess.run([sys.executable, "-c", "import matplotlib"], capture_output=True).returncode != 0,
    reason="matplotlib not installed",
)
def test_plot_function_loads_matplotlib_on_call():
    out = _run(
        "import sys, polars as pl\n"
        "from TACTICS.library_analysis.diagnostic_plots import plot_criticality_trajectory\n"
        "before = 'matplotlib' in sys.modules\n"
        "df = pl.DataFrame({'current_cycle': [0, 1], 'component_idx': [0, 0], 'criticality': [0.1, 0.2]})\n"
        "fig = plot_criticality_trajectory(df)\n"
        "print(before, 'matplotlib' in sys.modules, type(fig).__name__)\n"
    )
    assert out.split() == ["False", "True", "Figure"]
