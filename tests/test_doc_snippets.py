"""Every code example in the docs is a real file under docs/source/snippets/
and is executed here against the bundled thrombin data. A snippet that
stops running fails the suite -- that is the point.

Snippet contract:
- optional ``# requires: viz`` / ``# requires: openeye`` header line
  (``pytest.importorskip`` per tag); ``# requires: none`` or absent = always run
- code lives in ``main()`` guarded by ``if __name__ == "__main__"`` so
  ``processes > 1`` snippets survive multiprocessing's spawn re-import
- relative writes are contained: each snippet runs in a temp cwd
"""

from __future__ import annotations

import importlib.util
import re
import runpy
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SNIPPET_DIR = ROOT / "docs" / "source" / "snippets"
SNIPPETS = sorted(SNIPPET_DIR.glob("*.py"))
REQUIRES = {"viz": ("matplotlib", "altair"), "openeye": ("openeye",)}


def _requires(path: Path) -> list[str]:
    for line in path.read_text().splitlines()[:10]:
        if line.startswith("# requires:"):
            tags = [t.strip() for t in line.split(":", 1)[1].split(",")]
            return [t for t in tags if t and t != "none"]
    return []


def test_snippets_exist():
    assert SNIPPETS, f"no snippets found under {SNIPPET_DIR}"


@pytest.mark.parametrize("snippet", SNIPPETS, ids=[p.stem for p in SNIPPETS])
def test_snippet_runs(snippet: Path, tmp_path, monkeypatch):
    for tag in _requires(snippet):
        for module in REQUIRES[tag]:
            pytest.importorskip(module)
    if importlib.util.find_spec("matplotlib"):
        import matplotlib

        matplotlib.use("Agg")  # plt.show() becomes a no-op
    monkeypatch.chdir(tmp_path)  # any relative write lands here
    runpy.run_path(str(snippet), run_name="__main__")


def _region(path: Path, tag: str) -> str:
    text = path.read_text()
    start = text.index(f"# [start:{tag}]") + len(f"# [start:{tag}]\n")
    end = text.index(f"# [end:{tag}]")
    body = text[start:end]
    # dedent one level (the region sits inside main())
    return "\n".join(line[4:] if line.startswith("    ") else line for line in body.splitlines()).strip("\n")


def test_readme_quickstart_matches_snippet():
    """The README quickstart is the [start:search] region of search_minimal.py."""
    readme = (ROOT / "README.md").read_text()
    m = re.search(r"## Quickstart.*?```python\n(.*?)```", readme, re.S)
    assert m, "README has no python block under ## Quickstart"
    assert m.group(1).strip("\n") == _region(SNIPPET_DIR / "search_minimal.py", "search")
