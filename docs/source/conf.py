"""Sphinx configuration for the TACTICS documentation.

The build imports the installed package (``pip install -e ".[docs]"``) so the
API reference is generated from the code. Nothing is mocked: the package's
optional dependencies are either lazy (OpenEye) or included in the ``docs``
extra (matplotlib, altair).
"""

from __future__ import annotations

from importlib import metadata

try:
    import TACTICS  # noqa: F401  -- verifies the package is importable for autodoc
except ImportError as exc:  # pragma: no cover
    raise RuntimeError(
        "TACTICS must be installed to build the docs: pip install -e '.[docs]'"
    ) from exc

# -- Project information -----------------------------------------------------

project = "TACTICS"
copyright = "2024-2026, Aakankschit Nandkeolyar"
author = "Aakankschit Nandkeolyar"
release = metadata.version("chem-tactics")
version = ".".join(release.split(".")[:2])

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.githubpages",
    "sphinx.ext.graphviz",
    "sphinx_design",
    "sphinx_copybutton",
    "sphinx_togglebutton",
    "myst_nb",  # loads myst_parser itself; do not list both
]

exclude_patterns = ["_build", "snippets/README*"]

graphviz_output_format = "svg"

# MyST (theory pages are Markdown with $...$ / $$...$$ math)
nb_execution_mode = "off"
myst_enable_extensions = ["colon_fence", "deflist", "dollarmath", "amsmath"]
myst_heading_anchors = 4

# Autodoc
autodoc_mock_imports: list[str] = []  # the real package is installed; mocks would shadow it
autodoc_typehints = "description"
autodoc_typehints_description_target = "documented_params"
autodoc_typehints_format = "short"
# "mixed": signature on the class line; __init__ docstrings are merged into
# the class doc by napoleon_include_init_with_doc. "separated" would add
# __init__/__new__ entries that bypass autodoc-skip-member.
autodoc_class_signature = "mixed"
autodoc_member_order = "bysource"
# No default `members`: a bool default silently overwrites explicit
# `:members: a, b` lists in directives (sphinx.ext.autodoc.directive.
# process_documenter_options). Every directive names its members.
autodoc_default_options: dict = {}

# Napoleon (Google + NumPy docstrings)
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_use_param = True
napoleon_use_rtype = True
# Render 'Attributes:' sections as :ivar: fields, not `.. attribute::` object
# descriptions -- the latter collide with real @property members.
napoleon_use_ivar = True

# Copy button: strip prompts
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True

# Intersphinx
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "polars": ("https://docs.pola.rs/api/python/stable/", None),
}
intersphinx_timeout = 30

# Linkcheck (run as a separate, non-blocking CI job)
linkcheck_anchors_ignore_for_url = [r"https://github\.com/.*"]
linkcheck_timeout = 15
linkcheck_retries = 2

# -- HTML output -------------------------------------------------------------

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
html_title = f"TACTICS {release}"

html_theme_options = {
    "show_toc_level": 2,
    "navigation_depth": 3,
    "show_nav_level": 2,
    "navbar_align": "left",
    "navbar_start": ["navbar-logo"],
    "navbar_center": ["navbar-nav"],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "secondary_sidebar_items": ["page-toc", "edit-this-page"],
    "use_edit_page_button": True,
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version"],
    "pygments_light_style": "default",
    "pygments_dark_style": "monokai",
    "logo": {
        "image_light": "_static/images/TACTICS_logo.png",
        "image_dark": "_static/images/TACTICS_logo.png",
        "text": "TACTICS",
        "alt_text": "TACTICS - Home",
    },
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/aakankschit/TACTICS",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/chem-tactics/",
            "icon": "fa-brands fa-python",
        },
    ],
}

html_context = {
    "default_mode": "auto",
    "github_user": "aakankschit",
    "github_repo": "TACTICS",
    "github_version": "main",
    "doc_path": "docs/source",
}

# -- Pydantic models under autodoc --------------------------------------------
#
# Plain autodoc shows a Pydantic model's class docstring and signature but
# never reads Field(description=...). These two hooks render a "Fields" block
# from ``model_fields`` (so descriptions and defaults come from the code) and
# hide Pydantic's own machinery (``model_*``, validators, the generated
# ``__init__``) from the member list.

_DISCRIMINATORS = {"strategy_type", "warmup_type", "evaluator_type"}


def _is_pydantic_model(obj) -> bool:
    try:
        from pydantic import BaseModel
    except ImportError:  # pragma: no cover
        return False
    return isinstance(obj, type) and issubclass(obj, BaseModel) and obj is not BaseModel


def _describe_default(field) -> str:
    from pydantic_core import PydanticUndefined

    if field.is_required():
        return "required"
    if field.default_factory is not None:
        return "computed"
    if field.default is PydanticUndefined:
        return "required"
    return f"``{field.default!r}``"


def _pydantic_fields(app, what, name, obj, options, lines):
    if what != "class" or not _is_pydantic_model(obj):
        return
    from sphinx.util.typing import stringify_annotation

    # Napoleon has already turned any "Attributes:" section into :ivar:/:vartype:
    # fields; the generated Fields block below supersedes them. It has also
    # appended Pydantic's boilerplate __init__ docstring -- drop from there on.
    for i, ln in enumerate(lines):
        if ln.lstrip().startswith("Create a new model by parsing"):
            del lines[i:]
            break
    lines[:] = [ln for ln in lines if not ln.lstrip().startswith((":ivar ", ":vartype "))]

    lines += ["", ".. rubric:: Fields", ""]
    for fname, field in obj.model_fields.items():
        if fname in _DISCRIMINATORS:
            continue
        desc = (field.description or "").strip()
        if desc and not desc.endswith("."):
            desc += "."
        lines.append(f":param {fname}: {desc} Default: {_describe_default(field)}.")
        lines.append(f":type {fname}: {stringify_annotation(field.annotation, 'smart')}")
    if obj.model_config.get("extra") == "forbid":
        lines += ["", "Unknown keyword arguments raise :class:`pydantic.ValidationError`.", ""]


def _skip_pydantic_internals(app, what, name, obj, skip, options):
    if name.startswith("model_"):
        return True
    # Field/model validators are classmethods; their bound __self__ is the
    # model, whose __pydantic_decorators__ lists them by name.
    owner = getattr(obj, "__self__", None)
    if _is_pydantic_model(owner):
        decorators = owner.__pydantic_decorators__
        if name in decorators.field_validators or name in decorators.model_validators:
            return True
    return None


def setup(app):
    app.connect("autodoc-process-docstring", _pydantic_fields)
    app.connect("autodoc-skip-member", _skip_pydantic_internals)
